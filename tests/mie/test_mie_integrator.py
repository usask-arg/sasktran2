from __future__ import annotations

import numpy as np
import pytest
import xarray as xr
from sasktran2._core_rust import PyMie, PyMieIntegrator, WignerD
from sasktran2.mie.distribution import integrate_mie_cpp
from scipy.special import roots_legendre
from scipy.stats import lognorm


def _outputs(distributions, angles, orders, layout="C"):
    shapes = (
        [(distributions,)] * 2
        + [(distributions, angles)] * 4
        + [(distributions, orders)] * 6
    )
    if layout == "strided":
        return [np.zeros((*shape[:-1], shape[-1] * 2))[..., ::2] for shape in shapes]
    return [np.zeros(shape, order=layout) for shape in shapes]


@pytest.mark.parametrize("threads", [1, 2, 4])
@pytest.mark.parametrize("layout", ["C", "F", "strided"])
@pytest.mark.parametrize("num_distributions", [1, 3])
def test_mie_integrator_matches_direct_quadrature(threads, layout, num_distributions):
    cos_angles, angle_weights = roots_legendre(64)
    size = np.geomspace(0.001, 20, 48)
    size_weights = np.linspace(0.1, 1, len(size))
    pdf = np.array([np.exp(-size), np.exp(-size / 3), np.exp(-size / 8)], order="F")
    pdf = pdf[:num_distributions]
    wavelength = 550.0
    refractive_index = 1.5 - 0.01j
    orders = 24
    inputs = (wavelength, refractive_index, size, pdf, size_weights, angle_weights)
    result = _outputs(len(pdf), len(cos_angles), orders, layout)
    PyMieIntegrator(cos_angles, orders, threads).integrate(*inputs, *result)

    # Independently integrate the standalone Mie amplitudes in NumPy.
    single = PyMie().calculate(size, refractive_index, cos_angles, False)
    radius = size * wavelength / (2 * np.pi)
    weighted_pdf = pdf * size_weights
    xs_total = weighted_pdf @ (np.pi * radius**2 * single.Qext)
    xs_scattering = weighted_pdf @ (np.pi * radius**2 * single.Qsca)
    norm = 2 * np.pi / (2 * np.pi / wavelength) ** 2 / xs_scattering[:, None]
    s1, s2 = single.S1, single.S2
    p11 = norm * (weighted_pdf @ (np.abs(s1) ** 2 + np.abs(s2) ** 2))
    p12 = norm * (weighted_pdf @ (np.abs(s1) ** 2 - np.abs(s2) ** 2))
    p33 = norm * (weighted_pdf @ (2 * (s1 * s2.conj()).real))
    p34 = norm * (weighted_pdf @ (2 * (s1 * s2.conj()).imag))
    expected = [xs_total, xs_scattering, p11, p12, p33, p34]
    coefficients = [np.zeros((len(pdf), orders)) for _ in range(6)]
    calculators = [WignerD(m, n) for m, n in [(0, 0), (2, 2), (2, -2), (0, 2)]]
    for order in range(orders):
        d00, d22, d2m2, d02 = (
            w.d(np.arccos(cos_angles), order) * angle_weights * (2 * order + 1) / 2
            for w in calculators
        )
        coefficients[0][:, order] = p11 @ d00
        coefficients[3][:, order] = p33 @ d00
        coefficients[4][:, order] = p12 @ d02
        coefficients[5][:, order] = -(p34 @ d02)
        plus, minus = (p11 + p33) @ d22, (p11 - p33) @ d2m2
        coefficients[1][:, order] = (plus + minus) / 2
        coefficients[2][:, order] = (plus - minus) / 2
    for actual, reference in zip(result, expected + coefficients, strict=True):
        np.testing.assert_allclose(actual, reference, atol=3e-12, rtol=3e-12)

    # The private pool must not change results or leak settings into other instances.
    serial = _outputs(len(pdf), len(cos_angles), orders, layout)
    PyMieIntegrator(cos_angles, orders, 1).integrate(*inputs, *serial)
    for actual, reference in zip(result, serial, strict=True):
        np.testing.assert_array_equal(actual, reference)


def test_mie_integrator_rejects_mismatched_outputs():
    cos_angles, weights = roots_legendre(8)
    integrator = PyMieIntegrator(cos_angles, 4, 2)
    output = _outputs(1, 8, 4)
    output[-1] = np.zeros((1, 3))
    with pytest.raises(ValueError, match="coefficient outputs"):
        integrator.integrate(
            500.0,
            1.5 - 0.01j,
            np.ones(1),
            np.ones((1, 1)),
            np.ones(1),
            weights,
            *output,
        )


@pytest.mark.parametrize("threads", [0, 2])
def test_mie_distribution_threading(threads):
    distributions = [lognorm(s=0.3, scale=80), lognorm(s=0.4, scale=120)]
    inputs = (distributions, lambda _: 1.5 - 0.01j, np.array([500.0, 700.0]))
    serial = integrate_mie_cpp(*inputs, num_coeffs=24, num_quad=8, num_threads=1)
    threaded = integrate_mie_cpp(
        *inputs, num_coeffs=24, num_quad=8, num_threads=threads
    )
    xr.testing.assert_identical(serial, threaded)
