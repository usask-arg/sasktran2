from __future__ import annotations

import numpy as np
import pytest
from sasktran2._core_rust import WignerD
from sasktran2.legendre import compute_greek_coefficients
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre


@pytest.mark.parametrize("num_coeff", [1, 2, 32, 128])
@pytest.mark.parametrize("layout", ["C", "F", "strided"])
def test_greek_coefficients_match_original_transform(num_coeff, layout):
    angles = np.linspace(0, 180, 201)
    x = np.cos(np.deg2rad(angles))
    # Six independent elements, both signs, several wavelengths and a forward peak.
    rng = np.random.default_rng(123)
    phase = [
        rng.normal(size=(3, 1)) * (1 + x * x)
        + rng.normal(size=(3, 1)) * np.exp(30 * (x - 1))
        for _ in range(6)
    ]
    if layout == "F":
        phase = [np.asfortranarray(p) for p in phase]
    elif layout == "strided":
        phase = [np.repeat(p, 2, axis=1)[:, ::2] for p in phase]

    nodes, weights = roots_legendre(num_coeff)
    split = 0.995
    nodes = np.concatenate(
        [
            (split + 1) / 2 * nodes + (split - 1) / 2,
            (1 - split) / 2 * nodes + (1 + split) / 2,
        ]
    )
    weights = np.concatenate([(split + 1) / 2 * weights, (1 - split) / 2 * weights])
    interpolated = [
        PchipInterpolator(x[::-1], p[:, ::-1], axis=1)(nodes) for p in phase
    ]
    p11, p12, p22, p33, p34, p44 = interpolated
    expected = np.zeros((6, 3, num_coeff))
    calculators = [WignerD(m, n) for m, n in [(0, 0), (2, 2), (2, -2), (0, 2)]]
    for order in range(num_coeff):
        d00, d22, d2m2, d02 = (
            w.d(np.arccos(nodes), order) * weights * (2 * order + 1) / 2
            for w in calculators
        )
        expected[0, :, order] = p11 @ d00
        expected[3, :, order] = p44 @ d00
        expected[4, :, order] = p12 @ d02
        expected[5, :, order] = -(p34 @ d02)
        plus = (p22 + p33) @ d22
        minus = (p22 - p33) @ d2m2
        expected[1, :, order] = (plus + minus) / 2
        expected[2, :, order] = (plus - minus) / 2

    actual = compute_greek_coefficients(*phase, angles, num_coeff)
    np.testing.assert_allclose(actual, expected, atol=3e-12, rtol=1e-11)
