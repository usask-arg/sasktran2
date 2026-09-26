from __future__ import annotations

from decimal import Decimal, localcontext

import numpy as np
import pytest
import sasktran2 as sk


def _cell_integral_reference(optical_depth):
    # Evaluate the analytic slab solution independently at high precision,
    # including its extinction derivative arbitrarily close to vacuum.
    integral = np.empty_like(optical_depth)
    derivative = np.empty_like(optical_depth)
    with localcontext() as context:
        context.prec = 60
        for index in np.ndindex(optical_depth.shape):
            depth = Decimal(str(optical_depth[index]))
            if depth == 0:
                integral[index] = 1.0
                derivative[index] = -0.5
            else:
                transmission = (-depth).exp()
                integral[index] = float((1 - transmission) / depth)
                derivative[index] = float(((1 + depth) * transmission - 1) / depth**2)
    return integral, derivative


@pytest.mark.parametrize("num_layers", [1, 8])
@pytest.mark.parametrize(
    "geometry_type", [sk.GeometryType.PlaneParallel, sk.GeometryType.Spherical]
)
@pytest.mark.parametrize(
    ("num_stokes", "batch_size", "derivatives"),
    [(1, 1, True), (3, 4, True), (1, 4, False)],
)
def test_volume_emission_absorption_matches_uniform_slab(
    num_layers, geometry_type, num_stokes, batch_size, derivatives
):
    config = sk.Config()
    config.num_stokes = num_stokes
    config.num_threads = 2
    config.wavelength_batch_size = batch_size
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.emission_source = sk.EmissionSource.VolumeEmissionRate

    radius, height = 6_372_000.0, 1000.0
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        radius,
        np.linspace(0.0, height, num_layers + 1),
        sk.InterpolationMethod.LinearInterpolation,
        geometry_type,
    )
    cosines = np.array([1.0, 0.4])
    viewing = sk.ViewingGeometry()
    for cosine in cosines:
        viewing.add_ray(sk.GroundViewingSolar(0.6, 0.0, cosine, 200_000.0))

    # Span vacuum, cancellation-sensitive depths, the Taylor transition,
    # and opaque cells, with a partial final wavelength batch.
    vertical_depth = np.array([0.0, 1e-10, 1e-5, 9e-4, 1.1e-3, 0.1, 1.0, 10.0, 1000.0])
    emission = 1e-3 * np.linspace(0.7, 1.3, vertical_depth.size)
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        numwavel=vertical_depth.size,
        calculate_derivatives=derivatives,
        legendre_derivative=False,
    )
    atmosphere.storage.total_extinction[:] = vertical_depth / height
    atmosphere.storage.ssa[:] = 0.0
    atmosphere.storage.emission_source[:] = emission
    atmosphere.storage.solar_irradiance[:] = 0.0
    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    if geometry_type == sk.GeometryType.PlaneParallel:
        distance = height / cosines
    else:
        # Chord length between the ground and top of the spherical shell.
        radial_difference = height * (2 * radius + height)
        distance = radial_difference / (
            np.sqrt((radius * cosines) ** 2 + radial_difference) + radius * cosines
        )
    optical_depth = vertical_depth[:, None] * distance[None, :] / height
    integral, derivative = _cell_integral_reference(optical_depth)
    expected = emission[:, None] * distance * integral
    np.testing.assert_allclose(result.radiance[..., 0], expected, rtol=2e-10)
    np.testing.assert_allclose(result.radiance[..., 1:], 0.0, atol=1e-14)
    if derivatives:
        np.testing.assert_allclose(
            result.wf_emission.sum("altitude")[..., 0],
            distance * integral,
            rtol=2e-10,
        )
        np.testing.assert_allclose(
            result.wf_extinction.sum("altitude")[..., 0],
            emission[:, None] * distance**2 * derivative,
            rtol=2e-10,
        )
        np.testing.assert_allclose(result.wf_emission[..., 1:], 0.0, atol=1e-14)
        np.testing.assert_allclose(result.wf_extinction[..., 1:], 0.0, atol=1e-14)
