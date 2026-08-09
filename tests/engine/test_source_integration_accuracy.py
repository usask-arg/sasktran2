from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _grid_linear_single_scatter(
    spacing_m: float,
    solar_transmission: sk.SingleScatterSolarTransmission,
) -> np.ndarray:
    top_of_atmosphere_m = 60_000.0
    altitude_m = np.arange(0.0, top_of_atmosphere_m + spacing_m, spacing_m)
    config = sk.Config()
    config.num_threads = 1
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.single_scatter_source_quadrature = True
    config.single_scatter_solar_transmission = solar_transmission
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource

    geometry = sk.Geometry1D(
        cos_sza=0.2,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitude_m,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(12_345.0, -0.4, 200_000.0, 0.2))
    viewing.add_ray(sk.GroundViewingSolar(0.2, 0.35, 0.45, 200_000.0))
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        calculate_derivatives=False,
        numwavel=1,
    )
    atmosphere.storage.total_extinction[:, 0] = (
        (1.5 - 0.5 * altitude_m / top_of_atmosphere_m)
        * 3.0
        / (1.25 * top_of_atmosphere_m)
    )
    atmosphere.storage.ssa[:] = 0.88
    atmosphere.leg_coeff.a1[0, :, :] = 1.0
    atmosphere.leg_coeff.a1[1, :, :] = 0.08
    atmosphere.leg_coeff.a1[2, :, :] = 0.45
    atmosphere.surface.albedo[:] = 0.0

    return (
        sk.Engine(config, geometry, viewing)
        .calculate_radiance(atmosphere)
        .radiance.values.squeeze()
    )


def test_exact_single_scatter_spherical_limb_and_ground_integration_accuracy():
    coarse = _grid_linear_single_scatter(
        2_000.0, sk.SingleScatterSolarTransmission.Exact
    )
    dense = _grid_linear_single_scatter(125.0, sk.SingleScatterSolarTransmission.Exact)
    relative_error = np.abs((coarse - dense) / dense)

    assert relative_error[0] < 6e-3  # 12.345 km limb tangent
    assert relative_error[1] < 5e-4  # ground-looking ray


def test_ray_table_single_scatter_is_checked_against_exact_solar_tracing():
    ray_table = _grid_linear_single_scatter(
        500.0, sk.SingleScatterSolarTransmission.RayTable
    )
    exact = _grid_linear_single_scatter(500.0, sk.SingleScatterSolarTransmission.Exact)

    assert np.allclose(ray_table, exact, rtol=4e-4, atol=1e-12)
