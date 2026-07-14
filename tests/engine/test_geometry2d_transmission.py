from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.array([0.0, 10_000.0, 30_000.0])
HORIZONTAL_ANGLES = np.array([-0.5, 0.0, 0.5])
WAVELENGTHS_NM = np.array([500.0, 600.0])


def transmission_config(*, output_optical_depth: bool = False) -> sk.Config:
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    config.occultation_source = sk.OccultationSource.Standard
    config.output_los_optical_depth = output_optical_depth
    return config


def geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=HORIZONTAL_ANGLES,
    )


def tangent_ray(tangent_altitude_m: float = 15_000.0):
    return sk.TangentAltitudeSolar(
        tangent_altitude_m=tangent_altitude_m,
        relative_azimuth=0.0,
        observer_altitude_m=100_000.0,
        cos_sza=0.6,
    )


def test_constant_extinction_matches_analytic_chord_transmission():
    config = transmission_config(output_optical_depth=True)
    config.num_threads = 2
    geometry = geometry2d()
    viewing = sk.ViewingGeometry()
    tangent_altitude_m = 15_000.0
    viewing.add_ray(tangent_ray(tangent_altitude_m))
    extinction = np.array([1.0e-6, 3.0e-6])
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    atmosphere["extinction"] = sk.constituent.Manual(
        np.broadcast_to(extinction, (*geometry.shape, extinction.size)).copy(),
        np.zeros((*geometry.shape, extinction.size)),
    )

    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    top_radius = EARTH_RADIUS_M + ALTITUDES_M[-1]
    tangent_radius = EARTH_RADIUS_M + tangent_altitude_m
    path_length_m = 2.0 * np.sqrt(top_radius**2 - tangent_radius**2)
    expected_od = extinction * path_length_m
    np.testing.assert_allclose(
        result.radiance[:, 0, 0], np.exp(-expected_od), rtol=2.0e-12
    )
    np.testing.assert_allclose(
        result.los_optical_depth[:, 0], expected_od, rtol=2.0e-12
    )


def test_geometry_relative_tangent_horizontal_angle_selects_expected_profile():
    config = transmission_config(output_optical_depth=True)
    geometry = geometry2d()
    viewing = sk.ViewingGeometry()
    tangent_altitude_m = 15_000.0
    requested_angles = np.array([-0.2, 0.2])
    for horizontal_angle in requested_angles:
        viewing.add_ray(
            sk.TangentAltitude(
                tangent_altitude_m=tangent_altitude_m,
                observer_altitude_m=100_000.0,
                horizontal_angle_radians=horizontal_angle,
                viewing_azimuth_radians=np.pi / 2.0,
            )
        )

    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        calculate_derivatives=False,
    )
    horizontal_extinction = np.array([1.0e-6, 2.0e-6, 4.0e-6])
    atmosphere.storage.total_extinction[:] = np.broadcast_to(
        horizontal_extinction[:, np.newaxis, np.newaxis], (*geometry.shape, 1)
    ).reshape(geometry.shape[0] * geometry.shape[1], 1)
    atmosphere.storage.ssa[:] = 0.0

    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    top_radius = EARTH_RADIUS_M + ALTITUDES_M[-1]
    tangent_radius = EARTH_RADIUS_M + tangent_altitude_m
    path_length_m = 2.0 * np.sqrt(top_radius**2 - tangent_radius**2)
    sampled_extinction = np.interp(
        requested_angles, HORIZONTAL_ANGLES, horizontal_extinction
    )
    expected_od = sampled_extinction * path_length_m
    np.testing.assert_allclose(result.los_optical_depth[0], expected_od, rtol=2.0e-12)
    np.testing.assert_allclose(
        result.radiance[0, :, 0], np.exp(-expected_od), rtol=2.0e-12
    )


def test_geometry_relative_tangent_ray_does_not_move_with_solar_angles():
    config = transmission_config(output_optical_depth=True)
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.TangentAltitude(
            tangent_altitude_m=15_000.0,
            observer_altitude_m=100_000.0,
            horizontal_angle_radians=0.17,
            viewing_azimuth_radians=0.35,
        )
    )
    extinction = np.arange(1.0, 10.0).reshape(3, 3, 1) * 1.0e-7

    optical_depths = []
    for cos_sza, solar_azimuth in [(0.8, 0.0), (-0.2, 1.4)]:
        geometry = sk.Geometry2D(
            cos_sza=cos_sza,
            solar_azimuth=solar_azimuth,
            earth_radius_m=EARTH_RADIUS_M,
            altitude_grid_m=ALTITUDES_M,
            horizontal_angle_grid_radians=HORIZONTAL_ANGLES,
        )
        atmosphere = sk.Atmosphere(
            geometry,
            config,
            wavelengths_nm=np.array([500.0]),
            calculate_derivatives=False,
        )
        atmosphere.storage.total_extinction[:] = extinction.reshape(9, 1)
        atmosphere.storage.ssa[:] = 0.0
        result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
        optical_depths.append(result.los_optical_depth.copy())

    np.testing.assert_allclose(optical_depths[0], optical_depths[1], rtol=0.0, atol=0.0)


def test_transparent_limb_ray_is_unity_and_ground_ray_is_blocked():
    config = transmission_config()
    geometry = geometry2d()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(tangent_ray())
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.6,
            relative_azimuth=0.0,
            cos_viewing_zenith=0.5,
            observer_altitude_m=100_000.0,
        )
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    atmosphere["transparent"] = sk.constituent.Manual(
        np.zeros((*geometry.shape, WAVELENGTHS_NM.size)),
        np.zeros((*geometry.shape, WAVELENGTHS_NM.size)),
    )

    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    np.testing.assert_array_equal(result.radiance[:, 0, 0], 1.0)
    np.testing.assert_array_equal(result.radiance[:, 1, 0], 0.0)


def test_horizontally_uniform_2d_transmission_matches_1d():
    config = transmission_config()
    geometry_1d = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    geometry_2d = geometry2d()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(tangent_ray(12_000.0))
    viewing.add_ray(tangent_ray(22_000.0))
    extinction = np.array([[1.0e-5, 2.0e-5], [2.0e-5, 1.0e-5], [0.5e-5, 0.7e-5]])
    atmosphere_1d = sk.Atmosphere(
        geometry_1d,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    atmosphere_2d = sk.Atmosphere(
        geometry_2d,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    atmosphere_1d["extinction"] = sk.constituent.Manual(
        extinction, np.zeros_like(extinction)
    )
    atmosphere_2d["extinction"] = sk.constituent.Manual(
        extinction, np.zeros_like(extinction)
    )

    result_1d = sk.Engine(config, geometry_1d, viewing).calculate_radiance(
        atmosphere_1d
    )
    result_2d = sk.Engine(config, geometry_2d, viewing).calculate_radiance(
        atmosphere_2d
    )

    np.testing.assert_allclose(result_2d.radiance, result_1d.radiance, rtol=5.0e-12)


def test_native_2d_extinction_derivative_matches_finite_difference():
    config = transmission_config()
    geometry = geometry2d()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(tangent_ray())
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        legendre_derivative=False,
    )
    extinction = np.array(
        [
            [[0.8e-5, 1.2e-5], [1.0e-5, 1.4e-5], [0.7e-5, 0.9e-5]],
            [[1.1e-5, 1.3e-5], [1.4e-5, 1.0e-5], [0.6e-5, 0.8e-5]],
            [[1.5e-5, 0.7e-5], [1.2e-5, 0.9e-5], [0.5e-5, 1.1e-5]],
        ]
    )
    atmosphere.storage.total_extinction[:] = extinction.reshape(9, 2)
    atmosphere.storage.ssa[:] = 0.0
    engine = sk.Engine(config, geometry, viewing)

    base = engine.calculate_radiance(atmosphere)
    assert base.wf_extinction.dims == (
        "horizontal_angle",
        "altitude",
        "wavelength",
        "los",
        "stokes",
    )
    assert base.wf_extinction.shape == (3, 3, 2, 1, 1)

    horizontal_index = 1
    altitude_index = 1
    wavelength_index = 0
    location_index = geometry.location_index(altitude_index, horizontal_index)
    delta = 1.0e-10
    atmosphere.storage.total_extinction[location_index, wavelength_index] += delta
    perturbed = engine.calculate_radiance(atmosphere)

    numeric = (
        perturbed.radiance[wavelength_index, 0, 0]
        - base.radiance[wavelength_index, 0, 0]
    ) / delta
    analytic = base.wf_extinction[
        horizontal_index, altitude_index, wavelength_index, 0, 0
    ]
    np.testing.assert_allclose(numeric, analytic, rtol=4.0e-5)
    np.testing.assert_array_equal(perturbed.radiance[1, 0, 0], base.radiance[1, 0, 0])


def test_vector_transmission_has_only_intensity_component():
    config = transmission_config()
    config.num_stokes = 3
    geometry = geometry2d()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(tangent_ray())
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    atmosphere.storage.total_extinction[:] = 1.0e-6
    atmosphere.storage.ssa[:] = 0.0

    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    assert np.all(result.radiance[:, 0, 0] > 0.0)
    np.testing.assert_array_equal(result.radiance[:, 0, 1:], 0.0)


def test_geometry2d_engine_rejects_unsupported_sources_and_mismatched_geometry():
    geometry = geometry2d()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(tangent_ray())
    unsupported = sk.Config()
    unsupported.single_scatter_source = sk.SingleScatterSource.Table
    with pytest.raises(NotImplementedError, match="supports exact single scattering"):
        sk.Engine(unsupported, geometry, viewing)

    config = transmission_config()
    engine = sk.Engine(config, geometry, viewing)
    other_geometry = geometry2d()
    atmosphere = sk.Atmosphere(
        other_geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    with pytest.raises(ValueError, match="same Geometry2D object"):
        engine.calculate_radiance(atmosphere)


def test_geometry2d_engine_rejects_flux_and_unconfigured_refraction():
    geometry = geometry2d()
    config = transmission_config()
    refracted = transmission_config()
    refracted.los_refraction = True
    viewing = sk.ViewingGeometry()
    viewing.add_ray(tangent_ray())

    with pytest.raises(NotImplementedError, match="refractive-index profiles"):
        sk.Engine(refracted, geometry, viewing)

    viewing.add_flux_observer(
        sk.FluxObserverSolar(cos_sza=0.6, observer_altitude_m=10_000.0)
    )
    with pytest.raises(NotImplementedError, match="flux observers"):
        sk.Engine(config, geometry, viewing)
