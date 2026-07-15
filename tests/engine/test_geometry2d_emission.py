from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.array([0.0, 10_000.0, 30_000.0])
HORIZONTAL_ANGLES = np.array([-0.5, 0.0, 0.5])
WAVELENGTHS_NM = np.array([8_000.0, 10_000.0])


def emission_config(source: sk.EmissionSource) -> sk.Config:
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = source
    return config


def geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=HORIZONTAL_ANGLES,
    )


def geometry1d() -> sk.Geometry1D:
    return sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )


def tangent_ray(tangent_altitude_m: float = 15_000.0):
    return sk.TangentAltitudeSolar(
        tangent_altitude_m=tangent_altitude_m,
        relative_azimuth=0.0,
        observer_altitude_m=100_000.0,
        cos_sza=0.6,
    )


def tangent_viewing_geometry() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(tangent_ray())
    return viewing


def tangent_path_length_m(tangent_altitude_m: float = 15_000.0) -> float:
    top_radius = EARTH_RADIUS_M + ALTITUDES_M[-1]
    tangent_radius = EARTH_RADIUS_M + tangent_altitude_m
    return 2.0 * np.sqrt(top_radius**2 - tangent_radius**2)


def test_constant_volume_emission_matches_analytic_path_integral():
    config = emission_config(sk.EmissionSource.VolumeEmissionRate)
    config.num_threads = 2
    geometry = geometry2d()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    emission = np.array([2.0e-6, 3.0e-6])
    atmosphere.storage.total_extinction[:] = 0.0
    atmosphere.storage.ssa[:] = 0.0
    atmosphere.storage.emission_source[:] = emission

    result = sk.Engine(config, geometry, tangent_viewing_geometry()).calculate_radiance(
        atmosphere
    )

    np.testing.assert_allclose(
        result.radiance[:, 0, 0],
        emission * tangent_path_length_m(),
        rtol=2.0e-12,
    )


def test_standard_emission_and_native_derivatives_match_analytic_and_numeric():
    config = emission_config(sk.EmissionSource.Standard)
    geometry = geometry2d()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        legendre_derivative=False,
    )
    extinction = np.array([1.0e-6, 3.0e-6])
    emission = np.array([2.0, 3.0])
    atmosphere.storage.total_extinction[:] = extinction
    atmosphere.storage.ssa[:] = 0.0
    atmosphere.storage.emission_source[:] = emission
    engine = sk.Engine(config, geometry, tangent_viewing_geometry())

    base = engine.calculate_radiance(atmosphere)
    path_length = tangent_path_length_m()
    np.testing.assert_allclose(
        base.radiance[:, 0, 0],
        emission * (1.0 - np.exp(-extinction * path_length)),
        rtol=2.0e-12,
    )
    assert base.wf_extinction.shape == (3, 3, 2, 1, 1)
    assert base.wf_emission.shape == (3, 3, 2, 1, 1)

    horizontal_index = 1
    altitude_index = 1
    wavelength_index = 0
    location_index = geometry.location_index(altitude_index, horizontal_index)

    extinction_delta = 1.0e-10
    atmosphere.storage.total_extinction[
        location_index, wavelength_index
    ] += extinction_delta
    perturbed_extinction = engine.calculate_radiance(atmosphere)
    numeric_extinction = (
        perturbed_extinction.radiance[wavelength_index, 0, 0]
        - base.radiance[wavelength_index, 0, 0]
    ) / extinction_delta
    analytic_extinction = base.wf_extinction[
        horizontal_index, altitude_index, wavelength_index, 0, 0
    ]
    np.testing.assert_allclose(numeric_extinction, analytic_extinction, rtol=5.0e-5)

    atmosphere.storage.total_extinction[
        location_index, wavelength_index
    ] -= extinction_delta
    emission_delta = 1.0e-6
    atmosphere.storage.emission_source[
        location_index, wavelength_index
    ] += emission_delta
    perturbed_emission = engine.calculate_radiance(atmosphere)
    numeric_emission = (
        perturbed_emission.radiance[wavelength_index, 0, 0]
        - base.radiance[wavelength_index, 0, 0]
    ) / emission_delta
    analytic_emission = base.wf_emission[
        horizontal_index, altitude_index, wavelength_index, 0, 0
    ]
    np.testing.assert_allclose(numeric_emission, analytic_emission, rtol=2.0e-7)


def test_occultation_and_standard_emission_sources_add_consistently():
    config = emission_config(sk.EmissionSource.Standard)
    config.occultation_source = sk.OccultationSource.Standard
    geometry = geometry2d()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    extinction = np.array([1.0e-6, 3.0e-6])
    emission = np.array([2.0, 3.0])
    atmosphere.storage.total_extinction[:] = extinction
    atmosphere.storage.ssa[:] = 0.0
    atmosphere.storage.emission_source[:] = emission

    result = sk.Engine(config, geometry, tangent_viewing_geometry()).calculate_radiance(
        atmosphere
    )

    transmission = np.exp(-extinction * tangent_path_length_m())
    expected = transmission + emission * (1.0 - transmission)
    np.testing.assert_allclose(result.radiance[:, 0, 0], expected, rtol=2.0e-12)


@pytest.mark.parametrize(
    "source",
    [sk.EmissionSource.Standard, sk.EmissionSource.VolumeEmissionRate],
)
def test_horizontally_uniform_2d_emission_matches_1d(source):
    config = emission_config(source)
    geometry_1d = geometry1d()
    geometry_2d = geometry2d()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(tangent_ray(12_000.0))
    viewing.add_ray(tangent_ray(22_000.0))
    extinction = np.array([[1.0e-6, 2.0e-6], [2.0e-6, 1.0e-6], [0.5e-6, 0.7e-6]])
    emission = np.array([[1.0, 4.0], [2.0, 3.0], [5.0, 2.0]])
    ssa = (
        np.array([[0.1, 0.2], [0.3, 0.1], [0.2, 0.4]])
        if source == sk.EmissionSource.Standard
        else np.zeros_like(extinction)
    )
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
    atmosphere_1d.storage.total_extinction[:] = extinction
    atmosphere_1d.storage.ssa[:] = ssa
    atmosphere_1d.storage.emission_source[:] = emission
    atmosphere_2d.storage.total_extinction[:] = np.tile(extinction, (3, 1))
    atmosphere_2d.storage.ssa[:] = np.tile(ssa, (3, 1))
    atmosphere_2d.storage.emission_source[:] = np.tile(emission, (3, 1))

    result_1d = sk.Engine(config, geometry_1d, viewing).calculate_radiance(
        atmosphere_1d
    )
    result_2d = sk.Engine(config, geometry_2d, viewing).calculate_radiance(
        atmosphere_2d
    )

    np.testing.assert_allclose(result_2d.radiance, result_1d.radiance, rtol=5.0e-12)


def test_native_2d_thermal_emission_temperature_derivative():
    config = emission_config(sk.EmissionSource.Standard)
    geometry = geometry2d()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        legendre_derivative=False,
    )
    temperature = np.array(
        [[220.0, 230.0, 240.0], [225.0, 235.0, 245.0], [230.0, 240.0, 250.0]]
    )
    atmosphere.temperature_k = temperature
    extinction = np.full((*geometry.shape, WAVELENGTHS_NM.size), 1.0e-6)
    atmosphere["absorption"] = sk.constituent.Manual(
        extinction, np.zeros_like(extinction)
    )
    atmosphere["thermal"] = sk.constituent.ThermalEmission()
    engine = sk.Engine(config, geometry, tangent_viewing_geometry())

    base = engine.calculate_radiance(atmosphere)
    assert base.wf_temperature_k.shape == (3, 3, 2, 1, 1)
    horizontal_index = 1
    altitude_index = 1
    wavelength_index = 0
    delta = 1.0e-3
    atmosphere.temperature_k[horizontal_index, altitude_index] += delta
    perturbed = engine.calculate_radiance(atmosphere)

    numeric = (
        perturbed.radiance[wavelength_index, 0, 0]
        - base.radiance[wavelength_index, 0, 0]
    ) / delta
    analytic = base.wf_temperature_k[
        horizontal_index, altitude_index, wavelength_index, 0, 0
    ]
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-5)


def test_altitude_profile_volume_emission_constituent_matches_1d_and_numeric_wf():
    config = emission_config(sk.EmissionSource.VolumeEmissionRate)
    geometry_1d = geometry1d()
    geometry_2d = geometry2d()
    viewing = tangent_viewing_geometry()
    wavelength = np.array([557.7])
    ver = np.array([1.0, 3.0, 2.0])
    atmosphere_1d = sk.Atmosphere(geometry_1d, config, wavelengths_nm=wavelength)
    atmosphere_2d = sk.Atmosphere(geometry_2d, config, wavelengths_nm=wavelength)
    extinction = np.full((ALTITUDES_M.size, 1), 1.0e-6)
    atmosphere_1d["absorption"] = sk.constituent.Manual(
        extinction, np.zeros_like(extinction)
    )
    atmosphere_2d["absorption"] = sk.constituent.Manual(
        extinction, np.zeros_like(extinction)
    )
    atmosphere_1d["emitter"] = sk.constituent.MonochromaticVolumeEmissionRate(
        ALTITUDES_M, ver.copy(), wavelength[0], out_of_bounds_mode="extend"
    )
    atmosphere_2d["emitter"] = sk.constituent.MonochromaticVolumeEmissionRate(
        ALTITUDES_M, ver.copy(), wavelength[0], out_of_bounds_mode="extend"
    )
    engine_1d = sk.Engine(config, geometry_1d, viewing)
    engine_2d = sk.Engine(config, geometry_2d, viewing)

    result_1d = engine_1d.calculate_radiance(atmosphere_1d)
    result_2d = engine_2d.calculate_radiance(atmosphere_2d)
    np.testing.assert_allclose(result_2d.radiance, result_1d.radiance, rtol=5.0e-12)
    assert result_2d.wf_emitter_ver.dims[0] == "emitter_altitude"

    altitude_index = 1
    delta = 1.0e-5
    atmosphere_2d["emitter"].ver[altitude_index] += delta
    perturbed = engine_2d.calculate_radiance(atmosphere_2d)
    numeric = (perturbed.radiance[0, 0, 0] - result_2d.radiance[0, 0, 0]) / delta
    analytic = result_2d.wf_emitter_ver[altitude_index, 0, 0, 0]
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-6)


def test_ground_surface_emission_is_attenuated_and_matches_1d():
    config = emission_config(sk.EmissionSource.Standard)
    geometry_1d = geometry1d()
    geometry_2d = geometry2d()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.6,
            relative_azimuth=0.0,
            cos_viewing_zenith=0.5,
            observer_altitude_m=100_000.0,
        )
    )
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
    extinction = np.full((3, 2), 1.0e-6)
    surface_emission = np.array([2.0, 3.0])
    atmosphere_1d.storage.total_extinction[:] = extinction
    atmosphere_1d.storage.ssa[:] = 0.0
    atmosphere_1d.storage.emission_source[:] = 0.0
    atmosphere_1d.surface.emission[:] = surface_emission
    atmosphere_2d.storage.total_extinction[:] = np.tile(extinction, (3, 1))
    atmosphere_2d.storage.ssa[:] = 0.0
    atmosphere_2d.storage.emission_source[:] = 0.0
    atmosphere_2d.surface.emission[:] = surface_emission

    result_1d = sk.Engine(config, geometry_1d, viewing).calculate_radiance(
        atmosphere_1d
    )
    result_2d = sk.Engine(config, geometry_2d, viewing).calculate_radiance(
        atmosphere_2d
    )

    assert np.all(result_2d.radiance > 0.0)
    assert np.all(result_2d.radiance < surface_emission[:, np.newaxis, np.newaxis])
    np.testing.assert_allclose(result_2d.radiance, result_1d.radiance, rtol=5.0e-12)


def test_geometry2d_engine_rejects_unimplemented_emission_sources():
    config = emission_config(sk.EmissionSource.TwoStream)
    with pytest.raises(NotImplementedError, match="supports exact single scattering"):
        sk.Engine(config, geometry2d(), tangent_viewing_geometry())
