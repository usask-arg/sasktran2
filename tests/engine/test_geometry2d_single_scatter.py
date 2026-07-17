from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.array([0.0, 10_000.0, 30_000.0])
WAVELENGTHS_NM = np.array([500.0, 750.0])


def single_scatter_config(num_stokes: int = 1) -> sk.Config:
    config = sk.Config()
    config.num_stokes = num_stokes
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    config.num_singlescatter_moments = 6
    return config


def geometry2d(
    horizontal_angles: np.ndarray | None = None,
    interpolation_method: sk.InterpolationMethod = (
        sk.InterpolationMethod.LinearInterpolation
    ),
) -> sk.Geometry2D:
    if horizontal_angles is None:
        horizontal_angles = np.array([-1.5, 1.5])
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.2,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=horizontal_angles,
        interpolation_method=interpolation_method,
    )


def geometry1d(
    interpolation_method: sk.InterpolationMethod = (
        sk.InterpolationMethod.LinearInterpolation
    ),
) -> sk.Geometry1D:
    return sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.2,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=interpolation_method,
        geometry_type=sk.GeometryType.Spherical,
    )


def vertical_geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=1.0,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=np.array([-0.5, 0.5]),
    )


def vertical_geometry1d() -> sk.Geometry1D:
    return sk.Geometry1D(
        cos_sza=1.0,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )


def vertical_ground_viewing() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=1.0,
            relative_azimuth=0.0,
            cos_viewing_zenith=1.0,
            observer_altitude_m=100_000.0,
        )
    )
    return viewing


def parity_viewing() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.TangentAltitudeSolar(
            tangent_altitude_m=12_000.0,
            relative_azimuth=0.35,
            observer_altitude_m=100_000.0,
            cos_sza=0.6,
        )
    )
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.6,
            relative_azimuth=-0.25,
            cos_viewing_zenith=0.7,
            observer_altitude_m=100_000.0,
        )
    )
    return viewing


def ground_viewing() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.6,
            relative_azimuth=-0.25,
            cos_viewing_zenith=0.7,
            observer_altitude_m=100_000.0,
        )
    )
    return viewing


def set_phase(atmosphere: sk.Atmosphere) -> None:
    atmosphere.leg_coeff.a1[0] = 1.0
    atmosphere.leg_coeff.a1[2] = 0.5
    if atmosphere.nstokes == 3:
        atmosphere.leg_coeff.a2[2] = 3.0
        atmosphere.leg_coeff.b1[2] = np.sqrt(6.0) / 2.0


def test_vertical_isotropic_single_scatter_matches_analytic_solution():
    config = single_scatter_config()
    geometry = vertical_geometry2d()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    extinction = np.array([1.0e-5, 2.5e-5])
    ssa = np.array([0.4, 0.7])
    atmosphere.storage.total_extinction[:] = extinction
    atmosphere.storage.ssa[:] = ssa
    atmosphere.leg_coeff.a1[0] = 1.0

    result = sk.Engine(config, geometry, vertical_ground_viewing()).calculate_radiance(
        atmosphere
    )

    path_length = ALTITUDES_M[-1] - ALTITUDES_M[0]
    expected = ssa / (8.0 * np.pi) * (1.0 - np.exp(-2.0 * extinction * path_length))
    np.testing.assert_allclose(result.radiance[:, 0, 0], expected, rtol=2.0e-12)


def test_vertical_cell_corner_path_matches_analytic_solution():
    config = single_scatter_config()
    geometry = sk.Geometry2D(
        cos_sza=1.0,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=np.array([-0.5, 0.0, 0.5]),
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        calculate_derivatives=False,
    )
    extinction = 1.3e-5
    ssa = 0.65
    atmosphere.storage.total_extinction[:] = extinction
    atmosphere.storage.ssa[:] = ssa
    atmosphere.leg_coeff.a1[0] = 1.0

    result = sk.Engine(config, geometry, vertical_ground_viewing()).calculate_radiance(
        atmosphere
    )

    path_length = ALTITUDES_M[-1] - ALTITUDES_M[0]
    expected = ssa / (8.0 * np.pi) * (1.0 - np.exp(-2.0 * extinction * path_length))
    np.testing.assert_allclose(result.radiance.item(), expected, rtol=2.0e-12)


def test_nonzero_bottom_altitude_sets_shadow_and_surface_radius():
    config = single_scatter_config()
    altitudes = np.array([5_000.0, 15_000.0, 30_000.0])
    geometry = sk.Geometry2D(
        cos_sza=1.0,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=altitudes,
        horizontal_angle_grid_radians=np.array([-0.5, 0.5]),
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        calculate_derivatives=False,
    )
    extinction = 1.0e-5
    ssa = 0.5
    atmosphere.storage.total_extinction[:] = extinction
    atmosphere.storage.ssa[:] = ssa
    atmosphere.leg_coeff.a1[0] = 1.0

    result = sk.Engine(config, geometry, vertical_ground_viewing()).calculate_radiance(
        atmosphere
    )

    path_length = altitudes[-1] - altitudes[0]
    expected = ssa / (8.0 * np.pi) * (1.0 - np.exp(-2.0 * extinction * path_length))
    np.testing.assert_allclose(result.radiance.item(), expected, rtol=2.0e-12)


def test_outside_observer_looking_up_has_empty_finite_output():
    config = single_scatter_config()
    geometry = vertical_geometry2d()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.SolarAnglesObserverLocation(
            cos_sza=1.0,
            relative_azimuth=0.0,
            cos_viewing_zenith=1.0,
            observer_altitude_m=100_000.0,
        )
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        legendre_derivative=False,
    )
    atmosphere.storage.total_extinction[:] = 1.0e-5
    atmosphere.storage.ssa[:] = 0.5
    atmosphere.leg_coeff.a1[0] = 1.0

    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    assert result.radiance.item() == pytest.approx(0.0, abs=0.0)
    assert np.all(np.isfinite(result.wf_extinction))
    assert np.all(result.wf_extinction == 0.0)
    assert np.all(result.wf_ssa == 0.0)


@pytest.mark.parametrize("num_stokes", [1, 3])
@pytest.mark.parametrize(
    "interpolation_method",
    [
        sk.InterpolationMethod.LinearInterpolation,
        sk.InterpolationMethod.LowerInterpolation,
    ],
)
def test_horizontally_uniform_2d_single_scatter_matches_1d(
    num_stokes: int, interpolation_method: sk.InterpolationMethod
):
    config = single_scatter_config(num_stokes)
    geometry_1d = geometry1d(interpolation_method)
    geometry_2d = geometry2d(interpolation_method=interpolation_method)
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
    extinction = np.array([[1.0e-6, 1.4e-6], [2.0e-6, 2.5e-6], [0.4e-6, 0.7e-6]])
    ssa = np.array([[0.7, 0.5], [0.6, 0.8], [0.4, 0.3]])
    atmosphere_1d.storage.total_extinction[:] = extinction
    atmosphere_1d.storage.ssa[:] = ssa
    atmosphere_2d.storage.total_extinction[:] = np.tile(extinction, (2, 1))
    atmosphere_2d.storage.ssa[:] = np.tile(ssa, (2, 1))
    set_phase(atmosphere_1d)
    set_phase(atmosphere_2d)

    viewing = parity_viewing()
    result_1d = sk.Engine(config, geometry_1d, viewing).calculate_radiance(
        atmosphere_1d
    )
    result_2d = sk.Engine(config, geometry_2d, viewing).calculate_radiance(
        atmosphere_2d
    )

    np.testing.assert_allclose(result_2d.radiance, result_1d.radiance, rtol=2.0e-11)


def test_refining_same_continuous_2d_field_does_not_change_radiance():
    config = single_scatter_config()
    viewing = parity_viewing()

    def calculate(horizontal_angles: np.ndarray):
        geometry = geometry2d(horizontal_angles)
        atmosphere = sk.Atmosphere(
            geometry,
            config,
            wavelengths_nm=np.array([500.0]),
            calculate_derivatives=False,
        )
        horizontal, altitude = np.meshgrid(
            horizontal_angles, ALTITUDES_M, indexing="ij"
        )
        atmosphere.storage.total_extinction[:, 0] = (
            1.2e-6
            + 0.35e-6 * altitude.ravel() / ALTITUDES_M[-1]
            + 0.4e-6 * horizontal.ravel()
        )
        atmosphere.storage.ssa[:, 0] = (
            0.55 + 0.08 * horizontal.ravel() - 0.05 * altitude.ravel() / ALTITUDES_M[-1]
        )
        atmosphere.leg_coeff.a1[0] = 1.0
        atmosphere.leg_coeff.a1[2] = 0.3
        return sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    coarse = calculate(np.array([-0.6, 0.0, 0.6]))
    refined = calculate(np.array([-0.6, -0.3, 0.0, 0.3, 0.6]))

    np.testing.assert_allclose(refined.radiance, coarse.radiance, rtol=2.0e-13)


def test_broadcast_rayleigh_and_solar_irradiance_constituents_match_1d():
    config = single_scatter_config()
    wavelengths = np.array([350.0, 500.0])
    geometry_1d = geometry1d()
    geometry_2d = geometry2d()
    atmosphere_1d = sk.Atmosphere(
        geometry_1d,
        config,
        wavelengths_nm=wavelengths,
        calculate_derivatives=False,
    )
    atmosphere_2d = sk.Atmosphere(
        geometry_2d,
        config,
        wavelengths_nm=wavelengths,
        calculate_derivatives=False,
    )
    for atmosphere in (atmosphere_1d, atmosphere_2d):
        sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
        atmosphere["rayleigh"] = sk.constituent.Rayleigh()
        atmosphere["solar"] = sk.constituent.SolarIrradiance()

    viewing = parity_viewing()
    result_1d = sk.Engine(config, geometry_1d, viewing).calculate_radiance(
        atmosphere_1d
    )
    result_2d = sk.Engine(config, geometry_2d, viewing).calculate_radiance(
        atmosphere_2d
    )

    assert np.all(atmosphere_2d.storage.solar_irradiance != 1.0)
    np.testing.assert_allclose(result_2d.radiance, result_1d.radiance, rtol=2.0e-11)


def test_native_2d_manual_constituent_matches_raw_storage():
    config = single_scatter_config()
    geometry = geometry2d(np.array([-0.6, 0.0, 0.6]))
    wavelengths = np.array([500.0])
    horizontal, altitude = np.meshgrid(
        geometry.horizontal_angles(), geometry.altitudes(), indexing="ij"
    )
    extinction = (1.0e-6 + 0.2e-6 * horizontal + 0.3e-6 * altitude / ALTITUDES_M[-1])[
        ..., np.newaxis
    ]
    ssa = (0.55 + 0.1 * horizontal)[..., np.newaxis]
    legendre = np.zeros(
        (config.num_singlescatter_moments, *geometry.shape, wavelengths.size)
    )
    legendre[0] = 1.0
    legendre[2] = 0.3 + 0.1 * horizontal[..., np.newaxis]

    constituent_atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=wavelengths,
        calculate_derivatives=False,
    )
    constituent_atmosphere["scatterer"] = sk.constituent.Manual(
        extinction, ssa, legendre
    )
    raw_atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=wavelengths,
        calculate_derivatives=False,
    )
    raw_atmosphere.storage.total_extinction[:] = extinction.reshape(-1, 1)
    raw_atmosphere.storage.ssa[:] = ssa.reshape(-1, 1)
    raw_atmosphere.storage.leg_coeff[:] = legendre.reshape(legendre.shape[0], -1, 1)

    viewing = parity_viewing()
    constituent_result = sk.Engine(config, geometry, viewing).calculate_radiance(
        constituent_atmosphere
    )
    raw_result = sk.Engine(config, geometry, viewing).calculate_radiance(raw_atmosphere)

    np.testing.assert_allclose(
        constituent_result.radiance, raw_result.radiance, rtol=2.0e-13
    )


def test_native_2d_extinction_ssa_and_phase_derivatives_match_finite_difference():
    config = single_scatter_config()
    geometry = geometry2d(np.array([-0.6, 0.0, 0.6]))
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
    )
    location_factor = np.arange(geometry.shape[0] * geometry.shape[1]).reshape(
        geometry.shape
    )
    atmosphere.storage.total_extinction[:, 0] = (
        1.0e-6 + 0.15e-6 * location_factor.ravel()
    )
    atmosphere.storage.ssa[:, 0] = 0.45 + 0.02 * location_factor.ravel()
    atmosphere.leg_coeff.a1[0] = 1.0
    atmosphere.leg_coeff.a1[2] = 0.3
    viewing = parity_viewing()
    engine = sk.Engine(config, geometry, viewing)

    base = engine.calculate_radiance(atmosphere)
    assert base.wf_extinction.shape == (3, 3, 1, 2, 1)
    assert base.wf_ssa.shape == (3, 3, 1, 2, 1)
    assert base.wf_leg_coeff_2.shape == (3, 3, 1, 2, 1)

    horizontal_index = 1
    altitude_index = 1
    location_index = geometry.location_index(altitude_index, horizontal_index)

    extinction_delta = 1.0e-10
    atmosphere.storage.total_extinction[location_index, 0] += extinction_delta
    perturbed_above = engine.calculate_radiance(atmosphere)
    atmosphere.storage.total_extinction[location_index, 0] -= 2 * extinction_delta
    perturbed_below = engine.calculate_radiance(atmosphere)
    numeric = (perturbed_above.radiance - perturbed_below.radiance) / (
        2 * extinction_delta
    )
    analytic = base.wf_extinction[horizontal_index, altitude_index]
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-7, atol=1.0e-10)
    atmosphere.storage.total_extinction[location_index, 0] += extinction_delta

    ssa_delta = 1.0e-6
    atmosphere.storage.ssa[location_index, 0] += ssa_delta
    perturbed = engine.calculate_radiance(atmosphere)
    numeric = (perturbed.radiance - base.radiance) / ssa_delta
    analytic = base.wf_ssa[horizontal_index, altitude_index]
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-6, atol=1.0e-12)
    atmosphere.storage.ssa[location_index, 0] -= ssa_delta

    phase_delta = 1.0e-6
    atmosphere.leg_coeff.a1[2, location_index, 0] += phase_delta
    perturbed = engine.calculate_radiance(atmosphere)
    numeric = (perturbed.radiance - base.radiance) / phase_delta
    analytic = base.wf_leg_coeff_2[horizontal_index, altitude_index]
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-6, atol=1.0e-12)


def test_native_2d_polarized_phase_derivative_matches_finite_difference():
    config = single_scatter_config(num_stokes=3)
    geometry = geometry2d(np.array([-0.6, 0.0, 0.6]))
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
    )
    atmosphere.storage.total_extinction[:] = 1.5e-6
    atmosphere.storage.ssa[:] = 0.65
    set_phase(atmosphere)
    engine = sk.Engine(config, geometry, parity_viewing())
    base = engine.calculate_radiance(atmosphere)

    horizontal_index = 1
    altitude_index = 1
    location_index = geometry.location_index(altitude_index, horizontal_index)
    raw_b1_order_2_index = 4 * 2 + 3
    derivative_name = f"wf_leg_coeff_{raw_b1_order_2_index}"
    delta = 1.0e-6
    atmosphere.leg_coeff.b1[2, location_index, 0] += delta
    perturbed = engine.calculate_radiance(atmosphere)

    numeric = (perturbed.radiance - base.radiance) / delta
    analytic = base[derivative_name][horizontal_index, altitude_index]
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-6, atol=1.0e-12)
    assert np.any(np.abs(analytic[..., 1:]) > 0.0)


@pytest.mark.parametrize(
    "geometry_factory",
    [vertical_geometry1d, vertical_geometry2d],
    ids=["geometry1d", "geometry2d"],
)
@pytest.mark.parametrize("zero_quantity", ["extinction", "ssa"])
def test_zero_scattering_factors_have_finite_nonzero_boundary_derivatives(
    geometry_factory,
    zero_quantity: str,
):
    config = single_scatter_config()
    geometry = geometry_factory()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        legendre_derivative=False,
    )
    atmosphere.storage.total_extinction[:] = 1.0e-5
    atmosphere.storage.ssa[:] = 0.5
    atmosphere.leg_coeff.a1[0] = 1.0
    if zero_quantity == "extinction":
        atmosphere.storage.total_extinction[:] = 0.0
        derivative_name = "wf_extinction"
    else:
        atmosphere.storage.ssa[:] = 0.0
        derivative_name = "wf_ssa"

    result = sk.Engine(config, geometry, vertical_ground_viewing()).calculate_radiance(
        atmosphere
    )

    assert result.radiance.item() == pytest.approx(0.0, abs=1.0e-15)
    derivative = result[derivative_name].values
    assert np.all(np.isfinite(derivative))
    assert np.any(np.abs(derivative) > 0.0)


@pytest.mark.parametrize(
    "geometry_factory",
    [vertical_geometry1d, vertical_geometry2d],
    ids=["geometry1d", "geometry2d"],
)
def test_zero_extinction_and_ssa_have_finite_zero_derivatives(geometry_factory):
    config = single_scatter_config()
    geometry = geometry_factory()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        legendre_derivative=False,
    )
    atmosphere.storage.total_extinction[:] = 0.0
    atmosphere.storage.ssa[:] = 0.0
    atmosphere.leg_coeff.a1[0] = 1.0

    result = sk.Engine(config, geometry, vertical_ground_viewing()).calculate_radiance(
        atmosphere
    )

    assert result.radiance.item() == pytest.approx(0.0, abs=1.0e-15)
    for derivative_name in ("wf_extinction", "wf_ssa"):
        derivative = result[derivative_name].values
        assert np.all(np.isfinite(derivative))
        assert np.all(derivative == 0.0)


def test_ground_shadow_blocks_atmosphere_and_surface_single_scatter():
    config = single_scatter_config()
    geometry = sk.Geometry2D(
        cos_sza=-1.0,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=np.array([-0.5, 0.5]),
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=-1.0,
            relative_azimuth=0.0,
            cos_viewing_zenith=1.0,
            observer_altitude_m=100_000.0,
        )
    )
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([500.0]))
    atmosphere.storage.total_extinction[:] = 1.0e-5
    atmosphere.storage.ssa[:] = 0.8
    atmosphere.leg_coeff.a1[0] = 1.0
    atmosphere.surface.albedo[:] = 1.0

    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    assert result.radiance.item() == pytest.approx(0.0, abs=1.0e-15)
    assert np.all(np.isfinite(result.wf_extinction))
    assert np.all(np.isfinite(result.wf_ssa))
    assert result.wf_albedo.item() == pytest.approx(0.0, abs=1.0e-15)


def test_ground_single_scatter_and_surface_derivative_match_1d():
    config = single_scatter_config()
    geometry_1d = geometry1d()
    geometry_2d = geometry2d()
    atmosphere_1d = sk.Atmosphere(geometry_1d, config, wavelengths_nm=np.array([500.0]))
    atmosphere_2d = sk.Atmosphere(geometry_2d, config, wavelengths_nm=np.array([500.0]))
    extinction = np.array([[1.0e-6], [2.0e-6], [0.5e-6]])
    ssa = np.full_like(extinction, 0.6)
    atmosphere_1d.storage.total_extinction[:] = extinction
    atmosphere_1d.storage.ssa[:] = ssa
    atmosphere_2d.storage.total_extinction[:] = np.tile(extinction, (2, 1))
    atmosphere_2d.storage.ssa[:] = np.tile(ssa, (2, 1))
    atmosphere_1d.surface.albedo[:] = 0.35
    atmosphere_2d.surface.albedo[:] = 0.35
    set_phase(atmosphere_1d)
    set_phase(atmosphere_2d)
    viewing = ground_viewing()

    result_1d = sk.Engine(config, geometry_1d, viewing).calculate_radiance(
        atmosphere_1d
    )
    result_2d = sk.Engine(config, geometry_2d, viewing).calculate_radiance(
        atmosphere_2d
    )

    np.testing.assert_allclose(result_2d.radiance, result_1d.radiance, rtol=2.0e-12)
    np.testing.assert_allclose(result_2d.wf_albedo, result_1d.wf_albedo, rtol=2.0e-12)


def test_single_scatter_occultation_and_emission_sources_add_linearly():
    geometry = vertical_geometry2d()
    viewing = vertical_ground_viewing()
    results = {}

    source_settings = {
        "single": (
            sk.SingleScatterSource.Exact,
            sk.OccultationSource.NoSource,
            sk.EmissionSource.NoSource,
        ),
        "occultation": (
            sk.SingleScatterSource.NoSource,
            sk.OccultationSource.Standard,
            sk.EmissionSource.NoSource,
        ),
        "emission": (
            sk.SingleScatterSource.NoSource,
            sk.OccultationSource.NoSource,
            sk.EmissionSource.Standard,
        ),
        "combined": (
            sk.SingleScatterSource.Exact,
            sk.OccultationSource.Standard,
            sk.EmissionSource.Standard,
        ),
    }
    for name, (single, occultation, emission) in source_settings.items():
        config = single_scatter_config()
        config.single_scatter_source = single
        config.occultation_source = occultation
        config.emission_source = emission
        atmosphere = sk.Atmosphere(
            geometry,
            config,
            wavelengths_nm=np.array([500.0]),
            calculate_derivatives=False,
        )
        atmosphere.storage.total_extinction[:] = 1.0e-5
        atmosphere.storage.ssa[:] = 0.4
        atmosphere.storage.emission_source[:] = 0.7
        atmosphere.leg_coeff.a1[0] = 1.0
        results[name] = sk.Engine(config, geometry, viewing).calculate_radiance(
            atmosphere
        )

    expected = (
        results["single"].radiance
        + results["occultation"].radiance
        + results["emission"].radiance
    )
    np.testing.assert_allclose(results["combined"].radiance, expected, rtol=1.0e-14)


@pytest.mark.parametrize(
    "single_scatter_source",
    [
        sk.SingleScatterSource.Table,
        sk.SingleScatterSource.DiscreteOrdinates,
    ],
)
def test_geometry2d_rejects_nonexact_single_scatter_sources(
    single_scatter_source: sk.SingleScatterSource,
):
    config = single_scatter_config()
    config.single_scatter_source = single_scatter_source

    with pytest.raises(NotImplementedError, match="exact single scattering"):
        sk.Engine(config, geometry2d(), parity_viewing())


def test_geometry2d_wavelength_batch_matches_scalar():
    wavelengths = np.linspace(500.0, 800.0, 5)
    geometry = geometry2d(np.array([-0.6, 0.0, 0.6]))
    viewing = parity_viewing()

    def calculate(batch_size: int):
        config = single_scatter_config(num_stokes=3)
        config.wavelength_batch_size = batch_size
        atmosphere = sk.Atmosphere(
            geometry,
            config,
            wavelengths_nm=wavelengths,
            calculate_derivatives=True,
        )
        location = np.arange(np.prod(geometry.shape))[:, np.newaxis]
        spectral = np.linspace(0.8, 1.2, wavelengths.size)[np.newaxis, :]
        atmosphere.storage.total_extinction[:] = (1.0e-6 + location * 0.1e-6) * spectral
        atmosphere.storage.ssa[:] = 0.45 + 0.01 * location
        atmosphere.surface.albedo[:] = np.linspace(0.1, 0.3, wavelengths.size)
        set_phase(atmosphere)
        return sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    scalar = calculate(1)
    batched = calculate(3)
    assert scalar.data_vars.keys() == batched.data_vars.keys()
    for variable in scalar.data_vars:
        np.testing.assert_allclose(
            batched[variable], scalar[variable], rtol=5.0e-12, atol=2.0e-13
        )
