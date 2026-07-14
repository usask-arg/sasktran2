from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr

ALTITUDES_M = np.array([0.0, 10_000.0, 30_000.0])
HORIZONTAL_ANGLES = np.array([-0.3, 0.3])
WAVELENGTHS_NM = np.array([500.0, 600.0])
XS_M2 = np.array([2.0e-12, 5.0e-12])
SSA = np.array([0.8, 0.7])
ASYMMETRY = np.array([0.4, 0.6])


def _config() -> sk.Config:
    config = sk.Config()
    config.num_streams = 4
    config.delta_m_scaling = False
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    return config


def _single_scatter_config() -> sk.Config:
    config = _config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.num_singlescatter_moments = 16
    return config


def _viewing_geometry() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.TangentAltitude(
            tangent_altitude_m=5_000.0,
            observer_altitude_m=100_000.0,
            horizontal_angle_radians=-0.1,
            viewing_azimuth_radians=np.pi / 2.0,
        )
    )
    return viewing


def _geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=HORIZONTAL_ANGLES,
    )


def _geometry1d() -> sk.Geometry1D:
    return sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )


def _atmosphere(
    geometry: sk.Geometry1D | sk.Geometry2D, *, derivatives: bool = True
) -> sk.Atmosphere:
    return sk.Atmosphere(
        geometry,
        _config(),
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=derivatives,
        legendre_derivative=False,
    )


def _constant_optical_property() -> sk.optical.HenyeyGreenstein:
    return sk.optical.HenyeyGreenstein.from_parameters(
        wavelength_nm=WAVELENGTHS_NM,
        xs_total=XS_M2,
        ssa=XS_M2 * SSA,
        g=ASYMMETRY,
        max_num_moments=16,
    )


def _background_optical_property() -> sk.optical.HenyeyGreenstein:
    cross_section = np.array([1.0e-12, 3.0e-12])
    return sk.optical.HenyeyGreenstein.from_parameters(
        wavelength_nm=WAVELENGTHS_NM,
        xs_total=cross_section,
        ssa=cross_section * np.array([0.35, 0.45]),
        g=np.array([-0.2, -0.1]),
        max_num_moments=16,
    )


def _radius_dependent_optical_property() -> sk.optical.HenyeyGreenstein:
    radius = np.array([1.0, 3.0])
    database = xr.Dataset(
        {
            "xs_total": (
                ("radius", "wavelength_nm"),
                np.array([[2.0e-12, 4.0e-12], [6.0e-12, 18.0e-12]]),
            ),
            "ssa": (
                ("radius", "wavelength_nm"),
                np.array(
                    [[1.6e-12, 2.8e-12], [4.8e-12, 12.6e-12]]
                ),
            ),
            "asymmetry_parameter": (
                ("radius", "wavelength_nm"),
                np.array([[0.4, 0.5], [0.6, 0.7]]),
            ),
        },
        coords={"radius": radius, "wavelength_nm": WAVELENGTHS_NM},
    )
    return sk.optical.HenyeyGreenstein(db=database, max_num_moments=16)


def _number_density() -> np.ndarray:
    return np.array([[1.0e6, 2.0e6, 3.0e6], [4.0e6, 5.0e6, 6.0e6]])


def test_number_density_scatterer_2d_preserves_native_shape_and_mutability():
    number_density = _number_density()
    constituent = sk.constituent.NumberDensityScatterer2D(
        _constant_optical_property(), number_density, particle_size=2.0
    )

    assert constituent.volume_spatial_mode == "native_2d"
    np.testing.assert_array_equal(constituent.number_density, number_density)
    np.testing.assert_array_equal(
        constituent.particle_size, np.full(number_density.shape, 2.0)
    )

    constituent.number_density[1, 2] *= 2.0
    assert constituent.number_density[1, 2] == 2.0 * number_density[1, 2]

    replacement = number_density * 3.0
    constituent.number_density = replacement
    np.testing.assert_array_equal(constituent.number_density, replacement)


@pytest.mark.parametrize(
    "number_density", [np.ones(3), np.ones((1, 2, 3)), np.ones((0, 3))]
)
def test_number_density_scatterer_2d_rejects_non_native_shapes(number_density):
    with pytest.raises(ValueError, match="number_density"):
        sk.constituent.NumberDensityScatterer2D(
            _constant_optical_property(), number_density
        )


def test_number_density_scatterer_2d_validates_setters_aux_inputs_and_atmosphere():
    constituent = sk.constituent.NumberDensityScatterer2D(
        _constant_optical_property(), _number_density(), particle_size=2.0
    )
    with pytest.raises(ValueError, match="retain shape"):
        constituent.number_density = np.ones((1, 6))
    with pytest.raises(ValueError, match="particle_size"):
        constituent.particle_size = np.ones(3)

    atmosphere1d = _atmosphere(_geometry1d(), derivatives=False)
    atmosphere1d["aerosol"] = constituent
    with pytest.raises(TypeError, match="Geometry2D"):
        atmosphere1d.internal_object()

    atmosphere2d = _atmosphere(_geometry2d(), derivatives=False)
    atmosphere2d["aerosol"] = sk.constituent.NumberDensityScatterer2D(
        _constant_optical_property(), np.ones((3, 3))
    )
    with pytest.raises(ValueError, match="does not match"):
        atmosphere2d.internal_object()


def test_number_density_scatterer_2d_populates_native_optical_storage():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry, derivatives=False)
    number_density = _number_density()
    atmosphere["aerosol"] = sk.constituent.NumberDensityScatterer2D(
        _constant_optical_property(), number_density
    )

    atmosphere.internal_object()

    expected_extinction = number_density.reshape(-1, 1) * XS_M2
    np.testing.assert_allclose(atmosphere.storage.total_extinction, expected_extinction)
    np.testing.assert_allclose(
        atmosphere.storage.ssa, np.broadcast_to(SSA, expected_extinction.shape)
    )
    assert np.all(np.isfinite(atmosphere.storage.leg_coeff))


def test_uniform_number_density_scatterer_2d_matches_altitude_profile():
    profile = np.array([1.0e6, 2.0e6, 3.0e6])
    atmosphere1d = _atmosphere(_geometry1d(), derivatives=False)
    atmosphere1d["aerosol"] = sk.constituent.NumberDensityScatterer(
        _constant_optical_property(), ALTITUDES_M, profile
    )
    atmosphere2d = _atmosphere(_geometry2d(), derivatives=False)
    atmosphere2d["aerosol"] = sk.constituent.NumberDensityScatterer2D(
        _constant_optical_property(), np.tile(profile, (2, 1))
    )

    atmosphere1d.internal_object()
    atmosphere2d.internal_object()

    np.testing.assert_allclose(
        atmosphere2d.storage.total_extinction,
        np.tile(atmosphere1d.storage.total_extinction, (2, 1)),
    )
    np.testing.assert_allclose(
        atmosphere2d.storage.ssa, np.tile(atmosphere1d.storage.ssa, (2, 1))
    )
    legendre2d = atmosphere2d.storage.leg_coeff.reshape(
        atmosphere2d.storage.leg_coeff.shape[0], 2, 3, 2
    )
    np.testing.assert_allclose(legendre2d[:, 0], atmosphere1d.storage.leg_coeff)
    np.testing.assert_allclose(legendre2d[:, 1], atmosphere1d.storage.leg_coeff)


def test_native_number_density_derivative_is_local_and_has_no_dense_mapping():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry)
    constituent = sk.constituent.NumberDensityScatterer2D(
        _constant_optical_property(), _number_density()
    )
    atmosphere["aerosol"] = constituent
    atmosphere.internal_object()

    base_extinction = atmosphere.storage.total_extinction.copy()
    mapping_name = "wf_aerosol_number_density"
    mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    analytic_d_extinction = mapping.d_extinction.copy()
    assert mapping.interp_dim == "location"
    assert mapping.interpolator.shape == (0, 0)
    assert atmosphere.derivative_output_shape(mapping_name) == geometry.shape
    for field in ("d_extinction", "d_ssa", "d_leg_coeff", "scat_factor"):
        assert np.all(np.isfinite(getattr(mapping, field)))

    horizontal_index = 1
    altitude_index = 1
    location_index = geometry.location_index(altitude_index, horizontal_index)
    delta = 1.0
    constituent.number_density[horizontal_index, altitude_index] += delta
    atmosphere.internal_object()
    numeric = (atmosphere.storage.total_extinction - base_extinction) / delta

    expected = np.zeros_like(numeric)
    expected[location_index] = analytic_d_extinction[location_index]
    np.testing.assert_allclose(numeric, expected, rtol=2.0e-10, atol=1.0e-18)


def test_native_optical_input_derivative_is_local_and_matches_finite_difference():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry)
    radius = np.array([[1.2, 1.4, 1.6], [1.8, 2.0, 2.2]])
    constituent = sk.constituent.NumberDensityScatterer2D(
        _radius_dependent_optical_property(), _number_density(), radius=radius
    )
    atmosphere["aerosol"] = constituent
    atmosphere.internal_object()

    base_extinction = atmosphere.storage.total_extinction.copy()
    mapping_name = "wf_aerosol_radius"
    mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    analytic_d_extinction = mapping.d_extinction.copy()
    assert mapping.interp_dim == "location"
    assert mapping.interpolator.shape == (0, 0)
    assert atmosphere.derivative_output_shape(mapping_name) == geometry.shape

    horizontal_index = 0
    altitude_index = 1
    location_index = geometry.location_index(altitude_index, horizontal_index)
    delta = 1.0e-5
    constituent.radius[horizontal_index, altitude_index] += delta
    atmosphere.internal_object()
    numeric = (atmosphere.storage.total_extinction - base_extinction) / delta

    expected = np.zeros_like(numeric)
    expected[location_index] = analytic_d_extinction[location_index]
    np.testing.assert_allclose(numeric, expected, rtol=2.0e-9, atol=1.0e-18)


@pytest.mark.parametrize(
    ("profile_name", "mapping_name", "delta"),
    [
        ("number_density", "wf_aerosol_number_density", 100.0),
        ("radius", "wf_aerosol_radius", 1.0e-4),
    ],
)
def test_native_scatterer_single_scatter_wf_matches_finite_difference(
    profile_name, mapping_name, delta
):
    config = _single_scatter_config()
    geometry = _geometry2d()
    radius = np.array([[1.2, 1.4, 1.6], [1.8, 2.0, 2.2]])
    constituent = sk.constituent.NumberDensityScatterer2D(
        _radius_dependent_optical_property(), _number_density(), radius=radius
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        legendre_derivative=False,
    )
    atmosphere["background"] = sk.constituent.NumberDensityScatterer2D(
        _background_optical_property(), np.full(geometry.shape, 2.0e6)
    )
    atmosphere["aerosol"] = constituent
    engine = sk.Engine(config, geometry, _viewing_geometry())

    base = engine.calculate_radiance(atmosphere)
    horizontal_index = 0
    altitude_index = 1
    profile = getattr(constituent, profile_name)
    profile[horizontal_index, altitude_index] += delta
    above = engine.calculate_radiance(atmosphere)
    profile[horizontal_index, altitude_index] -= 2.0 * delta
    below = engine.calculate_radiance(atmosphere)
    numeric = (above.radiance - below.radiance) / (2.0 * delta)
    analytic = base[mapping_name][horizontal_index, altitude_index]

    assert np.any(np.abs(analytic) > 0.0)
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-5, atol=1.0e-18)


def test_extinction_scatterer_2d_normalizes_each_native_location():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry)
    extinction = np.array([[1.0e-6, 2.0e-6, 3.0e-6], [4.0e-6, 5.0e-6, 6.0e-6]])
    constituent = sk.constituent.ExtinctionScatterer2D(
        _constant_optical_property(), extinction, 500.0
    )
    atmosphere["aerosol"] = constituent
    atmosphere.internal_object()

    expected = extinction.reshape(-1, 1) * (XS_M2 / XS_M2[0])
    np.testing.assert_allclose(atmosphere.storage.total_extinction, expected)
    np.testing.assert_allclose(constituent.number_density, extinction / XS_M2[0])

    mapping_name = "wf_aerosol_extinction"
    mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    analytic_d_extinction = mapping.d_extinction.copy()
    assert mapping.interp_dim == "location"
    assert mapping.interpolator.shape == (0, 0)
    assert atmosphere.derivative_output_shape(mapping_name) == geometry.shape
    np.testing.assert_allclose(mapping.d_extinction[:, 0], 1.0)

    location_index = geometry.location_index(1, 1)
    base_extinction = atmosphere.storage.total_extinction.copy()
    delta = 1.0e-10
    constituent.extinction_per_m[1, 1] += delta
    atmosphere.internal_object()
    numeric = (atmosphere.storage.total_extinction - base_extinction) / delta
    expected_derivative = np.zeros_like(numeric)
    expected_derivative[location_index] = analytic_d_extinction[location_index]
    np.testing.assert_allclose(
        numeric, expected_derivative, rtol=2.0e-11, atol=1.0e-13
    )


def test_extinction_scatterer_2d_aux_derivative_includes_normalization_change():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry)
    extinction = np.full(geometry.shape, 2.0e-6)
    radius = np.array([[1.2, 1.4, 1.6], [1.8, 2.0, 2.2]])
    constituent = sk.constituent.ExtinctionScatterer2D(
        _radius_dependent_optical_property(), extinction, 500.0, radius=radius
    )
    atmosphere["aerosol"] = constituent
    atmosphere.internal_object()

    mapping = atmosphere.storage.get_derivative_mapping("wf_aerosol_radius")
    analytic_d_extinction = mapping.d_extinction.copy()
    np.testing.assert_allclose(analytic_d_extinction[:, 0], 0.0, atol=1.0e-20)

    horizontal_index = 1
    altitude_index = 0
    location_index = geometry.location_index(altitude_index, horizontal_index)
    delta = 1.0e-5
    constituent.radius[horizontal_index, altitude_index] += delta
    atmosphere.internal_object()
    above = atmosphere.storage.total_extinction.copy()
    constituent.radius[horizontal_index, altitude_index] -= 2.0 * delta
    atmosphere.internal_object()
    below = atmosphere.storage.total_extinction.copy()
    numeric = (above - below) / (2.0 * delta)

    expected = np.zeros_like(numeric)
    expected[location_index] = analytic_d_extinction[location_index]
    np.testing.assert_allclose(numeric, expected, rtol=3.0e-9, atol=1.0e-17)


@pytest.mark.parametrize(
    ("profile_name", "mapping_name", "delta"),
    [
        ("extinction_per_m", "wf_aerosol_extinction", 1.0e-10),
        ("radius", "wf_aerosol_radius", 1.0e-4),
    ],
)
def test_native_extinction_scatterer_single_scatter_wf_matches_finite_difference(
    profile_name, mapping_name, delta
):
    config = _single_scatter_config()
    geometry = _geometry2d()
    extinction = np.array(
        [[1.0e-6, 2.0e-6, 3.0e-6], [4.0e-6, 5.0e-6, 6.0e-6]]
    )
    radius = np.array([[1.2, 1.4, 1.6], [1.8, 2.0, 2.2]])
    constituent = sk.constituent.ExtinctionScatterer2D(
        _radius_dependent_optical_property(), extinction, 500.0, radius=radius
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        legendre_derivative=False,
    )
    atmosphere["background"] = sk.constituent.NumberDensityScatterer2D(
        _background_optical_property(), np.full(geometry.shape, 2.0e6)
    )
    atmosphere["aerosol"] = constituent
    engine = sk.Engine(config, geometry, _viewing_geometry())

    base = engine.calculate_radiance(atmosphere)
    horizontal_index = 0
    altitude_index = 1
    profile = getattr(constituent, profile_name)
    profile[horizontal_index, altitude_index] += delta
    above = engine.calculate_radiance(atmosphere)
    profile[horizontal_index, altitude_index] -= 2.0 * delta
    below = engine.calculate_radiance(atmosphere)
    numeric = (above.radiance - below.radiance) / (2.0 * delta)
    analytic = base[mapping_name][horizontal_index, altitude_index]

    assert np.any(np.abs(analytic) > 0.0)
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-5, atol=1.0e-13)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"extinction_per_m": np.ones(3), "extinction_wavelength_nm": 500.0},
        {
            "extinction_per_m": np.ones((2, 3)),
            "extinction_wavelength_nm": np.nan,
        },
        {
            "extinction_per_m": np.ones((2, 3)),
            "extinction_wavelength_nm": 0.0,
        },
    ],
)
def test_extinction_scatterer_2d_validates_constructor(kwargs):
    with pytest.raises(ValueError, match="extinction"):
        sk.constituent.ExtinctionScatterer2D(_constant_optical_property(), **kwargs)


def test_extinction_scatterer_2d_rejects_zero_reference_cross_section_cleanly():
    optical_property = sk.optical.HenyeyGreenstein.from_parameters(
        wavelength_nm=WAVELENGTHS_NM,
        xs_total=np.array([0.0, 1.0e-12]),
        ssa=np.array([0.0, 0.5e-12]),
        g=np.array([0.0, 0.2]),
        max_num_moments=16,
    )
    atmosphere = _atmosphere(_geometry2d(), derivatives=False)
    atmosphere["aerosol"] = sk.constituent.ExtinctionScatterer2D(
        optical_property, np.full(_geometry2d().shape, 1.0e-6), 500.0
    )

    with pytest.raises(ValueError, match="finite and positive"):
        atmosphere.internal_object()

    np.testing.assert_array_equal(atmosphere.storage.total_extinction, 0.0)
