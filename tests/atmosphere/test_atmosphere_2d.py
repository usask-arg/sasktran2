from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk

from sasktran2.constituent.base import Constituent


def geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=np.array([0.0, 10_000.0, 30_000.0]),
        horizontal_angle_grid_radians=np.array([-0.2, 0.3]),
    )


def atmosphere2d(
    *, calculate_derivatives: bool = True, num_stokes: int = 1
) -> sk.Atmosphere:
    config = sk.Config()
    config.num_stokes = num_stokes
    config.num_streams = 2
    config.delta_m_scaling = False
    return sk.Atmosphere(
        geometry2d(),
        config,
        wavelengths_nm=np.array([350.0, 600.0]),
        calculate_derivatives=calculate_derivatives,
    )


class _LegacyAltitudeConstituent(Constituent):
    def __init__(self) -> None:
        self.seen_altitudes = None
        self.add_atmosphere = None
        self.derivative_atmosphere = None
        self.add_pressure = None
        self.derivative_pressure = None

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        self.add_atmosphere = atmo
        self.add_pressure = atmo.pressure_pa
        self.seen_altitudes = atmo.model_geometry.altitudes().copy()
        atmo.storage.total_extinction[:] += self.seen_altitudes[:, np.newaxis] + 1.0

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        self.derivative_atmosphere = atmo
        self.derivative_pressure = atmo.pressure_pa
        mapping = atmo.storage.get_derivative_mapping(f"wf_{name}")
        mapping.d_extinction[:] = 1.0
        mapping.d_ssa[:] = 0.0
        mapping.interp_dim = f"{name}_altitude"
        num_altitudes = atmo.volume_shape[-1]
        mapping.interpolator = np.tile(
            np.eye(num_altitudes), (atmo.volume_shape[0], 1)
        )


class _FailingDerivativeConstituent(Constituent):
    def __init__(self) -> None:
        self.fail = True

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        atmo.storage.total_extinction[:] += 2.0

    def register_derivative(self, atmo: sk.Atmosphere, name: str):  # noqa: ARG002
        if self.fail:
            raise RuntimeError("intentional derivative failure")


class _InvalidSpatialModeConstituent(Constituent):
    @property
    def volume_spatial_mode(self) -> str:
        return "horizontal_profile"

    def add_to_atmosphere(self, atmo: sk.Atmosphere):  # noqa: ARG002
        pass

    def register_derivative(self, atmo: sk.Atmosphere, name: str):  # noqa: ARG002
        pass


class _NativeStateConstituent(Constituent):
    def __init__(self) -> None:
        self.add_pressure = None
        self.derivative_pressure = None

    @property
    def volume_spatial_mode(self) -> str:
        return "native_2d"

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        self.add_pressure = atmo._native_state("pressure_pa")

    def register_derivative(self, atmo: sk.Atmosphere, name: str):  # noqa: ARG002
        self.derivative_pressure = atmo._native_state("pressure_pa")


def test_atmosphere_2d_allocates_altitude_fastest_native_storage():
    atmosphere = atmosphere2d()

    assert atmosphere.volume_shape == (2, 3)
    assert atmosphere.num_locations == 6
    assert atmosphere.storage.total_extinction.shape == (6, 2)
    assert atmosphere.storage.ssa.shape == (6, 2)
    assert atmosphere.storage.leg_coeff.shape[1:] == (6, 2)

    native = np.arange(12).reshape(6, 2)
    structured = atmosphere.reshape_native(native)
    assert structured.shape == (2, 3, 2)
    np.testing.assert_array_equal(structured[0], native[:3])
    np.testing.assert_array_equal(structured[1], native[3:])


def test_atmosphere_2d_accepts_broadcast_and_native_state():
    atmosphere = atmosphere2d()

    pressure = np.array([100_000.0, 30_000.0, 1_000.0])
    temperature = np.array([[280.0, 240.0, 210.0], [285.0, 245.0, 215.0]])
    atmosphere.pressure_pa = pressure
    atmosphere.temperature_k = temperature

    assert atmosphere.pressure_pa.shape == (3,)
    assert atmosphere.temperature_k.shape == (2, 3)
    np.testing.assert_array_equal(
        atmosphere._native_state("pressure_pa"), np.tile(pressure, 2)
    )
    np.testing.assert_array_equal(
        atmosphere._native_state("temperature_k"), temperature.reshape(-1)
    )

    atmosphere.specific_humidity = np.full((2, 3), 0.01)
    atmosphere.specific_humidity = None
    assert atmosphere.specific_humidity is None


@pytest.mark.parametrize("shape", [(2,), (6,), (3, 2), (1, 3)])
def test_atmosphere_2d_rejects_ambiguous_state_shapes(shape):
    atmosphere = atmosphere2d()
    with pytest.raises(ValueError, match="temperature_k must have shape"):
        atmosphere.temperature_k = np.ones(shape)


def test_legacy_custom_constituent_is_broadcast_without_code_changes():
    atmosphere = atmosphere2d(calculate_derivatives=False)
    constituent = _LegacyAltitudeConstituent()
    atmosphere["legacy"] = constituent

    atmosphere.internal_object()

    expected_altitudes = np.tile(atmosphere.model_geometry.altitudes(), 2)
    np.testing.assert_array_equal(constituent.seen_altitudes, expected_altitudes)
    np.testing.assert_array_equal(
        atmosphere.storage.total_extinction,
        expected_altitudes[:, np.newaxis] + np.ones((1, 2)),
    )


def test_legacy_constituent_receives_one_cached_profile_view_per_build():
    atmosphere = atmosphere2d()
    atmosphere.pressure_pa = np.array([100_000.0, 30_000.0, 1_000.0])
    constituent = _LegacyAltitudeConstituent()
    atmosphere["legacy"] = constituent

    atmosphere.internal_object()

    assert constituent.add_atmosphere is constituent.derivative_atmosphere
    assert constituent.add_pressure is constituent.derivative_pressure


def test_native_constituent_receives_cached_flattened_state_per_build():
    atmosphere = atmosphere2d()
    pressure = np.array(
        [[100_000.0, 30_000.0, 1_000.0], [90_000.0, 25_000.0, 800.0]]
    )
    atmosphere.pressure_pa = pressure
    constituent = _NativeStateConstituent()
    atmosphere["native"] = constituent

    atmosphere.internal_object()

    assert constituent.add_pressure is constituent.derivative_pressure
    np.testing.assert_array_equal(constituent.add_pressure, pressure.reshape(-1))


def test_invalid_volume_spatial_mode_is_rejected():
    atmosphere = atmosphere2d()
    atmosphere["invalid"] = _InvalidSpatialModeConstituent()

    with pytest.raises(ValueError, match="Unsupported volume_spatial_mode"):
        atmosphere.internal_object()


def test_manual_altitude_profile_broadcasts_before_normalization():
    atmosphere = atmosphere2d(calculate_derivatives=False)
    extinction = np.array([[1.0, 2.0], [3.0, 5.0], [7.0, 11.0]])
    ssa = np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]])
    atmosphere["manual"] = sk.constituent.Manual(extinction, ssa)

    atmosphere.internal_object()

    np.testing.assert_array_equal(
        atmosphere.storage.total_extinction, np.tile(extinction, (2, 1))
    )
    np.testing.assert_allclose(atmosphere.storage.ssa, np.tile(ssa, (2, 1)))


def test_manual_accepts_native_2d_fields_and_preserves_public_shape():
    atmosphere = atmosphere2d(calculate_derivatives=False)
    extinction = np.arange(1.0, 13.0).reshape(2, 3, 2)
    ssa = np.linspace(0.1, 0.6, 12).reshape(2, 3, 2)
    constituent = sk.constituent.Manual(extinction, ssa)
    atmosphere["manual"] = constituent

    atmosphere.internal_object()

    assert constituent.extinction.shape == (2, 3, 2)
    assert constituent.ssa.shape == (2, 3, 2)
    np.testing.assert_array_equal(
        atmosphere.storage.total_extinction, extinction.reshape(6, 2)
    )
    np.testing.assert_allclose(atmosphere.storage.ssa, ssa.reshape(6, 2))


def test_manual_rejects_shape_mismatches_before_entering_rust_kernel():
    atmosphere = atmosphere2d(calculate_derivatives=False)

    atmosphere["wavelength"] = sk.constituent.Manual(
        np.ones((3, 1)), np.zeros((3, 1))
    )
    with pytest.raises(ValueError, match="wavelength dimension"):
        atmosphere.internal_object()

    atmosphere = atmosphere2d(calculate_derivatives=False)
    wrong_num_moments = atmosphere.storage.leg_coeff.shape[0] - 1
    atmosphere["moments"] = sk.constituent.Manual(
        np.ones((3, 2)),
        np.full((3, 2), 0.5),
        np.ones((wrong_num_moments, 3, 2)),
    )
    with pytest.raises(ValueError, match="Legendre moment dimension"):
        atmosphere.internal_object()

    with pytest.raises(ValueError, match="legendre_moments must have shape"):
        sk.constituent.Manual(
            np.ones((3, 2)), np.zeros((3, 2)), np.ones((3, 2))
        )

    with pytest.raises(ValueError, match="delta_scale"):
        sk.constituent.Manual(
            np.ones((3, 2)), np.zeros((3, 2)), delta_scale=True
        )


def test_broadcast_and_native_manual_constituents_mix_additively():
    atmosphere = atmosphere2d(calculate_derivatives=False)
    profile_extinction = np.array([[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]])
    profile_ssa = np.full((3, 2), 0.25)
    native_extinction = np.arange(1.0, 13.0).reshape(2, 3, 2) / 10.0
    native_ssa = np.full((2, 3, 2), 0.75)

    atmosphere["profile"] = sk.constituent.Manual(
        profile_extinction, profile_ssa
    )
    atmosphere["native"] = sk.constituent.Manual(native_extinction, native_ssa)
    atmosphere.internal_object()

    profile_extinction = np.tile(profile_extinction, (2, 1))
    native_extinction = native_extinction.reshape(6, 2)
    expected_extinction = profile_extinction + native_extinction
    expected_ssa = (
        profile_extinction * 0.25 + native_extinction * 0.75
    ) / expected_extinction
    np.testing.assert_allclose(
        atmosphere.storage.total_extinction, expected_extinction
    )
    np.testing.assert_allclose(atmosphere.storage.ssa, expected_ssa)


def test_rayleigh_broadcast_state_and_derivative_mapping():
    atmosphere = atmosphere2d()
    atmosphere.pressure_pa = np.array([100_000.0, 30_000.0, 1_000.0])
    atmosphere.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere.internal_object()

    extinction = atmosphere.reshape_native(atmosphere.storage.total_extinction)
    np.testing.assert_allclose(extinction[0], extinction[1])

    pressure_mapping = atmosphere.storage.get_derivative_mapping(
        "wf_rayleigh_pressure_pa"
    )
    temperature_mapping = atmosphere.storage.get_derivative_mapping(
        "wf_rayleigh_temperature_k"
    )
    assert pressure_mapping.interpolator.shape == (6, 3)
    assert temperature_mapping.interpolator.shape == (6, 3)
    assert pressure_mapping.interp_dim == "altitude"
    assert temperature_mapping.interp_dim == "altitude"


def test_broadcast_pressure_derivative_matches_finite_difference():
    atmosphere = atmosphere2d()
    atmosphere.pressure_pa = np.array([100_000.0, 30_000.0, 1_000.0])
    atmosphere.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere.internal_object()

    base_extinction = atmosphere.storage.total_extinction.copy()
    mapping = atmosphere.storage.get_derivative_mapping("wf_rayleigh_pressure_pa")
    altitude_index = 1
    analytic = (
        mapping.d_extinction
        * mapping.interpolator[:, altitude_index, np.newaxis]
    )

    delta = 0.1
    atmosphere.pressure_pa[altitude_index] += delta
    atmosphere.internal_object()
    numeric = (atmosphere.storage.total_extinction - base_extinction) / delta

    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-9, atol=1.0e-25)


def test_rayleigh_native_2d_state_keeps_independent_location_derivatives():
    atmosphere = atmosphere2d()
    atmosphere.pressure_pa = np.array(
        [[100_000.0, 30_000.0, 1_000.0], [90_000.0, 25_000.0, 800.0]]
    )
    atmosphere.temperature_k = np.array(
        [[280.0, 240.0, 210.0], [285.0, 245.0, 215.0]]
    )
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere.internal_object()

    extinction = atmosphere.reshape_native(atmosphere.storage.total_extinction)
    assert not np.allclose(extinction[0], extinction[1])

    mapping_name = "wf_rayleigh_pressure_pa"
    mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    assert mapping.interpolator.shape == (6, 6)
    assert mapping.interp_dim == "location"
    assert atmosphere.derivative_output_shape(mapping_name) == (2, 3)


def test_optional_specific_humidity_uses_broadcast_derivative_layout():
    atmosphere = atmosphere2d()
    atmosphere.pressure_pa = np.array([100_000.0, 30_000.0, 1_000.0])
    atmosphere.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        atmosphere.model_geometry.altitudes(),
        np.full(3, 1.0e-6),
    )

    atmosphere.internal_object()

    humidity_mappings = [
        atmosphere.storage.get_derivative_mapping(name)
        for name in atmosphere.storage.derivative_mapping_names()
        if atmosphere.storage.get_derivative_mapping(name).assign_name
        == "wf_specific_humidity"
    ]
    assert humidity_mappings
    for mapping in humidity_mappings:
        assert mapping.interpolator.shape == (6, 3)
        assert mapping.interp_dim == "altitude"


def test_altitude_profile_scatterer_broadcasts_phase_and_profile_derivatives():
    atmosphere = atmosphere2d()
    optical_property = sk.optical.HenyeyGreenstein.from_parameters(
        wavelength_nm=np.array([350.0, 600.0]),
        xs_total=np.array([1.0e-12, 2.0e-12]),
        ssa=np.array([0.8, 0.7]),
        g=np.array([0.6, 0.5]),
        max_num_moments=16,
    )
    altitudes = atmosphere.model_geometry.altitudes()
    atmosphere["aerosol"] = sk.constituent.NumberDensityScatterer(
        optical_property,
        altitudes,
        np.array([1.0e6, 2.0e6, 3.0e6]),
    )

    atmosphere.internal_object()

    extinction = atmosphere.reshape_native(atmosphere.storage.total_extinction)
    ssa = atmosphere.reshape_native(atmosphere.storage.ssa)
    legendre = atmosphere.storage.leg_coeff.reshape(
        atmosphere.storage.leg_coeff.shape[0], 2, 3, 2
    )
    np.testing.assert_allclose(extinction[0], extinction[1])
    np.testing.assert_allclose(ssa[0], ssa[1])
    np.testing.assert_allclose(legendre[:, 0], legendre[:, 1])

    mapping = atmosphere.storage.get_derivative_mapping(
        "wf_aerosol_number_density"
    )
    assert mapping.interpolator.shape == (6, 3)
    assert mapping.interp_dim == "aerosol_altitude"


def test_native_pressure_derivative_is_local_and_matches_finite_difference():
    atmosphere = atmosphere2d()
    atmosphere.pressure_pa = np.array(
        [[100_000.0, 30_000.0, 1_000.0], [90_000.0, 25_000.0, 800.0]]
    )
    atmosphere.temperature_k = np.array(
        [[280.0, 240.0, 210.0], [285.0, 245.0, 215.0]]
    )
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere.internal_object()

    base_extinction = atmosphere.storage.total_extinction.copy()
    mapping = atmosphere.storage.get_derivative_mapping("wf_rayleigh_pressure_pa")
    horizontal_index = 1
    altitude_index = 1
    location_index = atmosphere.model_geometry.location_index(
        altitude_index, horizontal_index
    )
    analytic = (
        mapping.d_extinction
        * mapping.interpolator[:, location_index, np.newaxis]
    )

    delta = 0.1
    atmosphere.pressure_pa[horizontal_index, altitude_index] += delta
    atmosphere.internal_object()
    numeric = (atmosphere.storage.total_extinction - base_extinction) / delta

    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-9, atol=1.0e-25)
    changed_locations = np.flatnonzero(np.any(numeric != 0.0, axis=1))
    np.testing.assert_array_equal(changed_locations, [location_index])


def test_thermal_emission_broadcasts_and_lifts_temperature_derivative():
    config = sk.Config()
    config.emission_source = sk.EmissionSource.Standard
    atmosphere = sk.Atmosphere(
        geometry2d(), config, wavelengths_nm=np.array([8_000.0, 10_000.0])
    )
    atmosphere.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere["thermal"] = sk.constituent.ThermalEmission()

    atmosphere.internal_object()

    emission = atmosphere.reshape_native(atmosphere.storage.emission_source)
    np.testing.assert_allclose(emission[0], emission[1])
    mapping = atmosphere.storage.get_derivative_mapping("wf_thermal_temperature_k")
    assert mapping.interpolator.shape == (6, 3)
    assert mapping.interp_dim == "altitude"


def test_state_derivative_layout_can_switch_between_broadcast_and_native():
    config = sk.Config()
    config.emission_source = sk.EmissionSource.Standard
    atmosphere = sk.Atmosphere(
        geometry2d(), config, wavelengths_nm=np.array([8_000.0, 10_000.0])
    )
    atmosphere.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere["thermal"] = sk.constituent.ThermalEmission()
    mapping_name = "wf_thermal_temperature_k"

    atmosphere.internal_object()
    first_mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    first_interpolator = first_mapping.interpolator.copy()
    atmosphere.internal_object()
    repeated_mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    np.testing.assert_array_equal(
        repeated_mapping.interpolator, first_interpolator
    )

    atmosphere.temperature_k = np.array(
        [[280.0, 240.0, 210.0], [285.0, 245.0, 215.0]]
    )
    atmosphere.internal_object()
    native_mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    assert native_mapping.interpolator.shape == (0, 0)
    assert native_mapping.interp_dim == "location"
    assert atmosphere.derivative_output_shape(mapping_name) == (2, 3)

    atmosphere.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere.internal_object()
    broadcast_mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    assert broadcast_mapping.interpolator.shape == (6, 3)
    assert broadcast_mapping.interp_dim == "altitude"


def test_legacy_constituents_rebuild_without_double_accumulation():
    atmosphere = atmosphere2d(calculate_derivatives=False)
    atmosphere.pressure_pa = np.array([100_000.0, 30_000.0, 1_000.0])
    atmosphere.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere.internal_object()
    first_extinction = atmosphere.storage.total_extinction.copy()
    first_ssa = atmosphere.storage.ssa.copy()
    first_legendre = atmosphere.storage.leg_coeff.copy()
    atmosphere.internal_object()

    np.testing.assert_array_equal(
        atmosphere.storage.total_extinction, first_extinction
    )
    np.testing.assert_array_equal(atmosphere.storage.ssa, first_ssa)
    np.testing.assert_array_equal(atmosphere.storage.leg_coeff, first_legendre)


def test_failed_constituent_build_is_transactional_and_retryable():
    atmosphere = atmosphere2d()
    constituent = _FailingDerivativeConstituent()
    atmosphere["failing"] = constituent

    with pytest.raises(RuntimeError, match="intentional derivative failure"):
        atmosphere.internal_object()

    np.testing.assert_array_equal(atmosphere.storage.total_extinction, 0.0)
    np.testing.assert_array_equal(atmosphere.storage.ssa, 0.0)
    np.testing.assert_array_equal(atmosphere.storage.leg_coeff, 0.0)

    constituent.fail = False
    atmosphere.internal_object()
    np.testing.assert_array_equal(atmosphere.storage.total_extinction, 2.0)
    np.testing.assert_array_equal(atmosphere.storage.ssa, 0.0)


def test_native_manual_supports_polarized_legendre_storage():
    atmosphere = atmosphere2d(calculate_derivatives=False, num_stokes=3)
    extinction = np.arange(1.0, 13.0).reshape(2, 3, 2)
    ssa = np.full((2, 3, 2), 0.5)
    legendre = np.zeros(
        (atmosphere.storage.leg_coeff.shape[0], 2, 3, 2), dtype=float
    )
    legendre[0] = 1.0
    constituent = sk.constituent.Manual(extinction, ssa, legendre)
    atmosphere["manual"] = constituent

    atmosphere.internal_object()

    assert constituent.leg_coeff.shape == legendre.shape
    np.testing.assert_allclose(
        atmosphere.storage.leg_coeff,
        legendre.reshape(legendre.shape[0], 6, 2),
    )


def test_delta_m_scaling_matches_1d_for_broadcast_manual_profile():
    config = sk.Config()
    config.num_streams = 2
    config.num_singlescatter_moments = 6
    config.delta_m_scaling = True
    altitudes = geometry2d().altitudes()
    geometry_1d = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitudes,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    atmosphere_1d = sk.Atmosphere(
        geometry_1d,
        config,
        wavelengths_nm=np.array([350.0, 600.0]),
        calculate_derivatives=False,
    )
    atmosphere_2d = sk.Atmosphere(
        geometry2d(),
        config,
        wavelengths_nm=np.array([350.0, 600.0]),
        calculate_derivatives=False,
    )
    extinction = np.array([[1.0, 2.0], [3.0, 5.0], [7.0, 11.0]])
    ssa = np.full((3, 2), 0.8)
    legendre = np.zeros_like(atmosphere_1d.storage.leg_coeff)
    legendre[0] = 1.0
    legendre[2] = 0.3

    atmosphere_1d["manual"] = sk.constituent.Manual(extinction, ssa, legendre)
    atmosphere_2d["manual"] = sk.constituent.Manual(extinction, ssa, legendre)
    atmosphere_1d.internal_object()
    atmosphere_2d.internal_object()

    for horizontal_index in range(2):
        location_slice = slice(horizontal_index * 3, (horizontal_index + 1) * 3)
        np.testing.assert_allclose(
            atmosphere_2d.storage.total_extinction[location_slice],
            atmosphere_1d.storage.total_extinction,
        )
        np.testing.assert_allclose(
            atmosphere_2d.storage.ssa[location_slice], atmosphere_1d.storage.ssa
        )
        np.testing.assert_allclose(
            atmosphere_2d.storage.leg_coeff[:, location_slice],
            atmosphere_1d.storage.leg_coeff,
        )


def test_surface_constituents_remain_global_in_2d():
    atmosphere = atmosphere2d()
    atmosphere["surface"] = sk.constituent.LambertianSurface(0.37)

    atmosphere.internal_object()

    np.testing.assert_allclose(atmosphere.surface.brdf_args[0], 0.37)
    mapping = atmosphere.surface.get_derivative_mapping("wf_surface_albedo")
    assert mapping.interpolator.shape == (2, 1)


def test_raw_2d_derivatives_are_marked_with_structured_output_shape():
    atmosphere = atmosphere2d()
    atmosphere.storage.total_extinction[:] = 1.0e-5
    atmosphere.storage.ssa[:] = 0.0

    atmosphere.internal_object()

    mapping = atmosphere.storage.get_derivative_mapping("wf_extinction")
    assert mapping.interp_dim == "location"
    assert atmosphere.derivative_output_shape("wf_extinction") == (2, 3)


def test_2d_air_mass_factor_is_explicitly_deferred():
    atmosphere = atmosphere2d()
    atmosphere["amf"] = sk.constituent.AirMassFactor()

    with pytest.raises(NotImplementedError, match="AirMassFactor"):
        atmosphere.internal_object()


def test_python_engine_rejects_atmosphere_geometry_dimension_mismatch():
    config = sk.Config()
    geometry_1d = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=geometry2d().altitudes(),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    engine = sk.Engine(config, geometry_1d, sk.ViewingGeometry())

    with pytest.raises(ValueError, match="geometry dimensions do not match"):
        engine.calculate_radiance(atmosphere2d())
