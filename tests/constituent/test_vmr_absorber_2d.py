from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk

from sasktran2.constants import K_BOLTZMANN
from sasktran2.optical.base import OpticalProperty, OpticalQuantities


ALTITUDES_M = np.array([0.0, 10_000.0, 30_000.0])
HORIZONTAL_ANGLES = np.array([-0.5, 0.0, 0.5])
WAVELENGTHS_NM = np.array([500.0, 600.0])
CROSS_SECTIONS_M2 = np.array([2.0e-25, 5.0e-25])


class _ConstantAbsorber(OpticalProperty):
    def __init__(self) -> None:
        self.seen_vmr: list[np.ndarray] = []

    def atmosphere_quantities(
        self, atmo: sk.Atmosphere, **kwargs
    ) -> OpticalQuantities:
        self.seen_vmr.append(np.asarray(kwargs["vmr"]).copy())
        cross_section = np.broadcast_to(
            CROSS_SECTIONS_M2, (atmo.num_locations, CROSS_SECTIONS_M2.size)
        ).copy()
        return OpticalQuantities(
            extinction=cross_section,
            ssa=np.zeros_like(cross_section),
        )


def _config(*, transmission: bool = False) -> sk.Config:
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    config.occultation_source = (
        sk.OccultationSource.Standard
        if transmission
        else sk.OccultationSource.NoSource
    )
    return config


def _geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=HORIZONTAL_ANGLES,
    )


def _atmosphere(
    geometry: sk.Geometry2D,
    *,
    calculate_derivatives: bool = True,
    transmission: bool = False,
) -> sk.Atmosphere:
    atmosphere = sk.Atmosphere(
        geometry,
        _config(transmission=transmission),
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=calculate_derivatives,
        legendre_derivative=False,
    )
    atmosphere.pressure_pa = np.array(
        [
            [100_000.0, 30_000.0, 1_000.0],
            [90_000.0, 25_000.0, 800.0],
            [80_000.0, 20_000.0, 600.0],
        ]
    )
    atmosphere.temperature_k = np.array(
        [[280.0, 240.0, 210.0], [285.0, 245.0, 215.0], [290.0, 250.0, 220.0]]
    )
    return atmosphere


def _vmr() -> np.ndarray:
    return np.array(
        [
            [1.0e-6, 2.0e-6, 3.0e-6],
            [1.5e-6, 2.5e-6, 3.5e-6],
            [2.0e-6, 3.0e-6, 4.0e-6],
        ]
    )


def test_construction_preserves_native_shape_and_mutability():
    vmr = _vmr()
    constituent = sk.constituent.VMRAbsorber2D(_ConstantAbsorber(), vmr)

    assert constituent.volume_spatial_mode == "native_2d"
    np.testing.assert_array_equal(constituent.vmr, vmr)

    constituent.vmr[1, 2] *= 2.0
    assert constituent.vmr[1, 2] == 2.0 * vmr[1, 2]

    replacement = vmr * 3.0
    constituent.vmr = replacement
    np.testing.assert_array_equal(constituent.vmr, replacement)


@pytest.mark.parametrize("vmr", [np.ones(3), np.ones((1, 2, 3)), np.ones((0, 3))])
def test_construction_rejects_non_native_2d_shapes(vmr):
    with pytest.raises(ValueError, match="vmr"):
        sk.constituent.VMRAbsorber2D(_ConstantAbsorber(), vmr)


def test_vmr_setter_rejects_shape_changes():
    constituent = sk.constituent.VMRAbsorber2D(_ConstantAbsorber(), _vmr())
    with pytest.raises(ValueError, match="retain shape"):
        constituent.vmr = np.ones((1, 9))


def test_requires_geometry2d_and_matching_native_shape():
    constituent = sk.constituent.VMRAbsorber2D(_ConstantAbsorber(), _vmr())
    geometry1d = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    atmosphere1d = sk.Atmosphere(
        geometry1d,
        _config(),
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    atmosphere1d.pressure_pa = np.array([100_000.0, 30_000.0, 1_000.0])
    atmosphere1d.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere1d["gas"] = constituent
    with pytest.raises(TypeError, match="Geometry2D"):
        atmosphere1d.internal_object()

    geometry2d = _geometry2d()
    atmosphere2d = _atmosphere(geometry2d, calculate_derivatives=False)
    atmosphere2d["gas"] = sk.constituent.VMRAbsorber2D(
        _ConstantAbsorber(), np.ones((2, 3))
    )
    with pytest.raises(ValueError, match="does not match"):
        atmosphere2d.internal_object()


def test_native_vmr_populates_expected_extinction_and_is_passed_to_optics():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry, calculate_derivatives=False)
    optical_property = _ConstantAbsorber()
    vmr = _vmr()
    atmosphere["gas"] = sk.constituent.VMRAbsorber2D(optical_property, vmr)

    atmosphere.internal_object()

    number_density = (
        atmosphere.pressure_pa / (K_BOLTZMANN * atmosphere.temperature_k)
    ).reshape(-1)
    expected = (
        number_density[:, np.newaxis]
        * vmr.reshape(-1, 1)
        * CROSS_SECTIONS_M2[np.newaxis, :]
    )
    np.testing.assert_allclose(atmosphere.storage.total_extinction, expected)
    np.testing.assert_array_equal(optical_property.seen_vmr[-1], vmr.reshape(-1))


def test_uniform_native_vmr_matches_altitude_profile_absorber():
    geometry2d = _geometry2d()
    geometry1d = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    profile = np.array([1.0e-6, 2.0e-6, 3.0e-6])
    pressure = np.array([100_000.0, 30_000.0, 1_000.0])
    temperature = np.array([280.0, 240.0, 210.0])

    atmosphere1d = sk.Atmosphere(
        geometry1d,
        _config(),
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    atmosphere1d.pressure_pa = pressure
    atmosphere1d.temperature_k = temperature
    atmosphere1d["gas"] = sk.constituent.VMRAltitudeAbsorber(
        _ConstantAbsorber(), ALTITUDES_M, profile
    )

    atmosphere2d = sk.Atmosphere(
        geometry2d,
        _config(),
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=False,
    )
    atmosphere2d.pressure_pa = pressure
    atmosphere2d.temperature_k = temperature
    atmosphere2d["gas"] = sk.constituent.VMRAbsorber2D(
        _ConstantAbsorber(), np.tile(profile, (3, 1))
    )

    atmosphere1d.internal_object()
    atmosphere2d.internal_object()
    np.testing.assert_allclose(
        atmosphere2d.storage.total_extinction,
        np.tile(atmosphere1d.storage.total_extinction, (3, 1)),
    )


def test_native_vmr_with_production_optical_database_matches_1d_columns():
    geometry2d = _geometry2d()
    wavelengths_nm = np.array([310.0, 330.0])
    atmosphere2d = sk.Atmosphere(
        geometry2d,
        _config(),
        wavelengths_nm=wavelengths_nm,
        calculate_derivatives=False,
    )
    pressure = np.array(
        [
            [100_000.0, 30_000.0, 1_000.0],
            [90_000.0, 25_000.0, 800.0],
            [80_000.0, 20_000.0, 600.0],
        ]
    )
    temperature = np.array(
        [[280.0, 240.0, 210.0], [285.0, 245.0, 215.0], [290.0, 250.0, 220.0]]
    )
    vmr = _vmr()
    atmosphere2d.pressure_pa = pressure
    atmosphere2d.temperature_k = temperature
    atmosphere2d["ozone"] = sk.constituent.VMRAbsorber2D(
        sk.optical.O3DBM(), vmr
    )
    atmosphere2d.internal_object()
    extinction2d = atmosphere2d.reshape_native(
        atmosphere2d.storage.total_extinction
    )

    geometry1d = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    for horizontal_index in range(geometry2d.shape[0]):
        atmosphere1d = sk.Atmosphere(
            geometry1d,
            _config(),
            wavelengths_nm=wavelengths_nm,
            calculate_derivatives=False,
        )
        atmosphere1d.pressure_pa = pressure[horizontal_index]
        atmosphere1d.temperature_k = temperature[horizontal_index]
        atmosphere1d["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), ALTITUDES_M, vmr[horizontal_index]
        )
        atmosphere1d.internal_object()
        np.testing.assert_allclose(
            extinction2d[horizontal_index],
            atmosphere1d.storage.total_extinction,
            rtol=2.0e-14,
        )


def test_native_vmr_derivative_is_local_finite_and_matches_storage_difference():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry)
    constituent = sk.constituent.VMRAbsorber2D(_ConstantAbsorber(), _vmr())
    atmosphere["gas"] = constituent
    atmosphere.internal_object()

    base_extinction = atmosphere.storage.total_extinction.copy()
    mapping_name = "wf_gas_vmr"
    mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
    assert mapping.interp_dim == "location"
    assert mapping.interpolator.shape == (0, 0)
    assert atmosphere.derivative_output_shape(mapping_name) == geometry.shape
    assert np.all(np.isfinite(mapping.d_extinction))
    assert np.all(np.isfinite(mapping.d_ssa))

    horizontal_index = 1
    altitude_index = 1
    location_index = geometry.location_index(altitude_index, horizontal_index)
    delta = 1.0e-9
    constituent.vmr[horizontal_index, altitude_index] += delta
    atmosphere.internal_object()
    numeric = (atmosphere.storage.total_extinction - base_extinction) / delta

    expected = np.zeros_like(numeric)
    expected[location_index] = mapping.d_extinction[location_index]
    np.testing.assert_allclose(numeric, expected, rtol=2.0e-12, atol=1.0e-18)


def test_native_pressure_derivative_matches_storage_difference():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry)
    atmosphere["gas"] = sk.constituent.VMRAbsorber2D(
        _ConstantAbsorber(), _vmr()
    )
    atmosphere.internal_object()

    base_extinction = atmosphere.storage.total_extinction.copy()
    mapping = atmosphere.storage.get_derivative_mapping("wf_gas_pressure_pa")
    horizontal_index = 1
    altitude_index = 1
    location_index = geometry.location_index(altitude_index, horizontal_index)
    analytic = (
        mapping.d_extinction
        * mapping.interpolator[:, location_index, np.newaxis]
    )

    delta = 0.1
    atmosphere.pressure_pa[horizontal_index, altitude_index] += delta
    atmosphere.internal_object()
    numeric = (atmosphere.storage.total_extinction - base_extinction) / delta
    np.testing.assert_allclose(numeric, analytic, rtol=2.0e-9, atol=1.0e-22)


def test_broadcast_pressure_derivative_retains_altitude_dimension():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry)
    atmosphere.pressure_pa = np.array([100_000.0, 30_000.0, 1_000.0])
    atmosphere.temperature_k = np.array([280.0, 240.0, 210.0])
    atmosphere["gas"] = sk.constituent.VMRAbsorber2D(
        _ConstantAbsorber(), _vmr()
    )
    atmosphere.internal_object()

    mapping = atmosphere.storage.get_derivative_mapping("wf_gas_pressure_pa")
    assert mapping.interp_dim == "altitude"
    assert mapping.interpolator.shape == (9, 3)


def test_zero_vmr_has_finite_derivative_mappings():
    geometry = _geometry2d()
    atmosphere = _atmosphere(geometry)
    atmosphere["gas"] = sk.constituent.VMRAbsorber2D(
        _ConstantAbsorber(), np.zeros(geometry.shape)
    )
    atmosphere.internal_object()

    mapping = atmosphere.storage.get_derivative_mapping("wf_gas_vmr")
    assert np.all(np.isfinite(mapping.d_extinction))
    assert np.all(np.isfinite(mapping.d_ssa))
    np.testing.assert_array_equal(mapping.d_ssa, 0.0)


def test_transmission_vmr_weighting_function_matches_finite_difference():
    geometry = _geometry2d()
    config = _config(transmission=True)
    atmosphere = _atmosphere(geometry, transmission=True)
    constituent = sk.constituent.VMRAbsorber2D(_ConstantAbsorber(), _vmr())
    atmosphere["gas"] = constituent
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.TangentAltitudeSolar(
            tangent_altitude_m=15_000.0,
            relative_azimuth=0.0,
            observer_altitude_m=100_000.0,
            cos_sza=0.6,
        )
    )
    engine = sk.Engine(config, geometry, viewing)

    base = engine.calculate_radiance(atmosphere)
    assert base.wf_gas_vmr.dims == (
        "horizontal_angle",
        "altitude",
        "wavelength",
        "los",
        "stokes",
    )
    assert base.wf_gas_vmr.shape == (*geometry.shape, 2, 1, 1)

    horizontal_index = 1
    altitude_index = 1
    delta = 1.0e-10
    constituent.vmr[horizontal_index, altitude_index] += delta
    perturbed_above = engine.calculate_radiance(atmosphere)
    constituent.vmr[horizontal_index, altitude_index] -= 2.0 * delta
    perturbed_below = engine.calculate_radiance(atmosphere)
    numeric = (perturbed_above.radiance - perturbed_below.radiance) / (2.0 * delta)
    analytic = base.wf_gas_vmr[horizontal_index, altitude_index]
    np.testing.assert_allclose(numeric, analytic, rtol=3.0e-7, atol=1.0e-10)
