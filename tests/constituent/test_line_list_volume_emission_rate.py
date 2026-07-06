from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _geometry():
    return sk.Geometry1D(
        cos_sza=-0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.array([0.0, 1000.0]),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )


def _atmosphere(wavelengths_nm, *, calculate_derivatives=False):
    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    geometry = _geometry()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.asarray(wavelengths_nm, dtype=np.float64),
        calculate_derivatives=calculate_derivatives,
    )
    atmosphere.temperature_k = np.array([200.0, 200.0])
    return atmosphere


def _integrated_emission_source(atmosphere):
    return (
        np.trapezoid(
            atmosphere.storage.emission_source,
            atmosphere.wavelengths_nm,
            axis=1,
        )
        * 4.0
        * np.pi
    )


def test_line_list_volume_emission_rate_adds_weighted_source():
    wavelengths_nm = np.arange(759.9, 760.10001, 0.0001)
    atmosphere = _atmosphere(wavelengths_nm)

    photon_ver = np.array([4.0 * np.pi, 8.0 * np.pi])
    constituent = sk.constituent.LineListVolumeEmissionRate(
        np.array([0.0, 1000.0]),
        photon_ver,
        np.array([760.0]),
        np.array([10.0]),
    )

    constituent.add_to_atmosphere(atmosphere)

    np.testing.assert_allclose(
        _integrated_emission_source(atmosphere),
        photon_ver,
        rtol=1.0e-4,
    )
    assert np.all(atmosphere.storage.emission_source.max(axis=1) > 0.0)
    np.testing.assert_allclose(constituent.weights, [[1.0], [1.0]])


def test_line_list_volume_emission_rate_accepts_altitude_dependent_weights():
    wavelengths_nm = np.arange(759.9, 761.10001, 0.0001)
    atmosphere = _atmosphere(wavelengths_nm)

    constituent = sk.constituent.LineListVolumeEmissionRate(
        np.array([0.0, 1000.0]),
        np.array([4.0 * np.pi, 8.0 * np.pi]),
        np.array([760.0, 761.0]),
        np.array([[1.0, 0.0], [0.0, 1.0]]),
    )

    constituent.add_to_atmosphere(atmosphere)

    peak_wavelengths = wavelengths_nm[
        np.argmax(atmosphere.storage.emission_source, axis=1)
    ]
    np.testing.assert_allclose(peak_wavelengths, [760.0, 761.0], atol=0.001)
    np.testing.assert_allclose(
        _integrated_emission_source(atmosphere),
        [4.0 * np.pi, 8.0 * np.pi],
        rtol=1.0e-4,
    )


def test_line_list_volume_emission_rate_registers_native_grid_weighted_derivative():
    wavelengths_nm = np.arange(759.9, 761.10001, 0.0001)
    atmosphere = _atmosphere(wavelengths_nm, calculate_derivatives=True)

    constituent = sk.constituent.LineListVolumeEmissionRate(
        np.array([0.0, 1000.0]),
        np.array([4.0 * np.pi, 8.0 * np.pi]),
        np.array([760.0, 761.0]),
        np.array([[1.0, 0.0], [0.0, 1.0]]),
    )

    constituent.add_to_atmosphere(atmosphere)
    constituent.register_derivative(atmosphere, "aband")

    mapping = atmosphere.storage.get_derivative_mapping("wf_aband_photon_ver")

    np.testing.assert_allclose(mapping.interpolator, np.eye(2))
    derivative_integrals = (
        np.trapezoid(mapping.d_emission, wavelengths_nm, axis=1) * 4.0 * np.pi
    )
    np.testing.assert_allclose(derivative_integrals, [1.0, 1.0], rtol=1.0e-4)
