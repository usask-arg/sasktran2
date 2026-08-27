from __future__ import annotations

import hashlib

import numpy as np
import sasktran2 as sk
from sasktran2.constituent.base import Constituent


class _VMRProfile(Constituent):
    def __init__(self, altitudes_m: np.ndarray, vmr: np.ndarray):
        self.altitudes_m = altitudes_m
        self.vmr = vmr

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        pass

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        pass


def _inputs():
    return (
        np.array([101325.0]),
        np.array([288.15]),
        np.array([0.01]),
        np.array([400.0e-6]),
        np.array([2.0e-6]),
    )


def test_first_data_download_prints_and_caches_aer_license(
    tmp_path, monkeypatch, capsys
):
    coefficient_data = b"coefficient data"
    license_data = b"AER license"
    files = {
        sk.continuum._MT_CKD_DATA_KEY: coefficient_data,
        sk.continuum._MT_CKD_LICENSE_KEY: license_data,
    }

    class _TestDatabase:
        def path(self, key):
            path = tmp_path.joinpath(key)
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                path.write_bytes(files[key])
            return path

    monkeypatch.setattr(sk.appconfig, "database_root", lambda: tmp_path)
    monkeypatch.setattr(sk.database, "StandardDatabase", _TestDatabase)
    monkeypatch.setattr(
        sk.continuum,
        "_MT_CKD_DATA_SHA256",
        hashlib.sha256(coefficient_data).hexdigest(),
    )
    monkeypatch.setattr(
        sk.continuum,
        "_MT_CKD_LICENSE_SHA256",
        hashlib.sha256(license_data).hexdigest(),
    )

    sk.continuum._mt_ckd_data_file.cache_clear()
    try:
        assert sk.continuum._mt_ckd_data_file().read_bytes() == coefficient_data
        output = capsys.readouterr().out
        assert "licensed separately from SASKTRAN2" in output
        assert sk.continuum.MT_CKD_AER_LICENSE in output
        assert tmp_path.joinpath(sk.continuum._MT_CKD_LICENSE_KEY).exists()

        sk.continuum._mt_ckd_data_file()
        assert capsys.readouterr().out == ""
    finally:
        sk.continuum._mt_ckd_data_file.cache_clear()


def test_mt_ckd_fixed_grid_interface():
    result = sk.mt_ckd(*_inputs())

    assert result.shape == (1, 1991)
    np.testing.assert_array_equal(
        sk.MT_CKD_WAVENUMBERS_CM_INV,
        np.arange(0.0, 19_901.0, 10.0),
    )
    np.testing.assert_allclose(
        result[0, [1, 100, 240, 800, 1300, 1990]],
        [
            7.904694654058591e-4,
            4.542272170306234e-5,
            2.2776831103477272e-4,
            4.661442799889581e-6,
            4.317920214252337e-6,
            2.3854890426228372e-5,
        ],
        rtol=5e-14,
        atol=1e-16,
    )


def test_mt_ckd_analytic_jacobians():
    inputs = _inputs()
    result = sk.mt_ckd_linearized(*inputs)
    perturbations = (1.0, 1.0e-3, 1.0e-7, 1.0e-8, 1.0e-8)

    assert len(result) == 6
    for input_index, (analytic, delta) in enumerate(
        zip(result[1:], perturbations, strict=True)
    ):
        above = [value.copy() for value in inputs]
        below = [value.copy() for value in inputs]
        above[input_index] += delta
        below[input_index] -= delta
        numeric = (sk.mt_ckd(*above) - sk.mt_ckd(*below)) / (2.0 * delta)
        np.testing.assert_allclose(
            analytic[:, 1:],
            numeric[:, 1:],
            rtol=3e-5,
            atol=1e-12,
        )


def test_mt_ckd_pure_water_branch_matches_upstream_and_is_linearized():
    inputs = list(_inputs())
    inputs[1] = np.array([296.0])
    inputs[2] = np.array([1.0])
    result = sk.mt_ckd_linearized(*inputs)

    assert all(np.all(np.isfinite(field)) for field in result)
    np.testing.assert_allclose(
        result[0][0, [1, 100, 240, 800, 1000, 1300, 1990]],
        [
            1.7199222942376495,
            0.3249611845717048,
            0.02554401487894435,
            0.000733663141266962,
            9.614781340354584e-05,
            5.429603526351588e-06,
            3.431332907137031e-05,
        ],
        rtol=5e-14,
        atol=1e-16,
    )

    delta = 1e-7
    above = [value.copy() for value in inputs]
    below = [value.copy() for value in inputs]
    above[2] += delta
    below[2] -= delta
    numeric = (sk.mt_ckd(*above) - sk.mt_ckd(*below)) / (2.0 * delta)
    np.testing.assert_allclose(
        result[3][:, 1:],
        numeric[:, 1:],
        rtol=3e-6,
        atol=1e-12,
    )


def test_full_spectra_match_mt_ckd_4_3_double_precision_reference():
    reference_path = sk.database.StandardDatabase().path(
        "continuum/mt_ckd_4_3_reference.bin"
    )
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == (
        "0fdc39229ae020f2979511acfc81600ddfaa0b9bc0bbb03bb86a5ce3d01a9760"
    )
    expected = np.fromfile(reference_path, dtype="<f8").reshape(4, 1991)
    scenarios = (
        (101_325.0, 288.15, 0.01, 400.0e-6, 2.0e-6, 100.0),
        (50_000.0, 250.0, 0.002, 420.0e-6, 1.0e-6, 100.0),
        (8_000.0, 220.0, 5.0e-5, 380.0e-6, 5.0e-6, 37.5),
        (90_000.0, 310.0, 0.2, 1.0e-3, 1.0e-5, 250.0),
    )

    for index, scenario in enumerate(scenarios):
        pressure, temperature, h2o, co2, o3, reference_path_length_cm = scenario
        actual = sk.mt_ckd(
            np.array([pressure]),
            np.array([temperature]),
            np.array([h2o]),
            np.array([co2]),
            np.array([o3]),
        )[0]
        expected_extinction = expected[index] * (100.0 / reference_path_length_cm)
        relative_error = np.abs(actual - expected_extinction) / np.maximum(
            np.abs(expected_extinction), 2.0e-15
        )
        assert np.max(relative_error) <= 5.0e-14


def test_compatibility_constructor_options_are_accepted():
    continuum = sk.continuum.MTCKDContinuum(
        h2o_name="water",
        co2_name="carbon_dioxide",
        o3_name="ozone",
        numeric_wf_fractional_change=0.2,
        numeric_wf_central_difference=False,
    )

    assert continuum._h2o_name == "water"


def test_constituent_adds_extinction_and_analytic_mappings():
    altitudes = np.array([0.0, 10_000.0, 30_000.0])
    wavenumbers = np.array([1000.0, 2400.0, 10_000.0])
    geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitudes,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        sk.Config(),
        wavenumber_cminv=wavenumbers,
        calculate_derivatives=True,
    )
    atmosphere.pressure_pa = np.array([101_325.0, 50_000.0, 8_000.0])
    atmosphere.temperature_k = np.array([288.15, 250.0, 220.0])
    atmosphere["H2O"] = _VMRProfile(altitudes, np.array([0.01, 0.002, 5e-5]))
    atmosphere["CO2"] = _VMRProfile(altitudes, np.array([400e-6, 420e-6, 380e-6]))
    atmosphere["O3"] = _VMRProfile(altitudes, np.array([2e-6, 1e-6, 5e-6]))
    atmosphere["continuum"] = sk.continuum.MTCKDContinuum()

    atmosphere.internal_object()

    low_level = sk.mt_ckd_linearized(
        atmosphere.pressure_pa,
        atmosphere.temperature_k,
        atmosphere["H2O"].vmr,
        atmosphere["CO2"].vmr,
        atmosphere["O3"].vmr,
        include_rayleigh=False,
    )
    indices = (wavenumbers / 10.0).astype(int)
    np.testing.assert_allclose(
        atmosphere.storage.total_extinction,
        low_level[0][:, indices],
        rtol=2e-14,
    )

    mapping_specs = (
        ("wf_continuum_pressure_pa", "wf_pressure_pa"),
        ("wf_continuum_temperature_k", "wf_temperature_k"),
        ("wf_continuum_H2O_vmr", "wf_H2O_vmr"),
        ("wf_continuum_CO2_vmr", "wf_CO2_vmr"),
        ("wf_continuum_O3_vmr", "wf_O3_vmr"),
    )
    for (mapping_name, assign_name), expected in zip(
        mapping_specs, low_level[1:], strict=True
    ):
        mapping = atmosphere.storage.get_derivative_mapping(mapping_name)
        np.testing.assert_allclose(
            mapping.d_extinction,
            expected[:, indices],
            rtol=2e-14,
            atol=1e-18,
        )
        assert mapping.assign_name == assign_name


def test_constituent_without_derivatives_uses_value_only_path(monkeypatch):
    altitudes = np.array([0.0, 10_000.0])
    geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitudes,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        sk.Config(),
        wavenumber_cminv=np.array([1000.0]),
        calculate_derivatives=False,
    )
    atmosphere.pressure_pa = np.array([101_325.0, 50_000.0])
    atmosphere.temperature_k = np.array([288.15, 250.0])
    atmosphere["H2O"] = _VMRProfile(altitudes, np.array([0.01, 0.002]))
    atmosphere["CO2"] = _VMRProfile(altitudes, np.full(2, 400e-6))
    atmosphere["O3"] = _VMRProfile(altitudes, np.full(2, 2e-6))

    value_calls = 0
    value_function = sk.continuum.mt_ckd

    def tracked_value(*args, **kwargs):
        nonlocal value_calls
        value_calls += 1
        return value_function(*args, **kwargs)

    def unexpected_linearized(*args, **kwargs):
        msg = "the linearized MT_CKD path should not be evaluated"
        raise AssertionError(msg)

    monkeypatch.setattr(sk.continuum, "mt_ckd", tracked_value)
    monkeypatch.setattr(sk.continuum, "mt_ckd_linearized", unexpected_linearized)
    atmosphere["continuum"] = sk.continuum.MTCKDContinuum()

    atmosphere.internal_object()

    assert value_calls == 1
    assert np.all(atmosphere.storage.total_extinction > 0.0)


def test_species_derivative_mapping_chains_to_profile_grid():
    altitudes = np.array([0.0, 10_000.0, 30_000.0])
    h2o_altitudes = np.array([0.0, 30_000.0])
    h2o_vmr = np.array([0.01, 5e-5])
    wavenumbers = np.array([1000.0, 2400.0])
    geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitudes,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        sk.Config(),
        wavenumber_cminv=wavenumbers,
        calculate_derivatives=True,
    )
    atmosphere.pressure_pa = np.array([101_325.0, 50_000.0, 8_000.0])
    atmosphere.temperature_k = np.array([288.15, 250.0, 220.0])
    atmosphere["H2O"] = _VMRProfile(h2o_altitudes, h2o_vmr)
    atmosphere["CO2"] = _VMRProfile(altitudes, np.full(3, 400e-6))
    atmosphere["O3"] = _VMRProfile(altitudes, np.full(3, 2e-6))
    atmosphere["continuum"] = sk.continuum.MTCKDContinuum()
    atmosphere.internal_object()

    expected_interpolator = np.array(
        [
            [1.0, 0.0],
            [2.0 / 3.0, 1.0 / 3.0],
            [0.0, 1.0],
        ]
    )
    mapping = atmosphere.storage.get_derivative_mapping("wf_continuum_H2O_vmr")
    np.testing.assert_allclose(mapping.interpolator, expected_interpolator)
    assert mapping.interp_dim == "H2O_altitude"
    assert mapping.assign_name == "wf_H2O_vmr"

    delta = 1e-7
    spectral_indices = (wavenumbers / 10.0).astype(int)
    for parameter_index in range(h2o_vmr.size):
        above = h2o_vmr.copy()
        below = h2o_vmr.copy()
        above[parameter_index] += delta
        below[parameter_index] -= delta
        numeric = (
            sk.mt_ckd(
                atmosphere.pressure_pa,
                atmosphere.temperature_k,
                expected_interpolator @ above,
                atmosphere["CO2"].vmr,
                atmosphere["O3"].vmr,
                include_rayleigh=False,
            )[:, spectral_indices]
            - sk.mt_ckd(
                atmosphere.pressure_pa,
                atmosphere.temperature_k,
                expected_interpolator @ below,
                atmosphere["CO2"].vmr,
                atmosphere["O3"].vmr,
                include_rayleigh=False,
            )[:, spectral_indices]
        ) / (2.0 * delta)
        analytic = (
            mapping.d_extinction * expected_interpolator[:, parameter_index, np.newaxis]
        )
        np.testing.assert_allclose(analytic, numeric, rtol=3e-5, atol=1e-12)


def test_constituent_excludes_legacy_rayleigh_term():
    altitudes = np.array([0.0, 10_000.0])
    wavenumber = np.array([10_000.0])
    geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitudes,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        sk.Config(),
        wavenumber_cminv=wavenumber,
        calculate_derivatives=False,
    )
    atmosphere.pressure_pa = np.array([101_325.0, 50_000.0])
    atmosphere.temperature_k = np.array([288.15, 250.0])
    atmosphere["H2O"] = _VMRProfile(altitudes, np.array([0.01, 0.002]))
    atmosphere["CO2"] = _VMRProfile(altitudes, np.array([400e-6, 420e-6]))
    atmosphere["O3"] = _VMRProfile(altitudes, np.array([2e-6, 1e-6]))
    atmosphere["continuum"] = sk.continuum.MTCKDContinuum()

    atmosphere.internal_object()

    inputs = _inputs()
    full = sk.mt_ckd(*inputs)[0, 1000]
    absorption = sk.mt_ckd(*inputs, include_rayleigh=False)[0, 1000]
    assert full > 5.0 * absorption
    np.testing.assert_allclose(
        atmosphere.storage.total_extinction[0, 0],
        absorption,
        rtol=2e-14,
    )


def test_constituent_supports_expanded_profiles_in_2d_atmosphere():
    altitudes = np.array([0.0, 10_000.0, 30_000.0])
    horizontal_angles = np.array([-0.2, 0.3])
    wavenumbers = np.array([1000.0, 10_000.0])
    geometry = sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitudes,
        horizontal_angle_grid_radians=horizontal_angles,
    )
    config = sk.Config()
    config.num_streams = 2
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavenumber_cminv=wavenumbers,
        calculate_derivatives=True,
    )
    pressure = np.array([101_325.0, 50_000.0, 8_000.0])
    temperature = np.array([288.15, 250.0, 220.0])
    h2o = np.array([0.01, 0.002, 5e-5])
    co2 = np.array([400e-6, 420e-6, 380e-6])
    o3 = np.array([2e-6, 1e-6, 5e-6])
    atmosphere.pressure_pa = pressure
    atmosphere.temperature_k = temperature
    atmosphere["H2O"] = _VMRProfile(altitudes, h2o)
    atmosphere["CO2"] = _VMRProfile(altitudes, co2)
    atmosphere["O3"] = _VMRProfile(altitudes, o3)
    atmosphere["continuum"] = sk.continuum.MTCKDContinuum()

    atmosphere.internal_object()

    repeats = len(horizontal_angles)
    low_level = sk.mt_ckd_linearized(
        np.tile(pressure, repeats),
        np.tile(temperature, repeats),
        np.tile(h2o, repeats),
        np.tile(co2, repeats),
        np.tile(o3, repeats),
        include_rayleigh=False,
    )
    spectral_indices = (wavenumbers / 10.0).astype(int)
    np.testing.assert_allclose(
        atmosphere.storage.total_extinction,
        low_level[0][:, spectral_indices],
        rtol=2e-14,
    )

    mapping = atmosphere.storage.get_derivative_mapping("wf_continuum_H2O_vmr")
    np.testing.assert_allclose(
        mapping.d_extinction,
        low_level[3][:, spectral_indices],
        rtol=2e-14,
        atol=1e-18,
    )
    np.testing.assert_allclose(
        mapping.interpolator,
        np.tile(np.eye(altitudes.size), (repeats, 1)),
    )
