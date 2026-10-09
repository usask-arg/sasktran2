from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.database.hitran_line import HITRANLineDatabase


def _has_local_o2_hitran_cache():
    db = HITRANLineDatabase()
    return (db._db_root / "O2.data").exists() and (db._db_root / "O2.header").exists()


def test_population_emission_rate_exposes_a_and_b_components():
    altitude = np.array([90_000.0, 95_000.0])
    state = xr.Dataset(
        {
            "temperature": (["altitude"], np.array([220.0, 230.0])),
            "O(1S)": (["altitude"], np.array([0.0, 0.0])),
            "O2(b)": (["altitude"], np.array([10.0, 20.0])),
            "O2(b, v=1)": (["altitude"], np.array([5.0, 5.0])),
        },
        coords={"altitude": altitude},
    )

    constituent = sk.constituent.PopulationEmissionRate(state)

    assert constituent.num_line_list_emissions == 2
    assert np.all(constituent.wavelengths_nm >= 759.0)
    assert np.all(constituent.wavelengths_nm <= 776.0)
    np.testing.assert_allclose(
        constituent.photon_ver,
        np.array([10.0, 20.0]) * 7.58e-2 + np.array([5.0, 5.0]) * 7.0e-2,
    )
    np.testing.assert_allclose(constituent.weights.sum(axis=1), 1.0)
    assert np.all(constituent.line_list_wavelengths_nm(1) >= 675.0)
    assert np.all(constituent.line_list_wavelengths_nm(1) <= 705.0)
    # B band (1-0): A = 7.34e-3 s^-1.
    np.testing.assert_allclose(constituent.line_list_photon_ver(1), 5.0 * 7.34e-3)


def test_population_emission_rate_source_integral_matches_photon_ver():
    altitude = np.array([90_000.0, 95_000.0])
    temperature = np.array([220.0, 230.0])
    state = xr.Dataset(
        {
            "temperature": (["altitude"], temperature),
            "O2(b)": (["altitude"], np.array([1.0e10, 2.0e10])),
            "O2(b, v=1)": (["altitude"], np.array([5.0e9, 5.0e9])),
        },
        coords={"altitude": altitude},
    )

    constituent = sk.constituent.PopulationEmissionRate(state)

    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    geometry = sk.Geometry1D(
        cos_sza=-0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitude,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    wavelengths_nm = np.arange(758.5, 776.5001, 0.001)
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=wavelengths_nm,
        calculate_derivatives=False,
    )
    atmosphere.temperature_k = temperature

    constituent.add_to_atmosphere(atmosphere)

    source_integral = np.trapezoid(
        atmosphere.storage.emission_source,
        wavelengths_nm,
        axis=1,
    )
    np.testing.assert_allclose(
        source_integral,
        constituent.photon_ver / (4.0 * np.pi),
        rtol=1.0e-3,
    )


def test_population_emission_rate_line_strength_fallback():
    altitude = np.array([90_000.0, 95_000.0])
    state = xr.Dataset(
        {
            "temperature": (["altitude"], np.array([220.0, 230.0])),
            "O2(b)": (["altitude"], np.array([10.0, 20.0])),
            "O2(b, v=1)": (["altitude"], np.array([5.0, 5.0])),
        },
        coords={"altitude": altitude},
    )

    constituent = sk.constituent.PopulationEmissionRate(
        state,
        line_weight_model="hitran_line_strength",
    )

    np.testing.assert_allclose(constituent.weights.sum(axis=1), 1.0)
    assert np.all(constituent.weights >= 0.0)
    assert constituent.num_line_list_emissions == 2
    np.testing.assert_allclose(constituent.line_list_weights(1).sum(axis=1), 1.0)


def test_oxygen_a_band_emission_absorption_engine_smoke():
    pytest.importorskip("hapi")

    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate

    altitude_grid_m = np.arange(0.0, 121_000.0, 5_000.0)
    geometry = sk.Geometry1D(
        cos_sza=-0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitude_grid_m,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()
    viewing_geo.add_ray(
        sk.TangentAltitudeSolar(
            tangent_altitude_m=95_000.0,
            relative_azimuth=0.0,
            observer_altitude_m=200_000.0,
            cos_sza=-0.6,
        )
    )

    wavelengths_nm = np.arange(759.0, 776.1, 0.1)
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=wavelengths_nm,
        calculate_derivatives=False,
    )
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    o2_b_population = 2.5e10 * np.exp(
        -0.5 * ((altitude_grid_m - 95_000.0) / 8_000.0) ** 2
    )
    state = xr.Dataset(
        {
            "temperature": ("altitude", atmosphere.temperature_k),
            "O2(b)": ("altitude", o2_b_population),
        },
        coords={"altitude": altitude_grid_m},
    )

    atmosphere["o2"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.HITRANAbsorber("O2"),
        altitude_grid_m,
        np.full_like(altitude_grid_m, 0.21),
        out_of_bounds_mode="extend",
    )
    atmosphere["photochemical_emission"] = sk.constituent.PopulationEmissionRate(state)

    engine = sk.Engine(config, geometry, viewing_geo)
    radiance = engine.calculate_radiance(atmosphere)

    radiance_i = radiance["radiance"].sel(stokes="I").to_numpy()
    assert np.isfinite(radiance_i).all()
    assert radiance_i.max() > 0.0
    assert atmosphere.storage.total_extinction.max() > 0.0
    assert atmosphere.storage.emission_source.max() > 0.0
