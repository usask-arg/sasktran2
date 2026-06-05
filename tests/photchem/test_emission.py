from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.database.hitran_line import HITRANLineDatabase
from sasktran2.photchem import Yankovsky


def _has_local_o2_hitran_cache():
    db = HITRANLineDatabase()
    return (db._db_root / "O2.data").exists() and (db._db_root / "O2.header").exists()


def test_yankovsky_emission_rates_from_state_profiles():
    altitude = np.array([90_000.0, 95_000.0])
    state = xr.Dataset(
        {
            "O(1S)": (["altitude"], np.array([1.0, 2.0])),
            "O2(b)": (["altitude"], np.array([10.0, 20.0])),
            "O2(b, v=1)": (["altitude"], np.array([1.0, 2.0])),
        },
        coords={"altitude": altitude},
    )

    emissions = Yankovsky().emissions(state=state)

    np.testing.assert_allclose(
        emissions["oxygen_green_5577_photon_ver"].to_numpy(),
        np.array([1.26, 2.52]),
    )
    np.testing.assert_allclose(
        emissions["oxygen_a_band_photon_ver"].to_numpy(),
        np.array([0.828, 1.656]),
    )
    np.testing.assert_allclose(
        emissions["oxygen_a_band_b0_photon_ver"].to_numpy(),
        np.array([0.758, 1.516]),
    )
    np.testing.assert_allclose(
        emissions["oxygen_a_band_b1_photon_ver"].to_numpy(),
        np.array([0.07, 0.14]),
    )
    np.testing.assert_allclose(emissions["altitude"].to_numpy(), altitude)
    assert emissions.attrs["emission_units"] == "photons m^-3 s^-1"


def test_oxygen_green_line_mcdade_parameterization():
    altitude = np.array([90_000.0])
    atmosphere_state = xr.Dataset(
        {
            "temperature": (["altitude"], np.array([200.0])),
            "o_density": (["altitude"], np.array([5.0e17])),
            "o2_density": (["altitude"], np.array([1.0e18])),
            "n2_density": (["altitude"], np.array([4.0e18])),
        },
        coords={"altitude": altitude},
    )

    green_line = Yankovsky().oxygen_green_line_mcdade(atmosphere_state)

    atomic_o_cm3 = 5.0e11
    o2_cm3 = 1.0e12
    n2_cm3 = 4.0e12
    k1 = 4.7e-33 * (300.0 / 200.0) ** 2
    three_k5 = 4.0e-12 * np.exp(-865.0 / 200.0)
    expected_cm3 = (
        k1
        * atomic_o_cm3**2
        * (n2_cm3 + o2_cm3)
        * atomic_o_cm3
        / (211.0 * atomic_o_cm3 + 15.0 * o2_cm3)
        * 1.18
        / (1.35 + three_k5 * o2_cm3)
    )

    np.testing.assert_allclose(
        green_line["oxygen_green_5577_photon_ver"].to_numpy(),
        [expected_cm3 * 1.0e6],
    )
    np.testing.assert_allclose(
        green_line["O(1S)"].to_numpy(), [expected_cm3 * 1.0e6 / 1.18]
    )
    assert green_line.attrs["model"].startswith("McDade/Murtagh")


def test_oxygen_green_line_mcdade_constituent_uses_parameterized_ver():
    altitude = np.array([90_000.0, 95_000.0])
    atmosphere_state = xr.Dataset(
        {
            "temperature": (["altitude"], np.array([200.0, 210.0])),
            "o_density": (["altitude"], np.array([5.0e17, 6.0e17])),
            "o2_density": (["altitude"], np.array([1.0e18, 0.8e18])),
            "n2_density": (["altitude"], np.array([4.0e18, 3.2e18])),
        },
        coords={"altitude": altitude},
    )

    model = Yankovsky()
    green_line = model.oxygen_green_line_mcdade(atmosphere_state)
    constituent = model.oxygen_green_line_mcdade_constituent(green_line=green_line)

    np.testing.assert_allclose(
        constituent.ver,
        green_line["oxygen_green_5577_photon_ver"].to_numpy(),
    )
    np.testing.assert_allclose(constituent.altitudes_m, altitude)


@pytest.mark.skipif(
    not _has_local_o2_hitran_cache(),
    reason="O2 HITRAN line cache is required for the O2 band constituent",
)
def test_yankovsky_emission_constituents_include_green_and_a_band():
    altitude = np.array([90_000.0, 95_000.0])
    state = xr.Dataset(
        {
            "temperature": (["altitude"], np.array([220.0, 230.0])),
            "O(1S)": (["altitude"], np.array([1.0, 2.0])),
            "O2(b)": (["altitude"], np.array([10.0, 20.0])),
            "O2(b, v=1)": (["altitude"], np.array([5.0, 5.0])),
        },
        coords={"altitude": altitude},
    )

    model = Yankovsky()
    constituents = model.emission_constituents(state=state)

    assert set(constituents) == {"oxygen_green", "oxygen_a_band"}
    np.testing.assert_allclose(constituents["oxygen_green"].ver, [1.26, 2.52])
    assert constituents["oxygen_a_band"].num_line_list_emissions == 2
    np.testing.assert_allclose(constituents["oxygen_a_band"].weights.sum(axis=1), 1.0)
    np.testing.assert_allclose(
        constituents["oxygen_a_band"].photon_ver,
        np.array([10.0, 20.0]) * 7.58e-2 + np.array([5.0, 5.0]) * 7.0e-2,
    )

    line_strength_constituent = model.oxygen_a_band_constituent(
        state=state,
        line_weight_model="hitran_line_strength",
    )
    np.testing.assert_allclose(
        line_strength_constituent.weights.sum(axis=1),
        1.0,
    )


@pytest.mark.skipif(
    not _has_local_o2_hitran_cache(),
    reason="O2 HITRAN line cache is required for the population emission constituent",
)
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
    assert np.all(constituent.wavelengths_nm <= 772.0)
    np.testing.assert_allclose(
        constituent.photon_ver,
        np.array([10.0, 20.0]) * 7.58e-2 + np.array([5.0, 5.0]) * 7.0e-2,
    )
    np.testing.assert_allclose(constituent.weights.sum(axis=1), 1.0)
    assert np.all(constituent.line_list_wavelengths_nm(1) >= 680.0)
    assert np.all(constituent.line_list_wavelengths_nm(1) <= 700.0)
    np.testing.assert_allclose(constituent.line_list_photon_ver(1), 5.0 * 7.0e-2)


@pytest.mark.skipif(
    not _has_local_o2_hitran_cache(),
    reason="O2 HITRAN line cache is required for the population emission constituent",
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


@pytest.mark.skipif(
    not _has_local_o2_hitran_cache(),
    reason="O2 HITRAN line cache is required for the A-band end-to-end smoke test",
)
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

    wavelengths_nm = np.arange(759.0, 772.1, 0.1)
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
