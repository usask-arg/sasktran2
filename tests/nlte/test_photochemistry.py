from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.database.hitran_line import HITRANLineDatabase

K_BOLTZMANN = 1.380649e-23


def _has_local_o2_hitran_cache():
    db = HITRANLineDatabase()
    return (db._db_root / "O2.data").exists() and (db._db_root / "O2.header").exists()


needs_o2_lines = pytest.mark.skipif(
    not _has_local_o2_hitran_cache(), reason="HITRAN O2 lines are not available"
)

ALTITUDES = np.arange(0.0, 120_001.0, 2_000.0)
COS_SZA = 0.6


def _state():
    z_km = ALTITUDES / 1e3
    temperature = 240.0 - 50.0 * np.exp(-(((z_km - 88.0) / 12.0) ** 2))
    pressure = 101325.0 * np.exp(-z_km / 7.0)
    return temperature, pressure, pressure / (K_BOLTZMANN * temperature)


def _o3_vmr():
    z_km = ALTITUDES / 1e3
    return 8e-6 * np.exp(-(((z_km - 32.0) / 10.0) ** 2)) + 1e-6 * np.exp(
        -(((z_km - 90.0) / 5.0) ** 2)
    )


def _atmosphere(emission_source=sk.EmissionSource.VolumeEmissionRate, absorbers=True):
    config = sk.Config()
    config.emission_source = emission_source
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    geometry = sk.Geometry1D(
        COS_SZA,
        0.0,
        6_372_000.0,
        ALTITUDES,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.arange(758.0, 772.0, 0.0005)
    )
    temperature, pressure, _ = _state()
    atmosphere.temperature_k = temperature
    atmosphere.pressure_pa = pressure
    if absorbers:
        atmosphere["O3"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), ALTITUDES, _o3_vmr()
        )
        atmosphere["O2"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O2UV(), ALTITUDES, np.full(ALTITUDES.size, 0.2095)
        )
    return atmosphere, geometry, config


def _background(absorbers=False):
    z_km = ALTITUDES / 1e3
    atomic_oxygen = 4e17 * np.exp(-(((z_km - 95.0) / 8.0) ** 2)) + 1e10
    variables = {"O": ("altitude", atomic_oxygen)}
    if absorbers:
        air = _state()[2]
        variables["O2"] = ("altitude", 0.2095 * air)
        variables["O3"] = ("altitude", _o3_vmr() * air)
    return xr.Dataset(variables, coords={"altitude": ALTITUDES})


def _rates():
    mechanism = sk.nlte.Mechanism.bundled("oxygen")
    z_km = ALTITUDES / 1e3
    shape = 1.0 / (1.0 + np.exp(-(z_km - 50.0) / 5.0))
    return xr.Dataset(
        {name: ("altitude", 1e-5 * shape) for name in mechanism.rate_inputs},
        coords={"altitude": ALTITUDES},
    )


@needs_o2_lines
def test_populations_match_the_kinetics():
    atmosphere, _, _ = _atmosphere()
    solution = sk.nlte.add_photochemical_species(
        atmosphere, ["o2(b)"], cos_sza=COS_SZA, background=_background(), rates=_rates()
    )
    assert atmosphere["O2(b) emission"] is not None

    # The same solve, with the background assembled by hand.
    temperature, pressure, air = _state()
    mechanism = sk.nlte.Mechanism.bundled("oxygen")
    chemistry = xr.Dataset(
        {
            "temperature_k": ("altitude", temperature),
            "pressure_pa": ("altitude", pressure),
            "O2": ("altitude", 0.2095 * air),
            "O3": ("altitude", atmosphere["O3"].vmr * air),
            "N2": ("altitude", 0.7808 * air),
            "CO2": ("altitude", 4.2e-4 * air),
            "O(3P)": ("altitude", _background()["O"].to_numpy()),
        },
        coords={"altitude": ALTITUDES},
    )
    expected = sk.nlte.solve(mechanism, chemistry, _rates())
    # Log-space interpolation of the inputs changes them at round-off, which
    # the high vibrational levels amplify; their densities are near zero.
    np.testing.assert_allclose(
        solution["density"].to_numpy(),
        expected["density"].to_numpy(),
        rtol=1e-6,
        atol=1.0,
    )
    assert "J_O3_O1D" in solution


@needs_o2_lines
def test_optically_thin_limb_radiance_is_the_ver_integral():
    # Absorbers only in the background, so the line of sight is transparent.
    atmosphere, geometry, config = _atmosphere(absorbers=False)
    solution = sk.nlte.add_photochemical_species(
        atmosphere,
        ["O2(b)"],
        cos_sza=COS_SZA,
        background=_background(absorbers=True),
        rates=_rates(),
    )
    viewing = sk.ViewingGeometry()
    tangent = 85_000.0
    viewing.add_ray(sk.TangentAltitudeSolar(tangent, 0.0, 600_000.0, COS_SZA))
    radiance = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
    band = float(radiance["radiance"].sum()) * 0.0005

    # The A band and its 1-1 hot band (O2(b, v=1) -> O2(X, v=1)).
    ver = (
        solution["photon_ver"]
        .sel(transition=["o2b0_a_band", "o2b1_x1"])
        .sum("transition")
        .to_numpy()
    )
    r_t = 6_372_000.0 + tangent
    s = np.linspace(-1_300e3, 1_300e3, 200_001)
    altitude = np.sqrt(r_t**2 + s**2) - 6_372_000.0
    expected = np.trapezoid(np.interp(altitude, ALTITUDES, ver, right=0.0), s)
    np.testing.assert_allclose(band, expected / (4.0 * np.pi), rtol=1e-3)


def test_missing_atomic_oxygen():
    atmosphere, _, _ = _atmosphere()
    with pytest.raises(ValueError, match=r"O\(3P\)"):
        sk.nlte.add_photochemical_species(
            atmosphere, ["O2(b)"], cos_sza=COS_SZA, rates=_rates()
        )


def test_unknown_and_pending_species():
    atmosphere, _, _ = _atmosphere()
    with pytest.raises(ValueError, match="available: O2\\(b\\)"):
        sk.nlte.add_photochemical_species(atmosphere, ["OH(v=8)"], cos_sza=COS_SZA)
    with pytest.raises(NotImplementedError, match="1.27 um"):
        sk.nlte.add_photochemical_species(atmosphere, ["O2(a)"], cos_sza=COS_SZA)


@needs_o2_lines
def test_warns_without_the_ver_emission_source():
    atmosphere, _, _ = _atmosphere(sk.EmissionSource.Standard)
    with pytest.warns(UserWarning, match="VolumeEmissionRate"):
        sk.nlte.add_photochemical_species(
            atmosphere,
            ["O2(b)"],
            cos_sza=COS_SZA,
            background=_background(),
            rates=_rates(),
        )


@needs_o2_lines
def test_rates_from_the_actinic_flux():
    atmosphere, _, _ = _atmosphere()
    coarse = np.arange(0.0, 120_001.0, 10_000.0)
    solution = sk.nlte.add_photochemical_species(
        atmosphere,
        ["O2(b)"],
        cos_sza=COS_SZA,
        background=_background(),
        actinic_flux=sk.photolysis.ActinicFlux(coarse),
    )
    np.testing.assert_allclose(solution["altitude"], coarse)
    o2b = solution["density"].sel(state="O2(b)")
    # Daytime O2(b): about 1e12-1e14 m^-3 through the mesosphere.
    assert np.all((o2b.sel(altitude=[50e3, 70e3, 90e3]) > 1e11).to_numpy())
    assert np.all((o2b.sel(altitude=[50e3, 70e3, 90e3]) < 1e15).to_numpy())
    assert float(solution["J_O2_EXC_B0"].sel(altitude=120e3)) > 0.0
