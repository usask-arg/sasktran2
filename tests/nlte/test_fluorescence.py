from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.database.hitran_line import HITRANLineDatabase
from sasktran2.nlte import fluorescence


def _has_oh_lines():
    return (HITRANLineDatabase()._db_root / "OH.data").exists()


pytestmark = pytest.mark.skipif(not _has_oh_lines(), reason="HITRAN OH not cached")


@pytest.fixture(scope="module")
def lines():
    return fluorescence.oh_ax_lines()


def test_line_list(lines):
    assert lines.wavenumber_cminv.size > 9000
    v00 = (lines.v_upper == 0) & (lines.v_lower == 0)
    assert 300 < lines.wavelength_nm[v00].min() < lines.wavelength_nm[v00].max() < 360
    # Radiative rate of low-J v'=0 levels: about 1.47e6 s^-1 (lifetime 0.68 us).
    total = np.bincount(lines.upper_level, weights=lines.einstein_a_s)
    low_j = np.unique(lines.upper_level[(lines.v_upper == 0) & (lines.g_upper <= 12)])
    np.testing.assert_allclose(total[low_j], 1.47e6, rtol=0.03)


def test_photons_absorbed_equal_photons_emitted(lines):
    flux = np.full((2, lines.wavenumber_cminv.size), 1.0e18)
    ver, weights = fluorescence.oh_ax_fluorescence(
        lines, [180.0, 250.0], [1.0e12, 2.0e12], flux
    )
    np.testing.assert_allclose(weights.sum(axis=1), 1.0)
    # Linear in OH density; the rotational distribution changes little.
    assert 1.8 < ver[1] / ver[0] < 2.2


def test_top_of_atmosphere_g_factor(lines):
    wavelength, photons = fluorescence._hsrs_photons()
    solar = np.interp(lines.wavelength_nm, wavelength, photons)[np.newaxis, :]
    g, weights = fluorescence.oh_ax_fluorescence(lines, [200.0], [1.0], solar)
    # Literature (0,0) g-factors are a few 1e-4 s^-1 per molecule.
    assert 5e-4 < g[0] < 9e-4
    v00 = (lines.v_upper == 0) & (lines.v_lower == 0)
    assert 0.8 < weights[0][v00].sum() < 0.95


def test_add_oh_fluorescence():
    z = np.arange(0.0, 100_001.0, 5_000.0)
    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_372_000.0,
        z,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([308.0]))
    atmosphere.temperature_k = np.full(z.size, 220.0)
    atmosphere.pressure_pa = 101325.0 * np.exp(-z / 7000.0)
    oh = 1.0e13 * np.exp(-(((z - 60e3) / 10e3) ** 2))
    background = xr.Dataset({"OH": ("altitude", oh)}, coords={"altitude": z})
    result = sk.nlte.add_photochemical_species(
        atmosphere,
        ["OH(A)"],
        cos_sza=0.6,
        background=background,
        actinic_flux=sk.photolysis.ActinicFlux(
            z, wavelengths_nm=np.arange(270.0, 360.0, 0.5)
        ),
    )
    assert atmosphere["OH(A) emission"] is not None
    ver = result["OH(A) photon_ver"].to_numpy()
    # No absorbers: the per-molecule rate is the top-of-atmosphere g (7e-4)
    # plus Rayleigh-scattered and reflected light.
    rate = ver[z == 60e3] / oh[z == 60e3]
    assert 7e-4 < rate[0] < 2.5 * 7e-4


def test_zero_oh_levels(lines):
    flux = np.full((3, lines.wavenumber_cminv.size), 1.0e18)
    ver, weights = fluorescence.oh_ax_fluorescence(
        lines, [200.0, 200.0, 200.0], [1.0e12, 0.0, 1.0e12], flux
    )
    assert ver[1] == 0.0
    # Every row must still sum to one for LineListVolumeEmissionRate.
    np.testing.assert_allclose(weights.sum(axis=1), 1.0)
    sk.constituent.LineListVolumeEmissionRate(
        np.array([0.0, 1.0, 2.0]), ver, lines.wavelength_nm, weights
    )
