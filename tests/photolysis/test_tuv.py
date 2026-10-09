from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.photolysis import (
    LinePhotolysis,
    Photolysis,
    TUVActinicFlux,
    TUVXQuantumYield,
    chabrillat_kockarts_cross_section,
    koppers_murtagh_cross_section,
    photolysis_rates,
    presets,
)


@pytest.fixture(scope="module")
def tables():
    try:
        return sk.photolysis.tuvx_v54_tables()
    except OSError:
        pytest.skip("TUV-x v5.4 tables are not available")


def _binned_flux(edges, actinic=1.0, cross_section=2.0e-24, temperature=(250.0,)):
    edges = np.asarray(edges, dtype=float)
    centres = 0.5 * (edges[1:] + edges[:-1])
    temperature = np.asarray(temperature, dtype=float)
    shape = (centres.size, temperature.size)
    return xr.Dataset(
        {
            "actinic_flux": (("wavelength", "altitude"), np.full(shape, actinic)),
            "solar_flux": (("wavelength",), np.full(centres.size, 2.0 * actinic)),
            "cross_section": (
                ("species", "wavelength", "altitude"),
                np.full((1, *shape), cross_section),
            ),
            "temperature_k": (("altitude",), temperature),
        },
        coords={
            "wavelength": centres,
            "wavelength_edge": ("wavelength_edge", edges),
            "altitude": np.arange(temperature.size) * 1.0e4,
            "species": ["O3"],
        },
    )


def test_binned_rates_sum_over_bins():
    flux = _binned_flux([300.0, 301.0, 303.0, 306.0, 310.0])
    rates = photolysis_rates(
        flux,
        [
            Photolysis("all", "O3", 0.5),
            # Bins whose centres are in the range: 302 and 304.5 nm.
            Photolysis("some", "O3", wavelength_range_nm=(301.5, 305.0)),
            Photolysis("linear", "O3", lambda w, t: w / 300.0 + 0.0 * t),
        ],
    )
    np.testing.assert_allclose(rates["all"], 0.5 * 2.0e-24 * 10.0)
    np.testing.assert_allclose(rates["some"], 2.0e-24 * 5.0)
    expected = 2.0e-24 * np.sum(
        np.array([300.5, 302.0, 304.5, 308.0]) / 300.0 * [1, 2, 3, 4]
    )
    np.testing.assert_allclose(rates["linear"], expected)


def test_binned_rates_check_bin_width_and_range():
    flux = _binned_flux([300.0, 301.0, 303.0])
    with pytest.raises(ValueError, match="at most 1.5 nm"):
        photolysis_rates(flux, [Photolysis("x", "O3", max_grid_spacing_nm=1.5)])
    with pytest.raises(ValueError, match="no bin centres"):
        photolysis_rates(
            flux, [Photolysis("x", "O3", wavelength_range_nm=(310.0, 320.0))]
        )


def test_binned_line_rate_uses_the_bin_containing_the_line():
    flux = _binned_flux([121.0, 121.4, 121.9, 122.3])
    flux["actinic_flux"][1] = 0.5  # transmission 0.25 in the line's bin
    rates = photolysis_rates(flux, [LinePhotolysis("lya", 121.567, 1.0e-24, 3.0e15)])
    np.testing.assert_allclose(rates["lya"], 1.0e-24 * 3.0e15 * 0.25)


def test_oxygen_photolysis_without_excitation():
    names = [r.name for r in presets.oxygen_photolysis(excitation=False)]
    assert names == [
        "J_O3_O1D",
        "J_O3_O3P",
        *(f"J_O3_O1D_A{v}" for v in range(6)),
        "J_O3_O1D_X",
        "J_O2_SRC",
        "J_O2_LYA",
    ]
    assert len(presets.oxygen_photolysis()) == 15


def test_chabrillat_kockarts_cross_section():
    sigma = chabrillat_kockarts_cross_section(np.array([0.0, 1.0e24, 1.0e30]))
    d = np.array([6.0073e-21, 4.28569e-21, 1.28059e-20])
    b = np.array([6.8431e-01, 2.29841e-01, 8.65412e-02])
    np.testing.assert_allclose(sigma[0], d.sum() / b.sum() * 1.0e-4)
    assert 0.5e-24 < sigma[1] < 1.1e-24
    # TUV's floor once the line is extinguished.
    assert sigma[2] == 1.0e-24


def _tuv_chebyshev(coefficients, x, lower, upper):
    """TUV-x's chebyshev_evaluation, line by line."""
    y = (2.0 * x - (lower + upper)) / (upper - lower)
    di = di1 = 0.0
    for c in coefficients[:0:-1]:
        di, di1 = 2.0 * y * di - di1 + c, di
    return y * di - di1 + 0.5 * coefficients[0]


def test_koppers_murtagh_matches_tuv_evaluation(tables):
    lower, upper = tables.attrs["srb_log_column_limits"]
    t0 = tables.attrs["srb_reference_temperature_k"]
    a = tables["srb_chebyshev_a"].to_numpy()
    b = tables["srb_chebyshev_b"].to_numpy()
    columns_cm2 = np.array([1.0e17, 1.0e19, 1.0e21, 1.0e23])
    temperature = np.array([190.0, 210.0, 240.0, 280.0])
    sigma = koppers_murtagh_cross_section(columns_cm2 * 1.0e4, temperature, tables)
    for k, (n, t) in enumerate(zip(columns_cm2, temperature, strict=True)):
        x = np.log(n)
        expected = [
            np.exp(
                _tuv_chebyshev(a[:, j], x, lower, upper) * (t - t0)
                + _tuv_chebyshev(b[:, j], x, lower, upper)
            )
            * 1.0e-4
            for j in range(a.shape[1])
        ]
        np.testing.assert_allclose(sigma[k], expected, rtol=1e-12)


def test_koppers_murtagh_outside_the_fitted_columns(tables):
    lower, upper = tables.attrs["srb_log_column_limits"]
    # Upward: two levels beyond the fit, two within, two above it.
    columns_m2 = 1.0e4 * np.exp(
        [upper + 2.0, upper + 1.0, upper - 1.0, lower + 1.0, lower - 1.0, lower - 2.0]
    )
    sigma = koppers_murtagh_cross_section(columns_m2, 220.0, tables)
    np.testing.assert_allclose(sigma[0], sigma[2])
    np.testing.assert_allclose(sigma[1], sigma[2])
    np.testing.assert_allclose(sigma[4], sigma[3])
    np.testing.assert_allclose(sigma[5], sigma[3])


def test_tuvx_quantum_yield_interpolates_in_temperature(tables):
    table = tables["o3_o1d_quantum_yield"]
    centre = float(tables["wavelength"].sel(wavelength=310.0, method="nearest"))
    j = int(np.argmin(abs(tables["wavelength"].to_numpy() - centre)))
    phi = TUVXQuantumYield("o3_o1d")(
        np.array([[centre]]), np.array([[230.0, 230.5, 400.0]])
    )
    np.testing.assert_allclose(phi[0, 0], table.sel(temperature_k=230.0)[j])
    np.testing.assert_allclose(
        phi[0, 1],
        0.5 * (table.sel(temperature_k=230.0)[j] + table.sel(temperature_k=231.0)[j]),
    )
    np.testing.assert_allclose(phi[0, 2], table.sel(temperature_k=300.0)[j])


def _atmosphere(altitudes, o2=True, pressure_scale=1.0):
    z_km = altitudes / 1e3
    temperature = 250.0 - 30.0 * np.exp(-(((z_km - 85.0) / 15.0) ** 2))
    pressure = pressure_scale * 101325.0 * np.exp(-z_km / 7.0)
    air = pressure / (1.380649e-23 * temperature)
    variables = {
        "temperature_k": ("altitude", temperature),
        "pressure_pa": ("altitude", pressure),
        "O3": ("altitude", 5.0e18 * np.exp(-(((z_km - 25.0) / 8.0) ** 2))),
    }
    if o2:
        variables["O2"] = ("altitude", 0.2095 * air)
    return xr.Dataset(variables, coords={"altitude": altitudes})


@pytest.mark.parametrize("data", ["sasktran2", "tuv-x"])
@pytest.mark.usefixtures("tables")
def test_overhead_direct_beam_in_homogeneous_layers(data):
    # Almost no air, so no scattering: the flux is the direct beam, attenuated
    # by O3 in homogeneous layers.
    altitudes = np.arange(0.0, 80.1e3, 2.0e3)
    atmosphere = _atmosphere(altitudes, o2=False, pressure_scale=1e-12)
    flux = TUVActinicFlux(altitudes, data=data, num_streams=2).calculate(
        atmosphere, cos_sza=1.0
    )
    assert flux.attrs["spectral_sampling"] == "bins"
    o3 = atmosphere["O3"].to_numpy()
    lo, hi = o3[:-1], o3[1:]
    ratio = np.log(lo / hi)
    safe = np.where(np.abs(ratio) > 1e-8, ratio, 1.0)
    columns = np.diff(altitudes) * np.where(
        np.abs(ratio) > 1e-8, (lo - hi) / safe, 0.5 * (lo + hi)
    )
    sigma_layer = flux["cross_section"].sel(species="O3").to_numpy()
    # O3 cross sections barely depend on temperature here; use the layer's
    # lower level (TUV evaluates them at the layer midpoint temperature).
    tau = np.cumsum((sigma_layer[:, :-1] * columns)[:, ::-1], axis=1)[:, ::-1]
    expected = flux["solar_flux"].to_numpy()[:, np.newaxis] * np.exp(-tau)
    window = ((flux["wavelength"] > 280.0) & (flux["wavelength"] < 340.0)).to_numpy()
    thin = window[:, np.newaxis] & (tau < 5.0)
    np.testing.assert_allclose(
        flux["actinic_flux"].to_numpy()[:, :-1][thin], expected[thin], rtol=0.01
    )


@pytest.mark.parametrize("data", ["sasktran2", "tuv-x"])
@pytest.mark.usefixtures("tables")
def test_tuv_mode_oxygen_rates(data):
    altitudes = np.arange(0.0, 120.1e3, 2.0e3)
    atmosphere = _atmosphere(altitudes)
    flux = TUVActinicFlux(altitudes, data=data).calculate(atmosphere, cos_sza=0.6)
    assert flux["wavelength"].size == 156
    rates = photolysis_rates(flux, presets.tuvx_v54_photolysis())
    for name in rates.data_vars:
        assert np.all(np.isfinite(rates[name]))
        assert np.all(rates[name] >= 0.0)
    # O2 photolysis grows with altitude; at the top it is a few 1e-6 s^-1.
    j_o2 = rates["J_O2"].to_numpy()
    assert np.all(np.diff(j_o2[altitudes > 40e3]) > 0)
    assert 1.0e-7 < j_o2[-1] < 1.0e-5
    # O3 -> O(1D) near the top is about 8e-3 s^-1 (TUV-x sun) or less (HSRS).
    assert 5.0e-3 < float(rates["J_O3_O1D"][-1]) < 1.0e-2


def test_tuv_mode_rejects_unknown_data():
    with pytest.raises(ValueError, match="data must be"):
        TUVActinicFlux([0.0, 1.0], data="tuv")
