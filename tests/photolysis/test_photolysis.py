from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.photolysis import (
    LinePhotolysis,
    LymanAlphaPhotolysis,
    Photolysis,
    lyman_alpha_o2_rate_per_photon_m2,
    lyman_alpha_reduction_factor,
    o3_o1d_matsumi2002,
    o3_o3p_matsumi2002,
    photolysis_rates,
    presets,
    slant_columns,
)


def _flux(
    wavelength,
    actinic=1.0,
    cross_sections=None,
    temperature=(250.0, 220.0),
    slant_column=0.0,
):
    """A synthetic ActinicFlux.calculate result with constant spectra."""
    wavelength = np.asarray(wavelength, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    cross_sections = cross_sections or {"O3": 3.0e-24}
    shape = (wavelength.size, temperature.size)
    columns = np.broadcast_to(slant_column, (len(cross_sections), temperature.size))
    return xr.Dataset(
        {
            "actinic_flux": (("wavelength", "altitude"), np.full(shape, actinic)),
            "solar_flux": (("wavelength",), np.full(wavelength.size, 2.0 * actinic)),
            "cross_section": (
                ("species", "wavelength", "altitude"),
                np.stack([np.full(shape, s) for s in cross_sections.values()]),
            ),
            "slant_column": (("species", "altitude"), columns),
            "temperature_k": (("altitude",), temperature),
        },
        coords={
            "wavelength": wavelength,
            "altitude": np.arange(temperature.size) * 1.0e4,
            "species": list(cross_sections),
        },
    )


def test_o1d_quantum_yield_regimes():
    assert o3_o1d_matsumi2002(250.0, 250.0) == 0.90
    assert o3_o1d_matsumi2002(330.0, 250.0) == 0.08
    assert o3_o1d_matsumi2002(345.0, 250.0) == 0.0
    # Regression value of the Matsumi et al. (2002) parameterisation.
    np.testing.assert_allclose(o3_o1d_matsumi2002(308.0, 298.0), 0.7934, atol=1e-4)
    # The fall-off grows with temperature (vibrationally excited O3).
    yields = o3_o1d_matsumi2002(315.0, np.array([200.0, 250.0, 300.0]))
    assert np.all(np.diff(yields) > 0)


def test_quantum_yields_broadcast_and_complement():
    wavelength = np.linspace(290.0, 350.0, 7)[:, np.newaxis]
    temperature = np.array([[200.0, 300.0]])
    o1d = o3_o1d_matsumi2002(wavelength, temperature)
    assert o1d.shape == (7, 2)
    assert np.all((o1d >= 0) & (o1d <= 1))
    np.testing.assert_allclose(o1d + o3_o3p_matsumi2002(wavelength, temperature), 1.0)


def test_continuum_rate_integrates_over_range():
    flux = _flux(np.arange(300.0, 601.0, 1.0), actinic=2.0)
    rates = photolysis_rates(
        flux,
        [
            Photolysis("constant", "O3", 0.5, wavelength_range_nm=(400.0, 500.0)),
            Photolysis(
                "linear",
                "O3",
                lambda w, t: w / 1000.0 + 0.0 * t,
                wavelength_range_nm=(400.0, 500.0),
            ),
        ],
    )
    f_sigma = 2.0 * 3.0e-24
    np.testing.assert_allclose(rates["constant"], f_sigma * 0.5 * 100.0, rtol=1e-12)
    np.testing.assert_allclose(
        rates["linear"], f_sigma * (500.0**2 - 400.0**2) / 2000.0, rtol=1e-12
    )
    assert rates["constant"].attrs["units"] == "s^-1"


def test_quantum_yield_sees_each_altitude_temperature():
    flux = _flux(np.arange(300.0, 401.0, 1.0), temperature=(200.0, 300.0))
    rates = photolysis_rates(
        flux, [Photolysis("t", "O3", lambda w, t: t / 300.0 + 0.0 * w)]
    )
    np.testing.assert_allclose(rates["t"][1] / rates["t"][0], 1.5)


def test_line_rate_scales_toa_flux_by_transmission():
    flux = _flux([121.0, 121.567, 122.0])
    line = LinePhotolysis("lya", 121.567, 1.0e-24, 3.0e15, quantum_yield=0.5)
    rates = photolysis_rates(flux, [line])
    # solar_flux is twice the actinic flux in _flux, so transmission is 0.5.
    np.testing.assert_allclose(rates["lya"], 0.5 * 1.0e-24 * 3.0e15 * 0.5)


def test_line_rate_scales_with_earth_sun_distance():
    flux = _flux([121.0, 121.567, 122.0])
    line = LinePhotolysis("lya", 121.567, 1.0e-24, 3.0e15)
    at_1au = photolysis_rates(flux, [line])["lya"]
    flux.attrs["earth_sun_distance_au"] = 2.0
    np.testing.assert_allclose(photolysis_rates(flux, [line])["lya"], at_1au / 4.0)


def test_lyman_alpha_rate_at_the_top_of_the_atmosphere():
    flux = _flux([121.0, 121.567, 122.0], cross_sections={"O2": 1.0e-24})
    rates = photolysis_rates(
        flux, [LymanAlphaPhotolysis("lya", 3.0e15, quantum_yield=0.5)]
    )
    # Unattenuated: the sum of the d_i of Chabrillat and Kockarts (1997).
    np.testing.assert_allclose(rates["lya"], 0.5 * 3.0e15 * 2.30989e-24, rtol=1e-5)
    np.testing.assert_allclose(lyman_alpha_reduction_factor(0.0), 1.0, atol=1e-3)
    flux.attrs["earth_sun_distance_au"] = 2.0
    np.testing.assert_allclose(
        photolysis_rates(flux, [LymanAlphaPhotolysis("lya", 3.0e15, 0.5)])["lya"],
        rates["lya"] / 4.0,
    )


def test_lyman_alpha_effective_cross_section_falls_with_depth():
    column = np.array([0.0, 1.0e20, 1.0e22, 1.0e24, np.inf])
    reduction = lyman_alpha_reduction_factor(column)
    effective = lyman_alpha_o2_rate_per_photon_m2(column[:-1]) / reduction[:-1]
    assert np.all(np.diff(reduction) < 0)
    assert np.all(np.diff(effective) < 0)
    # The line wings go first; deep in the atmosphere the core cross section
    # of about 1e-20 cm^2 remains.
    assert 2.0e-24 < effective[0] < 2.5e-24
    assert 0.6e-24 < effective[-1] < 1.1e-24
    assert reduction[-1] == 0.0

    flux = _flux(
        [121.0, 121.567, 122.0],
        cross_sections={"O2": 1.0e-24},
        temperature=np.full(4, 200.0),
        slant_column=column[1:],
    )
    rates = photolysis_rates(flux, [LymanAlphaPhotolysis("lya", 3.0e15)])
    np.testing.assert_allclose(
        rates["lya"], 3.0e15 * lyman_alpha_o2_rate_per_photon_m2(column[1:])
    )


def test_lyman_alpha_needs_an_o2_column():
    with pytest.raises(ValueError, match="no O2 slant column"):
        photolysis_rates(
            _flux([121.0, 121.567, 122.0]), [LymanAlphaPhotolysis("lya", 3.0e15)]
        )


def test_slant_columns_overhead_and_oblique():
    altitude = np.arange(0.0, 200.1e3, 1.0e3)
    scale_height = 7.0e3
    density = 1.0e25 * np.exp(-altitude / scale_height)
    vertical = (
        1.0e25
        * scale_height
        * (np.exp(-altitude / scale_height) - np.exp(-altitude[-1] / scale_height))
    )
    np.testing.assert_allclose(
        slant_columns(altitude, density, 1.0)[0, :-20], vertical[:-20], rtol=1e-5
    )
    # At 60 degrees, curvature lowers the plane-parallel column by 0.325%
    # (from an independent quadrature).
    np.testing.assert_allclose(
        slant_columns(altitude, density, 0.5)[0, [0, 30]],
        vertical[[0, 30]] / 0.5 * (1.0 - np.array([3.2536e-3, 3.2386e-3])),
        rtol=1e-5,
    )


def test_slant_columns_horizontal_and_shadowed():
    altitude = np.arange(0.0, 300.1e3, 1.0e3)
    scale_height = 7.0e3
    density = 1.0e25 * np.exp(-altitude / scale_height)
    columns = slant_columns(altitude, density, 0.0)[0]
    # Chapman function at 90 degrees, sqrt(pi x / 2) with x = r / H.
    for i in (20, 50, 80):
        r = sk.photolysis.flux.EARTH_RADIUS_M + altitude[i]
        chapman = np.sqrt(np.pi * r / (2.0 * scale_height))
        np.testing.assert_allclose(
            columns[i], density[i] * scale_height * chapman, rtol=2e-3
        )
    # 30 degrees below the horizon: blocked below about 1 / cos(30) - 1 of R.
    below = slant_columns(altitude, density, -0.5)[0]
    shadow_km = sk.photolysis.flux.EARTH_RADIUS_M * (2.0 / np.sqrt(3.0) - 1.0) / 1e3
    assert np.all(np.isinf(below[altitude / 1e3 < shadow_km]))
    assert np.all(np.isfinite(below[altitude / 1e3 > shadow_km]))


def test_negative_flux_noise_is_clipped():
    flux = _flux([121.0, 121.567, 122.0, 123.0])
    flux["actinic_flux"][:, 1] = -1.0e-20
    rates = photolysis_rates(
        flux,
        [
            LinePhotolysis("line", 121.567, 1.0e-24, 3.0e15),
            Photolysis("continuum", "O3"),
        ],
    )
    assert float(rates["line"][1]) == 0.0
    assert float(rates["continuum"][1]) == 0.0


@pytest.mark.parametrize(
    ("reaction", "message"),
    [
        (Photolysis("x", "NO2"), "no cross section for NO2"),
        (Photolysis("x", "O3", wavelength_range_nm=(1.0, 2.0)), "fewer than two"),
        (Photolysis("x", "O3", max_grid_spacing_nm=0.1), "at most 0.1 nm"),
        (LinePhotolysis("x", 900.0, 1.0, 1.0), "outside the wavelength grid"),
    ],
)
def test_invalid_reactions(reaction, message):
    with pytest.raises(ValueError, match=message):
        photolysis_rates(_flux(np.arange(300.0, 601.0, 1.0)), [reaction])


def test_duplicate_rate_names():
    with pytest.raises(ValueError, match="twice"):
        photolysis_rates(
            _flux(np.arange(300.0, 601.0, 1.0)), [Photolysis("x", "O3")] * 2
        )


def test_oxygen_yankovsky_rates_match_the_mechanism():
    flux = _flux(
        sk.photolysis.airglow_wavelength_grid(),
        cross_sections={"O3": 1.0e-23, "O2": 1.0e-26},
    )
    channels = photolysis_rates(flux, presets.oxygen_photolysis())
    rates = presets.oxygen_yankovsky_rates(flux)

    mechanism = sk.nlte.Mechanism.bundled("oxygen_yankovsky")
    assert sorted(rates.data_vars) == sorted(mechanism.rate_inputs)
    np.testing.assert_allclose(
        sum(rates[f"J_O3_A{v}"] for v in range(6)), channels["J_O3_O1D"], rtol=1e-12
    )
    np.testing.assert_allclose(
        sum(rates[f"J_O3_X{v}"] for v in range(1, 36)),
        channels["J_O3_O3P"],
        rtol=1e-12,
    )
    # The Schumann-Runge continuum rate covers 130-175 nm only.
    np.testing.assert_allclose(channels["J_O2_SRC"], 1.0e-26 * 45.0, rtol=1e-9)


def test_actinic_flux_ozone_only_atmosphere():
    altitudes = np.arange(0.0, 100.1e3, 10.0e3)
    z_km = altitudes / 1e3
    pressure = 101325.0 * np.exp(-z_km / 7.0)
    temperature = np.full(altitudes.size, 250.0)
    atmosphere = xr.Dataset(
        {
            "temperature_k": ("altitude", temperature),
            "pressure_pa": ("altitude", pressure),
            "O3": ("altitude", 5.0e18 * np.exp(-(((z_km - 25.0) / 8.0) ** 2))),
            "O(3P)": ("altitude", np.zeros(altitudes.size)),
        },
        coords={"altitude": altitudes},
    )
    calculator = sk.photolysis.ActinicFlux(
        altitudes, wavelengths_nm=np.arange(200.0, 700.1, 1.0)
    )
    flux = calculator.calculate(atmosphere, cos_sza=0.8, albedo=0.1)

    assert list(flux["species"].values) == ["O3"]
    assert flux["actinic_flux"].shape == (501, altitudes.size)
    np.testing.assert_allclose(
        flux["slant_column"].sel(species="O3")[-2:],
        slant_columns(altitudes, atmosphere["O3"].to_numpy()[np.newaxis], 0.8)[0, -2:],
    )

    rates = photolysis_rates(
        flux,
        [
            Photolysis("J_O3_O1D", "O3", o3_o1d_matsumi2002),
            Photolysis("J_O3_O3P", "O3", o3_o3p_matsumi2002),
        ],
    )
    total = rates["J_O3_O1D"] + rates["J_O3_O3P"]
    # Top-of-atmosphere O3 photolysis is about 8e-3 s^-1, mostly O(1D).
    assert 6.0e-3 < float(total[-1]) < 1.0e-2
    assert 0.75 < float(rates["J_O3_O1D"][-1] / total[-1]) < 0.9
    # The ozone layer removes most O(1D)-producing UV below it.
    assert float(rates["J_O3_O1D"][0]) < 0.05 * float(rates["J_O3_O1D"][-1])
