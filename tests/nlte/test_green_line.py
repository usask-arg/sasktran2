from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.photolysis.quantum_yields import o2_o1s_yield


def _column(o2a=0.0):
    z = np.array([85e3, 95e3])
    temperature = np.array([190.0, 200.0])
    n = np.array([5e19, 1e19])
    return (
        z,
        temperature,
        n,
        xr.Dataset(
            {
                "temperature_k": ("altitude", temperature),
                "O2": ("altitude", 0.21 * n),
                "N2": ("altitude", 0.78 * n),
                "O(3P)": ("altitude", np.array([3e17, 5e17])),
                "M": ("altitude", n),
                "O2(a)": ("altitude", np.full(2, o2a)),
            },
            coords={"altitude": z},
        ),
    )


def test_barth_matches_mcdade():
    _, temperature, n, column = _column()
    mechanism = sk.nlte.Mechanism.bundled("oxygen_green")
    solution = sk.nlte.solve(
        mechanism, column, {"J_O2_O1S": np.zeros(2), "P_O1S_ION": np.zeros(2)}
    )
    cm = 1e-6
    o, o2 = column["O(3P)"].to_numpy() * cm, column["O2"].to_numpy() * cm
    production = (
        4.7e-33 * (300 / temperature) ** 2 * o**3 * n * cm / (15 * o2 + 211 * o)
    )
    loss = 1.26 + 7.54e-2 + 2.42e-4 + 4.0e-12 * np.exp(-865 / temperature) * o2
    np.testing.assert_allclose(
        solution["photon_ver"].sel(transition="o1s_green_line") * cm,
        1.26 * production / loss,
        rtol=1e-6,
    )


def test_o2a_quenches_o1s():
    mechanism = sk.nlte.Mechanism.bundled("oxygen_green")
    rates = {"J_O2_O1S": np.full(2, 1e-9), "P_O1S_ION": np.zeros(2)}
    without = sk.nlte.solve(mechanism, _column()[3], rates)
    with_o2a = sk.nlte.solve(mechanism, _column(o2a=1e16)[3], rates)
    ratio = (with_o2a["density"] / without["density"]).sel(state="O(1S)").to_numpy()
    # 1.7e-10 cm3 s-1 x 1e10 cm-3 = 1.7 s-1 against ~1.34 s-1 radiative.
    assert np.all((ratio > 0.3) & (ratio < 0.6))


def test_o2_o1s_yield_bins():
    np.testing.assert_allclose(
        o2_o1s_yield(np.array([80.5, 83.0, 88.0, 95.0, 105.0, 112.0, 118.0, 121.5])),
        [0.0, 0.01, 0.03, 0.07, 0.10, 0.03, 0.01, 0.0],
    )


def test_vuv_tables():
    try:
        o2 = sk.optical.VUVAbsorber("O2")
    except OSError:
        pytest.skip("VUV tables are not available")
    assert o2 is not None
    with pytest.raises(ValueError, match="available"):
        sk.optical.VUVAbsorber("CO2")


def test_add_green_line_with_given_rates():
    z = np.arange(80e3, 130_001.0, 5e3)
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
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([557.9]))
    atmosphere.temperature_k = np.full(z.size, 200.0)
    atmosphere.pressure_pa = 101325.0 * np.exp(-z / 7000.0)
    air = atmosphere.pressure_pa / (1.380649e-23 * 200.0)
    background = xr.Dataset(
        {
            "O2": ("altitude", 0.21 * air),
            "O3": ("altitude", 1e-7 * air),
            "O": ("altitude", np.full(z.size, 3e17)),
        },
        coords={"altitude": z},
    )
    oxygen = sk.nlte.Mechanism.bundled("oxygen")
    rates = xr.Dataset(
        {
            **{
                name: ("altitude", np.full(z.size, 1e-9)) for name in oxygen.rate_inputs
            },
            "J_O2_O1S": ("altitude", np.full(z.size, 1e-9)),
        },
        coords={"altitude": z},
    )
    production = np.full(z.size, 1e5)
    result = sk.nlte.add_photochemical_species(
        atmosphere,
        ["O(1S)"],
        cos_sza=0.6,
        background=background,
        rates=rates,
        ionospheric_o1s_production=production,
    )
    assert atmosphere["O(1S) emission"] is not None
    assert "O2(a)" in result["state"].values
    assert "O(1S)" in result["state"].values
    budget = sk.nlte.budget(sk.nlte.Mechanism.bundled("oxygen_green"), result, "O(1S)")
    np.testing.assert_allclose(
        budget.sel(process="ionospheric_o1s"), production, rtol=1e-6
    )
