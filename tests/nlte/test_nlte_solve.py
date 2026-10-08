from __future__ import annotations

import numpy as np
import pytest
import xarray as xr
from sasktran2 import nlte

CHAIN = """
[mechanism]
name = "chain"
version = "1"
background = ["O2", "O(3P)"]

[references]
ref = "test"

[[state]]
id = "O(1D)"

[[state]]
id = "O2(b, v=0)"

[[photo_reaction]]
id = "src"
reactant = "O2"
products = ["O(1D)", "O(3P)"]
rate_input = "J"
reference = "ref"

[[reaction]]
id = "o1d_o2"
reactants = ["O(1D)", "O2"]
rate = { law = "constant", value = 1.0e-6, units = "cm3 s-1" }
channels = [{ yield = 0.25, products = ["O2(b)", "O(3P)"] }]
reference = "ref"

[[transition]]
id = "a_band"
upper = "O2(b)"
lower = "O2"
einstein_a_s = 2.0
wavelength_nm = 762.0
reference = "ref"
"""


def _atmosphere():
    altitude = np.array([80.0e3, 90.0e3, 100.0e3])
    return xr.Dataset(
        {
            "temperature_k": ("altitude", np.array([200.0, 190.0, 210.0])),
            "O2": ("altitude", np.array([1.0e20, 1.0e19, 1.0e18])),
            "O(3P)": ("altitude", np.zeros(3)),
            "unused": ("altitude", np.ones(3)),
        },
        coords={"altitude": altitude},
    )


def test_chain_matches_closed_form():
    mechanism = nlte.Mechanism.from_toml(CHAIN)
    atmosphere = _atmosphere()
    j = np.array([1.0e-6, 2.0e-6, 3.0e-6])

    solution = nlte.solve(mechanism, atmosphere, {"J": j})

    o2 = atmosphere["O2"].to_numpy()
    k = 1.0e-12
    o1d = j * o2 / (k * o2)
    o2b = 0.25 * k * o2 * o1d / 2.0
    np.testing.assert_allclose(solution["density"].sel(state="O(1D)"), o1d, rtol=1e-12)
    np.testing.assert_allclose(solution["density"].sel(state="O2(b)"), o2b, rtol=1e-12)
    np.testing.assert_allclose(
        solution["photon_ver"].sel(transition="a_band"), 2.0 * o2b, rtol=1e-12
    )
    assert solution["transition_wavelength_nm"].sel(transition="a_band") == 762.0
    np.testing.assert_array_equal(solution["altitude"], atmosphere["altitude"])
    assert (solution["relative_residual"] < 1e-12).all()


def test_budget_contributions_balance():
    mechanism = nlte.Mechanism.from_toml(CHAIN)
    solution = nlte.solve(mechanism, _atmosphere(), {"J": np.full(3, 1.0e-6)})

    o2b = nlte.budget(mechanism, solution, "O2(b)")

    assert list(o2b["process"].values) == ["o1d_o2", "a_band"]
    np.testing.assert_allclose(o2b.sum("process"), 0.0, atol=1e-9 * o2b.max())
    with pytest.raises(ValueError, match="not a state"):
        nlte.budget(mechanism, solution, "O2")


def test_missing_inputs_raise_value_error():
    mechanism = nlte.Mechanism.from_toml(CHAIN)
    atmosphere = _atmosphere().drop_vars("O2")

    with pytest.raises(ValueError, match="O2"):
        nlte.solve(mechanism, atmosphere, {"J": np.ones(3)})
    with pytest.raises(ValueError, match="J"):
        nlte.solve(mechanism, _atmosphere(), {})


def test_invalid_mechanism_raises_value_error():
    with pytest.raises(ValueError, match="'O3'"):
        nlte.Mechanism.from_toml(CHAIN.replace('["O(1D)", "O2"]', '["O(1D)", "O3"]'))


def test_bundled_oxygen_mechanism():
    assert "oxygen_yankovsky" in nlte.Mechanism.bundled_names()
    mechanism = nlte.Mechanism.bundled("oxygen_yankovsky")

    assert len(mechanism.states) == 46
    assert mechanism.background == ["O2", "O3", "N2", "CO2", "O(3P)"]
    assert "J_O3_A0" in mechanism.rate_inputs
    assert "legacy_photchem" in mechanism.references
    assert "46 states" in repr(mechanism)
