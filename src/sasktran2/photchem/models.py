from __future__ import annotations

import numpy as np
import xarray as xr

from sasktran2 import nlte
from sasktran2._core_rust import PyYankovsky
from sasktran2.constituent import (
    MonochromaticVolumeEmissionRate,
    PopulationEmissionRate,
)
from sasktran2.photchem import actinic_flux

OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING = "einstein_a_branching"
OXYGEN_A_BAND_WEIGHT_MODEL_HITRAN_LINE_STRENGTH = "hitran_line_strength"


def _atmosphere_from_flux(flux: xr.Dataset) -> xr.Dataset:
    """Background densities for the ``oxygen_yankovsky`` mechanism.

    Reads the variables of :func:`sasktran2.photchem.actinic_flux` output.
    CO2 and O(3P) default to zero when absent, as in the earlier solver.
    """
    altitude = flux["altitude"]
    zeros = np.zeros(altitude.size)

    def density(*names: str) -> np.ndarray:
        for name in names:
            if name in flux:
                return flux[name].to_numpy()
        return zeros

    return xr.Dataset(
        {
            "temperature_k": ("altitude", flux["temperature"].to_numpy()),
            "O2": ("altitude", flux["o2_density"].to_numpy()),
            "O3": ("altitude", flux["o3_density"].to_numpy()),
            "N2": ("altitude", flux["n2_density"].to_numpy()),
            "CO2": ("altitude", density("co2_density")),
            "O(3P)": ("altitude", density("o3p_density", "o_density")),
        },
        coords={"altitude": altitude},
    )


class Yankovsky:
    """The ``oxygen_yankovsky`` mechanism driven by an actinic-flux spectrum.

    Photolysis rates come from the earlier photochem calculation (integrated
    rates scaled to fixed top-of-atmosphere values); the kinetics are solved
    with :mod:`sasktran2.nlte`.
    """

    _model: PyYankovsky

    def __init__(self):
        self._model = PyYankovsky()
        self.mechanism = nlte.Mechanism.bundled("oxygen_yankovsky")

    def run(self):
        return self.solve(actinic_flux())

    def photolysis_rates(self, flux) -> dict[str, np.ndarray]:
        """The mechanism's rate inputs [s^-1] derived from an actinic-flux dataset."""
        return self._model.photolysis_rates(flux)

    def solve_full(self, flux) -> xr.Dataset:
        """The :func:`sasktran2.nlte.solve` result, including budgets."""
        return nlte.solve(
            self.mechanism, _atmosphere_from_flux(flux), self.photolysis_rates(flux)
        )

    def solve(self, flux):
        """State populations [m^-3], one variable per state."""
        solution = self.solve_full(flux)
        return xr.Dataset(
            {
                state: ("altitude", solution["density"].sel(state=state).to_numpy())
                for state in sorted(self.mechanism.states)
            },
            coords={"altitude": flux["altitude"]},
        )

    def emissions(self, state=None, flux=None):
        if state is None:
            state = self.solve(flux if flux is not None else actinic_flux())

        return self._model.emission_rates(state)

    def oxygen_green_line_mcdade(self, atmosphere_state=None, flux=None):
        if atmosphere_state is None:
            atmosphere_state = flux if flux is not None else actinic_flux()

        return self._model.oxygen_green_line_mcdade(atmosphere_state)

    def oxygen_green_line_mcdade_constituent(
        self, green_line=None, atmosphere_state=None, flux=None
    ):
        if green_line is None:
            green_line = self.oxygen_green_line_mcdade(
                atmosphere_state=atmosphere_state,
                flux=flux,
            )

        return self.oxygen_green_line_constituent(emissions=green_line)

    def oxygen_green_line_constituent(self, emissions=None, state=None, flux=None):
        if emissions is None:
            emissions = self.emissions(state=state, flux=flux)

        return MonochromaticVolumeEmissionRate(
            emissions["altitude"].to_numpy(),
            emissions["oxygen_green_5577_photon_ver"].to_numpy(),
            emissions.attrs["oxygen_green_wavelength_nm"],
        )

    def oxygen_a_band_constituent(
        self,
        emissions=None,
        state=None,
        flux=None,
        line_data=None,
        temperature_k=None,
        line_weight_model=OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING,
    ):
        _ = (line_data, temperature_k)

        if state is None:
            if emissions is not None:
                msg = (
                    "oxygen_a_band_constituent now requires population state. "
                    "Use PopulationEmissionRate directly for O2 band emission."
                )
                raise ValueError(msg)
            state = self.solve(flux if flux is not None else actinic_flux())

        return PopulationEmissionRate(
            state,
            species=["O2"],
            line_weight_model=line_weight_model,
        )

    def emission_constituents(
        self,
        emissions=None,
        state=None,
        flux=None,
        line_data=None,
        temperature_k=None,
        line_weight_model=OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING,
        include_oxygen_green_line=True,
        include_oxygen_a_band=True,
    ):
        _ = (line_data, temperature_k)

        if state is None and include_oxygen_a_band:
            state = self.solve(flux if flux is not None else actinic_flux())
        if emissions is None and include_oxygen_green_line:
            emissions = self.emissions(state=state, flux=flux)

        constituents = {}
        if include_oxygen_green_line:
            constituents["oxygen_green"] = self.oxygen_green_line_constituent(
                emissions=emissions
            )

        if include_oxygen_a_band:
            constituents["oxygen_a_band"] = self.oxygen_a_band_constituent(
                state=state,
                line_weight_model=line_weight_model,
            )

        return constituents

    def add_emissions_to_atmosphere(
        self,
        atmosphere,
        emissions=None,
        state=None,
        flux=None,
        line_data=None,
        temperature_k=None,
        line_weight_model=OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING,
        include_oxygen_green_line=True,
        include_oxygen_a_band=True,
    ):
        _ = (line_data, temperature_k)

        constituents = self.emission_constituents(
            emissions=emissions,
            state=state,
            flux=flux,
            line_weight_model=line_weight_model,
            include_oxygen_green_line=include_oxygen_green_line,
            include_oxygen_a_band=include_oxygen_a_band,
        )

        for name, constituent in constituents.items():
            atmosphere[name] = constituent

        return constituents
