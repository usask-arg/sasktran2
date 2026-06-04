from __future__ import annotations

import numpy as np

from sasktran2._core_rust import PyYankovsky
from sasktran2.constituent import (
    LineListVolumeEmissionRate,
    MonochromaticVolumeEmissionRate,
)
from sasktran2.database.hitran_line import HITRANLineDatabase
from sasktran2.photchem import actinic_flux

C2_K_CM = 1.4387769
OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING = "einstein_a_branching"
OXYGEN_A_BAND_WEIGHT_MODEL_HITRAN_LINE_STRENGTH = "hitran_line_strength"


class Yankovsky:
    _model: PyYankovsky

    def __init__(self):
        self._model = PyYankovsky()

    def run(self):
        flux = actinic_flux()

        return self._model.solve(flux)

    def solve(self, flux):
        return self._model.solve(flux)

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

    def _temperature_on_emission_grid(self, emissions, temperature_k):
        if temperature_k is None:
            return None

        emission_altitude = emissions["altitude"].to_numpy()
        temperature_k = np.atleast_1d(temperature_k).astype(np.float64)

        if len(temperature_k) == 1:
            return np.full_like(emission_altitude, temperature_k.item())

        if len(temperature_k) != len(emission_altitude):
            msg = (
                "temperature_k must be scalar or have the same length as the "
                "emission altitude grid"
            )
            raise ValueError(msg)

        return temperature_k

    def oxygen_green_line_constituent(self, emissions=None, state=None, flux=None):
        if emissions is None:
            emissions = self.emissions(state=state, flux=flux)

        return MonochromaticVolumeEmissionRate(
            emissions["altitude"].to_numpy(),
            emissions["oxygen_green_5577_photon_ver"].to_numpy(),
            emissions.attrs["oxygen_green_wavelength_nm"],
        )

    def oxygen_a_band_line_data(self, db: HITRANLineDatabase | None = None):
        db = db or HITRANLineDatabase()
        db_path = db.path("O2")

        return self._model.oxygen_a_band_line_data(db_path.as_posix())

    def oxygen_a_band_lte_line_weights(self, line_data, temperature_k):
        temperature_k = np.atleast_1d(temperature_k).astype(np.float64)
        upper_state_id = line_data["upper_state_id"].to_numpy()
        upper_vibrational_state = line_data["upper_vibrational_state"].to_numpy()
        upper_energy = line_data["upper_energy_cminv"].to_numpy()
        upper_g = line_data["upper_statistical_weight"].to_numpy()
        branching = line_data["upper_branching_ratio"].to_numpy()
        isotope_abundance = line_data["isotope_abundance"].to_numpy()

        weights = np.zeros((len(temperature_k), line_data.sizes["line"]))

        for vibrational_state in np.unique(upper_vibrational_state):
            band_mask = upper_vibrational_state == vibrational_state
            band_indices = np.flatnonzero(band_mask)
            unique_upper, inverse = np.unique(
                upper_state_id[band_mask],
                return_inverse=True,
            )
            upper_energy_by_level = np.zeros(len(unique_upper))
            upper_g_by_level = np.zeros(len(unique_upper))
            isotope_abundance_by_level = np.zeros(len(unique_upper))

            for i in range(len(unique_upper)):
                idx = band_indices[np.flatnonzero(inverse == i)[0]]
                upper_energy_by_level[i] = upper_energy[idx]
                upper_g_by_level[i] = upper_g[idx]
                isotope_abundance_by_level[i] = isotope_abundance[idx]

            upper_population = (
                isotope_abundance_by_level[None, :]
                * upper_g_by_level[None, :]
                * np.exp(
                    -C2_K_CM * upper_energy_by_level[None, :] / temperature_k[:, None]
                )
            )
            band_weights = upper_population[:, inverse] * branching[None, band_mask]
            band_weights /= band_weights.sum(axis=1, keepdims=True)
            weights[:, band_mask] = band_weights

        return weights

    def oxygen_a_band_lte_line_strength_weights(self, line_data, temperature_k):
        temperature_k = np.atleast_1d(temperature_k).astype(np.float64)
        line_intensity_296 = line_data["line_intensity_296"].to_numpy()
        lower_energy = line_data["lower_energy_cminv"].to_numpy()
        wavenumber = line_data["wavenumber_cminv"].to_numpy()
        upper_vibrational_state = line_data["upper_vibrational_state"].to_numpy()

        temperature = temperature_k[:, None]
        log_line_strength = (
            np.log(line_intensity_296[None, :])
            + np.log(296.0 / temperature)
            + C2_K_CM
            * lower_energy[None, :]
            * (temperature - 296.0)
            / (temperature * 296.0)
            + 2.0 * np.log(wavenumber[None, :])
            - C2_K_CM * wavenumber[None, :] / temperature
        )

        weights = np.zeros_like(log_line_strength)
        for vibrational_state in np.unique(upper_vibrational_state):
            band_mask = upper_vibrational_state == vibrational_state
            band_log_strength = log_line_strength[:, band_mask]
            band_log_strength -= band_log_strength.max(axis=1, keepdims=True)
            band_strength = np.exp(band_log_strength)
            band_sum = band_strength.sum(axis=1, keepdims=True)
            band_weights = np.divide(
                band_strength,
                band_sum,
                out=np.zeros_like(band_strength),
                where=band_sum > 0.0,
            )
            weights[:, band_mask] = band_weights

        return weights

    def _oxygen_a_band_line_weights(self, line_data, temperature_k, line_weight_model):
        if temperature_k is None:
            if line_weight_model != OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING:
                msg = f"{line_weight_model} requires temperature_k"
                raise ValueError(msg)
            return line_data["weight"].to_numpy()

        if line_weight_model == OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING:
            return self.oxygen_a_band_lte_line_weights(line_data, temperature_k)
        if line_weight_model == OXYGEN_A_BAND_WEIGHT_MODEL_HITRAN_LINE_STRENGTH:
            return self.oxygen_a_band_lte_line_strength_weights(line_data, temperature_k)

        msg = (
            "line_weight_model must be one of "
            f"{OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING!r}, "
            f"{OXYGEN_A_BAND_WEIGHT_MODEL_HITRAN_LINE_STRENGTH!r}"
        )
        raise ValueError(msg)

    def _oxygen_a_band_vibrational_photon_ver(self, emissions, line_data):
        altitude = emissions["altitude"].to_numpy()
        upper_vibrational_state = line_data["upper_vibrational_state"].to_numpy()
        photon_ver_by_state = {}

        state_to_var = {
            "O2(b)": "oxygen_a_band_b0_photon_ver",
            "O2(b, v=1)": "oxygen_a_band_b1_photon_ver",
        }

        for state_name in np.unique(upper_vibrational_state):
            var_name = state_to_var.get(state_name)
            if var_name is not None and var_name in emissions:
                photon_ver_by_state[state_name] = emissions[var_name].to_numpy()
            elif state_name == "O2(b)" and "oxygen_a_band_photon_ver" in emissions:
                photon_ver_by_state[state_name] = emissions[
                    "oxygen_a_band_photon_ver"
                ].to_numpy()
            else:
                photon_ver_by_state[state_name] = np.zeros_like(altitude)

        return photon_ver_by_state

    def oxygen_a_band_constituent(
        self,
        emissions=None,
        state=None,
        flux=None,
        line_data=None,
        temperature_k=None,
        line_weight_model=OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING,
    ):
        if emissions is None:
            emissions = self.emissions(state=state, flux=flux)
        if line_data is None:
            line_data = self.oxygen_a_band_line_data()
        if temperature_k is None:
            if line_weight_model != OXYGEN_A_BAND_WEIGHT_MODEL_EINSTEIN_BRANCHING:
                msg = f"{line_weight_model} requires temperature_k"
                raise ValueError(msg)
            weights = line_data["weight"].to_numpy()
        else:
            temperature_k = self._temperature_on_emission_grid(emissions, temperature_k)
            weights = self._oxygen_a_band_line_weights(
                line_data,
                temperature_k,
                line_weight_model,
            )

        if weights.ndim == 1:
            weights = np.broadcast_to(weights[None, :], (emissions.sizes["altitude"], len(weights)))

        upper_vibrational_state = line_data["upper_vibrational_state"].to_numpy()
        photon_ver_by_state = self._oxygen_a_band_vibrational_photon_ver(emissions, line_data)
        line_photon_ver = np.zeros_like(weights)

        for state_name, photon_ver in photon_ver_by_state.items():
            band_mask = upper_vibrational_state == state_name
            line_photon_ver[:, band_mask] = photon_ver[:, None] * weights[:, band_mask]

        total_photon_ver = line_photon_ver.sum(axis=1)
        combined_weights = np.zeros_like(line_photon_ver)
        np.divide(
            line_photon_ver,
            total_photon_ver[:, None],
            out=combined_weights,
            where=total_photon_ver[:, None] > 0.0,
        )

        zero_rows = total_photon_ver <= 0.0
        if np.any(zero_rows):
            fallback_weights = line_data["weight"].to_numpy()
            fallback_weights = fallback_weights / fallback_weights.sum()
            combined_weights[zero_rows, :] = fallback_weights[None, :]

        return LineListVolumeEmissionRate(
            emissions["altitude"].to_numpy(),
            total_photon_ver,
            line_data["wavelength_nm"].to_numpy(),
            combined_weights,
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
        if emissions is None:
            emissions = self.emissions(state=state, flux=flux)

        constituents = {}
        if include_oxygen_green_line:
            constituents["oxygen_green"] = self.oxygen_green_line_constituent(
                emissions=emissions
            )

        if include_oxygen_a_band:
            constituents["oxygen_a_band"] = self.oxygen_a_band_constituent(
                emissions=emissions,
                line_data=line_data,
                temperature_k=temperature_k,
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
        if emissions is None:
            emissions = self.emissions(state=state, flux=flux)

        if include_oxygen_a_band and temperature_k is None:
            temperature_k = np.interp(
                emissions["altitude"].to_numpy(),
                atmosphere.model_geometry.altitudes(),
                atmosphere.temperature_k,
            )

        constituents = self.emission_constituents(
            emissions=emissions,
            line_data=line_data,
            temperature_k=temperature_k,
            line_weight_model=line_weight_model,
            include_oxygen_green_line=include_oxygen_green_line,
            include_oxygen_a_band=include_oxygen_a_band,
        )

        for name, constituent in constituents.items():
            atmosphere[name] = constituent

        return constituents


if __name__ == "__main__":
    test = Yankovsky()

    state = test.run()

    import matplotlib.pyplot as plt

    plt.subplot(1, 4, 1)
    (state["O2(a)"] / 1e6).plot(y="altitude", label="O2(a)")
    (state["O2(a, v=1)"] / 1e6).plot(y="altitude", label="O2(a, v=1)")
    (state["O2(a, v=2)"] / 1e6).plot(y="altitude", label="O2(a, v=2)")
    (state["O2(a, v=3)"] / 1e6).plot(y="altitude", label="O2(a, v=3)")
    (state["O2(a, v=4)"] / 1e6).plot(y="altitude", label="O2(a, v=4)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e11)

    plt.ylim(60000, 120000)

    plt.legend()
    plt.subplot(1, 4, 2)

    (state["O2(b)"] / 1e6).plot(y="altitude", label="O2(b)")
    (state["O2(b, v=1)"] / 1e6).plot(y="altitude", label="O2(b, v=1)")
    (state["O2(b, v=2)"] / 1e6).plot(y="altitude", label="O2(b, v=2)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e11)

    plt.ylim(60000, 120000)

    plt.legend()

    plt.xlabel("Density (cm$^{-3}$)")
    plt.subplot(1, 4, 3)

    (state["O2(X, v=1)"] / 1e6).plot(y="altitude", label="O2(X, v=1)")
    (state["O2(X, v=2)"] / 1e6).plot(y="altitude", label="O2(X, v=2)")
    (state["O2(X, v=3)"] / 1e6).plot(y="altitude", label="O2(X, v=3)")
    (state["O2(X, v=4)"] / 1e6).plot(y="altitude", label="O2(X, v=4)")
    (state["O2(X, v=5)"] / 1e6).plot(y="altitude", label="O2(X, v=5)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e11)

    plt.ylim(60000, 120000)

    plt.legend()

    plt.xlabel("Density (cm$^{-3}$)")

    plt.subplot(1, 4, 4)

    (state["O(1D)"] / 1e6).plot(y="altitude", label="O(1D)")

    plt.xscale("log")
    plt.xlim(1e-3, 1e5)

    plt.ylim(60000, 120000)

    plt.legend()

    plt.xlabel("Density (cm$^{-3}$)")

    plt.show()
