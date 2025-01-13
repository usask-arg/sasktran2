from __future__ import annotations

import numpy as np

from sasktran2.constants import K_BOLTZMANN


class EquationOfState:
    def __init__(self, molar_mass_dry_air: float = 28.96e-3) -> None:
        """
        The equation of state relates pressure and temperature to number density and its derivatives.
        In addition we include specific humidity to convert between wet air and dry air

        Parameters
        ----------
        molar_mass_dry_air : float, optional
            Molar mass of dry air, by default 28.96e-3
        """
        self._pressure_pa = None
        self._temperature_k = None
        self._specific_humidity = None
        self._M = molar_mass_dry_air

    @property
    def pressure_pa(self) -> np.array:
        """
        Returns
        -------
        np.array
            Pressure in Pa
        """
        return self._pressure_pa

    @pressure_pa.setter
    def pressure_pa(self, pressure: np.array):
        self._pressure_pa = pressure

    @property
    def temperature_k(self) -> np.array:
        """

        Returns
        -------
        np.array
            Temperature in K
        """
        return self._temperature_k

    @temperature_k.setter
    def temperature_k(self, temperature: np.array):
        self._temperature_k = temperature

    @property
    def specific_humidity(self) -> np.array:
        """

        Returns
        -------
        np.array
            Specific humidity
        """
        return self._specific_humidity

    @specific_humidity.setter
    def specific_humidity(self, specific_humidity: np.array):
        self._specific_humidity = specific_humidity

    @property
    def air_numberdensity(self) -> np.array:
        """
        Converts pressure and temperature to number density using the ideal gas law at the grid sample points

        Returns
        -------
        dict with keys
            "N": Number density in [molecules / m^3]
            "dN_dP": Derivative of number density with respect to pressure in [molecules / m^3 / Pa]
            "dN_dT": Derivative of number density with respect to temperature in [molecules / m^3 / K]
        """
        # Ideal gas law is PV = N k T
        # Or (N/V) = P / (k T)

        return {
            "N": self.pressure_pa / (K_BOLTZMANN * self.temperature_k),
            "dN_dP": 1 / (K_BOLTZMANN * self.temperature_k),
            "dN_dT": -self.pressure_pa / (K_BOLTZMANN * self.temperature_k**2),
        }

    @property
    def dry_air_numberdensity(self) -> np.array:
        """
        Converts pressure and temperature to dry air number density using the ideal gas law at the grid sample points

        Returns
        -------
        dict with keys
            "N": Number density in [molecules / m^3]
            "dN_dP": Derivative of number density with respect to pressure in [molecules / m^3 / Pa]
            "dN_dT": Derivative of number density with respect to temperature in [molecules / m^3 / K]
            "dN_dsh": Derivative of number density with respect to specific humidity in [molecules / m^3 / kg/kg]
        """
        wet = self.air_numberdensity

        q = 0 if self.specific_humidity is None else self.specific_humidity

        return {
            "N": wet["N"] * (1 - q),
            "dN_dP": wet["dN_dP"] * (1 - q),
            "dN_dT": wet["dN_dT"] * (1 - q),
            "dN_dsh": -wet["N"],
        }
