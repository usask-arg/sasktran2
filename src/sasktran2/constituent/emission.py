import numpy as np

import sasktran2 as sk
from sasktran2.constants import K_BOLTZMANN, PLANCK, SPEED_OF_LIGHT

from .base import Constituent


class ThermalEmission(Constituent):
    def __init__(self):
        """
        An implementation of thermal emissions calculated from the Planck function. The emission is
        calculated with units of [W / (m^2 nm sr)].

        This Constituent requires that the atmosphere object have `temperature_k` and
        `wavelength_nm` defined inside the :py:class:`sasktran2.Atmosphere` object.
        """
        return

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere

        :meta private:
        """
        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the ThermalEmission constituent"
            raise ValueError(msg)

        if atmo.temperature_k is None:
            msg = "It is required to set the temperature_k property in the Atmosphere object to use the ThermalEmission Constituent"
            raise ValueError(msg)

        # calculate radiance in W / (m^2 nm sr)
        wavelengths_m = atmo.wavelengths_nm * 1e-9
        exponent = PLANCK * SPEED_OF_LIGHT / (wavelengths_m * K_BOLTZMANN * atmo.temperature_k[:, np.newaxis])
        blackbodyradiance = (2 * PLANCK * SPEED_OF_LIGHT ** 2 / wavelengths_m ** 5) / (np.exp(exponent) - 1) * 1e-9

        atmo.storage.emission_source += blackbodyradiance

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        return super().register_derivative(atmo=atmo, name=name)
