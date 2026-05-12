from __future__ import annotations

from scipy import constants

from sasktran2.atmosphere import Atmosphere
from sasktran2.solar import SolarModel

from .base import Constituent


class SolarIrradiance(Constituent):
    def __init__(
        self,
        photon_units: bool = False,
        **kwargs,
    ) -> None:
        """
        A constituent that calculates the solar irradiance at TOA using a reference solar irradiance
        spectrum.

        Parameters
        ----------
        photon_units: bool
            If True, the solar irradiance will be given in photon units (photons / m^2 / s / nm). If False, the solar irradiance will be given in energy units (W / m^2 / nm).
        """
        self._photon_units = photon_units
        self._model = SolarModel(**kwargs)

    def add_to_atmosphere(self, atmo: Atmosphere):
        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the SolarIrradiance constituent"
            raise ValueError(msg)

        atmo.storage.solar_irradiance[:] = self._model.irradiance(atmo.wavelengths_nm)

        if self._photon_units:
            # Photon energy in J
            photon_energy = constants.h * constants.c / (atmo.wavelengths_nm * 1e-9)
            atmo.storage.solar_irradiance[:] /= photon_energy

    def register_derivative(self, atmo, name):
        pass
