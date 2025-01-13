from __future__ import annotations

from sasktran2.atmosphere import Atmosphere
from sasktran2.solar import SolarModel

from .base import Constituent


class SolarIrradiance(Constituent):
    def __init__(
        self,
        **kwargs,
    ) -> None:
        """
        A constituent that calculates the solar irradiance at TOA using a reference solar irradiance
        spectrum.
        """
        self._model = SolarModel(**kwargs)

    def add_to_atmosphere(self, atmo: Atmosphere):
        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the SolarIrradiance constituent"
            raise ValueError(msg)

        atmo.storage.solar_irradiance = self._model.irradiance(atmo.wavelengths_nm)

    def register_derivative(self, atmo, name):
        pass
