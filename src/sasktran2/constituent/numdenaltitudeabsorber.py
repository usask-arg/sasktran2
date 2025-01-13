from __future__ import annotations

from sasktran2.atmosphere import Atmosphere
from sasktran2.optical.base import OpticalProperty

from .base import Constituent


class AltitudeAbsorber(Constituent):
    def __init__(self, optical_property: OpticalProperty) -> None:
        super().__init__()

        self._optical_property = optical_property

    def add_to_atmosphere(self, atmo: Atmosphere):
        return super().add_to_atmosphere(atmo)

    def register_derivative(self, atmo: Atmosphere, name: str):  # noqa: ARG002
        return super().register_derivative(atmo)
