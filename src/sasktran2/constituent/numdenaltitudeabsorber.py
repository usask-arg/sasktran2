from sasktran2 import Atmosphere
from .base import Constituent
from sasktran2.optical.base import OpticalProperty


class AltitudeAbsorber(Constituent):
    def __init__(self, optical_property: OpticalProperty) -> None:
        super().__init__()

        self._optical_property = optical_property
    
    def name(self) -> str:
        return super().name()
    
    def add_to_atmosphere(self, atmo: Atmosphere):
        optical_quants = self._optical_property.atmosphere_quantities(atmo)

        return super().add_to_atmosphere(atmo)
    
    def register_derivative(self, atmo: Atmosphere):
        return super().register_derivative(atmo)
    
