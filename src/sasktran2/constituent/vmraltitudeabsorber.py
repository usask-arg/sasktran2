from sasktran2 import Atmosphere
from sasktran2.atmosphere import NativeGridDerivative, DerivativeMapping
import numpy as np
from . import Constituent
from sasktran2.optical.base import OpticalProperty
from sasktran2.optical import pressure_temperature_to_numberdensity


class VMRAltitudeAbsorber(Constituent):
    def __init__(self, 
                 optical_property: OpticalProperty,
                 altitudes_m: np.array,
                 vmr: np.array,
                 out_of_bounds_mode: str = 'zero'
                 ) -> None:
        super().__init__()

        self._out_of_bounds_mode = out_of_bounds_mode
        self._altitudes_m = altitudes_m
        self._vmr = vmr
        self._optical_property = optical_property

    @property
    def vmr(self):
        return self._vmr

    @vmr.setter
    def vmr(self, vmr: np.array):
        self._vmr = vmr

    def name(self) -> str:
        return super().name()

    def add_to_atmosphere(self, atmo: Atmosphere):
        self._optical_quants = self._optical_property.atmosphere_quantities(atmo)

        if atmo.pressure_pa is None or atmo.temperature_k is None:
            msg = 'Both pressure_pa and temperature_k have to be specified in the atmosphere to use VMRAltitudeAbsorber'
            raise ValueError(msg)

        number_density = pressure_temperature_to_numberdensity(atmo.pressure_pa, atmo.temperature_k)

        if self._out_of_bounds_mode.lower() == 'zero':
            interp_vmr = np.interp(atmo.model_geometry.altitudes(), self._altitudes_m, self._vmr, left=0, right=0)
        elif self._out_of_bounds_mode.lower() == 'extend':
            interp_vmr = np.interp(atmo.model_geometry.altitudes(), self._altitudes_m, self._vmr, left=self._vmr[0], right=self._vmr[-1])

        atmo.storage.total_extinction += self._optical_quants.extinction * (number_density * interp_vmr)[:, np.newaxis]

    def register_derivative(self, atmo: Atmosphere):
        number_density = pressure_temperature_to_numberdensity(atmo.pressure_pa, atmo.temperature_k)
        derivs = {}

        derivs['vmr'] = DerivativeMapping(NativeGridDerivative(
            d_extinction = self._optical_quants.extinction * number_density[:, np.newaxis],
            d_ssa = self._optical_quants.extinction * (self._optical_quants.ssa - atmo.storage.ssa) / atmo.storage.total_extinction * number_density[:, np.newaxis]
        ))

        return derivs

