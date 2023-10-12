from sasktran2 import Atmosphere
from sasktran2.atmosphere import NativeGridDerivative, DerivativeMapping
import numpy as np
from .base import Constituent
from sasktran2.optical.base import OpticalProperty
from sasktran2.optical import pressure_temperature_to_numberdensity


class AirMassFactor(Constituent):
    def __init__(self,
                 ) -> None:
        """
        A dummy atmospheric constituent that does not add any terms to the atmosphere, but rather
        enables the calculation of Air Mass Factor derivatives.
        """
        super().__init__()

    def name(self) -> str:
        return 'air_mass_factor'

    def add_to_atmosphere(self, atmo: Atmosphere):
        pass

    def register_derivative(self, atmo: Atmosphere):
        altitudes = atmo.model_geometry.altitudes()

        alt_factors = -1 / np.gradient(altitudes)

        derivs = {}

        derivs['amf'] = DerivativeMapping(NativeGridDerivative(
            d_extinction = alt_factors[:, np.newaxis],
            d_ssa = alt_factors[:, np.newaxis] * (0 - atmo.storage.ssa) / atmo.storage.total_extinction
        ), summable=True, log_radiance_space=True)

        return derivs

