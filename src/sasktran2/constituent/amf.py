import numpy as np

from sasktran2 import Atmosphere
from sasktran2.atmosphere import DerivativeMapping, NativeGridDerivative

from .base import Constituent


class AirMassFactor(Constituent):
    def __init__(
        self,
    ) -> None:
        """
        A dummy atmospheric constituent that does not add any terms to the atmosphere, but rather
        enables the calculation of Air Mass Factor derivatives.
        """
        super().__init__()

    def add_to_atmosphere(self, atmo: Atmosphere):
        pass

    def register_derivative(self, atmo: Atmosphere, name: str):  # noqa: ARG002
        altitudes = atmo.model_geometry.altitudes()

        alt_factors = -1 / np.gradient(altitudes)

        # Need to adjust bottom of atmosphere and top to match AMF definition
        alt_factors[0] *= 2
        alt_factors[-1] *= 2

        derivs = {}

        derivs["air_mass_factor"] = DerivativeMapping(
            NativeGridDerivative(
                d_extinction=alt_factors[:, np.newaxis],
                d_ssa=alt_factors[:, np.newaxis]
                * (0 - atmo.storage.ssa)
                / atmo.storage.total_extinction,
            ),
            summable=True,
            log_radiance_space=True,
            name_prefix="",
        )

        return derivs
