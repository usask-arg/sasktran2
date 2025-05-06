from __future__ import annotations

import numpy as np

from sasktran2.atmosphere import Atmosphere

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

        deriv_mapping = atmo.storage.get_derivative_mapping("air_mass_factor")
        deriv_mapping.d_extinction[:] += alt_factors[:, np.newaxis]
        deriv_mapping.d_ssa[:] += (
            alt_factors[:, np.newaxis]
            * (0 - atmo.storage.ssa)
            / atmo.storage.total_extinction
        )
        deriv_mapping.log_radiance_space = True
        deriv_mapping.interp_dim = "altitude"

        return derivs
