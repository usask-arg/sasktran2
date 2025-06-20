from __future__ import annotations

from typing import Any

import numpy as np

from sasktran2.atmosphere import Atmosphere
from sasktran2.optical.base import OpticalProperty
from sasktran2.util.interpolation import linear_interpolating_matrix

from .base import Constituent


class GaussianHeight(Constituent):
    def __init__(
        self,
        optical_property: OpticalProperty,
        altitudes_m: np.array,
        cloud_height_m: float,
        cloud_width_fwhm_m: float,
        vertical_optical_depth: float,
        vertical_optical_depth_wavel_nm: float,
        **kwargs,
    ) -> None:
        """
        A constituent that is defined by a gaussian-shaped extinction profile,
        such as a cloud.

        Parameters
        ----------
        TODO
        """
        super().__init__()

    def add_to_atmosphere(self, atmo):
        return super().add_to_atmosphere(atmo)

    def register_derivative(self, atmo, name):
        return super().register_derivative(atmo, name)
