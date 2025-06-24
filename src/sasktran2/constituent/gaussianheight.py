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
        altitudes_m: np.array,  # TODO: is it even necessary to define altitudes_m (maybe if there are altitude-dependent kwargs)
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

        self._altitudes_m = altitudes_m
        self._optical_property = optical_property

        # inputs with derivatives
        self._cloud_height_m = cloud_height_m
        self._cloud_width_fwhm_m = cloud_width_fwhm_m
        self._vertical_optical_depth = vertical_optical_depth

        self._vertical_optical_depth_wavel_nm = vertical_optical_depth_wavel_nm

        self._kwargs = kwargs

    # TODO: not sure if getattr and setattr are needed...
    def __getattr__(self, __name: str) -> Any:
        if __name in self.__dict__.get("_kwargs", {}):
            return self._kwargs[__name]
        return None

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name in self.__dict__.get("_kwargs", {}):
            self._kwargs[__name] = __value
        else:
            super().__setattr__(__name, __value)

    @property
    def cloud_height_m(self):
        return self._cloud_height_m

    @cloud_height_m.setter
    def cloud_height_m(self, cloud_height_m: float):
        self._cloud_height_m = cloud_height_m

    @property
    def cloud_width_fwhm_m(self):
        return self._cloud_width_fwhm_m

    @cloud_width_fwhm_m.setter
    def cloud_width_fwhm_m(self, cloud_width_fwhm_m: float):
        self._cloud_width_fwhm_m = cloud_width_fwhm_m

    @property
    def vertical_optical_depth(self):
        return self._vertical_optical_depth

    @vertical_optical_depth.setter
    def vertical_optical_depth(self, vertical_optical_depth: float):
        self._vertical_optical_depth = vertical_optical_depth

    def add_to_atmosphere(self, atmo):
        alt_interp_matrix = linear_interpolating_matrix(
            self._altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )
        wavel_interp_matrix = linear_interpolating_matrix(
            atmo.wavelengths_nm,
            [self._vertical_optical_depth_wavel_nm],
            "zero",
        )

        interped_kwargs = {k: alt_interp_matrix @ v for k, v in self._kwargs.items()}

        # Cloud top is defined as cloud height + 2 sigma
        cloud_sigma = self._cloud_width_fwhm_m / (2 * np.sqrt(2 * np.log(2)))

        # Unnormalized gaussian since we will normalize to vertical optical depth anyways
        cloud_numden = np.exp(-(atmo.model_geometry.altitudes() - self._cloud_height_m)**2 / (2 * cloud_sigma**2))

        # optical_quants_at_wavel = self._optical_property.atmosphere_quantities(
        #     atmo, wavelengths_nm=[self._vertical_optical_depth_wavel_nm], **interped_kwargs
        # )

        self._optical_quants = self._optical_property.atmosphere_quantities(
            atmo, **interped_kwargs
        )

        ext_at_wavel = self._optical_quants.extinction @ wavel_interp_matrix.T
        current_vert_optical_depth = np.trapezoid(cloud_numden * ext_at_wavel.flatten(), atmo.model_geometry.altitudes())
        cloud_numden *= self._vertical_optical_depth / current_vert_optical_depth

        atmo.storage.total_extinction[:] += (
            self._optical_quants.extinction * (cloud_numden)[:, np.newaxis]
        )

        atmo.storage.ssa[:] += self._optical_quants.ssa
        # atmo.storage.leg_coeff[:] += self._optical_quants.leg_coeff

    def register_derivative(self, atmo, name):
        return super().register_derivative(atmo, name)

# class SpeciesGaussianHeight(Species):
#     def __init__(self, optical_property: sk.OpticalProperty, altitude_grid: np.array, name: str, cloud_height_m: float,
#                  cloud_width_fwhm_m: float, vertical_optical_depth: float,
#                  vertical_optical_depth_wavel_nm: float):
#         cloud_sigma = cloud_width_fwhm_m / (2 * np.sqrt(2 * np.log(2)))

#         # Cloud top is defined as middle + 2 sigma
#         cloud_middle_m = cloud_height_m
#         cloud_heights = altitude_grid

#         # Unnormalized gaussian since we will normalize to vertical optical depth anyways
#         cloud_numden = np.exp(-(cloud_heights - cloud_middle_m)**2 / (2 * cloud_sigma**2))

#         xs_at_wavel = optical_property.calculate_cross_sections(sk.MSIS90(), 0, 0, 0, 54372,
#                                                                 vertical_optical_depth_wavel_nm).total

#         current_vert_optical_depth = np.trapz(cloud_numden * xs_at_wavel, altitude_grid)

#         cloud_numden *= vertical_optical_depth / current_vert_optical_depth

#         climatology = sk.ClimatologyUserDefined(altitude_grid, {name: cloud_numden},
#                                                 out_of_bounds_value=0.0)

#         super().__init__(optical_property, climatology)
