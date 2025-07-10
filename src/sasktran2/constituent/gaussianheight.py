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
        self._cloud_height_m = np.array(cloud_height_m, dtype=float)
        self._cloud_width_fwhm_m = np.array(cloud_width_fwhm_m, dtype=float)
        self._vertical_optical_depth = np.array(vertical_optical_depth, dtype=float)

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
    def cloud_height_m(self, cloud_height_m: np.ndarray):
        self._cloud_height_m = cloud_height_m

    @property
    def cloud_width_fwhm_m(self):
        return self._cloud_width_fwhm_m

    @cloud_width_fwhm_m.setter
    def cloud_width_fwhm_m(self, cloud_width_fwhm_m: np.ndarray):
        self._cloud_width_fwhm_m = cloud_width_fwhm_m

    @property
    def vertical_optical_depth(self):
        return self._vertical_optical_depth

    @vertical_optical_depth.setter
    def vertical_optical_depth(self, vertical_optical_depth: np.ndarray):
        self._vertical_optical_depth = vertical_optical_depth

    def add_to_atmosphere(self, atmo: Atmosphere):
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
        # Unnormalized gaussian since we will normalize to vertical optical depth anyways
        # g
        self._gaussian = np.exp(-4*np.log(2)*(atmo.model_geometry.altitudes() - self._cloud_height_m)**2 / self._cloud_width_fwhm_m**2)

        # k_o
        self._optical_quants = self._optical_property.atmosphere_quantities(
            atmo, **interped_kwargs
        )

        # k_ol
        self._ext_at_wavel = (self._optical_quants.extinction @ wavel_interp_matrix.T).flatten()
        # g * k_ol
        gaussian_extinction = self._gaussian * self._ext_at_wavel.flatten()
        # Integrate[g * k_ol]
        self._current_vert_optical_depth = np.trapezoid(gaussian_extinction, atmo.model_geometry.altitudes())
        # g * k_o * Integrate[g * k_ol]
        self._extinction = self._optical_quants.extinction * (self._gaussian * self._vertical_optical_depth / self._current_vert_optical_depth)[:, np.newaxis]

        atmo.storage.total_extinction[:] += (
            self._extinction
        )

        atmo.storage.ssa[:] += self._optical_quants.ssa
        # TODO: what do?
        # atmo.storage.leg_coeff[:] += self._optical_quants.leg_coeff

    def register_derivative(self, atmo: Atmosphere, name):
        alt_interp_matrix = linear_interpolating_matrix(
            self._altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )
        interped_kwargs = {k: alt_interp_matrix @ v for k, v in self._kwargs.items()}

        # TODO: matrix to interp from alt grid to scalar (single-value array), then set mapping dimension to "scalar"

        derivs = {}  # necessary?

        # common terms
        d_gaussian_d_height = self._gaussian * 8 * np.log(2) * (atmo.model_geometry.altitudes() - self._cloud_height_m) / self._cloud_width_fwhm_m**2
        d_gaussian_d_width = self._gaussian * 8 * np.log(2) * (atmo.model_geometry.altitudes() - self._cloud_height_m)**2 / self._cloud_width_fwhm_m**3
        ext_over_tau = self._optical_quants.extinction / self._vertical_optical_depth

        # height derivative
        h_deriv_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_cloud_height_m"
        )
        h1 = self._gaussian * np.trapezoid(self._ext_at_wavel * d_gaussian_d_height, atmo.model_geometry.altitudes())
        h2 = d_gaussian_d_height * self._current_vert_optical_depth
        h_deriv_mapping.d_extinction[:] += ext_over_tau * (h1 + h2)[:, np.newaxis]
        h_deriv_mapping.d_ssa[:] += (
            self._optical_quants.extinction
            * (self._optical_quants.ssa - atmo.storage.ssa)
            / atmo.storage.total_extinction
        )
        h_deriv_mapping.d_leg_coeff[:] += (
            self._optical_quants.leg_coeff - atmo.storage.leg_coeff
        )
        # h_deriv_mapping.interp_dim = 'altitude'
        h_deriv_mapping.interp_dim = 'cloud_height_m'
        h_deriv_mapping.interpolator = np.ones_like(atmo.model_geometry.altitudes())[:, np.newaxis]

        # width derivative
        w_deriv_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_cloud_width_fwhm_m"
        )
        w1 = self._gaussian * np.trapezoid(self._ext_at_wavel * d_gaussian_d_width, atmo.model_geometry.altitudes())
        w2 = d_gaussian_d_width * self._current_vert_optical_depth
        w_deriv_mapping.d_extinction[:] += ext_over_tau * (w1 + w2)[:, np.newaxis]
        w_deriv_mapping.d_ssa[:] += (
            self._optical_quants.extinction
            * (self._optical_quants.ssa - atmo.storage.ssa)
            / atmo.storage.total_extinction
        )
        w_deriv_mapping.d_leg_coeff[:] += (
            self._optical_quants.leg_coeff - atmo.storage.leg_coeff
        )
        # w_deriv_mapping.interp_dim = 'altitude'
        w_deriv_mapping.interp_dim = 'cloud_width_fwhm'  # TODO: change to cloud_width_fwhm_m
        w_deriv_mapping.interpolator = np.ones_like(atmo.model_geometry.altitudes())[:, np.newaxis]

        # optical depth derivative
        tau_deriv_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_vertical_optical_depth"
        )
        tau_deriv_mapping.d_extinction[:] += -self._extinction / self._vertical_optical_depth
        tau_deriv_mapping.d_ssa[:] += (
            self._optical_quants.extinction
            * (self._optical_quants.ssa - atmo.storage.ssa)
            / atmo.storage.total_extinction
        )
        tau_deriv_mapping.d_leg_coeff[:] += (
            self._optical_quants.leg_coeff - atmo.storage.leg_coeff
        )
        # tau_deriv_mapping.interp_dim = 'altitude'
        tau_deriv_mapping.interp_dim = 'vertical_optical_depth'
        # tau_interp_matrix = linear_interpolating_matrix(
        #     -self._extinction / self._vertical_optical_depth
        # )
        tau_deriv_mapping.interpolator = np.ones_like(atmo.model_geometry.altitudes())[:, np.newaxis]

        return derivs

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
