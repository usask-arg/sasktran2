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
        TODO: documentation
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

        # Unnormalized gaussian since we will normalize to vertical optical depth anyways
        self._gaussian = np.exp(
            -4
            * np.log(2)
            * (atmo.model_geometry.altitudes() - self._cloud_height_m) ** 2
            / self._cloud_width_fwhm_m**2
        )

        self._optical_quants = self._optical_property.atmosphere_quantities(
            atmo, **interped_kwargs
        )

        # TODO: should this be done by calling optical property with wavelength instead of interpolation?
        self._xs_at_wavel = (
            self._optical_quants.extinction @ wavel_interp_matrix.T
        ).flatten()
        gaussian_xs = self._gaussian * self._xs_at_wavel.flatten()
        self._gaussian_od = np.trapezoid(gaussian_xs, atmo.model_geometry.altitudes())
        number_density = (
            self._gaussian * self._vertical_optical_depth / self._gaussian_od
        )

        atmo.storage.total_extinction[:] += (
            self._optical_quants.extinction * (number_density)[:, np.newaxis]
        )

        atmo.storage.ssa[:] += (
            self._optical_quants.ssa * (number_density)[:, np.newaxis]
        )

        atmo.storage.leg_coeff[:] += (
            self._optical_quants.ssa[np.newaxis, :, :]
            * (number_density)[np.newaxis, :, np.newaxis]
            * self._optical_quants.leg_coeff
        )

        # Convert back to SSA for ease of use later in the derivatives
        self._optical_quants.ssa[:] /= self._optical_quants.extinction
        self._optical_quants.ssa[~np.isfinite(self._optical_quants.ssa)] = 1

    def register_derivative(self, atmo: Atmosphere, name):
        alt_interp_matrix = linear_interpolating_matrix(
            self._altitudes_m,
            atmo.model_geometry.altitudes(),
            "zero",
        )
        interped_kwargs = {k: alt_interp_matrix @ v for k, v in self._kwargs.items()}

        derivs = {}  # necessary?

        # common terms
        d_gaussian_d_height = (
            self._gaussian
            * 8
            * np.log(2)
            * (atmo.model_geometry.altitudes() - self._cloud_height_m)
            / self._cloud_width_fwhm_m**2
        )
        d_gaussian_d_width = (
            self._gaussian
            * 8
            * np.log(2)
            * (atmo.model_geometry.altitudes() - self._cloud_height_m) ** 2
            / self._cloud_width_fwhm_m**3
        )
        outer_term = self._vertical_optical_depth / self._gaussian_od

        d_extinction = self._optical_quants.extinction
        d_ssa = (
            self._optical_quants.extinction
            * (self._optical_quants.ssa - atmo.storage.ssa)
            / atmo.storage.total_extinction
        )
        d_leg_coeff = self._optical_quants.leg_coeff - atmo.storage.leg_coeff
        scat_factor = (self._optical_quants.ssa * self._optical_quants.extinction) / (
            atmo.storage.ssa * atmo.storage.total_extinction
        )

        # height derivative
        h_deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_height_m")
        h1 = d_gaussian_d_height
        h2 = (
            self._gaussian
            * np.trapezoid(
                self._xs_at_wavel * d_gaussian_d_height, atmo.model_geometry.altitudes()
            )
            / self._gaussian_od
        )
        h_deriv_mapping.d_extinction[:] += d_extinction
        h_deriv_mapping.d_ssa[:] += d_ssa
        h_deriv_mapping.d_leg_coeff[:] += d_leg_coeff
        h_deriv_mapping.scat_factor[:] += scat_factor
        h_deriv_mapping.interp_dim = "cloud_height_m"
        h_deriv_mapping.interpolator = outer_term * (h1 - h2)[:, np.newaxis]

        # width derivative
        w_deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_width_fwhm_m")
        w1 = d_gaussian_d_width
        w2 = (
            self._gaussian
            * np.trapezoid(
                self._xs_at_wavel * d_gaussian_d_width, atmo.model_geometry.altitudes()
            )
            / self._gaussian_od
        )
        w_deriv_mapping.d_extinction[:] += d_extinction
        w_deriv_mapping.d_ssa[:] += d_ssa
        w_deriv_mapping.d_leg_coeff[:] += d_leg_coeff
        w_deriv_mapping.scat_factor[:] += scat_factor
        w_deriv_mapping.interp_dim = (
            "cloud_width_fwhm"  # TODO: change to cloud_width_fwhm_m
        )
        w_deriv_mapping.interpolator = outer_term * (w1 - w2)[:, np.newaxis]

        # optical depth derivative
        tau_deriv_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_vertical_optical_depth"
        )
        tau_deriv_mapping.d_extinction[:] += d_extinction
        tau_deriv_mapping.d_ssa[:] += d_ssa
        tau_deriv_mapping.d_leg_coeff[:] += d_leg_coeff
        tau_deriv_mapping.scat_factor[:] += scat_factor
        tau_deriv_mapping.interp_dim = "vertical_optical_depth"
        # tau_interp_matrix = linear_interpolating_matrix(
        #     -self._extinction / self._vertical_optical_depth
        # )
        tau_deriv_mapping.interpolator = (
            self._gaussian[:, np.newaxis] / self._gaussian_od
        )

        # # TODO: not sure...
        # optical_derivs = self._optical_property.optical_derivatives(
        #     atmo, **interped_kwargs
        # )

        # for key, val in optical_derivs.items():
        #     deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_{key}")

        #     deriv_mapping.d_extinction[:] += val.d_extinction
        #     # First, the optical property returns back d_scattering extinction in the d_ssa container,
        #     # convert this to d_ssa
        #     deriv_mapping.d_ssa[:] += (
        #         val.d_ssa - val.d_extinction * self._optical_quants.ssa
        #     ) / self._optical_quants.extinction
        #     deriv_mapping.d_leg_coeff[:] += val.d_leg_coeff

        #     # Start with leg_coeff
        #     deriv_mapping.d_leg_coeff[:] += (
        #         self._optical_quants.leg_coeff - atmo.storage.leg_coeff
        #     ) * (
        #         1 / self._optical_quants.ssa * deriv_mapping.d_ssa
        #         + 1 / self._optical_quants.extinction * deriv_mapping.d_extinction
        #     )[
        #         np.newaxis, :, :
        #     ]

        #     # Then adjust d_ssa
        #     deriv_mapping.d_ssa[:] *= self._optical_quants.extinction
        #     deriv_mapping.d_ssa[:] += deriv_mapping.d_extinction * (
        #         self._optical_quants.ssa - atmo.storage.ssa
        #     )
        #     deriv_mapping.d_ssa[:] /= atmo.storage.total_extinction

        #     # TODO: The model should probably handle this
        #     norm_factor = deriv_mapping.d_leg_coeff.max(axis=0)
        #     norm_factor[norm_factor == 0] = 1

        #     deriv_mapping.scat_factor[:] = (
        #         self._optical_quants.ssa * self._optical_quants.extinction
        #     ) / (atmo.storage.ssa * atmo.storage.total_extinction)

        #     deriv_mapping.d_leg_coeff[:] /= norm_factor[np.newaxis, :, :]
        #     deriv_mapping.scat_factor[:] *= norm_factor

        #     deriv_mapping.interpolator = (
        #         alt_interp_matrix * self._number_density[np.newaxis, :]
        #     )
        #     deriv_mapping.interp_dim = f"{name}_altitude"

        return derivs
