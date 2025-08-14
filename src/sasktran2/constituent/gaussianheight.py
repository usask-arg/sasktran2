from __future__ import annotations

from typing import Any

import numpy as np

from sasktran2.atmosphere import Atmosphere
from sasktran2.optical.base import OpticalProperty
from sasktran2.util.interpolation import linear_interpolating_matrix

from .base import Constituent


class GaussianHeightExtinction(Constituent):
    def __init__(
        self,
        optical_property: OpticalProperty,
        cloud_height_m: float,
        cloud_width_fwhm_m: float,
        vertical_optical_depth: float,
        vertical_optical_depth_wavel_nm: float,
        altitudes_m: np.array,
        out_of_bounds_mode: str = "zero",
        **kwargs,
    ) -> None:
        """
        A constituent that is defined by a gaussian-shaped extinction profile.

        Parameters
        ----------
        optical_property : OpticalProperty
            The optical property defining the scattering information
        cloud_height_m : float
            Height of the centre of the gaussian extinction profile in [m]
        cloud_width_fwhm_m : float
            FWHM of the gaussian extinction profile in [m]
        vertical_optical_depth : float
            Vertical optical depth
        vertical_optical_depth_wavel_nm
            Wavelength that the vertical optical depth is specified at
        altitudes_m : np.array
            The altitude grid in [m] over which the optical depth is calculated, as well as the grid for any optical property arguments passed in through kwargs.
        out_of_bounds_mode : str, optional
            Interpolation mode for outside of the boundaries of the altitude grid, "extend" and "zero" are supported, by default "zero"
        kwargs : dict
            Additional arguments to pass to the optical property.
        """
        super().__init__()

        self._out_of_bounds_mode = out_of_bounds_mode
        self._altitudes_m = altitudes_m
        self._optical_property = optical_property

        # Save inputs as array datatype so that they are mutable
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
        interp_matrix = linear_interpolating_matrix(
            self._altitudes_m,
            atmo.model_geometry.altitudes(),
            self._out_of_bounds_mode,
        )

        interped_kwargs = {k: interp_matrix @ v for k, v in self._kwargs.items()}

        self._xs_at_wavel = self._optical_property.cross_sections(
            np.array([self._vertical_optical_depth_wavel_nm]),
            altitudes_m=self._altitudes_m,
            **self._kwargs,
        ).extinction.flatten()

        # Unnormalized gaussian since we will normalize to vertical optical depth anyways
        self._gaussian = np.exp(
            -4
            * np.log(2)
            * (self._altitudes_m - self._cloud_height_m) ** 2
            / self._cloud_width_fwhm_m**2
        )

        self._optical_quants = self._optical_property.atmosphere_quantities(
            atmo, **interped_kwargs
        )

        self._gaussian_od = np.trapezoid(self._gaussian, self._altitudes_m)
        self._number_density = (
            self._gaussian
            * self._vertical_optical_depth
            / self._gaussian_od
            / self._xs_at_wavel
        )
        self._interp_numden = interp_matrix @ self._number_density

        atmo.storage.total_extinction[:] += (
            self._optical_quants.extinction * (self._interp_numden)[:, np.newaxis]
        )

        atmo.storage.ssa[:] += (
            self._optical_quants.ssa * (self._interp_numden)[:, np.newaxis]
        )

        atmo.storage.leg_coeff[:] += (
            self._optical_quants.ssa[np.newaxis, :, :]
            * (self._interp_numden)[np.newaxis, :, np.newaxis]
            * self._optical_quants.leg_coeff
        )

        # Convert back to SSA for ease of use later in the derivatives
        self._optical_quants.ssa[:] /= self._optical_quants.extinction
        self._optical_quants.ssa[~np.isfinite(self._optical_quants.ssa)] = 1

    def register_derivative(self, atmo: Atmosphere, name):
        interp_matrix = linear_interpolating_matrix(
            self._altitudes_m,
            atmo.model_geometry.altitudes(),
            self._out_of_bounds_mode,
        )
        interped_kwargs = {k: interp_matrix @ v for k, v in self._kwargs.items()}

        derivs = {}

        # common terms
        d_gaussian_d_height = (
            self._gaussian
            * 8
            * np.log(2)
            * (self._altitudes_m - self._cloud_height_m)
            / self._cloud_width_fwhm_m**2
        )
        d_gaussian_d_width = (
            self._gaussian
            * 8
            * np.log(2)
            * (self._altitudes_m - self._cloud_height_m) ** 2
            / self._cloud_width_fwhm_m**3
        )
        outer_term = (
            self._vertical_optical_depth / self._gaussian_od / self._xs_at_wavel
        )

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
            * np.trapezoid(d_gaussian_d_height, self._altitudes_m)
            / self._gaussian_od
        )
        h_deriv_mapping.d_extinction[:] += d_extinction
        h_deriv_mapping.d_ssa[:] += d_ssa
        h_deriv_mapping.d_leg_coeff[:] += d_leg_coeff
        h_deriv_mapping.scat_factor[:] += scat_factor
        h_deriv_mapping.interp_dim = "cloud_height_m"
        h_deriv_mapping.interpolator = interp_matrix @ (
            outer_term * (h1 - h2)[:, np.newaxis]
        )

        # width derivative
        w_deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_width_fwhm_m")
        w1 = d_gaussian_d_width
        w2 = (
            self._gaussian
            * np.trapezoid(d_gaussian_d_width, self._altitudes_m)
            / self._gaussian_od
        )
        w_deriv_mapping.d_extinction[:] += d_extinction
        w_deriv_mapping.d_ssa[:] += d_ssa
        w_deriv_mapping.d_leg_coeff[:] += d_leg_coeff
        w_deriv_mapping.scat_factor[:] += scat_factor
        w_deriv_mapping.interp_dim = "cloud_width_fwhm_m"
        w_deriv_mapping.interpolator = interp_matrix @ (
            outer_term * (w1 - w2)[:, np.newaxis]
        )

        # optical depth derivative
        tau_deriv_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_vertical_optical_depth"
        )
        tau_deriv_mapping.d_extinction[:] += d_extinction
        tau_deriv_mapping.d_ssa[:] += d_ssa
        tau_deriv_mapping.d_leg_coeff[:] += d_leg_coeff
        tau_deriv_mapping.scat_factor[:] += scat_factor
        tau_deriv_mapping.interp_dim = "vertical_optical_depth"
        tau_deriv_mapping.interpolator = interp_matrix @ (
            (self._gaussian / self._gaussian_od / self._xs_at_wavel)[:, np.newaxis]
        )

        optical_derivs = self._optical_property.optical_derivatives(
            atmo, **interped_kwargs
        )
        vertical_deriv_factor = 1 / self._xs_at_wavel
        d_vertical_deriv_factor = self._optical_property.cross_section_derivatives(
            np.array([self._vertical_optical_depth_wavel_nm]),
            altitudes_m=self._altitudes_m,
            **self._kwargs,
        )
        for _, val in d_vertical_deriv_factor.items():
            # convert from derivative of x to derivative of 1/x
            val *= -1 * vertical_deriv_factor**2  # noqa: PLW2901

        for key, val in optical_derivs.items():
            # Code copied from numdenscatterer.py
            deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_{key}")

            deriv_mapping.d_extinction[:] += val.d_extinction
            d_extinction_scat = val.d_ssa

            # First, the optical property returns back d_scattering extinction in the d_ssa container,
            # convert this to d_ssa
            deriv_mapping.d_ssa[:] += (
                d_extinction_scat - val.d_extinction * self._optical_quants.ssa
            ) / self._optical_quants.extinction
            deriv_mapping.d_leg_coeff[:] += val.d_leg_coeff

            if key in d_vertical_deriv_factor:
                # Have to make some adjustments

                # The change in extinction is adjusted
                deriv_mapping.d_extinction[:] += (
                    self._optical_quants.extinction
                    / (interp_matrix @ vertical_deriv_factor)[:, np.newaxis]
                    * (interp_matrix @ d_vertical_deriv_factor[key])[:, np.newaxis]
                )

                # Change in single scatter albedo should be invariant whether or not we are
                # in extinction space or number density space

            # Start with leg_coeff
            deriv_mapping.d_leg_coeff[:] += (
                self._optical_quants.leg_coeff - atmo.storage.leg_coeff
            ) * (
                1 / self._optical_quants.ssa * deriv_mapping.d_ssa
                + 1 / self._optical_quants.extinction * deriv_mapping.d_extinction
            )[
                np.newaxis, :, :
            ]

            # Then adjust d_ssa
            deriv_mapping.d_ssa[:] *= self._optical_quants.extinction
            deriv_mapping.d_ssa[:] += deriv_mapping.d_extinction * (
                self._optical_quants.ssa - atmo.storage.ssa
            )
            deriv_mapping.d_ssa[:] /= atmo.storage.total_extinction

            # TODO: The model should probably handle this
            norm_factor = deriv_mapping.d_leg_coeff.max(axis=0)
            norm_factor[norm_factor == 0] = 1

            deriv_mapping.scat_factor[:] = (
                self._optical_quants.ssa * self._optical_quants.extinction
            ) / (atmo.storage.ssa * atmo.storage.total_extinction)

            deriv_mapping.d_leg_coeff[:] /= norm_factor[np.newaxis, :, :]
            deriv_mapping.scat_factor[:] *= norm_factor

            deriv_mapping.interpolator = (
                interp_matrix * self._number_density[np.newaxis, :]
            )
            deriv_mapping.interp_dim = f"{name}_altitude"

        return derivs
