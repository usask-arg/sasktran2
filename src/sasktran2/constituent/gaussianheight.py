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
        cloud_height_m: float,
        cloud_width_fwhm_m: float,
        vertical_optical_depth: float,
        vertical_optical_depth_wavel_nm: float,
        altitudes_m: np.array = None,
        out_of_bounds_mode: str = "zero",
        fractional_change=0.01,
        numerical_debug=False,
        **kwargs,
    ) -> None:
        """
        A constituent that is defined by a gaussian-shaped extinction profile,
        such as a cloud.

        Parameters
        ----------
        optical_property : OpticalProperty
            The optical property defining the scattering information
        cloud_height_m : float
            Height of the centre of the gaussian profile in [m]
        cloud_width_fwhm_m : float
            FWHM of the gaussian profile in [m]
        vertical_optical_depth : float
            Vertical optical depth
        vertical_optical_depth_wavel_nm
            Wavelength that the vertical optical depth is specified at
        altitudes_m : np.array, optional
            The altitude grid in [m] for optical property arguments passed in through kwargs. If not specified,
            the atmosphere altitude grid is used.
        out_of_bounds_mode : str, optional
            Interpolation mode for outside of the boundaries, "extend" and "zero" are supported, by default "zero"
        kwargs : dict
            Additional arguments to pass to the optical property.
        """
        super().__init__()

        self._out_of_bounds_mode = out_of_bounds_mode
        self._altitudes_m = altitudes_m
        self._optical_property = optical_property

        # inputs with derivatives
        self._cloud_height_m = np.array(cloud_height_m, dtype=float)
        self._cloud_width_fwhm_m = np.array(cloud_width_fwhm_m, dtype=float)
        self._vertical_optical_depth = np.array(vertical_optical_depth, dtype=float)

        self._vertical_optical_depth_wavel_nm = vertical_optical_depth_wavel_nm

        self._kwargs = kwargs

        self._fractional_change = fractional_change
        self._numerical_debug = numerical_debug

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
        if self._altitudes_m is None:
            self._altitudes_m = atmo.model_geometry.altitudes()
        alt_interp_matrix = linear_interpolating_matrix(
            self._altitudes_m,
            atmo.model_geometry.altitudes(),
            self._out_of_bounds_mode,
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

        # TODO: should we account for the possibility of no leg_coeff?
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
            self._out_of_bounds_mode,
        )
        wavel_interp_matrix = linear_interpolating_matrix(
            atmo.wavelengths_nm,
            [self._vertical_optical_depth_wavel_nm],
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
            "cloud_width_fwhm_m"
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

        optical_derivs = self._optical_property.optical_derivatives(
            atmo, **interped_kwargs
        )

        number_density = (
            self._gaussian * self._vertical_optical_depth / self._gaussian_od
        )

        #####
        extinction_to_numden_factors = self._optical_property.cross_sections(
            np.array([self._vertical_optical_depth_wavel_nm]),
            altitudes_m=self._altitudes_m,
            **self._kwargs,
        ).extinction.flatten()
        vertical_deriv_factor = 1 / extinction_to_numden_factors
        d_vertical_deriv_factor = (
            self._optical_property.cross_section_derivatives(
                np.array([self._vertical_optical_depth_wavel_nm]),
                altitudes_m=self._altitudes_m,
                **self._kwargs,
            )
        )
        for _, val in d_vertical_deriv_factor.items():
            # convert from derivative of x to derivative of 1/x
            val *= -1 * vertical_deriv_factor**2  # noqa: PLW2901
        #####

        for (key, val) in optical_derivs.items():
            if self._numerical_debug:
                # only works if this is the only constituent
                deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_{key}")

                deriv_mapping.scat_factor[:] = np.ones_like(deriv_mapping.d_extinction)
                deriv_mapping.interpolator = (
                    alt_interp_matrix
                )
                deriv_mapping.interp_dim = f"{name}_altitude"

                d_xs_at_wavel = (
                    val.d_extinction @ wavel_interp_matrix.T
                ).flatten()

                d_n = np.trapezoid(d_xs_at_wavel * self._gaussian, atmo.model_geometry.altitudes()) / self._gaussian_od

                #######################

                x_base = interped_kwargs[key]
                # deriv_mapping.interpolator = np.zeros([len(x_base), len(atmo.model_geometry.altitudes())])
                for i in range(len(x_base)):
                    x_down = x_base.copy()
                    x_up = x_base.copy()

                    x_down[i] = x_down[i] * (1 - self._fractional_change)
                    x_up[i] = x_up[i] * (1 + self._fractional_change)

                    kwargs_down = {key: x_down}
                    kwargs_up = {key: x_up}

                    optical_quants_down = self._optical_property.atmosphere_quantities(
                        atmo, **kwargs_down
                    )
                    optical_quants_up = self._optical_property.atmosphere_quantities(
                        atmo, **kwargs_up
                    )

                    xs_at_wavel_down = (optical_quants_down.extinction @ wavel_interp_matrix.T).flatten()
                    xs_at_wavel_up = (optical_quants_up.extinction @ wavel_interp_matrix.T).flatten()

                    gaussian_xs_down = self._gaussian * xs_at_wavel_down
                    gaussian_xs_up = self._gaussian * xs_at_wavel_up

                    gaussian_od_down = np.trapezoid(gaussian_xs_down, atmo.model_geometry.altitudes())
                    gaussian_od_up = np.trapezoid(gaussian_xs_up, atmo.model_geometry.altitudes())

                    number_density_down = self._gaussian * self._vertical_optical_depth / gaussian_od_down
                    number_density_up = self._gaussian * self._vertical_optical_depth / gaussian_od_up

                    total_extinction_down = optical_quants_down.extinction * (number_density_down)[:, np.newaxis]
                    total_extinction_up = optical_quants_up.extinction * (number_density_up)[:, np.newaxis]

                    ssa_down = optical_quants_down.ssa * (number_density_down)[:, np.newaxis]
                    ssa_up = optical_quants_up.ssa * (number_density_up)[:, np.newaxis]

                    leg_coeff_down = optical_quants_down.leg_coeff
                    leg_coeff_up = optical_quants_up.leg_coeff

                    ssa_down = optical_quants_down.ssa / optical_quants_down.extinction
                    ssa_up = optical_quants_up.ssa / optical_quants_up.extinction

                    ssa_down[~np.isfinite(optical_quants_down.ssa)] = 1
                    ssa_up[~np.isfinite(optical_quants_up.ssa)] = 1

                    # deriv_mapping.d_extinction[i] += (
                    #     (total_extinction_up[i] - total_extinction_down[i]) / (x_up[i] - x_down[i])
                    # )
                    # deriv_mapping.d_ssa[i] += (
                    #     (ssa_up[i] - ssa_down[i]) / (x_up[i] - x_down[i])
                    # )
                    # deriv_mapping.d_leg_coeff[:, i] += (
                    #     (leg_coeff_up[:, i] - leg_coeff_down[:, i]) / (x_up[i] - x_down[i])
                    # )
                    deriv_mapping.d_extinction[:] += (
                        (total_extinction_up - total_extinction_down) / (x_up[i] - x_down[i])
                    )
                    deriv_mapping.d_ssa[:] += (
                        (ssa_up - ssa_down) / (x_up[i] - x_down[i])
                    )
                    deriv_mapping.d_leg_coeff[:] += (
                        (leg_coeff_up - leg_coeff_down) / (x_up[i] - x_down[i])
                    )

                    # deriv_mapping.d_leg_coeff[i] += (
                        # (self._leg_coeff_up / self._ssa_up[np.newaxis, :, :] - self._leg_coeff_down / self._ssa_down[np.newaxis, :, :]) / (x_up - x_down)[np.newaxis, :, np.newaxis]
                    # )
                    # deriv_mapping.d_leg_coeff[i] += (
                    #     (self._leg_coeff_up / (self._ssa_up * self._total_extinction_up)[np.newaxis, :, :] - self._leg_coeff_down / (self._ssa_down * self._total_extinction_down)[np.newaxis, :, :]) / (x_up - x_down)[np.newaxis, :, np.newaxis]
                    # )
                    # deriv_mapping.d_leg_coeff[:, i] += (
                        # (leg_coeff_up[:, i] / (ssa_up[np.newaxis, i] * total_extinction_up[np.newaxis, i]) - leg_coeff_down[:, i] / (ssa_down[np.newaxis, i] * total_extinction_down[np.newaxis, i])) / (x_up[i] - x_down[i])
                    # )

                    deriv_mapping.interpolator[i] = (number_density / number_density) * (1 - x_base[i] * d_n)

            else:
                deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_{key}")

                d_xs_at_wavel = (
                    val.d_extinction @ wavel_interp_matrix.T
                ).flatten()

                d_n = np.trapezoid(d_xs_at_wavel * self._gaussian, atmo.model_geometry.altitudes()) / self._gaussian_od

                # """
                #####
                deriv_mapping.d_extinction[:] += self._optical_quants.extinction * self._gaussian[:, np.newaxis] * self._vertical_optical_depth
                # d_extinction_scat = val.d_ssa - self._optical_quants.ssa * self._optical_quants.extinction * d_n[np.newaxis, np.newaxis]
                # deriv_mapping.d_extinction[:] += val.d_extinction
                d_extinction_scat = val.d_ssa
                #####

                # First, the optical property returns back d_scattering extinction in the d_ssa container,
                # convert this to d_ssa
                deriv_mapping.d_ssa[:] += (
                    d_extinction_scat - val.d_extinction * self._optical_quants.ssa
                ) / self._optical_quants.extinction
                deriv_mapping.d_leg_coeff[:] += val.d_leg_coeff

                #####
                if key in d_vertical_deriv_factor:
                    # Have to make some adjustments

                    # The change in extinction is adjusted
                    # deriv_mapping.d_extinction[:] += (
                    #     self._optical_quants.extinction
                    #     / (alt_interp_matrix @ vertical_deriv_factor)[:, np.newaxis]
                    #     * (alt_interp_matrix @ d_vertical_deriv_factor[key])[
                    #         :, np.newaxis
                    #     ]
                    # )
                    # deriv_mapping.d_extinction[:] += (
                    #     - self._optical_quants.extinction * d_n[np.newaxis, np.newaxis]
                    # )
                    pass

                    # Change in single scatter albedo should be invariant whether or not we are
                    # in extinction space or number density space
                #####

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

                # deriv_mapping.interpolator = -self._vertical_optical_depth * self._gaussian[np.newaxis, :] * np.trapezoid(d_xs_at_wavel * self._gaussian, atmo.model_geometry.altitudes()) / self._gaussian_od**2
                # deriv_mapping.interpolator = (
                #     alt_interp_matrix * number_density[np.newaxis, :]
                # )
                deriv_mapping.interp_dim = f"{name}_altitude"
                # """

                # construct interpolator
                interp_top = -np.trapezoid(np.outer(d_xs_at_wavel, self._gaussian), atmo.model_geometry.altitudes())
                interp_bot = np.trapezoid(np.outer(self._xs_at_wavel, self._gaussian), atmo.model_geometry.altitudes())**2
                deriv_mapping.interpolator = (
                    np.outer(np.ones_like(atmo.model_geometry.altitudes()), (interp_top / interp_bot))
                )

                # deriv_mapping.d_extinction[:] += d_extinction
                # deriv_mapping.d_ssa[:] += d_ssa
                # deriv_mapping.d_leg_coeff[:] += d_leg_coeff
                # deriv_mapping.scat_factor[:] += scat_factor
                # deriv_mapping.interp_dim = f"{name}_altitude"
                # deriv_mapping.interpolator = alt_interp_matrix * (-number_density * d_n)[np.newaxis, :]

        return derivs


class GaussianHeightExtinction(Constituent):
    def __init__(
        self,
        optical_property: OpticalProperty,
        cloud_height_m: float,
        cloud_width_fwhm_m: float,
        vertical_optical_depth: float,
        vertical_optical_depth_wavel_nm: float,
        altitudes_m: np.array = None,
        out_of_bounds_mode: str = "zero",
        fractional_change=0.01,
        numerical_debug=False,
        **kwargs,
    ) -> None:
        """
        A constituent that is defined by a gaussian-shaped extinction profile,
        such as a cloud.

        Parameters
        ----------
        optical_property : OpticalProperty
            The optical property defining the scattering information
        cloud_height_m : float
            Height of the centre of the gaussian profile in [m]
        cloud_width_fwhm_m : float
            FWHM of the gaussian profile in [m]
        vertical_optical_depth : float
            Vertical optical depth
        vertical_optical_depth_wavel_nm
            Wavelength that the vertical optical depth is specified at
        altitudes_m : np.array, optional
            The altitude grid in [m] for optical property arguments passed in through kwargs. If not specified,
            the atmosphere altitude grid is used.
        out_of_bounds_mode : str, optional
            Interpolation mode for outside of the boundaries, "extend" and "zero" are supported, by default "zero"
        kwargs : dict
            Additional arguments to pass to the optical property.
        """
        super().__init__()

        self._out_of_bounds_mode = out_of_bounds_mode
        self._altitudes_m = altitudes_m
        self._optical_property = optical_property

        # inputs with derivatives
        self._cloud_height_m = np.array(cloud_height_m, dtype=float)
        self._cloud_width_fwhm_m = np.array(cloud_width_fwhm_m, dtype=float)
        self._vertical_optical_depth = np.array(vertical_optical_depth, dtype=float)

        self._vertical_optical_depth_wavel_nm = vertical_optical_depth_wavel_nm

        self._kwargs = kwargs

        self._fractional_change = fractional_change
        self._numerical_debug = numerical_debug

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
        if self._altitudes_m is None:
            self._altitudes_m = atmo.model_geometry.altitudes()
        alt_interp_matrix = linear_interpolating_matrix(
            self._altitudes_m,
            atmo.model_geometry.altitudes(),
            self._out_of_bounds_mode,
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

        self._xs_at_wavel = (
            self._optical_quants.extinction @ wavel_interp_matrix.T
        ).flatten()
        self._gaussian_od = np.trapezoid(self._gaussian, atmo.model_geometry.altitudes())
        number_density = (
            self._gaussian * self._vertical_optical_depth / self._gaussian_od / self._xs_at_wavel
        )

        atmo.storage.total_extinction[:] += (
            self._optical_quants.extinction * (number_density)[:, np.newaxis]
        )

        atmo.storage.ssa[:] += (
            self._optical_quants.ssa * (number_density)[:, np.newaxis]
        )

        # TODO: should we account for the possibility of no leg_coeff?
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
            self._out_of_bounds_mode,
        )
        wavel_interp_matrix = linear_interpolating_matrix(
            atmo.wavelengths_nm,
            [self._vertical_optical_depth_wavel_nm],
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
        outer_term = self._vertical_optical_depth / self._gaussian_od / self._xs_at_wavel

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
                d_gaussian_d_height, atmo.model_geometry.altitudes()
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
                d_gaussian_d_width, atmo.model_geometry.altitudes()
            )
            / self._gaussian_od
        )
        w_deriv_mapping.d_extinction[:] += d_extinction
        w_deriv_mapping.d_ssa[:] += d_ssa
        w_deriv_mapping.d_leg_coeff[:] += d_leg_coeff
        w_deriv_mapping.scat_factor[:] += scat_factor
        w_deriv_mapping.interp_dim = (
            "cloud_width_fwhm_m"
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
            (self._gaussian / self._gaussian_od / self._xs_at_wavel)[:, np.newaxis]
        )

        optical_derivs = self._optical_property.optical_derivatives(
            atmo, **interped_kwargs
        )

        number_density = (
            self._gaussian * self._vertical_optical_depth / self._gaussian_od / self._xs_at_wavel
        )

        #####
        extinction_to_numden_factors = self._optical_property.cross_sections(
            np.array([self._vertical_optical_depth_wavel_nm]),
            altitudes_m=self._altitudes_m,
            **self._kwargs,
        ).extinction.flatten()
        vertical_deriv_factor = 1 / extinction_to_numden_factors
        d_vertical_deriv_factor = (
            self._optical_property.cross_section_derivatives(
                np.array([self._vertical_optical_depth_wavel_nm]),
                altitudes_m=self._altitudes_m,
                **self._kwargs,
            )
        )
        for _, val in d_vertical_deriv_factor.items():
            # convert from derivative of x to derivative of 1/x
            val *= -1 * vertical_deriv_factor**2  # noqa: PLW2901
        #####

        for (key, val) in optical_derivs.items():
            if self._numerical_debug:
                # only works if this is the only constituent
                deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_{key}")

                deriv_mapping.scat_factor[:] = np.ones_like(deriv_mapping.d_extinction)
                deriv_mapping.interpolator = (
                    alt_interp_matrix
                )
                deriv_mapping.interp_dim = f"{name}_altitude"

                #######################

                x_base = interped_kwargs[key]
                for i in range(len(x_base)):
                    x_down = x_base.copy()
                    x_up = x_base.copy()

                    x_down[i] = x_down[i] * (1 - self._fractional_change)
                    x_up[i] = x_up[i] * (1 + self._fractional_change)

                    kwargs_down = {key: x_down}
                    kwargs_up = {key: x_up}

                    optical_quants_down = self._optical_property.atmosphere_quantities(
                        atmo, **kwargs_down
                    )
                    optical_quants_up = self._optical_property.atmosphere_quantities(
                        atmo, **kwargs_up
                    )

                    xs_at_wavel_down = (optical_quants_down.extinction @ wavel_interp_matrix.T).flatten()
                    xs_at_wavel_up = (optical_quants_up.extinction @ wavel_interp_matrix.T).flatten()

                    gaussian_xs_down = self._gaussian * xs_at_wavel_down
                    gaussian_xs_up = self._gaussian * xs_at_wavel_up

                    gaussian_od_down = np.trapezoid(gaussian_xs_down, atmo.model_geometry.altitudes())
                    gaussian_od_up = np.trapezoid(gaussian_xs_up, atmo.model_geometry.altitudes())

                    number_density_down = self._gaussian * self._vertical_optical_depth / gaussian_od_down
                    number_density_up = self._gaussian * self._vertical_optical_depth / gaussian_od_up

                    total_extinction_down = optical_quants_down.extinction * (number_density_down)[:, np.newaxis]
                    total_extinction_up = optical_quants_up.extinction * (number_density_up)[:, np.newaxis]

                    ssa_down = optical_quants_down.ssa * (number_density_down)[:, np.newaxis]
                    ssa_up = optical_quants_up.ssa * (number_density_up)[:, np.newaxis]

                    leg_coeff_down = optical_quants_down.leg_coeff
                    leg_coeff_up = optical_quants_up.leg_coeff

                    ssa_down = optical_quants_down.ssa / optical_quants_down.extinction
                    ssa_up = optical_quants_up.ssa / optical_quants_up.extinction

                    ssa_down[~np.isfinite(optical_quants_down.ssa)] = 1
                    ssa_up[~np.isfinite(optical_quants_up.ssa)] = 1

                    # deriv_mapping.d_extinction[i] += (
                    #     (total_extinction_up[i] - total_extinction_down[i]) / (x_up[i] - x_down[i])
                    # )
                    # deriv_mapping.d_ssa[i] += (
                    #     (ssa_up[i] - ssa_down[i]) / (x_up[i] - x_down[i])
                    # )
                    # deriv_mapping.d_leg_coeff[:, i] += (
                    #     (leg_coeff_up[:, i] - leg_coeff_down[:, i]) / (x_up[i] - x_down[i])
                    # )
                    deriv_mapping.d_extinction[:] += (
                        (total_extinction_up - total_extinction_down) / (x_up[i] - x_down[i])
                    )
                    deriv_mapping.d_ssa[:] += (
                        (ssa_up - ssa_down) / (x_up[i] - x_down[i])
                    )
                    deriv_mapping.d_leg_coeff[:] += (
                        (leg_coeff_up - leg_coeff_down) / (x_up[i] - x_down[i])
                    )

                    # deriv_mapping.d_leg_coeff[i] += (
                        # (self._leg_coeff_up / self._ssa_up[np.newaxis, :, :] - self._leg_coeff_down / self._ssa_down[np.newaxis, :, :]) / (x_up - x_down)[np.newaxis, :, np.newaxis]
                    # )
                    # deriv_mapping.d_leg_coeff[i] += (
                    #     (self._leg_coeff_up / (self._ssa_up * self._total_extinction_up)[np.newaxis, :, :] - self._leg_coeff_down / (self._ssa_down * self._total_extinction_down)[np.newaxis, :, :]) / (x_up - x_down)[np.newaxis, :, np.newaxis]
                    # )
                    # deriv_mapping.d_leg_coeff[:, i] += (
                        # (leg_coeff_up[:, i] / (ssa_up[np.newaxis, i] * total_extinction_up[np.newaxis, i]) - leg_coeff_down[:, i] / (ssa_down[np.newaxis, i] * total_extinction_down[np.newaxis, i])) / (x_up[i] - x_down[i])
                    # )

            else:
                deriv_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_{key}")

                d_xs_at_wavel = (
                    val.d_extinction @ wavel_interp_matrix.T
                ).flatten()

                d_n = np.trapezoid(d_xs_at_wavel * self._gaussian, atmo.model_geometry.altitudes()) / self._gaussian_od

                # """
                #####
                # deriv_mapping.d_extinction[:] += (val.d_extinction - self._optical_quants.extinction * d_n[np.newaxis, np.newaxis])
                # d_extinction_scat = val.d_ssa - self._optical_quants.ssa * self._optical_quants.extinction * d_n[np.newaxis, np.newaxis]
                deriv_mapping.d_extinction[:] += val.d_extinction
                d_extinction_scat = val.d_ssa
                #####

                # First, the optical property returns back d_scattering extinction in the d_ssa container,
                # convert this to d_ssa
                deriv_mapping.d_ssa[:] += (
                    d_extinction_scat - val.d_extinction * self._optical_quants.ssa
                ) / self._optical_quants.extinction
                deriv_mapping.d_leg_coeff[:] += val.d_leg_coeff

                #####
                if key in d_vertical_deriv_factor:
                    # Have to make some adjustments

                    # The change in extinction is adjusted
                    deriv_mapping.d_extinction[:] += (
                        self._optical_quants.extinction
                        / (alt_interp_matrix @ vertical_deriv_factor)[:, np.newaxis]
                        * (alt_interp_matrix @ d_vertical_deriv_factor[key])[
                            :, np.newaxis
                        ]
                    )
                    # deriv_mapping.d_extinction[:] += (
                    #     - self._optical_quants.extinction * d_n[np.newaxis, np.newaxis]
                    # )
                    # pass

                    # Change in single scatter albedo should be invariant whether or not we are
                    # in extinction space or number density space
                #####

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

                ############################
                # TODO: testing
                # deriv_mapping.d_extinction[:] = np.zeros_like(deriv_mapping.d_extinction)
                # deriv_mapping.d_ssa[:] = (val.d_ssa - self._optical_quants.ssa * val.d_extinction) / atmo.storage.total_extinction
                ############################

                # deriv_mapping.interpolator = -self._vertical_optical_depth * self._gaussian[np.newaxis, :] * np.trapezoid(d_xs_at_wavel * self._gaussian, atmo.model_geometry.altitudes()) / self._gaussian_od**2
                deriv_mapping.interpolator = (
                    alt_interp_matrix * number_density[np.newaxis, :]
                )
                # deriv_mapping.interpolator = (
                #     alt_interp_matrix * number_density[np.newaxis, :] * 26897416669715.535 / 2. / 1.0527448224656957  # 12774898577376.762
                # )
                deriv_mapping.interp_dim = f"{name}_altitude"
                # """

                # deriv_mapping.d_extinction[:] += d_extinction
                # deriv_mapping.d_ssa[:] += d_ssa
                # deriv_mapping.d_leg_coeff[:] += d_leg_coeff
                # deriv_mapping.scat_factor[:] += scat_factor
                # deriv_mapping.interp_dim = f"{name}_altitude"
                # deriv_mapping.interpolator = alt_interp_matrix * (-number_density * d_n)[np.newaxis, :]

        return derivs