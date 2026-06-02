from __future__ import annotations

from typing import Any

import numpy as np

from sasktran2._core_rust import PyNumberDensityScatterer
from sasktran2.atmosphere import Atmosphere
from sasktran2.optical.base import OpticalProperty

from .base import Constituent


class NumberDensityScatterer(Constituent):
    _constituent: PyNumberDensityScatterer

    def __init__(
        self,
        optical_property: OpticalProperty,
        altitudes_m: np.array,
        number_density: np.array,
        out_of_bounds_mode: str = "zero",
        **kwargs,
    ) -> None:
        """
        A scattering constituent that is defined by a number density on an altitude grid and an optical property

        Parameters
        ----------
        optical_property : OpticalProperty
            The optical property defining the scattering information
        altitudes_m : np.array
            The altitude grid in [m]
        number_density : np.array
            Number density in [m^-3]
        out_of_bounds_mode : str, optional
            Interpolation mode outside of the boundaries, "extend" and "zero" are supported, by default "zero"
        kwargs : dict
            Additional arguments to pass to the optical property.
        """
        super().__init__()

        if number_density is None:
            number_density = np.zeros_like(altitudes_m, dtype=np.float64)

        self._out_of_bounds_mode = out_of_bounds_mode
        self._altitudes_m = altitudes_m
        self._optical_property = optical_property

        # Extra factor to apply to the vertical derivatives, used by the derived Extinction class
        self._vertical_deriv_factor = np.ones_like(number_density)

        # Optical derivatives can also have derivatives to this factor
        self._d_vertical_deriv_factor = {}

        self._wf_name = "number_density"

        self._kwargs = kwargs

        self._constituent = PyNumberDensityScatterer(
            optical_property,
            altitudes_m.astype(np.float64),
            number_density.astype(np.float64),
            out_of_bounds_mode,
            **{k: np.asarray(v, dtype=np.float64) for k, v in kwargs.items()},
        )
        self._sync_constituent_state()

    def _sync_constituent_state(self):
        vertical_deriv_factor = np.asarray(
            self._vertical_deriv_factor, dtype=np.float64
        )
        if vertical_deriv_factor.size == 1:
            vertical_deriv_factor = np.full_like(
                self._altitudes_m, vertical_deriv_factor.item(), dtype=np.float64
            )

        self._constituent.aux_inputs = {
            k: np.asarray(v, dtype=np.float64) for k, v in self._kwargs.items()
        }
        self._constituent.vertical_deriv_factor = vertical_deriv_factor
        self._constituent.d_vertical_deriv_factor = {
            k: (
                np.full_like(self._altitudes_m, np.asarray(v, dtype=np.float64).item())
                if np.asarray(v, dtype=np.float64).size == 1
                else np.asarray(v, dtype=np.float64)
            )
            for k, v in self._d_vertical_deriv_factor.items()
        }
        self._constituent.wf_name = self._wf_name

    def __getattr__(self, __name: str) -> Any:
        if __name in self.__dict__.get("_kwargs", {}):
            return self._kwargs[__name]
        return None

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name in self.__dict__.get("_kwargs", {}):
            self._kwargs[__name] = __value
            self._sync_constituent_state()
        else:
            super().__setattr__(__name, __value)

    @property
    def number_density(self):
        return self._constituent.number_density

    @number_density.setter
    def number_density(self, number_density: np.array):
        self._constituent.number_density = np.asarray(number_density, dtype=np.float64)

    @property
    def altitudes_m(self):
        return self._constituent.altitudes_m

    def add_to_atmosphere(self, atmo: Atmosphere):
        self._sync_constituent_state()
        self._constituent.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: Atmosphere, name: str):
        self._sync_constituent_state()
        self._constituent.register_derivative(atmo, name)


class ExtinctionScatterer(NumberDensityScatterer):
    def __init__(
        self,
        optical_property: OpticalProperty,
        altitudes_m: np.array,
        extinction_per_m: np.array,
        extinction_wavelength_nm: float,
        out_of_bounds_mode: str = "zero",
        **kwargs,
    ) -> None:
        """
        A scattering constituent that is defined by a number density on an altitude grid and an optical property

        Parameters
        ----------
        optical_property : OpticalProperty
            The optical property defining the scattering information
        altitudes_m : np.array
            The altitude grid in [m]
        extinction_per_m : np.array
            Extinction in [m^-1]
        extinction_wavelength_nm : float
            Wavelength that the extinction profile is specified at
        out_of_bounds_mode : str, optional
            Interpolation mode outside of the boundaries, "extend" and "zero" are supported, by default "zero"
        kwargs : dict
            Additional arguments passed to the optical property
        """
        self._extinction_per_m = extinction_per_m
        self._extinction_wavelength_nm = extinction_wavelength_nm

        super().__init__(
            optical_property, altitudes_m, None, out_of_bounds_mode, **kwargs
        )
        self._extinction_to_numden_factors = None
        self._update_numberdensity()
        self._wf_name = "extinction"

    def _update_numberdensity(self):
        self._extinction_to_numden_factors = self._optical_property.cross_sections(
            np.array([self._extinction_wavelength_nm]),
            altitudes_m=self._altitudes_m,
            **self._kwargs,
        ).extinction.flatten()
        self._vertical_deriv_factor = 1 / self._extinction_to_numden_factors

        self.number_density = (
            self._extinction_per_m / self._extinction_to_numden_factors
        )

        self._d_vertical_deriv_factor = (
            self._optical_property.cross_section_derivatives(
                np.array([self._extinction_wavelength_nm]),
                altitudes_m=self._altitudes_m,
                **self._kwargs,
            )
        )

        for _, val in self._d_vertical_deriv_factor.items():
            # convert from derivative of x to derivative of 1/x
            val *= -1 * self._vertical_deriv_factor**2  # noqa: PLW2901

    @property
    def extinction_per_m(self):
        return self._extinction_per_m

    @extinction_per_m.setter
    def extinction_per_m(self, extinction: np.array):
        self._extinction_per_m = extinction

    def add_to_atmosphere(self, atmo: Atmosphere):
        self._update_numberdensity()
        super().add_to_atmosphere(atmo)
