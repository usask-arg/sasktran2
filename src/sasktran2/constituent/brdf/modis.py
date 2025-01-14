from __future__ import annotations

import numpy as np

from sasktran2.atmosphere import Atmosphere

from ..base import Constituent
from . import (
    MODISStokes_1,
    MODISStokes_3,
    WavelengthInterpolatorMixin,
)


class MODIS(Constituent, WavelengthInterpolatorMixin):
    def __init__(
        self,
        isotropic: np.array,
        volumetric: np.array = 0.0,
        geometric: np.array = 0.0,
        wavelengths_nm: np.array = None,
        out_of_bounds_mode="zero",
    ) -> None:
        """
        Parameters
        ----------
        isotropic : np.array
            Isotropic component (contribution weight of Lambertian).
        volumetric : np.array, optional
            Volumetric component (contribution weight of RossThick kernel), by default 0.
        geometric : np.array, optional
            Geometric component (contribution weight of LiSparse-R kernel), by default 0.
        wavelengths_nm : np.array, optional
            Wavelengths in [nm] that the parameters isotropic, volumetric, geometric is specified at, by default None indicating that these parameters are scalar.
        out_of_bounds_mode : str, optional
            One of ["extend" or "zero"], "extend" will extend the last/first value if we are
            interpolating outside the grid. "zero" will set the albedo to 0 outside of the
            grid boundaries, by default "zero"
        """
        Constituent.__init__(self)
        WavelengthInterpolatorMixin.__init__(
            self,
            wavelengths_nm=wavelengths_nm,
            wavenumbers_cminv=None,
            out_of_bounds_mode=out_of_bounds_mode,
            param_length=len(np.atleast_1d(isotropic)),
        )
        self._iso = np.atleast_1d(isotropic)
        self._vol = np.atleast_1d(volumetric)
        self._geo = np.atleast_1d(geometric)

    @property
    def isotropic(self) -> np.array:
        return self._iso

    @isotropic.setter
    def isotropic(self, iso: np.array):
        self._iso = np.atleast_1d(iso)

    @property
    def volumetric(self) -> np.array:
        return self._vol

    @volumetric.setter
    def volumetric(self, vol: np.array):
        self._vol = np.atleast_1d(vol)

    @property
    def geometric(self) -> np.array:
        return self._geo

    @geometric.setter
    def geometric(self, geo: np.array):
        self._geo = np.atleast_1d(geo)

    def add_to_atmosphere(self, atmo: Atmosphere):
        if atmo.wavelengths_nm is None:
            msg = "Atmosphere must have wavelengths defined before using MODIS"
            raise ValueError(msg)

        atmo.surface.brdf = MODISStokes_1() if atmo.nstokes == 1 else MODISStokes_3()

        interp_matrix = self._interpolating_matrix(atmo)
        atmo.surface.brdf_args[0, :] = interp_matrix @ self._iso
        atmo.surface.brdf_args[1, :] = interp_matrix @ self._vol
        atmo.surface.brdf_args[2, :] = interp_matrix @ self._geo

    def register_derivative(self, atmo: Atmosphere, name: str):
        # TODO update once C++ derivatives are implemented
        # Start by constructing the interpolation matrix
        interp_matrix = self._interpolating_matrix(atmo)

        derivs = {}

        iso_deriv = np.zeros((atmo.num_wavel, 3))
        iso_deriv[:, 0] = 1

        deriv_mapping = atmo.surface.get_derivative_mapping(f"wf_{name}_isotropic")
        deriv_mapping.d_brdf[:] = iso_deriv
        deriv_mapping.interpolator = interp_matrix
        deriv_mapping.interp_dim = f"{name}_wavelength"

        vol_deriv = np.zeros((atmo.num_wavel, 3))
        vol_deriv[:, 1] = 1

        deriv_mapping = atmo.surface.get_derivative_mapping(f"wf_{name}_volumetric")
        deriv_mapping.d_brdf[:] = vol_deriv
        deriv_mapping.interpolator = interp_matrix
        deriv_mapping.interp_dim = f"{name}_wavelength"

        geo_deriv = np.zeros((atmo.num_wavel, 3))
        geo_deriv[:, 2] = 1

        deriv_mapping = atmo.surface.get_derivative_mapping(f"wf_{name}_geometric")
        deriv_mapping.d_brdf[:] = geo_deriv
        deriv_mapping.interpolator = interp_matrix
        deriv_mapping.interp_dim = f"{name}_wavelength"

        return derivs
