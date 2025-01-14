from __future__ import annotations

import numpy as np

from sasktran2.atmosphere import Atmosphere

from ..base import Constituent
from . import LambertianStokes_1, LambertianStokes_3, WavelengthInterpolatorMixin


class LambertianSurface(Constituent, WavelengthInterpolatorMixin):
    def __init__(
        self,
        albedo: np.array,
        wavelengths_nm: np.array = None,
        wavenumbers_cminv: np.array = None,
        out_of_bounds_mode="zero",
    ) -> None:
        """
        A Lambertian surface that is defined by albedo at discrete grid points.

        This can either operate in a "scalar" mode where the albedo is constant in wavelength,
        a "native" mode where the albedo is defined on the same grid as the atmosphere, or
        an "interpolated" mode where the albedo is interpolated either in wavenumber or wavelength

        Parameters
        ----------
        albedo : np.array
            Surface albedo.  Can be a scalar to indicate it is constant in wavelength.  If set to an
            array it must either match the atmosphere wavelength grid dimension, or one of
            wavelengths_nm or wavenumbers_cminv must be set.
        wavelengths_nm : np.array, optional
            Wavelengths in [nm] that the albedo is specified at, by default None
        wavenumbers_cminv : np.array, optional
            Wavenumbers in inverse cm that the albedo is specified at, by default None
        out_of_bounds_mode : str, optional
            One of ["extend" or "zero"], "extend" will extend the last/first value if we are
            interpolating outside the grid. "zero" will set the albedo to 0 outside of the
            grid boundaries, by default "zero"
        """
        Constituent.__init__(self)
        WavelengthInterpolatorMixin.__init__(
            self,
            wavelengths_nm=wavelengths_nm,
            wavenumbers_cminv=wavenumbers_cminv,
            out_of_bounds_mode=out_of_bounds_mode,
            param_length=len(np.atleast_1d(albedo)),
        )
        self._albedo = np.atleast_1d(albedo)

    @property
    def albedo(self) -> np.array:
        return self._albedo

    @albedo.setter
    def albedo(self, albedo: np.array):
        self._albedo = np.atleast_1d(albedo)

    def add_to_atmosphere(self, atmo: Atmosphere):
        atmo.surface.brdf = (
            LambertianStokes_1() if atmo.nstokes == 1 else LambertianStokes_3()
        )

        interp_matrix = self._interpolating_matrix(atmo)

        atmo.surface.brdf_args[0, :] = interp_matrix @ self._albedo
        atmo.surface.d_brdf_args[0][0, :] = 1

    def register_derivative(self, atmo: Atmosphere, name: str):
        # Start by constructing the interpolation matrix
        interp_matrix = self._interpolating_matrix(atmo)

        derivs = {}

        deriv_mapping = atmo.surface.get_derivative_mapping(f"wf_{name}_albedo")
        deriv_mapping.d_brdf[:] = np.ones((atmo.num_wavel, 1))
        deriv_mapping.interpolator = interp_matrix
        deriv_mapping.interp_dim = f"{name}_wavelength"

        return derivs
