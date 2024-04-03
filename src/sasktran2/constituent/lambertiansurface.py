import numpy as np

from sasktran2 import Atmosphere, LambertianStokes_1, LambertianStokes_3
from sasktran2.atmosphere import (
    NativeGridDerivative,
    SurfaceDerivativeMapping,
)
from sasktran2.util.interpolation import linear_interpolating_matrix

from .base import Constituent


class LambertianSurface(Constituent):
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
        super().__init__()

        self._out_of_bounds_mode = out_of_bounds_mode
        self._albedo = np.atleast_1d(albedo)

        if wavelengths_nm is not None:
            self._x = np.atleast_1d(wavelengths_nm)
            self._interp_var = "wavelengths_nm"
        elif wavenumbers_cminv is not None:
            self._x = np.atleast_1d(wavenumbers_cminv)
            self._interp_var = "wavenumbers_cminv"
        else:
            if len(self._albedo) == 1:
                self._interp_var = "constant"
            else:
                self._interp_var = "native"

    @property
    def albedo(self) -> np.array:
        return self._albedo

    @albedo.setter
    def albedo(self, albedo: np.array):
        self._albedo = np.atleast_1d(albedo)

    def _interpolating_matrix(self, atmo: Atmosphere):
        if self._interp_var == "constant":
            # Don't have to interpolate, just assign the constant value
            return np.ones(atmo.num_wavel).reshape(-1, 1)
        if self._interp_var == "native":
            # Also can just assign the user value, but first make sure that it is the correct length
            if len(self._albedo) != atmo.num_wavel:
                msg = "The number of albedo values must match the number of wavelengths in the atmosphere"
                raise ValueError(msg)
            return np.eye(atmo.num_wavel)
        # Now we have to interpolate
        grid_x = getattr(atmo, self._interp_var)

        if grid_x is None:
            msg = f"The atmosphere does not have {self._interp_var} defined, cannot interpolate the user albedo"
            raise ValueError(msg)
        return linear_interpolating_matrix(self._x, grid_x, self._out_of_bounds_mode)

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

        derivs["albedo"] = SurfaceDerivativeMapping(
            NativeGridDerivative(d_albedo=np.ones(atmo.num_wavel)),
            interpolating_matrix=interp_matrix,
            interp_dim="wavelength",
            result_dim=f"{name}_wavelength",
        )

        return derivs
