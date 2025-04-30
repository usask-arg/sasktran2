from __future__ import annotations

import numpy as np

from sasktran2.atmosphere import Atmosphere
from sasktran2.util.interpolation import linear_interpolating_matrix

from ..._core_rust import PyKokhanovsky, PyLambertian, PyMODIS  # noqa: F401


class WavelengthInterpolatorMixin:
    def __init__(
        self,
        wavelengths_nm: np.array = None,
        wavenumbers_cminv: np.array = None,
        out_of_bounds_mode="zero",
        param_length: int = 1,
    ):
        """
        A MixIn class that provides interpolation functionality over wavelength or wavenumber

        Parameters
        ----------
        wavelengths_nm : np.array, optional
            Wavelengths in [nm] that the parameters are specified at, by default None
        wavenumbers_cminv : np.array, optional
            Wavenumbers in inverse cm that the parameters are specified at, by default None
        out_of_bounds_mode : str, optional
            One of ["extend" or "zero"], "extend" will extend the last/first value if we are
            interpolating outside the grid. "zero" will set the parameters to 0 outside of the
            grid boundaries, by default "zero"
        """
        self._out_of_bounds_mode = out_of_bounds_mode

        if wavelengths_nm is not None:
            self._x = np.atleast_1d(wavelengths_nm)
            self._interp_var = "wavelengths_nm"
        elif wavenumbers_cminv is not None:
            self._x = np.atleast_1d(wavenumbers_cminv)
            self._interp_var = "wavenumbers_cminv"
        else:
            if param_length == 1:
                self._interp_var = "constant"
            else:
                self._interp_var = "native"

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
