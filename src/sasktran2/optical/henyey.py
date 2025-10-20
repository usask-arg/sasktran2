from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from sasktran2._core_rust import (
    PyScatteringDatabaseDim1,
    PyScatteringDatabaseDim2,
    PyScatteringDatabaseDim3,
)
from sasktran2.atmosphere import Atmosphere
from sasktran2.optical.database import OpticalDatabase
from sasktran2.optical.quantities import OpticalQuantities

from .quantities import OpticalQuantities as RustOpticalQuantities


class HenyeyGreenstein(OpticalDatabase):
    @classmethod
    def from_parameters(
        cls,
        wavelength_nm: np.array,
        xs_total: np.array,
        ssa: np.array,
        g: np.array,
        max_num_moments: int = 128,
    ) -> HenyeyGreenstein:
        """
        Create a Henyey-Greenstein optical database from parameter arrays.  Note that
        these parameters are thus geometry independent, i.e., do not depend on altitude or location.

        Parameters
        ----------
        wavelength_nm : np.array
            Wavelength in [nm]
        xs_total : np.array
            Total extinction cross-section in [m^2]
        ssa : np.array
            Single scattering albedo (unitless)
        g : np.array
            Asymmetry parameter (unitless)
        max_num_moments : int, optional
            Maximum number of moments to use, should be set higher than the number of streams
             used in the calculation. by default 128
        """

        ds = xr.Dataset(
            {
                "xs_total": (("wavelength_nm",), xs_total),
                "ssa": (("wavelength_nm",), ssa),
                "asymmetry_parameter": (("wavelength_nm",), g),
            },
            coords={"wavelength_nm": wavelength_nm},
        )

        return cls(db=ds, max_num_moments=max_num_moments)

    def __init__(
        self,
        db_filepath: Path | None = None,
        db: xr.Dataset | None = None,
        max_num_moments: int = 128,
    ) -> None:
        """
        An optical property that uses a Henyey-Greenstein phase function. I.e., the phase function is defined
        from a single asymmetry parameter g.

        The input should either be a path to a netCDF file containing the database, or an xarray Dataset.

        The dataset must contain the following variables:

        - xs_total: Total extinction cross-section in [m^2]
        - ssa: Single scattering albedo (unitless)
        - asymmetry_parameter: Asymmetry parameter (unitless)

        with Coordinates:
        - wavelength_nm: Wavelength in [nm]

        The variables may have additional dimensions corresponding to different parameters (e.g., temperature).

        Parameters
        ----------
        db_filepath : Path | None, optional
            Path to a netCDF file containing the database, by default None
        db : xr.Dataset | None, optional
            An xarray Dataset containing the database, by default None
        """
        super().__init__(db_filepath, db)

        self._validate_db()

        # Reorient the dimensions
        dims = list(self._database["xs_total"].isel(wavelength_nm=0).dims)
        db = self._database.transpose(*dims, "wavelength_nm", ...)

        # construct internal object
        xs = db["xs_total"].to_numpy()
        ssa = db["ssa"].to_numpy()
        g = db["asymmetry_parameter"].to_numpy()

        wvnum = 1e7 / db["wavelength_nm"].to_numpy()
        sidx = np.argsort(wvnum)

        if len(xs.shape) == 1:
            self._db = PyScatteringDatabaseDim1.from_asymmetry_parameter(
                xs[sidx], ssa[sidx], g[sidx], max_num_moments, wvnum[sidx]
            )
        elif len(xs.shape) == 2:
            param_names = list(db["xs_total"].dims)[:-1]
            param0 = db[param_names[0]].to_numpy()
            self._db = PyScatteringDatabaseDim2.from_asymmetry_parameter(
                xs[:, sidx],
                ssa[:, sidx],
                g[:, sidx],
                max_num_moments,
                wvnum[sidx],
                np.atleast_1d(param0).astype(np.float64),
                param_names,
            )
        elif len(xs.shape) == 3:
            param_names = list(db["xs_total"].dims)[:-1]
            param0 = db[param_names[0]].to_numpy()
            param1 = db[param_names[1]].to_numpy()
            self._db = PyScatteringDatabaseDim3.from_asymmetry_parameter(
                xs[:, :, sidx],
                ssa[:, :, sidx],
                g[:, :, sidx],
                max_num_moments,
                wvnum[sidx],
                np.atleast_1d(param0).astype(np.float64),
                np.atleast_1d(param1).astype(np.float64),
                param_names,
            )

    def _validate_db(self):
        pass

    def cross_sections(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs
    ) -> OpticalQuantities:
        return self._db.cross_sections(
            np.atleast_1d(wavelengths_nm).astype(float),
            np.atleast_1d(altitudes_m).astype(float),
            **kwargs,
        )

    def atmosphere_quantities(self, atmo: Atmosphere, **kwargs) -> OpticalQuantities:
        return RustOpticalQuantities(self._db.atmosphere_quantities(atmo, **kwargs))

    def optical_derivatives(self, atmo: Atmosphere, **kwargs) -> dict:
        result = self._db.optical_derivatives(atmo, **kwargs)

        return {k: RustOpticalQuantities(v) for k, v in result.items()}

    def cross_section_derivatives(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs
    ) -> dict:
        result = self._db.cross_section_derivatives(
            np.atleast_1d(wavelengths_nm).astype(float),
            np.atleast_1d(altitudes_m).astype(float),
            **kwargs,
        )
        return {k: v.extinction.flatten() for k, v in result.items()}
