from __future__ import annotations

import numpy as np
import xarray as xr
from scipy.integrate import cumulative_trapezoid

import sasktran2 as sk


class SolarModel:
    def __init__(
        self,
        source="solar_irradiance_hsrs_2022_11_30_extended",
        ds=None,
        mode="sample",
        resolution=None,
        resolution_in_wavelength=True,
    ):
        """
        Loads in a solar model for use inside SASKTRAN2.  By default the HSRS extended solar model is loaded in
        that has a useful wavelength range covering 115 nm to 200 microns and intgrates to 1362.8 W/m^2.

        Coddington, O. M., Richard, E. C., Harber, D., Pilewskie, P., Woods, T. N., Snow, M., et al. (2023). Version 2 of the TSIS-1 Hybrid Solar Reference Spectrum and Extension to the Full Spectrum. Earth and Space Science, 10, e2022EA002637. https://doi.org/10.1029/2022EA002637

        Optionally a manually specified dataset can be used, which must contain the variables "wavelength" and "irradiance". It is assumed
        that "wavelength" is in units of nm and "irradiance" has spectral units of /nm.

        Parameters
        ----------
        source: str, optional
            sasktran2 database source identifier.  Currently only "solar_irradiance_hsrs_2022_11_30_extended" is supported and is the default
        ds: xr.Dataset, optional
            If set, then the dataset is used instead of the source.  The dataset must contain the variables "wavelength" and "irradiance"
            default is None
        mode: str, optional
            The mode of the solar model.  Options are "integrate", "average", "sample".  Default is "sample".  If "integrate" is selected
            then the solar spectrum is integrated over the wavelength range of interest.  If "average" is selected then the solar spectrum
            is averaged over the wavelength range of interest.  If "sample" is selected then the solar spectrum is sampled at the wavelengths
            of interest.  The "integrate" and "average" modes are useful for calculating the solar irradiance at a specific wavelength range
            and the "sample" mode is useful for calculating the solar irradiance at specific wavelengths.
        resolution: float, optional
            The resolution to integrate the solar spectrum at.  If None, then the resolution is determined by the wavelength grid.  Default is None
        resolution_in_wavelength: bool, optional
            If True, then the resolution is in wavelength units.  If False, then the resolution is in wavenumber units.  Default is True
        """
        if ds is not None:
            self._ds = ds
        else:
            self._ds = sk.database.StandardDatabase().path(f"solar/{source}.nc")
            self._ds = xr.open_dataset(self._ds)

        self._wv = self._ds["wavelength"].to_numpy()
        self._irrad = self._ds["irradiance"].to_numpy()
        self._ds.close()

        self._mode = mode.lower()

        if self._mode == "integrate" or self._mode == "average":
            self._integral = cumulative_trapezoid(
                self._ds["irradiance"].values, self._ds["wavelength"].values
            )
        elif self._mode == "sample":
            self._integral = None
        else:
            msg = "Invalid mode"
            raise ValueError(msg)

        self._resolution = resolution
        self._resolution_in_wavelength = resolution_in_wavelength

    def irradiance(
        self, wavelengths: np.ndarray, solardistance: float | None = None
    ) -> np.array:
        """
        Calculates the solar irradiance at the specified wavelengths in nm.  Default units
        are in W/m^2/nm, but can be different if a manual database is used instead.  If the mode is set to "integrate"
        then the units will be W/m^2.

        The solar distance factor can be used to scale the solar irradiance to a different distance.  The default is None
        indicating that the solar irradiance is at 1 AU.

        Parameters
        ----------
        wavelengths : np.ndarray
            Wavelengths in [nm] to calculate the solar irradiance at
        solardistance : float, optional
            Solar distance to scale by in AU., by default None

        Returns
        -------
        np.array
            The solar irradiance at the specified wavelengths
        """
        solar_distance_factor = (
            1 / (solardistance**2) if solardistance is not None else 1
        )

        if self._mode == "sample":
            return np.interp(wavelengths, self._wv, self._irrad) * solar_distance_factor
        # Have to calculate the wavelength intervals
        # and integrate the spectrum over the intervals

        if self._resolution is not None:
            if self._resolution_in_wavelength:
                left_intervals = wavelengths - self._resolution / 2
                right_intervals = wavelengths + self._resolution / 2
            else:
                left_intervals = 1e7 / (1e7 / wavelengths + self._resolution / 2)
                right_intervals = 1e7 / (1e7 / wavelengths - self._resolution / 2)
        else:
            left_intervals = np.concat(
                ([wavelengths[0]], np.diff(wavelengths) / 2 + wavelengths[:-1])
            )
            right_intervals = np.concat(
                (np.diff(wavelengths) / 2 + wavelengths[:-1], [wavelengths[-1]])
            )

        lidx = np.searchsorted(self._ds["wavelength"].values, left_intervals) - 1
        ridx = np.searchsorted(self._ds["wavelength"].values, right_intervals) - 1

        left_vals = np.interp(left_intervals, self._wv, self._irrad)
        right_vals = np.interp(right_intervals, self._wv, self._irrad)

        # Cumlative integral to left side
        I_left = self._integral[lidx] + 0.5 * (self._irrad[lidx] + left_vals) * (
            left_intervals - self._wv[lidx]
        )
        I_right = self._integral[ridx] + 0.5 * (self._irrad[ridx] + right_vals) * (
            right_intervals - self._wv[ridx]
        )

        if self._mode == "integrate":
            return (I_right - I_left) * solar_distance_factor
        return (
            (I_right - I_left)
            / (right_intervals - left_intervals)
            * solar_distance_factor
        )

        # TODO: Implement doppler shift?
