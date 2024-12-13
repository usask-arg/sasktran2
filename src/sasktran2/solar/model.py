from __future__ import annotations

import numpy as np
import xarray as xr

import sasktran2 as sk


class SolarModel:
    def __init__(self, source="solar_irradiance_hsrs_2022_11_30_extended"):
        """

        Parameters
        ----------
        source: str, optional
        """

        self._ds = sk.database.StandardDatabase().path(f"solar/{source}.nc")

        self._ds = xr.open_dataset(self._ds)

    def irradiance(self, wavelengths, solardistance=None) -> np.array:  # noqa: ARG002
        return (
            self._ds["irradiance"]
            .swap_dims({"spectral": "wavelength"})
            .interp(wavelength=wavelengths)
            .to_numpy()
        )

        # TODO: Correct for solar distance
        # TODO: Implement doppler shift
