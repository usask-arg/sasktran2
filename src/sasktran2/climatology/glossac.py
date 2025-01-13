from __future__ import annotations

import numpy as np
import xarray as xr

import sasktran2 as sk


def load_glossac_raw_data(version: str = "2.21") -> xr.Dataset:
    """
    Load's in the GloSSAC raw data for a specific version.  First we check if the file exists within the database directory,
    if it doesn't we load it in through opendap.

    Parameters
    ----------
    version : str, optional
        GloSSAC File version to look for, by default "2.21"

    Returns
    -------
    xr.Dataset
        GloSSAC Dataset
    """
    db_root = sk.appconfig.database_root()

    if db_root is None:
        return xr.open_dataset(
            f"https://opendap.larc.nasa.gov/opendap/GloSSAC/GloSSAC_{version}/GloSSAC_V{version}.nc"
        )

    glossac_file = sk.appconfig.database_root().joinpath(
        f"climatology/glossac/GloSSAC_V{version}.nc"
    )

    if glossac_file.exists():
        return xr.open_dataset(glossac_file)

    return xr.open_dataset(
        f"https://opendap.larc.nasa.gov/opendap/GloSSAC/GloSSAC_{version}/GloSSAC_V{version}.nc"
    )


def stratospheric_background(
    month: int, lat: float, alts: np.array, wavelength_nm: float, version: str = "2.21"
):
    ds = load_glossac_raw_data(version)

    wavel_idx = np.argmin(np.abs(ds["wavelengths_glossac"].to_numpy() - wavelength_nm))
    lat_idx = np.argmin(np.abs(ds["lat"].to_numpy() - lat))

    if np.abs(ds["wavelengths_glossac"].to_numpy()[wavel_idx] - wavelength_nm) > 1e-4:
        msg = f"Could not find wavelength {wavelength_nm} in GloSSAC data.  Valid wavelengths are {ds['wavelengths_glossac'].to_numpy()}"
        raise ValueError(msg)

    data = ds.isel(wavelengths_glossac=wavel_idx)

    background = (
        data["Stratospheric_Background"]
        .sel(month=month)
        .isel(lat=lat_idx)
        .interp(alt=alts / 1000, method="linear")
        .to_numpy()
    )

    background[np.isnan(background)] = 0

    return background / 1e3  # Convert to /m
