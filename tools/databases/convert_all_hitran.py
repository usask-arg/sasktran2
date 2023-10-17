import pandas as pd
from pathlib import Path
import numpy as np
import xarray as xr
import sasktran2 as sk

HITRAN_FOLDER = Path(
    "/Users/dannyz/OneDrive - University of Saskatchewan/sasktran_databases/input/hitran_uv/"
)
OUT_DIR = Path(
    "/Users/dannyz/OneDrive - University of Saskatchewan/sasktran_databases/dist/cross_sections/"
)


for folder in HITRAN_FOLDER.iterdir():
    name = folder.stem

    all_wv = []
    all_xs = []
    all_T = []
    for file in folder.glob("*.xsc"):
        header = pd.read_csv(file, header=None, nrows=1, sep="\s+")

        wv_num_low = header.values[0, 1]
        wv_num_high = header.values[0, 2]

        npoints = header.values[0, 3]
        temperature = header.values[0, 4]

        wv_num = 1e7 / np.linspace(wv_num_low, wv_num_high, npoints)

        if name == "no2":
            wv_num = sk.optical.air_wavelength_to_vacuum_wavelength(wv_num)

        xs_cm2 = pd.read_csv(file, header=None, skiprows=1, sep="\s+").values.flatten()[
            :npoints
        ]

        sort_idx = np.argsort(wv_num)

        all_wv.append(wv_num[sort_idx])
        all_xs.append(xs_cm2[sort_idx])
        all_T.append(temperature)

    if len(all_wv) == 0:
        continue

    combined_wv = np.sort(np.unique(np.hstack(all_wv)))

    idx_min_T = np.argmin(all_T)
    idx_max_T = np.argmax(all_T)

    all_T.append(0)
    all_xs.append(all_xs[idx_min_T])
    all_wv.append(all_wv[idx_min_T])
    all_T.append(1000)
    all_xs.append(all_xs[idx_max_T])
    all_wv.append(all_wv[idx_max_T])

    T_sort = np.argsort(all_T)

    combined_xs = np.zeros((len(all_T), len(combined_wv)))

    for i in range(len(all_T)):
        combined_xs[i] = np.interp(
            combined_wv, all_wv[i], all_xs[i], left=np.nan, right=np.nan
        )

    all_T = np.array(all_T)[T_sort]
    combined_xs = combined_xs[T_sort]

    # now we have to go through and remove all of the NaN values by extending the other temperature measurements
    for i in range(len(combined_wv)):
        isna = np.isnan(combined_xs[:, i])
        if sum(isna) > 0:
            combined_xs[isna, i] = np.interp(
                all_T[isna],
                all_T[~isna],
                combined_xs[~isna, i],
                left=combined_xs[~isna, i][0],
                right=combined_xs[~isna, i][-1],
            )

    ds = xr.Dataset(
        {"xs": (["temperature", "wavelength_nm"], combined_xs / 1e4)},
        coords={"temperature": all_T, "wavelength_nm": combined_wv},
    )

    species_dir = OUT_DIR.joinpath(name)

    if not species_dir.exists():
        species_dir.mkdir()

    ds.to_netcdf(species_dir.joinpath("hitran2022.nc").as_posix())
