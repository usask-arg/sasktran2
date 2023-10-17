from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import sasktran2 as sk

V_DIR = "/Users/dannyz/OneDrive - University of Saskatchewan/sasktran_databases/input/cross_sections/no2/vandaele/"
OUT_DIR = "/Users/dannyz/OneDrive - University of Saskatchewan/sasktran_databases/dist/cross_sections/no2/"


all_wv = []
all_xs = []
all_T = []
for no2_file in Path(V_DIR).iterdir():
    data = pd.read_csv(no2_file.as_posix(), header=0, sep="\s+")

    wavelength_nm = 1e7 / data.values[:, 1].astype("float")
    xs_cm2 = data.values[:, 2].astype("float")

    sort_idx = np.argsort(wavelength_nm)

    if "c" in no2_file.stem:
        temp_k = 220
    else:
        temp_k = 294

    all_wv.append(wavelength_nm[sort_idx])
    all_xs.append(xs_cm2[sort_idx])
    all_T.append(temp_k)

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

ds.to_netcdf(Path(OUT_DIR).joinpath("vandaele.nc").as_posix())
