from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import sasktran2 as sk

DBM_DIR = "/Users/dannyz/OneDrive - University of Saskatchewan/sasktran_databases/input/cross_sections/o3/dbm/"
OUT_DIR = "/Users/dannyz/OneDrive - University of Saskatchewan/sasktran_databases/dist/cross_sections/o3/"


all_wv = []
all_xs = []
all_T = []
for dbm_file in Path(DBM_DIR).iterdir():
    data = pd.read_csv(
        dbm_file,
        header=3,
        converters={
            0: lambda x: (
                float(x.replace("{", "").replace("}", "").replace(" ", ""))
                if x is not None
                else 0
            ),
            1: lambda x: (
                float(x.replace("{", "").replace("}", "").replace(" ", ""))
                if x is not None
                else 0
            ),
        },
        skipfooter=1,
    )

    wavelength_nm = data.values[:, 0].astype("float")
    xs_cm2 = data.values[:, 1].astype("float")

    temp_k = float(dbm_file.stem[7:10])

    all_wv.append(wavelength_nm)
    all_xs.append(xs_cm2)
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
    coords={
        "temperature": all_T,
        "wavelength_nm": sk.optical.air_wavelength_to_vacuum_wavelength(combined_wv),
    },
)

ds.to_netcdf(Path(OUT_DIR).joinpath("dbm.nc").as_posix())
