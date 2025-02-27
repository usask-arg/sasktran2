from errno import ENEEDAUTH
import xarray as xr
from pathlib import Path
import numpy as np
import os
import sasktran2 as sk

from sasktran2.climatology.us76 import _PRESSURE, _ALTS

wvnum_range = np.arange(1e7 / 100000, 1e7 / 1000 + 0.01, 1.0)

executable = Path("/Users/djz828/dev/MT_CKD_H2O/mt_ckd_h2o_4.3_OS_X_gnu_dbl")

PRESSURE_GRID = (
    np.array(
        [
            1.03181655e03,
            9.99910226e02,
            9.72231997e02,
            9.44960640e02,
            9.18098202e02,
            8.91646636e02,
            8.65607893e02,
            8.39983721e02,
            8.14775620e02,
            7.89985094e02,
            7.65613488e02,
            7.41661804e02,
            7.18131094e02,
            6.95022257e02,
            6.72335891e02,
            6.50072549e02,
            6.28232578e02,
            6.06816233e02,
            5.85823558e02,
            5.65254408e02,
            5.45108580e02,
            5.25385676e02,
            5.06085044e02,
            4.87205987e02,
            4.68747652e02,
            4.50708890e02,
            4.33088502e02,
            4.15885138e02,
            3.99097246e02,
            3.82723078e02,
            3.66760783e02,
            3.51208311e02,
            3.36063513e02,
            3.21323986e02,
            3.06987184e02,
            2.93050454e02,
            2.79510897e02,
            2.66365512e02,
            2.53611100e02,
            2.41244361e02,
            2.29261844e02,
            2.17659797e02,
            2.06434425e02,
            1.95581772e02,
            1.85097641e02,
            1.74977731e02,
            1.65217644e02,
            1.55812775e02,
            1.46758379e02,
            1.38049500e02,
            1.29681042e02,
            1.21647854e02,
            1.13944585e02,
            1.06565734e02,
            9.95056250e01,
            9.27584798e01,
            8.63183864e01,
            8.01792901e01,
            7.43350047e01,
            6.87792048e01,
            6.35054499e01,
            5.85071744e01,
            5.37776926e01,
            4.93102089e01,
            4.50978073e01,
            4.11334673e01,
            3.74100678e01,
            3.39203935e01,
            3.06571276e01,
            2.76128649e01,
            2.47801240e01,
            2.21513493e01,
            1.97189141e01,
            1.74751173e01,
            1.54122078e01,
            1.35223946e01,
            1.17978408e01,
            1.02306666e01,
            8.81297618e00,
            7.53686238e00,
            6.39440891e00,
            5.37770790e00,
            4.47887025e00,
            3.69003722e00,
            3.00339331e00,
            2.41118145e00,
            1.90571877e00,
            1.47941371e00,
            1.12478403e00,
            8.34477191e-01,
            6.01292929e-01,
            4.18208541e-01,
            2.78406802e-01,
            1.75308171e-01,
            1.02608084e-01,
            5.43204621e-02,
        ]
    )
    * 100
)
TEMP_GRID = np.arange(190, 311, 10)

all_xs = np.zeros((len(PRESSURE_GRID), len(TEMP_GRID), len(wvnum_range)))


def gen_config(
    pressure: float,
    temperature: float,
    h2o_vmr: float,
    low_wv: float,
    high_wv: float,
    wv_spacing: float,
):
    content = f"""&mt_ckd_input
        p_atm={pressure:.2f}
        t_atm={temperature:.2f}
        h2o_frac={h2o_vmr:.2f}
        wv1={low_wv:.2f}
        wv2={high_wv:.2f}
        dwv={wv_spacing:.2f}
    /
    """

    return content


def run_mt_ckd(config: str):
    os.chdir(executable.parent / "run_example")
    with open("mt_ckd.config", "w") as f:
        f.write(config)

    r = os.system(f"../mt_ckd_h2o_4.3_OS_X_gnu_dbl < mt_ckd.config")

    ds = xr.open_dataset("mt_ckd_h2o_output.nc").load()
    os.unlink("mt_ckd.config")
    os.unlink("mt_ckd_h2o_output.nc")

    return ds


mip_clim = sk.climatology.mipas.constituent("H2O", sk.optical.HITRANAbsorber("H2O"))
mip_vmr = mip_clim._vmr
mip_alts = mip_clim._altitudes_m

for i, pres in enumerate(PRESSURE_GRID):
    print(i)
    for j, temp in enumerate(TEMP_GRID):
        alt = np.interp(pres, _PRESSURE[::-1] * 1e4, _ALTS[::-1], left=0, right=80000)
        vmr = np.interp(alt, mip_alts, mip_vmr, left=mip_vmr[0], right=mip_vmr[-1])

        cfg = gen_config(pres / 100, temp, vmr, wvnum_range[0], wvnum_range[-1], 1.0)

        ds = run_mt_ckd(cfg)

        all_xs[i, j, :] = (
            ds["self_absorption"].to_numpy() + ds["frgn_absorption"].to_numpy()
        )

        pass

ds = xr.Dataset(
    {"xs": (["pressure_pa", "temperature_k", "wavenumber_cminv"], all_xs / 1e4)},
    coords={
        "pressure_pa": PRESSURE_GRID,
        "temperature_k": TEMP_GRID,
        "wavenumber_cminv": wvnum_range,
    },
)

ds.to_netcdf("mt_ckd_h2o.nc")
