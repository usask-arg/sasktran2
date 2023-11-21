import sasktran as sk_old
import sasktran2 as sk
import numpy as np
from tqdm import tqdm
import xarray as xr
from pathlib import Path


OUT_DIR = Path(
    "/Users/dannyz/OneDrive - University of Saskatchewan/sasktran_databases/dist/cross_sections/"
)


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
HIRES_WAVEL_GRID = np.arange(200.0, 1200.0, 0.001)

NUM_COMBINE = 100

MAX_WAVEL_RES = 0.1
SIGNIFICANT_VOD_THRESHOLD = 0.01


def typical_vod(species):
    const = sk.climatology.mipas.constituent(
        species, None, dataset="fascode", climatology="std"
    )

    # Hacks for P/T
    temperature_k = (
        sk.climatology.mipas.constituent(
            "TEM", None, dataset="fascode", climatology="std"
        ).vmr
        * 1e6
    )
    pressure_pa = (
        sk.climatology.mipas.constituent(
            "PRE", None, dataset="fascode", climatology="std"
        ).vmr
        * 1e6
        * 1e2
    )
    alts_m = (
        sk.climatology.mipas.constituent("HGT", None, "fascode", "std").vmr * 1e6 * 1e3
    )

    nden = sk.optical.pressure_temperature_to_numberdensity(pressure_pa, temperature_k)

    integrated_nden = np.trapz(nden, alts_m)

    return integrated_nden


def reduced_wv():
    all_wv = []

    i = 0
    while i + NUM_COMBINE < len(HIRES_WAVEL_GRID):
        all_wv.append(np.nanmean(HIRES_WAVEL_GRID[i : i + NUM_COMBINE]))

        i += NUM_COMBINE

    all_wv = np.array(all_wv)

    return all_wv


def reduce_xs(vod, xs):
    all_wv = []
    all_xs = []

    i = 0
    while i + NUM_COMBINE < len(HIRES_WAVEL_GRID):
        all_xs.append(np.nanmean(xs[i : i + NUM_COMBINE]))

        i += NUM_COMBINE

    all_xs = np.array(all_xs)

    return all_xs


def gen_db(species, ltol=1e-9):
    vod = typical_vod(species)

    r_wv = reduced_wv()
    reduced_xs = np.zeros((len(PRESSURE_GRID), len(TEMP_GRID), len(r_wv)))

    for idx, pres in tqdm(enumerate(PRESSURE_GRID)):
        for idy, temp in enumerate(TEMP_GRID):
            xs = (
                sk_old.HITRANChemical(species, line_tolerance=ltol)
                .calculate_cross_sections(
                    sk_old.ClimatologyUserDefined(
                        altitudes=[0, 100000],
                        values={
                            "SKCLIMATOLOGY_PRESSURE_PA": [pres, pres],
                            "SKCLIMATOLOGY_TEMPERATURE_K": [temp, temp],
                        },
                    ),
                    0,
                    0,
                    10000,
                    54372,
                    HIRES_WAVEL_GRID,
                )
                .total
                / 1e4
            )

            reduced_xs[idx, idy] = reduce_xs(vod, xs)

            # xs = sk_old.HITRANChemical(species, line_tolerance=0).calculate_cross_sections(sk_old.ClimatologyUserDefined(altitudes=[0, 100000], values={"SKCLIMATOLOGY_PRESSURE_PA": [pres, pres], "SKCLIMATOLOGY_TEMPERATURE_K": [temp, temp]}), 0, 0, 10000, 54372, HIRES_WAVEL_GRID).total / 1e4
            #
            # reduced_wv2, reduced_xs2 = reduce_xs(vod, xs)

    ds = xr.Dataset(
        {"xs": (["pressure", "temperature", "wavelength_nm"], reduced_xs)},
        coords={
            "pressure": PRESSURE_GRID,
            "temperature": TEMP_GRID,
            "wavelength_nm": r_wv,
        },
    )

    out_dir = OUT_DIR.joinpath(species.lower())
    if not out_dir.exists():
        out_dir.mkdir()

    ds.to_netcdf(out_dir.joinpath("hitran_01nm_res.nc"))


gen_db("O2")
# gen_db("H2O")
