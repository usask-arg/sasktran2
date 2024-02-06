import numpy as np
import xarray as xr
import sasktran2 as sk
import matplotlib.pyplot as plt


# https://hitran.org/data/CIA/O2-O2_2018b.cia

fname_in = "tools/databases/O2-O2_2018b.cia"
fname_out = R"C:\Users\luf542\AppData\Local\usask-arg\sasktran2\database\cross_sections\o2o2\hitran_cia.nc"

# header line numbers
groups = [
    dict(header_line_number=ln)
    for ln in [
        1,
        4004,
        8007,
        12010,
        16013,
        20016,
        24019,
        28022,
        32025,
        36028,
        40031,
        44034,
        48037,
        52040,
        56043,
        60046,
        64301,
        64808,
        65526,
        68099,
        68793,
        70169,
        71545,
        72919,
        74295,
        87387,
        100366,
        113414,
        126550,
    ]
]

# read header data
with open(fname_in, "r") as f:
    lines = f.readlines()
    for g in groups:
        line = lines[g["header_line_number"] - 1]
        tokens = line.split()
        g["size"] = int(tokens[3])
        g["temperature"] = float(tokens[4])
        g["wavenumber_range"] = (float(tokens[1]), float(tokens[2]))

# read wavenumber and cross section data
for g in groups:
    data = np.loadtxt(
        fname=fname_in, skiprows=g["header_line_number"], max_rows=g["size"]
    )
    g["wavenumber"] = data[:, 0]  # cm5
    g["cross_section"] = data[:, 1]

# construct a wavenumber grid that spans all the given ranges

# 0-14: 0.2 spacing (but with 0.1 spacing for 2 data points?)
# 15: ~0.245 spacing
# 16: 1.0 spacing
# 17: 1.0 spacing
# 18: ~0.482 spacing
# 19: 1.0 spacing
# 20-23: 1.0 spacing
# 24-28: 1.0 spacing

# 0-14 have identical grids
# 20-23 have slightly differing bounds, but this is not a problem because all points are integer wavenumbers (just span the integers from overall min to overall max)
# 24-28 have slightly differing bounds, but ^
# some overlap between 20-23 and 24-28, but ^


# pad regions with an extra data point so we return 0 between bands instead of trying to interpolate between
def pad(v: np.ndarray):
    return np.unique([np.min(v) - 0.01, *v, np.max(v) + 0.01])


wn1 = pad(groups[0]["wavenumber"])
wn2 = pad(groups[15]["wavenumber"])
wn3 = pad(groups[16]["wavenumber"])
wn4 = pad(groups[17]["wavenumber"])
wn5 = pad(groups[18]["wavenumber"])
wn6 = pad(groups[19]["wavenumber"])
wn7_min = np.min([np.min(groups[i]["wavenumber"]) for i in range(21, 29)])
wn7_max = np.max([np.max(groups[i]["wavenumber"]) for i in range(21, 29)])
wn7 = pad(np.arange(wn7_min, wn7_max + 0.1, 1.0))
all_wn = np.concatenate((wn1, wn2, wn3, wn4, wn5, wn6, wn7))

# combine all temperatures into common grid
all_t = np.unique([g["temperature"] for g in groups])

# construct single cross section array
all_xs = np.zeros((len(all_t), len(all_wn)))
regions = [
    list(range(15)),
    [15],
    [16],
    [17],
    [18],
    [19],
    list(range(20, 24)),
    list(range(24, 29)),
]
for region in regions:
    t0 = np.array(
        [groups[k]["temperature"] for k in region]
    )  # temperatures of given data
    for i in range(len(all_t)):
        t = all_t[i]
        if t in t0:  # insert data where all_t matches temperatures of given data
            k = np.searchsorted(t0, t, side="left") + region[0]
            wn = groups[k]["wavenumber"]
            xs = groups[k]["cross_section"]
            j = np.searchsorted(all_wn, wn, side="left")
            all_xs[i, j] = xs
        elif t > np.min(t0) and t < np.max(
            t0
        ):  # interpolate where all_t does not match but is bounded by the temperatures of the given data
            # bound indices
            k2 = np.searchsorted(t0, t, side="right")  # right: subscript 2
            k1 = k2 - 1  # left: subscript 1
            # weights
            w2 = (t - t0[k1]) / (t0[k2] - t0[k1])
            w1 = 1.0 - w2

            k1, k2 = k1 + region[0], k2 + region[0]
            wn1 = groups[k1]["wavenumber"]
            wn2 = groups[k2]["wavenumber"]
            wn = np.union1d(wn1, wn2)

            # pad interpolation data for where grids don't overlap
            xs1 = np.zeros_like(wn)
            j1 = np.searchsorted(wn, wn1, side="left")
            xs1[j1] = groups[k1]["cross_section"]
            xs2 = np.zeros_like(wn)
            j2 = np.searchsorted(wn, wn2, side="left")
            xs2[j2] = groups[k2]["cross_section"]

            # interpolate
            j = np.searchsorted(all_wn, wn, side="left")
            all_xs[i, j] = w1 * xs1 + w2 * xs2
        else:  # extend where all_t does not match and is not bounded by the temperatures of the given data
            k = region[0] + (
                len(t0) - 1 if t > np.max(t0) else 0
            )  # copy min temp data if t < min temp, max temp data if t > max temp
            wn = groups[k]["wavenumber"]
            xs = groups[k]["cross_section"]
            j = np.searchsorted(all_wn, wn, side="left")
            all_xs[i, j] = xs

all_xs_wl = all_xs[:, ::-1]
all_wl = 1e7 / all_wn[::-1]
ds = xr.Dataset(
    {"xs": (["temperature", "wavelength_nm"], all_xs_wl / 1e10)},  # cm^5 -> m^5
    coords={
        "temperature": all_t,
        "wavelength_nm": all_wl,  # assuming vacuum wavelength (see note below)
        # "wavelength_nm": sk.optical.air_wavelength_to_vacuum_wavelength(all_wl),
    },
)

ds.to_netcdf(fname_out)

# note on vacuum vs air wavelength:
# the wavelengths for the 477 nm absorption band are taken from Thalman and Volkamer (2013)
# they calibrated their detector using "Solar Flux Atlas from 296 to 1300 nm" (Kurucz, 1984)
# if the atlas is in vacuum wavelengths, then the o2o2 cross sections are in vacuum wavelength
# I didn't find a copy of the atlas, but I am assuming it is in vacuum wavelength
