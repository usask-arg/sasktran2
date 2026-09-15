---
file_format: mystnb
mystnb:
  execution_raise_on_error: true
---

(_example_ers)=
# CAIRT Extended Reference Scenarios

The {py:mod}`sasktran2.climatology.ers` module creates reference atmospheres from
[CAIRT ERS v7](https://doi.org/10.5281/zenodo.10022129), assembled by Quentin Errera
from WACCM, BASCOE, and ACE-FTS sources. It provides 31 gas species, pressure,
temperature, and climatological standard deviations on an altitude grid from
0 to 200 km at 1 km intervals.

ERS uses multiannual averages. **CH4, N2O, CO2, SF6, CCl4, HCFC22, CFC11, and
CFC12 are scaled to expected 2030 abundances.** Selecting a month does not select
a calendar year. Source models have different altitude coverage; numerical
values throughout 0–200 km do not establish equal scientific support for all
species over that range.

## Select a scenario

This implementation selects the supplied scenarios exactly. It does not
interpolate between seasons, latitude bands, or local times.

| Argument | Choices |
| --- | --- |
| `month` | `1`, `4`, `7`, `10` |
| `latitude_degrees` | `-80`, `-45`, `0`, `45`, `80` |
| `local_time_hours` | `9.5`, `21.5` |
| `solar_activity` | `"minimum"` (default), `"maximum"` |
| `volcanic_activity` | `"background"` (default), `"enhanced"` |

Latitude coordinates represent bands: 90–70°S, 55–35°S, 20°S–20°N, 35–55°N,
and 70–90°N. A request for 50°N raises an error; explicitly select the 45°N band
for that reference scenario. Local times are overpass times, so they should not
be interpreted as universal day/night categories, especially near the poles.

```{code-cell}
import numpy as np
import sasktran2 as sk
import matplotlib.pyplot as plt

scenario = {
    "month": 7,
    "latitude_degrees": 45,
    "local_time_hours": 9.5,
    "solar_activity": "minimum",
    "volcanic_activity": "background",
}

profiles = sk.climatology.ers.profile(species=["O3", "NO2"], **scenario)
profiles[["o3_mean", "o3_std", "temperature_k", "pressure_pa"]]
```

The first call downloads only the 4.9 MB NetCDF into
`<database_root>/climatology/ers/v07/CAIRT_ERS_v07.nc`. Later calls verify its
checksum and use the cache. The plot archive is not downloaded, and the optional
`zenodo-get` dependency is not needed. Transient server and connection failures
are retried up to twice. Persistent failures propagate, and incomplete downloads
are not retained as valid cache files.

The returned xarray dataset is loaded into memory with no open file handle.
Its altitude coordinate is `altitude_m` in metres. Gas means and standard
deviations retain lower-case names such as `o3_mean` and `o3_std`, in mol/mol.
Temperature is in K and pressure in Pa. The selected band bounds, source
checksum, attribution, and scenario choices are included.

Only variables with a given scenario dimension depend on that choice. For
example, CFC11 has monthly profiles without a local-time dimension, whereas
NO additionally depends on solar activity and SO2 depends on volcanism.

## Plot the selected gas profiles

We can plot the July northern-midlatitude profiles over the altitude range used
in the radiance example below. Shading shows one climatological standard
deviation about the mean, with the lower edge limited to zero for display.
This spread describes variability, not uncertainty in the mean.

```{code-cell}
lower_atmosphere = profiles.sel(altitude_m=slice(0, 80000))
altitude_km = lower_atmosphere.altitude_m / 1000

fig, axes = plt.subplots(1, 2, figsize=(9, 5), sharey=True)
for ax, gas, scale, label in zip(
    axes,
    ["o3", "no2"],
    [1e6, 1e9],
    ["O$_3$ VMR [ppmv]", "NO$_2$ VMR [ppbv]"],
):
    mean = lower_atmosphere[f"{gas}_mean"] * scale
    spread = lower_atmosphere[f"{gas}_std"] * scale
    ax.plot(mean, altitude_km, label="Mean")
    ax.fill_betweenx(
        altitude_km,
        np.maximum(mean - spread, 0),
        mean + spread,
        alpha=0.2,
        label="Mean ± one standard deviation",
    )
    ax.set_xlabel(label)
    ax.set_xlim(left=0)
    ax.grid(alpha=0.3)

axes[0].set_ylabel("Altitude [km]")
axes[0].set_ylim(0, 80)
axes[0].legend(fontsize="small")
fig.suptitle("July, 35–55°N, 09:30 local time")
fig.tight_layout()
plt.show()
```

## Compare seasons and local times

Each curve below is an exact supplied scenario. On the left we vary the month
at 09:30; on the right we compare the two local times in July. The latitude band,
solar activity, and volcanic state stay fixed.

```{code-cell}
fig, axes = plt.subplots(1, 2, figsize=(9, 5), sharey=True)

for month, label in [(1, "January"), (4, "April"), (7, "July"), (10, "October")]:
    seasonal = sk.climatology.ers.profile(
        species="O3", **{**scenario, "month": month}
    ).sel(altitude_m=slice(0, 60000))
    axes[0].plot(
        seasonal.o3_mean * 1e6, seasonal.altitude_m / 1000, label=label
    )

for hour, label in [(9.5, "09:30"), (21.5, "21:30")]:
    overpass = sk.climatology.ers.profile(
        species="NO2", **{**scenario, "local_time_hours": hour}
    ).sel(altitude_m=slice(0, 60000))
    axes[1].plot(
        overpass.no2_mean * 1e9, overpass.altitude_m / 1000, label=label
    )

axes[0].set_title("Seasonal O$_3$ at 09:30")
axes[0].set_xlabel("O$_3$ VMR [ppmv]")
axes[0].set_ylabel("Altitude [km]")
axes[1].set_title("July NO$_2$ at two local times")
axes[1].set_xlabel("NO$_2$ VMR [ppbv]")
axes[0].set_ylim(0, 60)
for ax in axes:
    ax.set_xlim(left=0)
    ax.legend()
    ax.grid(alpha=0.3)
fig.suptitle("ERS reference scenarios, 35–55°N")
fig.tight_layout()
plt.show()
```

## Create an atmosphere

The same scenario can initialize an atmosphere for a spectrum calculation.

```{code-cell}
config = sk.Config()
config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
config.num_streams = 4

geometry = sk.Geometry1D(
    cos_sza=0.6,
    solar_azimuth=0,
    earth_radius_m=6372000,
    altitude_grid_m=np.arange(0, 80001, 1000.0),
    interpolation_method=sk.InterpolationMethod.LinearInterpolation,
    geometry_type=sk.GeometryType.Spherical,
)
atmosphere = sk.Atmosphere(
    geometry, config, wavelengths_nm=np.arange(300.0, 801.0, 5.0)
)

sk.climatology.ers.add_to_atmosphere(
    atmosphere,
    {"O3": sk.optical.O3DBM(), "NO2": sk.optical.NO2Vandaele()},
    **scenario,
)
atmosphere["rayleigh"] = sk.constituent.Rayleigh()

viewing = sk.ViewingGeometry()
viewing.add_ray(sk.GroundViewingSolar(0.6, 0, 0.8, 200000))
engine = sk.Engine(config, geometry, viewing)
radiance = engine.calculate_radiance(atmosphere)
```

The following spectrum includes ERS O3 and NO2 absorption and Rayleigh
scattering. SASKTRAN2 assumes unit incident solar irradiance here, so the plotted
radiance is normalized by that irradiance.

```{code-cell}
fig, ax = plt.subplots(figsize=(8, 4))
radiance["radiance"].isel(los=0, stokes=0).plot(ax=ax, x="wavelength")
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Solar-normalized radiance [sr$^{-1}$]")
ax.set_title("ERS July, 35–55°N: O$_3$, NO$_2$, and Rayleigh scattering")
ax.grid(alpha=0.3)
fig.tight_layout()
plt.show()
```

Geometry and solar angles remain caller choices; selecting ERS local time does
not derive those angles. Optical properties are also supplied by the caller.
The ERS helper initializes standard {py:class}`sasktran2.constituent.VMRAltitudeAbsorber`
objects, retaining their VMR derivatives.

The helper loads the source once for all requested gases and validates the
inputs before updating the atmosphere. Temperature and VMR are interpolated
linearly in altitude; pressure is interpolated logarithmically. It retains
source levels that bracket the model grid. Pass `set_pressure_temperature=False`
to preserve existing pressure and temperature, or an empty species mapping to
set pressure and temperature only.

`add_to_atmosphere` rejects model altitudes outside the source grid by default.
Explicit `out_of_bounds_mode="extend"` holds boundary gas VMR, temperature, and
pressure constant; it does not infer a physical extrapolation.

A single gas can also be constructed on the full native ERS grid:

```{code-cell}
ozone = sk.climatology.ers.constituent("O3", sk.optical.O3DBM(), **scenario)
```

This lower-level helper uses the standard constituent boundary policy: zero VMR
outside the native grid by default, or `out_of_bounds_mode="extend"` to hold
boundary VMR constant. It validates the entire native profile.

Species names are case insensitive. Aliases include `HNO4` for `HO2NO2`,
`F11`/`CFCl3` for `CFC11`, `F12`/`CF2Cl2` for `CFC12`, `F22`/`CHClF2` for
`HCFC22`, and `F14` for `CF4`. HDO requires isotope-specific optical properties;
it should not be paired with a bulk H2O absorber without accounting for isotope
abundances and possible double counting.

## Offline files and data quality

Every helper accepts `path="/path/to/CAIRT_ERS_v07.nc"` to use a local source,
or `db_root="/path/to/cache"` to override the configured download cache.
An explicit local file bypasses the pinned checksum check and records its
actual checksum, so subsets and modified files can be used deliberately.
The only supported schema/version is `version="v07"`.

To demonstrate local-file access, reuse the cache populated above. Selecting
the July southern-polar scenario also exposes two known data-quality flags.

```{code-cell}
from sasktran2.database.ers import ERSDatabase

cached_path = ERSDatabase().path()
raw = sk.climatology.ers.load_dataset(path=cached_path)
polar_profiles = sk.climatology.ers.profile(
    path=cached_path,
    species=["H2O", "O"],
    **{**scenario, "latitude_degrees": -80},
)
for name in ["h2o_mean", "o_mean"]:
    print(f"{name}: {polar_profiles[name].attrs['quality_flags']}")
```

Raw access preserves the source coordinates and values. Selected profile access
converts altitude and names the state fields, but also preserves gas values and
adds per-variable `quality_flags`. Neither path silently repairs data.

- **Negative H2O:** v7 contains 14 negative values in polar scenarios between
  191 and 200 km. Constituent helpers reject negative VMRs by default. Explicit
  `negative_vmr="clip"` replaces them with zero and emits a warning. The source
  file is unchanged. `add_to_atmosphere` checks only levels needed by its model
  grid, including interpolation brackets.
- **Zero atomic oxygen:** some polar O/O1D scenarios are entirely zero and their
  physical interpretation is unresolved. Constituent helpers reject these
  profiles. They remain accessible through `profile` for inspection.
- **Humidity:** the source wet/dry convention is not explicit. Gas values are
  passed through as supplied, assuming total-air mole fractions with unset/zero
  SASKTRAN2 specific humidity. `add_to_atmosphere` rejects nonzero existing
  humidity when adding gases. It does not derive humidity from H2O; callers using
  `constituent` directly must preserve this assumption themselves.
- **Pressure:** `surface_pressure_pa` is separate from `pressure_pa` at zero
  altitude. The vertical pressure profile is used as supplied; the surface value
  is not substituted or used to reconstruct it.
- **Air molar mass:** this source field has a suspect altitude dependence and is
  exposed with a quality flag, but is not used in the helpers. Number density
  follows from supplied pressure and temperature without needing molar mass.

With `species=None`, `profile` also exposes `h2so4m_c_mean`/`h2so4m_c_std` and
`airmolmass` in their original units. Condensed sulfuric acid is in µg/m³ and
cannot be treated as a gas VMR or extinction coefficient. An aerosol model would
also need composition, particle size, density, and refractive index assumptions.

Standard deviations describe climatological variability. No vertical or
cross-species covariance is supplied, and the helpers do not turn the standard
deviations into a retrieval covariance or uncertainty in the mean.

## Attribution

Errera, Quentin. *Extended reference scenarios (ERS) version 7*.
[Zenodo, DOI: 10.5281/zenodo.10022129](https://doi.org/10.5281/zenodo.10022129).
Distributed under CC BY 4.0. See the record's `CAIRT_ERS_v07_readme.pdf` for the
underlying model and observational sources. Cached v7 files are pinned to
published MD5 `71a5ed74d7056538cfd6a99d20ca3599`.
