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

```python
import numpy as np
import sasktran2 as sk

scenario = {
    "month": 7,
    "latitude_degrees": 45,
    "local_time_hours": 9.5,
    "solar_activity": "minimum",
    "volcanic_activity": "background",
}

profiles = sk.climatology.ers.profile(species=["O3", "NO2"], **scenario)
print(profiles[["o3_mean", "o3_std", "temperature_k", "pressure_pa"]])
```

The first call downloads only the 4.9 MB NetCDF into
`<database_root>/climatology/ers/v07/CAIRT_ERS_v07.nc`. Later calls verify its
checksum and use the cache. The plot archive is not downloaded, and the optional
`zenodo-get` dependency is not needed. Download failures propagate, and incomplete
downloads are not retained as valid cache files.

The returned xarray dataset is loaded into memory with no open file handle.
Its altitude coordinate is `altitude_m` in metres. Gas means and standard
deviations retain lower-case names such as `o3_mean` and `o3_std`, in mol/mol.
Temperature is in K and pressure in Pa. The selected band bounds, source
checksum, attribution, and scenario choices are included.

Only variables with a given scenario dimension depend on that choice. For
example, CFC11 has monthly profiles without a local-time dimension, whereas
NO additionally depends on solar activity and SO2 depends on volcanism.

## Create an atmosphere

```python
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
    geometry, config, wavelengths_nm=np.array([350.0, 500.0, 600.0])
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

```python
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

```python
raw = sk.climatology.ers.load_dataset(path="/path/to/CAIRT_ERS_v07.nc")
selected = sk.climatology.ers.profile(
    path="/path/to/CAIRT_ERS_v07.nc", **scenario
)
print(selected["h2o_mean"].attrs["quality_flags"])
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
