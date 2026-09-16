---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

(_example_stratospheric_aerosol)=
# Add stratospheric aerosol

Use `sk.climatology.stratospheric_aerosol` to add a sulfate aerosol reference
scenario to an existing atmosphere:

```python
import sasktran2 as sk

aerosol = sk.climatology.stratospheric_aerosol
atmosphere["stratospheric_aerosol"] = aerosol.constituent(
    "tropical_typical",
    altitudes_m=atmosphere.model_geometry.altitudes(),
)
```

The helper supplies extinction, particle size, and sulfate scattering properties
from a scenario derived from SAGE III–ISS observations. By default it smooths the
extinction and extends the profile below and above the observations. Pass your
atmosphere's altitude grid as shown to cover the full calculation domain.

You can add this constituent alongside your chosen gases and Rayleigh scattering,
including atmospheres using ERS gas scenarios. Add a separate aerosol constituent
if you also want tropospheric aerosol.

## Choose a scenario

Combine a latitude band with a loading level, for example `nh_midlat_elevated`.
There are twelve scenarios:

| Latitude band | Name prefix | Available loading levels |
| --- | --- | --- |
| Southern midlatitudes, 55–35°S | `sh_midlat` | `low`, `typical`, `elevated`, `extreme` |
| Tropics, 20°S–20°N | `tropical` | `low`, `typical`, `elevated`, `extreme` |
| Northern midlatitudes, 35–55°N | `nh_midlat` | `low`, `typical`, `elevated`, `extreme` |

Start with `typical`, then use the other levels to explore the effect of aerosol
loading. The levels represent approximately the 10th, 50th, 90th, and 99th
percentiles of 756 nm optical depth between 18 and 30 km within each latitude
band. These are fixed reference cases, so they do not select conditions for a
particular date or season.

You can list the scenario names in Python:

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt
import sasktran2 as sk

aerosol = sk.climatology.stratospheric_aerosol
aerosol.scenarios()[["latitude_band", "loading"]]
```

## Compare the profiles

Use `profile()` to inspect extinction or particle size before adding a scenario.
The extinction is at 756 nm, in m⁻¹; altitude is in metres and median radius is in
nanometres.

```{code-cell}
fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
bands = {
    "sh_midlat": "Southern midlatitudes",
    "tropical": "Tropics",
    "nh_midlat": "Northern midlatitudes",
}
for ax, (band, title) in zip(axes, bands.items()):
    for loading in ["low", "typical", "elevated", "extreme"]:
        profile = aerosol.profile(f"{band}_{loading}")
        ax.semilogx(
            profile.extinction_per_m, profile.altitude_m / 1000, label=loading
        )
    ax.set(title=title, xlabel="Extinction at 756 nm [m$^{-1}$]",
           xlim=(1e-11, 1e-4), ylim=(0, 55))
    ax.grid(alpha=0.2)
axes[0].set_ylabel("Altitude [km]")
axes[-1].legend()
fig.tight_layout()
```

The curves include the smooth extensions outside the measured altitude range.

## Adjust smoothing and altitude extensions

The defaults work without any additional arguments:

| Option | Default | Effect |
| --- | --- | --- |
| `smoothing_fwhm_m` | `1500.` | Smooths extinction over a 1.5 km width while preserving optical depth over the observed interval. Set to `0.` to disable smoothing. |
| `lower_scale_height_m` | `2000.` | Controls the downward exponential taper, which reaches zero at the ground. A larger value gives a broader lower tail. |
| `ground_altitude_m` | `0.` | Sets the altitude at and below which this aerosol is zero. |
| `upper_scale_height_m` | `"reference"` | Uses 2.8 km for either midlatitude band and 3.6 km for the tropics. A larger value gives a more slowly decreasing upper tail. |

All four loading levels in a latitude band share the same default upper scale
height, based on regular aerosol conditions. Choosing `extreme` increases loading
without introducing a different upper decay rate. Particle size is held constant
outside the observed interval.

Pass overrides when adding the aerosol:

```python
atmosphere["stratospheric_aerosol"] = aerosol.constituent(
    "tropical_extreme",
    altitudes_m=atmosphere.model_geometry.altitudes(),
    smoothing_fwhm_m=2000.,
    lower_scale_height_m=1500.,
    upper_scale_height_m=3200.,
)
```

The same options work with `profile()` so you can plot a customised profile first.
To remove either extension, set `lower_extension="zero"` or
`upper_extension="zero"`. Extinction is then zero outside that end of the observed
interval, with an abrupt cutoff.

## Calculate radiance with aerosol

This complete example creates an atmosphere, adds a typical tropical aerosol
scenario, and calculates a limb radiance at three wavelengths:

```{code-cell}
config = sk.Config()
model_geometry = sk.Geometry1D(
    cos_sza=0.6,
    solar_azimuth=0.,
    earth_radius_m=6372000.,
    altitude_grid_m=np.arange(0., 65001., 1000.),
    interpolation_method=sk.InterpolationMethod.LinearInterpolation,
    geometry_type=sk.GeometryType.Spherical,
)
atmosphere = sk.Atmosphere(
    model_geometry, config, wavelengths_nm=np.array([525., 756., 1021.])
)
sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
atmosphere["rayleigh"] = sk.constituent.Rayleigh()
atmosphere["stratospheric_aerosol"] = aerosol.constituent(
    "tropical_typical", altitudes_m=model_geometry.altitudes()
)

viewing_geometry = sk.ViewingGeometry()
viewing_geometry.add_ray(sk.TangentAltitudeSolar(
    tangent_altitude_m=20000.,
    relative_azimuth=0.,
    observer_altitude_m=200000.,
    cos_sza=0.6,
))
engine = sk.Engine(config, model_geometry, viewing_geometry)
result = engine.calculate_radiance(atmosphere)
result.radiance
```

The first use may download the sulfate refractive-index data. Scattering is
calculated as needed; for repeated calculations across many wavelengths, you can
pass a compatible {py:class}`~sasktran2.database.MieDatabase` as the
`optical_property` argument to `constituent()`. Use sulfate with lognormal width
1.6, and include 756 nm as well as your calculation wavelengths in the table.
