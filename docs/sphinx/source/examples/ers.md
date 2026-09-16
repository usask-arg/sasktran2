---
file_format: mystnb
mystnb:
  execution_raise_on_error: true
---

(_example_ers)=
# CAIRT Extended Reference Scenarios

Use {py:mod}`sasktran2.climatology.ers` to explore seasonal gas profiles and
initialize a radiative-transfer calculation with
[CAIRT ERS v7](https://doi.org/10.5281/zenodo.10022129). ERS provides 31 gases,
pressure, and temperature for a choice of latitude band, month, local time,
solar activity, and volcanic activity. These are multiannual reference
conditions; several long-lived gases are scaled to projected 2030 abundances.

## Select a scenario

Start with July in the northern midlatitudes at 09:30 local time. Choose from
the supplied scenarios:

| Argument | Choices |
| --- | --- |
| `month` | `1`, `4`, `7`, `10` |
| `latitude_degrees` | `-80`, `-45`, `0`, `45`, `80` |
| `local_time_hours` | `9.5`, `21.5` |
| `solar_activity` | `"minimum"` (default), `"maximum"` |
| `volcanic_activity` | `"background"` (default), `"enhanced"` |

The latitude values select the bands 90–70°S, 55–35°S, 20°S–20°N, 35–55°N,
and 70–90°N, respectively. For example, use `45` for the 35–55°N band.
Selections must match the table; there is no interpolation between scenarios.
The data downloads automatically on first use and is cached for later calls.

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
```

The returned xarray dataset uses `altitude_m` in metres and gas volume mixing
ratios (VMRs) in mol/mol. Multiply by `1e6` for parts per million (ppmv) or
`1e9` for parts per billion (ppbv).

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

## Compare solar and volcanic activity

Solar activity is most visible in upper-atmosphere temperature. Volcanic activity
changes the SO2 profile. Here we compare both choices for the same July
midlatitude scenario, using an altitude range suited to each quantity.
`species=[]` selects pressure and temperature without any gases.

```{code-cell}
fig, axes = plt.subplots(1, 2, figsize=(9, 5))

for activity in ["minimum", "maximum"]:
    solar = sk.climatology.ers.profile(
        species=[], **{**scenario, "solar_activity": activity}
    ).sel(altitude_m=slice(100000, 200000))
    axes[0].plot(
        solar.temperature_k, solar.altitude_m / 1000, label=activity.capitalize()
    )

for activity in ["background", "enhanced"]:
    volcanic = sk.climatology.ers.profile(
        species="SO2", **{**scenario, "volcanic_activity": activity}
    ).sel(altitude_m=slice(10000, 40000))
    axes[1].plot(
        volcanic.so2_mean * 1e9,
        volcanic.altitude_m / 1000,
        label=activity.capitalize(),
    )

axes[0].set_title("Solar activity")
axes[0].set_xlabel("Temperature [K]")
axes[0].set_ylim(100, 200)
axes[1].set_title("Volcanic activity")
axes[1].set_xlabel("SO$_2$ VMR [ppbv]")
axes[1].set_xlim(left=0)
axes[1].set_ylim(10, 40)
for ax in axes:
    ax.set_ylabel("Altitude [km]")
    ax.legend()
    ax.grid(alpha=0.3)
fig.suptitle("July, 35–55°N, 09:30 local time")
fig.tight_layout()
plt.show()
```

## Create an atmosphere

Use `add_to_atmosphere` to set pressure, temperature, and gas abundances from
ERS on the model's altitude grid. Supply an optical property for each gas to
calculate absorption. Here we add O3, NO2, and Rayleigh scattering for a
satellite view toward the ground.

Set the solar and viewing angles for your calculation explicitly; ERS local
time only selects the reference profiles.

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

optical_properties = {"O3": sk.optical.O3DBM(), "NO2": sk.optical.NO2Vandaele()}
sk.climatology.ers.add_to_atmosphere(atmosphere, optical_properties, **scenario)
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

## Change a gas abundance

ERS constituents can be adjusted like any other
{py:class}`sasktran2.constituent.VMRAltitudeAbsorber`. Use `constituent` to create
an individual gas, then change its `vmr` before adding it to the atmosphere.
For example, halve the ozone abundance and compare the spectrum with the
original calculation. Pressure, temperature, NO2, and viewing geometry stay
fixed, so the difference shows the effect of changing ozone alone.

```{code-cell}
ozone = sk.climatology.ers.constituent("O3", optical_properties["O3"], **scenario)
ozone.vmr *= 0.5
atmosphere["O3"] = ozone
reduced_ozone_radiance = engine.calculate_radiance(atmosphere)

fig, ax = plt.subplots(figsize=(8, 4))
for result, label in [
    (radiance, "ERS ozone"),
    (reduced_ozone_radiance, "50% of ERS ozone"),
]:
    result["radiance"].isel(los=0, stokes=0).plot(
        ax=ax, x="wavelength", label=label
    )
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Solar-normalized radiance [sr$^{-1}$]")
ax.set_title("Effect of reducing ozone")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
plt.show()
```

## Keep an existing pressure and temperature profile

If your atmosphere already has pressure and temperature from another source,
pass `set_pressure_temperature=False` to add only ERS gases. For this example,
use the US76 standard atmosphere and compare its state with ERS.

```{code-cell}
custom_atmosphere = sk.Atmosphere(
    geometry, config, wavelengths_nm=atmosphere.wavelengths_nm
)
sk.climatology.us76.add_us76_standard_atmosphere(custom_atmosphere)
sk.climatology.ers.add_to_atmosphere(
    custom_atmosphere,
    optical_properties,
    set_pressure_temperature=False,
    **scenario,
)

altitude_km = geometry.altitudes() / 1000
fig, axes = plt.subplots(1, 2, figsize=(9, 5), sharey=True)
for model, label in [(atmosphere, "ERS"), (custom_atmosphere, "US76")]:
    axes[0].plot(model.temperature_k, altitude_km, label=label)
    axes[1].semilogx(model.pressure_pa / 100, altitude_km, label=label)

axes[0].set_xlabel("Temperature [K]")
axes[0].set_ylabel("Altitude [km]")
axes[0].set_ylim(0, 80)
axes[1].set_xlabel("Pressure [hPa]")
for ax in axes:
    ax.legend()
    ax.grid(alpha=0.3)
fig.suptitle("Choosing pressure and temperature independently of ERS gases")
fig.tight_layout()
plt.show()
```

To initialize only pressure and temperature from ERS, pass an empty gas
mapping: `sk.climatology.ers.add_to_atmosphere(atmosphere, {}, **scenario)`.
See the {ref}`climatology API <api_climatology>` for the full set of options.

## Attribution

Errera, Quentin. *Extended reference scenarios (ERS) version 7*.
[Zenodo, DOI: 10.5281/zenodo.10022129](https://doi.org/10.5281/zenodo.10022129).
Distributed under CC BY 4.0.
