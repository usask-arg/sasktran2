---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
  language: python
mystnb:
  output_stderr: remove
---

(_users_photochemical_emission)=
# Photochemical Emission

Use a volume emission rate (VER) profile to specify how much light the atmosphere
emits at each altitude. This example calculates O2 A-band limb radiance from a
VER profile, includes absorption by O2, and shows how to obtain weighting
functions for VER and temperature.

The examples use photon VER in photons m^-3 s^-1, summed over all directions.
SASKTRAN2 applies the conversion to emission per steradian for you.

## Calculate O2 A-band radiance

First, set up a limb view with a tangent altitude of 90 km. We enable volume
emission and turn off solar scattering for this example.

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import sasktran2 as sk

config = sk.Config()
config.emission_source = sk.EmissionSource.VolumeEmissionRate
config.single_scatter_source = sk.SingleScatterSource.NoSource
config.multiple_scatter_source = sk.MultipleScatterSource.NoSource

altitude_m = np.arange(0.0, 130_001.0, 2_000.0)
model_geometry = sk.Geometry1D(
    cos_sza=0.6,
    solar_azimuth=0.0,
    earth_radius_m=6_372_000.0,
    altitude_grid_m=altitude_m,
    interpolation_method=sk.InterpolationMethod.LinearInterpolation,
    geometry_type=sk.GeometryType.Spherical,
)
viewing_geo = sk.ViewingGeometry()
viewing_geo.add_ray(
    sk.TangentAltitudeSolar(
        tangent_altitude_m=90_000.0,
        relative_azimuth=0.0,
        observer_altitude_m=600_000.0,
        cos_sza=0.6,
    )
)
```

Choose a small wavelength interval within the A band to keep the example fast.
The 0.0002 nm spacing resolves the individual lines. For a wider band, keep a
fine enough calculation grid to resolve the emission and absorption before
averaging the radiance to your instrument's spectral resolution.

```{code-cell} ipython3
wavelength_nm = np.linspace(760.8, 761.2, 2001)
atmosphere = sk.Atmosphere(
    model_geometry,
    config,
    wavelengths_nm=wavelength_nm,
    temperature_derivative=True,
    pressure_derivative=False,
    specific_humidity_derivative=False,
)
sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
```

Next, supply the total VER of the 0-0 band as an altitude profile. Here we use an
illustrative Gaussian profile peaking at 94 km; replace it with your own values.
The VER grid can differ from the atmosphere's temperature grid. By default,
emission is zero outside the supplied VER altitude range.

Add an O2 absorber to account for self-absorption along the line of sight. The
O2 line data are downloaded and cached on first use.

```{code-cell} ipython3
:tags: [remove-output]

ver_altitude_m = np.arange(60_000.0, 130_001.0, 2_000.0)
ver_00 = 2.0e11 * np.exp(-0.5 * ((ver_altitude_m - 94_000.0) / 7_000.0) ** 2)

atmosphere["o2_00"] = sk.constituent.O2BandEmissionRate(
    ver_altitude_m, ver_00, band="0-0"
)
atmosphere["o2"] = sk.constituent.VMRAltitudeAbsorber(
    sk.optical.HITRANAbsorber("O2"),
    altitude_m,
    np.full_like(altitude_m, 0.2095),
)
```

Calculate and plot the spectrum. With photon VER and a wavelength grid in nm,
the radiance units are photons m^-2 s^-1 sr^-1 nm^-1.

```{code-cell} ipython3
engine = sk.Engine(config, model_geometry, viewing_geo)
result = engine.calculate_radiance(atmosphere)

result.radiance.isel(los=0, stokes=0).plot()
plt.xlabel("Wavelength [nm]")
plt.ylabel("Radiance [photons m$^{-2}$ s$^{-1}$ sr$^{-1}$ nm$^{-1}$]")
plt.title("O2 0-0 emission, 90 km tangent altitude")
plt.show()
```

You can use the same interface for other supported O2 bands:

| Band | Constructor argument |
| --- | --- |
| A-band 0-0 | `band="0-0"` |
| A-band 1-1 | `band="1-1"` |
| B-band 1-0 | `band="1-0"` |

Add each band under a different constituent name to give it an independent VER
profile. To also include scattered sunlight, enable a scattering source and add
`Rayleigh()` and `SolarIrradiance(photon_units=True)` when setting up the
atmosphere. See {ref}`_users_solar_irradiance` for absolute radiance units.

## Use VER and temperature weighting functions

The output contains `wf_o2_00_photon_ver` for changes to the input VER profile and
`wf_temperature_k` for changes to the atmospheric temperature profile. The
constituent name `o2_00` determines the name of its VER weighting function.

Attach the input altitude coordinates so you can select a particular height.
The following plots show the response to changing one profile value at 94 km.

```{code-cell} ipython3
result = result.assign_coords(
    altitude=altitude_m,
    o2_00_altitude=ver_altitude_m,
)
dI_dVER = result.wf_o2_00_photon_ver.isel(los=0, stokes=0)
dI_dT = result.wf_temperature_k.isel(los=0, stokes=0)

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
dI_dVER.sel(o2_00_altitude=94_000.0).plot(ax=axes[0])
dI_dT.sel(altitude=94_000.0).plot(ax=axes[1])
axes[0].set(ylabel="dI / dVER", xlabel="", title="VER change at 94 km")
axes[1].set(
    ylabel="dI / dT [radiance / K]",
    xlabel="Wavelength [nm]",
    title="Temperature change at 94 km",
)
fig.tight_layout()
plt.show()
```

The temperature weighting function holds band VER, pressure, and gas VMRs fixed.
It includes both the change in the emitted spectrum and the change in O2
self-absorption in this example. When comparing with an instrument, apply the
same spectral averaging to the radiance and its weighting functions.

If you only need radiance, construct the atmosphere with
`calculate_derivatives=False`. To keep VER weighting functions but omit
temperature weighting functions, use `temperature_derivative=False`.

## Update the profiles

Assign new values to `photon_ver` or `temperature_k` and calculate again. You can
reuse the constituent, atmosphere, and engine. Here we increase VER by 10% and
temperature by 5 K.

```{code-cell} ipython3
initial_temperature_k = atmosphere.temperature_k.copy()
atmosphere["o2_00"].photon_ver = 1.1 * ver_00
atmosphere.temperature_k = initial_temperature_k + 5.0
updated = engine.calculate_radiance(atmosphere)

result.radiance.isel(los=0, stokes=0).plot(label="Original profiles")
updated.radiance.isel(los=0, stokes=0).plot(label="VER +10%, temperature +5 K")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Radiance [photons m$^{-2}$ s$^{-1}$ sr$^{-1}$ nm$^{-1}$]")
plt.title("Updating VER and temperature")
plt.margins(y=0.2)
plt.legend()
plt.show()
```

## Start from excited-state populations

If you have an O2 excited-state population profile, use
{py:class}`sasktran2.constituent.PopulationEmissionRate` to convert it to band VER.
The dataset needs altitude in meters, temperature in kelvin, and `O2(b)` in
molecules m^-3. This example creates an illustrative population dataset using
the same altitude grid as above.

```{code-cell} ipython3
populations = xr.Dataset(
    {
        "temperature": (
            "altitude",
            np.interp(ver_altitude_m, altitude_m, initial_temperature_k),
        ),
        "O2(b)": (
            "altitude",
            3.0e12 * np.exp(-0.5 * ((ver_altitude_m - 94_000.0) / 7_000.0) ** 2),
        ),
    },
    coords={"altitude": ver_altitude_m},
)
population_emission = sk.constituent.PopulationEmissionRate(populations)
bands = population_emission.to_band_emissions()
atmosphere["o2_00"] = bands["0-0"]
atmosphere.temperature_k = initial_temperature_k
```

The conversion multiplies the population by the band's Einstein-A coefficient.
The returned band can then be updated through `photon_ver` just like the direct
VER example. Its spectrum follows `atmosphere.temperature_k` during each
calculation. Add the converted band instead of also adding `population_emission`,
which would count the emission twice.

## Calculate a single emission line

For an isolated line such as the oxygen green line at 557.7 nm, use
{py:class}`sasktran2.constituent.MonochromaticVolumeEmissionRate`. We reuse the
viewing geometry and define a separate atmosphere with an illustrative green-line
VER profile, again in photons m^-3 s^-1.

```{code-cell} ipython3
green_atmosphere = sk.Atmosphere(
    model_geometry,
    config,
    wavelengths_nm=np.linspace(557.5, 557.9, 41),
    calculate_derivatives=False,
)
green_ver = 3.0e9 * np.exp(-0.5 * ((ver_altitude_m - 97_000.0) / 5_000.0) ** 2)
green_atmosphere["green_line"] = sk.constituent.MonochromaticVolumeEmissionRate(
    ver_altitude_m, green_ver, 557.7
)
green_result = engine.calculate_radiance(green_atmosphere)

green_result.radiance.isel(los=0, stokes=0).plot()
plt.xlabel("Wavelength [nm]")
plt.ylabel("Radiance [photons m$^{-2}$ s$^{-1}$ sr$^{-1}$ nm$^{-1}$]")
plt.title("Oxygen green line, 90 km tangent altitude")
plt.show()
```

By default, this line is represented on the calculation grid so its integrated
radiance is preserved. Its plotted width therefore depends on the grid spacing;
use the integral over wavelength when comparing total line brightness.
