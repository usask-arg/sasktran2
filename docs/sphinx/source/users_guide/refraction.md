---
file_format: mystnb
---

(_users_refraction)=
# Including Refraction

```{code-cell}
import sasktran2 as sk
```

SASKTRAN2 supports refraction through three different configuration parameters:

- {py:attr}`sasktran2.Config.los_refraction` which controls refraction for the observer (line of sight) rays
- {py:attr}`sasktran2.Config.solar_refraction` which controls refraction for the solar rays.
- {py:attr}`sasktran2.Config.multiple_scatter_refraction` which controls refraction for rays in the multiple scatter source calculation


```{code-cell}
config = sk.Config()

print(f"los_refraction: {config.los_refraction}, solar_refraction: {config.solar_refraction}, multiple_scatter_refraction: {config.multiple_scatter_refraction}")
```
By default all of these parameters are set to False, indicating that refraction is disabled.

All of these settings have their own quirks and may only work with specific configuration options, see the sections below for
each option.

## Setting the Refractive Index
Refractive effects in the atmosphere require knowledge of the refractive index of air.
The refractive index is set through the {py:attr}`sasktran2.Geometry1D.refractive_index` property, which is
the air refractive index on an altitude grid.
In SASKTRAN2 refraction is considered a geometry (wavelength independent effect) for computational efficiency reasons.
If wavelength dependent refractive index calculations are important for your application, the model can be reset and executed
for a single wavelength at a time.

SASKTRAN2 includes an implementation of the refractive index calculations described by Ciddor (1996) in "Optical constants of Air"
where the refractive index can be set through profiles of pressure, temperature, humidity, and CO2 concentration.  For most applications
the humidity and CO2 concentration are very small perturbations of the refractive index.  An example of how to set the refractive
index is

```python
atmosphere = ...

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

model_geometry.refractive_index = sk.optical.refraction.ciddor_index_of_refraction(
    atmosphere.temperature_k, atmosphere.pressure_pa, 0.0, 400, 600
)
```

which sets the refractive index from the US76 standard atmosphere parameters, ignoring humidity effects, with a CO2
concentration of 400 ppm, and at 600 nm.

## Line of Sight Refraction
Line of sight refraction refers to refractive bending of the observer lines of sight, and is the dominant form of refraction
when viewing in the atmospheric limb or at high viewing zenith angles.  It is controlled through the
{py:attr}`sasktran2.Config.los_refraction` option.  It is generally an important effect when viewing in the limb with tangent altitudes
below approximately 20 km, and causes a bending of the ray downwards towards lower tangent altitudes, usually increasing the observed signal.

```{code-cell}
import sasktran2 as sk
import numpy as np
import matplotlib.pyplot as plt


tan_alts = np.arange(1000, 65000, 1000)


def run(include_refraction):
    config = sk.Config()
    config.los_refraction = include_refraction

    model_geometry = sk.Geometry1D(cos_sza=0.6,
                                    solar_azimuth=0,
                                    earth_radius_m=6372000,
                                    altitude_grid_m=np.arange(0, 100001, 1000),
                                    interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                    geometry_type=sk.GeometryType.Spherical)

    viewing_geo = sk.ViewingGeometry()


    for alt in tan_alts:
        ray = sk.TangentAltitudeSolar(tangent_altitude_m=alt,
                                        relative_azimuth=0,
                                        observer_altitude_m=200000,
                                        cos_sza=0.6)
        viewing_geo.add_ray(ray)

    wavel = np.array([600])
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    model_geometry.refractive_index = sk.optical.refraction.ciddor_index_of_refraction(atmosphere.temperature_k, atmosphere.pressure_pa, 0.0, 200, 600)

    atmosphere['rayleigh'] = sk.constituent.Rayleigh()
    #atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
    #atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())

    engine = sk.Engine(config, model_geometry, viewing_geo)
    return engine.calculate_radiance(atmosphere)

refracted = run(True)
unrefracted = run(False)

plt.plot((refracted["radiance"] / unrefracted["radiance"]).to_numpy().flatten(), tan_alts)
```

## Solar Refraction
Solar refraction refers to refraction of the incoming solar beam, it is most important for twilight conditions where the sun is low in the sky.
Typically solar refraction is unimportant until the solar zenith angle is above 85 degrees.  The option is controlled by
{py:attr}`sasktran2.Config.solar_refraction`.  Note that currently this option is considered experimental, and only works when
both the single scatter source and multiple scatter source are set to DiscreteOrdinates, and thus only works for nadir viewing
applications.

## Multiple Scatter Refraction
Refraction effects during the multiple scatter calculation are usually negligible.  The option {py:attr}`sasktran2.Config.multiple_scatter_refraction`
is included to estimate them.  The option only works when the multiple scatter source is set to SuccessiveOrders.
