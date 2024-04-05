---
file_format: mystnb
---

(_users_brdf)=
# Setting the Surface BRDF

```{code-cell}
import sasktran2 as sk
```

SASKTRAN2 has the capability to include surface reflectance in the form of a Bi-directional Reflectance Distribution Function (BRDF).
The simplest BRDF is a Lambertian surface where the outgoing radiation is the same in all directions.

Let's start by setting up a calculation where we view the ground at a variety of different viewing zenith angles.

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt

config = sk.Config()
#config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
config.num_streams = 4

model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=np.arange(0, 65001, 1000),
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

viewing_geo = sk.ViewingGeometry()

cos_viewing_angle = np.arange(0.01, 0.99, 0.05)

for mu in cos_viewing_angle:
    ray = sk.GroundViewingSolar(0.6, np.pi, mu, 100000)
    viewing_geo.add_ray(ray)

for mu in cos_viewing_angle[::-1]:
    ray = sk.GroundViewingSolar(0.6, 0, mu, 100000)
    viewing_geo.add_ray(ray)

wavel = np.array([600])

atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere['rayleigh'] = sk.constituent.Rayleigh()
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())

engine = sk.Engine(config, model_geometry, viewing_geo)
```

BRDFs are added to the atmosphere the same way any constituent is.  Let's add a Lambertian surface
with a constant surface albedo of 1.0

```{code-cell}
atmosphere['brdf'] = sk.constituent.LambertianSurface(1.0)
```

And calculate the radiance

```{code-cell}
radiance = engine.calculate_radiance(atmosphere)

radiance["radiance"].isel(stokes=0).sel(wavelength=600).plot()
```

We can change the surface to a snow surface and compare the results,

```{code-cell}
atmosphere['brdf'] = sk.constituent.SnowKokhanovsky()

radiance_snow = engine.calculate_radiance(atmosphere)

radiance["radiance"].isel(stokes=0).sel(wavelength=600).plot()
radiance_snow["radiance"].isel(stokes=0).sel(wavelength=600).plot()
```

And we see a shifting of radiance values towards the later lines of sight,
which have relative azimuth 0 and are in the glint region.

# Surface BRDF Linearization
Each BRDF provides linearizations with respect to it's input parameters.
If we look at the radiance calculated with the Lambertian surface,

```{code-cell}
print(radiance)
```
We see there is a field `wf_brdf_albedo` that represents the derivative of
the observations with respect to our input albedo value.

Looking at the snow calculation,
```{code-cell}
print(radiance_snow)
```
We see two parameters, `wf_brdf_L` and `wf_brdf_M`. These are the derivatives
with respect to the input parameters `L` and `M` of the Kokhanovsky model.

# Available BRDFs
A full list of available surface parameterizations can be found at [BRDFs](../api/constituents.rst#brdfs).
