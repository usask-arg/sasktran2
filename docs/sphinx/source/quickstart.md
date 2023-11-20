---
file_format: mystnb
---

(_quickstart)=
# Quick Start
```{code-cell}
import sasktran2 as sk
```

SASKTRAN2 requires four things in order to perform the radiative transfer calculation

 * A configuration object which stores user options, {py:class}`sasktran2.Config`
 * A specification of the model geometry, and a global grid where atmospheric parameters are specified on.  This is handled by the {py:class}`sasktran2.Geometry1D` object
 * A {py:class}`sasktran2.ViewingGeometry` object that specifies the viewing conditions we want output at.
 * A {py:class}`sasktran2.Atmosphere` class that defines the atmospheric state

## The Configuration object
The {py:class}`sasktran2.Config` object stores all of the user configuration options that are necessary for the calculation,

```{code-cell}
config = sk.Config()
```

## The Model Geometry
The model geometry is the grid that the radiative transfer calculation is actually performed on.
Usually this involves specifying the coordinates (spherical, plane parallel, etc.) the number of dimensions the
atmosphere is allowed to vary in (1, 2, 3), as well as the actual grid values themselves.  Sometimes an internal
definition of a `reference point` is also required to provide context to the viewing geometry policies.

Similar to other SASKTRAN2 components, there are multiple ways to construct the model geometry.
The standard one, and the one you probably want, is {py:class}`sasktran2.Geometry1D`,

```{code-cell}
import numpy as np

model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=np.arange(0, 65001, 1000),
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)
```

## The Viewing Geometry
We specify radiance output through a list of single positions and look vectors.  These are stored inside the
{py:class}`sasktran2.ViewingGeometry` object.

```{code-cell}
viewing_geo = sk.ViewingGeometry()
```

We can add a set of four limb viewing lines of sight through,


```{code-cell}
for alt in [10000, 20000, 30000, 40000]:
    ray = sk.TangentAltitudeSolar(tangent_altitude_m=alt,
                                    relative_azimuth=0,
                                    observer_altitude_m=200000,
                                    cos_sza=0.6)
    viewing_geo.add_ray(ray)
```

## The Atmospheric State
We start by constructing our {py:class}`sasktran2.Atmosphere` object,

```{code-cell}
wavel = np.arange(280.0, 800.0, 0.1)
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)
```

```{code-cell}
sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
```

```{code-cell}
atmosphere['rayleigh'] = sk.constituent.Rayleigh()
```

```{code-cell}
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())
```

## Performing the Calculation
The radiative transfer calculation is done through the {py:class}`sasktran2.Engine` object,

```{code-cell}
engine = sk.Engine(config, model_geometry, viewing_geo)
```

Upon construction, all of the geometry information is calculated and cached. In order to do the actual
calculation, we pass in the {py:class}`sasktran2.Atmosphere` object

```{code-cell}
output = engine.calculate_radiance(atmosphere)
```

```{code-cell}
output['radiance'].isel(los=0).plot()
```
