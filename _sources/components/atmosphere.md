---
file_format: mystnb
---

(_atmosphere)=
# The Atmospheric State

The atmospheric state at it's fundamental core is:

 - Extinction specified on the model geometry grid points.
 - Single scatter albedo specified on the model geometry grid points.
 - (optional) Scattering parameters (Legendre coefficients, phase functions) specified at the model geometry grid points
 - A representation of surface scattering (e.g. Lambertian surface)

It is possible to directly specify these quantities, however the easier way to use the model
is through what we refer to as constituents.

```{note}
SASKTRAN2 is designed to perform efficient calculations for a static geometry and
multiple atmospheric states.  For this reason all atmospheric optical quantities contain
an explicit `wavelength` dimension.  Instead of extinction being a 1-D array at the geometry
grid points, it is a 2-D array (wavelength, geometry).
We refer to this dimension as the wavelength dimension
because that is most commonly what it is, but in reality it can be any dimension where the
atmospheric state varies.
```

## Setting the Atmosphere from Constituents
Most of the time the easiest way to construct the atmosphere is using
something that we call the `constituent interface`.  All this means is
that the atmosphere is constructed in pieces through objects called constituents rather than all at once.
An example constituent could be Rayleigh scattering, or ozone absorption, or a Lambertian surface.

We start by setting up our {py:class}`sasktran2.Atmosphere` object

```{code-cell}
import numpy as np
import sasktran2 as sk

model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=np.arange(0, 65001, 1000),
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

config = sk.Config()

wavel = np.arange(280.0, 800.0, 0.1)
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)
```

Note that here we have explicity passed the `wavelength_nm` parameter to the atmospheric state.
Most atmospheric constituents require that the wavelength (or wavenumber) is known in order to look up
things like cross sections from databases.  On the same token, most constituents require that atmospheric
temperature and/or pressure is known for similar reasons, we can add that to the atmosphere


```{code-cell}
sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere.pressure_pa

atmosphere.temperature_k
```

In this case we used one of SASKTRAN2's built in climatologies to set temperature/pressure. Alternatively
the properties {py:attr}`sasktran2.Atmosphere.pressure_pa` and {py:attr}`sasktran2.Atmosphere.temperature_k`
can be directly accessed and set to any user desired values.

Now that our base atmosphere is constructed, we can start adding constituents

```{code-cell}
atmosphere["rayleigh"] = sk.constituent.Rayleigh()
```

will add Rayleigh scattering to the atmosphere.  If we want to add ozone absorption with a constant
VMR profile of 1ppm, we can do,

```{code-cell}
alt_grid = np.arange(0, 100001, 10000)

atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
    optical_property=sk.optical.O3DBM(),
    altitudes_m=alt_grid,
    vmr=np.ones_like(alt_grid) * 1e-6
)
```

There are a few things to pay attention to here. The first is that the altitude grid is different
than our global geometry grid. Constituents aren't limited to being specified on the global geometry grid,
interpolation to this grid will be internally performed at some stage though. The second thing is that we provided
something known as an `optical property`.  An optical property is, in essence, a database of particle
cross sections and possibly scattering properties.

Once a constituent is added to the atmosphere it is possible to access it directly and modify it

```{code-cell}
atmosphere["ozone"].vmr *= 2

print(atmosphere["ozone"].vmr)
```

## Available Constituents

```{eval-rst}
.. autosummary::
    sasktran2.constituent.Rayleigh
    sasktran2.constituent.VMRAltitudeAbsorber
    sasktran2.constituent.ExtinctionScatterer
    sasktran2.constituent.NumberDensityScatterer
    sasktran2.constituent.LambertianSurface
```
