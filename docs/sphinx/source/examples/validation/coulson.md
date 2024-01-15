---
file_format: mystnb
---

(_example_coulson)=
# Reproducing the Coulson Tables
The Coulson tables are benchmark results for polarized radiative transfer concerning a single homogenous
Rayleigh scattering layer with a Lambertian surface underneath.  They were originally published in

Coulson K. L., Dave J. V. and Sekera Z. 1960 Tables Related to Radiation Emerging From a Planetary Atmosphere With Rayleigh Scattering (Berkeley, CA: Univ. California Press)

However they have since been refined to be more accurate in the following publication

Natraj, Vijay, King-Fai Li, and Yuk L. Yung. "Rayleigh scattering in planetary atmospheres: corrected tables through accurate computation of X and Y functions." The Astrophysical Journal 691.2 (2009): 1909.

Here we calculate some of the table values for cosine solar zenith angle of 0.2 with SASKTRAN2.

```{code-cell}
cos_sza = 0.2
```

## Model Configuration
We start by setting up the model

```{code-cell}
import sasktran2 as sk
import numpy as np

config = sk.Config()
config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates

config.num_streams = 40
config.num_singlescatter_moments = 40
config.num_stokes = 3
```

Here we have:

- Set both the single scatter source and the multiple scatter source to be the Discrete Ordinates method.
  the discrete ordinates method is the most accurate for the homogenous, plane-parallel slab that the Coulson
  table uses
- Set the number of Stokes parameters to 3 to enable polarization
- Set the number of full streams to 40.  This is sufficient to get a few decimal places of accuracy with the tables.

## Geometry configuration
Next we set up the model geometry,

```{code-cell}
model_geometry = sk.Geometry1D(cos_sza=cos_sza,
                               solar_azimuth=0,
                               earth_radius_m=6372000,
                               altitude_grid_m=np.array([0, 1.0]),
                               interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                               geometry_type=sk.GeometryType.PlaneParallel)

```

There are a couple of things to note here:

- We have put in two entries into the altitude grid. This is because in SASKTRAN2 we specify quantities on the edges of the layers.
  Two altitude grid points corresponds to a single layer
- Since we are in PlaneParallel mode, the earth radius parameter is not used.
- We have set the top altitude to 1m, this is is also arbitrary since we are in plane parallel mode

## Setting up the viewing geometry
While the tables cover a wide range of different viewing geometries, as an example we will only calculate two:

- A viewing cosine zenith of 0.02 with a relative azimuth of 30 degrees
- A viewing cosine zenith of 0.92 with a relative azimuth of 60 degrees

```{code-cell}
viewing_geo = sk.ViewingGeometry()

viewing_geo.add_ray(sk.GroundViewingSolar(cos_sza, np.deg2rad(30), 0.02, 2.0))
viewing_geo.add_ray(sk.GroundViewingSolar(cos_sza, np.deg2rad(60), 0.92, 2.0))
```

Here we have set the observer altitude to be 2.0m so that it is outside of the atmosphere.  Note that in SASKTRAN2
a positive cosine zenith corresponds to downward viewing since it is measured with respect to the propagation direction
(negative of the viewing direction).

## Setting up the atmosphere
Lastly we need to set up the atmosphere parameters, we use the low level interface to the atmosphere to directly
specify the optical quantities rather than using constituents.

```{code-cell}
atmosphere = sk.Atmosphere(model_geometry, config, calculate_derivatives=False, numwavel=1)


atmosphere.storage.total_extinction[:] = 0.5 / 1.0
atmosphere.storage.ssa[:] = 1

atmosphere.leg_coeff.a1[0, :, 0] = 1
atmosphere.leg_coeff.a1[2, :, 0] = 0.5

atmosphere.leg_coeff.a2[2] = 3
atmosphere.leg_coeff.b1[2] = -np.sqrt(6.0) / 2

atmosphere.surface.albedo[:] = 0.0
```

Here we have

- Set the extection to be 0.5 / m, since our atmosphere layer is 1m tall this is an optical depth of 0.5
- Set the Legendre coefficient expansion to be that of Rayleigh scattering.
- Set the surface albedo to 0.0

## Run the calculation
Next we run the calculation and multiply the results by `pi` since the original calculations assume a solar flux of `pi`

```{code-cell}
engine = sk.Engine(config, model_geometry, viewing_geo)
radiance = engine.calculate_radiance(atmosphere)

radiance["radiance"] *= np.pi

print(radiance["radiance"])
```

We can add in the values from the table,

```{code-cell}
radiance["true_radiance"] = (["wavelength", "los", "stokes"], np.array(
    [[[0.39444956, -0.06485313, 0.04390364],
      [0.05643322, -0.01979730, 0.03822653]]]
))
```

And print the percent differences,

```{code-cell}

print((radiance["radiance"] - radiance["true_radiance"]) / radiance["true_radiance"] * 100)
```
We see agreement at the 0.001% level.  Agreement is not perfect due to the number of streams used
inside the calculation, as well as slight adjustments that have to be made in the discrete ordinate
calculation to account for purely conservative scattering layers.
