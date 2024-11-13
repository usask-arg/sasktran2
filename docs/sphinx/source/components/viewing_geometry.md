---
file_format: mystnb
---
(_viewing_geometry)=
# The Viewing Geometry
The viewing geometry can be thought of as specifying what output quantities we need from the model.
These can be things like radiances at the top of the atmosphere, ground, or anywhere inbetween.  Or they can
be fluxes that are automatically integrated by the model.

## Radiance Output
We specify radiance output through a list of single positions and look vectors.  These are stored inside the
{py:class}`sasktran2.ViewingGeometry` object.


```{code-cell}
import sasktran2 as sk

viewing_geo = sk.ViewingGeometry()
```

Rather than specifying look directions and positions directly, they are specified through objects that we
colloquially refer to as Viewing Policies.  One example of a viewing policy is the {py:class}`sasktran2.TangentAltitudeSolar`
object, which defines a position and viewing direction based off an observer altitude, as well as various parameters
that are defined at the tangent point of the measurement.  For example,

```{code-cell}
ray = sk.TangentAltitudeSolar(tangent_altitude_m=10000,
                                relative_azimuth=0,
                                observer_altitude_m=200000,
                                cos_sza=0.6)
```

To add it to our overall container, we can do,

```{code-cell}
viewing_geo.add_ray(ray)
```

And we can look at all of the rays in our container

```{code-cell}
viewing_geo.observer_rays
```


## Viewing Policies
```{eval-rst}
.. autosummary::

    sasktran2.TangentAltitudeSolar
    sasktran2.GroundViewingSolar
```

## Working with Real Measurements

### The Solar Geometry Handler
`sasktran2` contains utilities to compute the solar angles automatically for a given location and time
using the `astropy` package.  To use this functionality you must have `astropy` installed.

For example, to calculate the (solar zenith angle, solar azimuth angle) at a given location we can do,

```{code-cell}
import pandas as pd


solar = sk.solar.SolarGeometryHandlerAstropy()

print(solar.target_solar_angles(
    latitude=20,
    longitude=-100,
    altitude=0,
    time=pd.Timestamp("2024-11-12 20:00:00")
))
```

Solar azimuth angles are always measured from true north, with 90 degrees pointing in the east direction.

```{eval-rst}
.. note::
    The solar handler specifies angles in degrees, whereas most other aspects of SASKTRAN2 use radians.
```

### Converting from ECEF Coordinates
In the original `sasktran1` package, all rays were specified as a triplet of
(observer position, local look vector, time) in Earth Centered Earth Fixed coordinates.

Here we have pre-computed the position and look vector of a satellite directly above
Saskatoon, Canada, and we will specify a local time of noon (18 UTC).

```{code-cell}
import numpy as np

observer = np.array([-1158730.59368676, -3875262.18406142,  5170772.5134034 ])
look_vector = np.array([ 0.17579194,  0.58791911, -0.78958743])

time = pd.Timestamp("2024-11-12 18:00:00")
```

We can then convert the ray to one `sasktran2` recognizes,

```{code-cell}
ray = sk.viewinggeo.ecef_to_sasktran2_ray(
    observer=observer,
    look_vector=look_vector,
    time=time,
    solar_handler=solar
)
print(ray)
```

Note here we re-used our solar handler from the previous section.


## Flux Output
Flux output is planned for a future version.
