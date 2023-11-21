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

## Flux Output
Flux output is planned for a future version.
