---
file_format: mystnb
---

(_model_geometry)=
# The Model Geometry
The model geometry is the grid that the radiative transfer calculation is actually performed on.
Usually this involves specifying the coordinates (spherical, plane parallel, etc.) the number of dimensions the
atmosphere is allowed to vary in (1, 2, 3), as well as the actual grid values themselves.  Sometimes an internal
definition of a `reference point` is also required to provide context to the viewing geometry policies.

Similar to other SASKTRAN2 components, there are multiple ways to construct the model geometry.
The standard one, and the one you probably want, is {py:class}`sasktran2.Geometry1D`,

```{code-cell}
import sasktran2 as sk
import numpy as np

model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=np.arange(0, 100001, 1000),
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

```

The Reference Point
-------------------
Sometimes we refer to a location known as the `reference point`, which can be somewhat abstract and hard to understand
the significance of.  In the above construction, the reference point is a location with a cosine solar zenith angle
of 0.6, a solar azimuth angle of 0, and a spherical Earth radius of 6372 km.

As a user you should think of the reference point as being the location where you want the radiative transfer calculation
to be the most accurate.  For nadir viewing lines of sight, typically this means you should set the reference point to be
the ground viewing point.  For limb viewing lines of sight, the reference point is often best represented as the mean
tangent point of all the lines of sight.

The Location Grid
-----------------
The location grid is an abstract array of locations inside the atmosphere representative of the specific model geometry.
In the above construction it is simply an altitude grid above the surface.

The important thing to remember is that internally atmospheric quantities are sampled, and specified at, this location grid.
Therefore, it's resolution can have huge impacts on the accuracy and time complexity of the radiative transfer calculation.

Interpolation Method
--------------------
The interpolation method defines how we treat atmospheric properties inbetween the location grid points.
The method used above, `LinearInterpolation`, says that we should treat quantities such as extinction as varying
linearly between atmospheric grid points.  This usually comes into play when evaluating integrals of these quantities
in solving the radiative transfer equation.

Note that in some ways specifying the interpolation method is more of a suggestion than a rule.  Some radiative transfer
solution techniques, notably the discrete ordinates method, have to implicitly assume homogeneous constant shells, and are
thus incompatible with linearly varying optical properties.


Geometry Type
-------------
The geometry type is the global geometry for the radiative transfer calculation.  Think spherical, or plane parallel.

Available Geometries
--------------------
```{eval-rst}
.. autosummary::

    sasktran2.Geometry1D
```
