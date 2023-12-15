---
file_format: mystnb
---

(_users_geodetic)=
# Geodesy Tools
```{code-cell}
import sasktran2 as sk
```

SASKTRAN2 contains a simple implementation of ellipsoidal implementations.
These can be created with the generic {py:class}`sasktran2.Geodetic` object
which allows the creation of ellipsoids with arbitrary axis lengths and/or
flattening factors.  Typically it is recommended to use one of the built in
specializations,

```{eval-rst}
.. autosummary::
    sasktran2.WGS84
    sasktran2.SphericalGeoid
```

In this guide we will use the {py:class}`sasktran2.WGS84` object.

```{code-cell}
geodetic = sk.WGS84()
```

## Usage
In the beginning, the geodetic object is in an unitialized state,
```{code-cell}
print(geodetic)
```
to initialize the object, we use one of the available `from_*` methods,

```{code-cell}
geodetic.from_lat_lon_alt(latitude=60, longitude=120, altitude=10000)
print(geodetic)
```

Then we can do things such as convert to geocentric cartesian coordinates,

```{code-cell}
print(geodetic.location)
```

or access the local right handed coordinate system pointing (up, south, west),

```{code-cell}
print(f"up: {geodetic.local_up}, south: {geodetic.local_south}, west: {geodetic.local_west}")
```

## Initialization Methods
```{eval-rst}
.. autosummary::
    sasktran2.Geodetic.from_xyz
    sasktran2.Geodetic.from_lat_lon_alt
    sasktran2.Geodetic.from_tangent_point
    sasktran2.Geodetic.from_tangent_altitude
```
