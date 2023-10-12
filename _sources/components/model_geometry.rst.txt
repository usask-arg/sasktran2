.. model_geometry:

The Model Geometry
==================
The model geometry is the grid that the radiative transfer calculation is actually performed on.
Usually this involves specifying the coordinates (spherical, plane parallel, etc.) the number of dimensions the
atmosphere is allowed to vary in (1, 2, 3), as well as the actual grid values themselves.  Sometimes an internal
definition of a `reference point` is also required to provide context to the viewing geometry policies.

Similar to other SASKTRAN2 components, there are multiple ways to construct the model geometry.
The standard one, and the one you probably want, is :py:class:`sasktran2.Geometry1D`,

.. ipython:: python

    import sasktran2 as sk
    import numpy as np

    model_geometry = sk.Geometry1D(cos_sza=0.6,
                                   solar_azimuth=0,
                                   earth_radius_m=6372000,
                                   altitude_grid_m=np.arange(0, 100001, 1000),
                                   interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                   geometry_type=sk.GeometryType.Spherical)

