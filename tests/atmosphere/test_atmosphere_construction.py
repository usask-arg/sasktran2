from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_atmosphere_construction():
    config = sk.Config()

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    _ = sk.Atmosphere(geometry, config, numwavel=17)


def test_atmosphere_vector_construction():
    config = sk.Config()
    config.num_stokes = 3

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(geometry, config, numwavel=17)

    _ = atmosphere.storage
