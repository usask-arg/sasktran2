import sasktran2

import sasktran2 as sk
import numpy as np


def test_atmosphere_construction():
    config = sk.Config()

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical)
    atmosphere = sk.Atmosphere(1, geometry, config, 1)

    storage = atmosphere.storage


def test_atmosphere_vector_construction():
    config = sk.Config()

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical)
    atmosphere = sk.Atmosphere(1, geometry, config, 3)

    storage = atmosphere.storage