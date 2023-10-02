import sasktran2

import sasktran2 as sk
import numpy as np


def test_atmosphere_construction():
    config = sk.Config()

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sasktran2.LinearInterpolation, sasktran2.Spherical)

    atmosphere = sk.AtmosphereScalar(1, geometry, config, False)

    storage = atmosphere.storage


def test_atmosphere_vector_construction():
    config = sk.Config()

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sasktran2.LinearInterpolation, sasktran2.Spherical)

    atmosphere = sk.AtmosphereVector(1, geometry, config, False)

    storage = atmosphere.storage