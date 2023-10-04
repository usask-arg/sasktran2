import sasktran2 as sk
import numpy as np


def test_scalar_full_chain():
    config = sk.Config()

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical)

    viewing_geo = sk.ViewingGeometry()

    ray = sk.TangentAltitudeSolar(10000, 0, 600000, 0.6)

    viewing_geo.add_ray(ray)

    engine = sk.Engine(config, geometry, viewing_geo, 1)

    atmosphere = sk.Atmosphere(1, geometry, config, 1)

    output = engine.calculate_radiance(atmosphere)