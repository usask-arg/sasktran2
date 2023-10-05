import sasktran2 as sk
import numpy as np
from sasktran2.climatology.us76 import add_us76_standard_atmosphere
from sasktran2.constituent.rayleigh import Rayleigh


def test_scalar_full_chain():
    config = sk.Config()

    altitude_grid = np.arange(0, 65001, 1000)

    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical)

    viewing_geo = sk.ViewingGeometry()

    ray = sk.TangentAltitudeSolar(10000, 0, 600000, 0.6)

    viewing_geo.add_ray(ray)

    engine = sk.Engine(config, geometry, viewing_geo, 1)

    wavel = np.arange(280.0, 350.0, 0.001)

    atmosphere = sk.Atmosphere(len(wavel), geometry, config, 1)
    atmosphere.wavelengths_nm = wavel

    add_us76_standard_atmosphere(atmosphere)

    atmosphere['rayleigh'] = Rayleigh()

    output = engine.calculate_radiance(atmosphere)

    pass