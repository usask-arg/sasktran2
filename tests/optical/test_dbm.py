import sasktran2 as sk
import numpy as np
from sasktran2.climatology.us76 import add_us76_standard_atmosphere



def test_dbm():
    dbm = sk.optical.O3DBM()

    config = sk.Config()

    altitude_grid = np.arange(0, 65001, 1000)

    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical)

    viewing_geo = sk.ViewingGeometry()

    ray = sk.TangentAltitudeSolar(10000, 0, 600000, 0.6)

    viewing_geo.add_ray(ray)

    wavel = np.arange(280.0, 350.0, 0.001)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    dbm.atmosphere_quantities(atmo=atmosphere)
