import sasktran2 as sk
import numpy as np
from sasktran2.climatology.us76 import add_us76_standard_atmosphere
from sasktran2.constituent.rayleigh import Rayleigh
from sasktran2.constituent.vmraltitudeabsorber import VMRAltitudeAbsorber
from sasktran2.optical import O3DBM


def test_wf_basic():
    config = sk.Config()
    config.num_threads = 8

    altitude_grid = np.arange(0, 65001, 1000.0)

    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical)

    viewing_geo = sk.ViewingGeometry()

    for tan_alt in np.arange(10000, 60000, 2000):
        viewing_geo.add_ray(sk.TangentAltitudeSolar(tan_alt, 0, 600000, 0.6))

    engine = sk.Engine(config, geometry, viewing_geo, 1)

    wavel = np.arange(280.0, 800.0, 1)

    atmosphere = sk.Atmosphere(len(wavel), geometry, config, 1)
    atmosphere.wavelengths_nm = wavel

    add_us76_standard_atmosphere(atmosphere)

    atmosphere['rayleigh'] = Rayleigh()

    atmosphere['ozone'] = VMRAltitudeAbsorber(O3DBM(),
                                              altitude_grid,
                                              np.ones_like(altitude_grid)*1e-6
                                              )

    radiance_base = engine.calculate_radiance(atmosphere)

    #dvmr = 0.001e-6
    #atmosphere['ozone'].vmr[30] += dvmr

    #radiance_perturb = engine.calculate_radiance(atmosphere)

    #wf_numeric = (radiance_perturb['radiance'] - radiance_base['radiance']) / dvmr

    #pass