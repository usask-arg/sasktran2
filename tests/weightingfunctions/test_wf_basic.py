from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_wf_basic():
    """
    Runs a full calculation across multiple lines of sight and wavelengths, we don't verify the results, just check the chain
    """
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 2

    altitude_grid = np.arange(0, 65001, 1000.0)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    for tan_alt in np.arange(10000, 60000, 2000):
        viewing_geo.add_ray(sk.TangentAltitudeSolar(tan_alt, 0, 600000, 0.6))

    engine = sk.Engine(config, geometry, viewing_geo)

    wavel = np.arange(280.0, 800.0, 1)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
    )

    engine.calculate_radiance(atmosphere)
