from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_amf_basic():
    """
    Verifies that the AMF constituent runs without issue.  We do not check the results here
    """
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates

    altitude_grid = np.arange(0, 65001, 1000.0)
    geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0,
        earth_radius_m=6327000,
        altitude_grid_m=altitude_grid,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()
    viewing_geo.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.6,
            relative_azimuth=0,
            cos_viewing_zenith=-0.8,
            observer_altitude_m=200000,
        )
    )

    engine = sk.Engine(config, geometry, viewing_geo)

    wavel = np.array([330])
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["amf"] = sk.constituent.AirMassFactor()
    atmosphere.surface.albedo[:] = 0.3

    _ = engine.calculate_radiance(atmosphere)
