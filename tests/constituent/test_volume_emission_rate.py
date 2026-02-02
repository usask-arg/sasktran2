from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _test_scenarios():
    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate

    geometry = sk.Geometry1D(
        cos_sza=-0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.arange(0, 120001, 1000),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    for alt in [95000]:
        ray = sk.TangentAltitudeSolar(
            tangent_altitude_m=alt,
            relative_azimuth=0,
            observer_altitude_m=200000,
            cos_sza=-0.6,
        )
        viewing_geo.add_ray(ray)

    for ele in [0, 10, 20, 30, 40, 50, 60]:
        ray = sk.SolarAnglesObserverLocation(
            -0.6, 0.0, np.cos(np.radians(90 - ele)), 30000
        )

        viewing_geo.add_ray(ray)

    wavel = np.arange(280.0, 800.0, 0.5)
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=wavel, finite_resolution_mode=False
    )

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    altitude = (
        np.array([140, 135, 130, 125, 120, 115, 110, 105, 100, 95, 90, 85]).astype(
            np.float64
        )[::-1]
        * 1000
    )

    # VER in ph cm^-3 s^-1 (red dots)
    VER = np.array(
        [400, 600, 800, 1100, 1500, 2000, 2600, 3200, 3800, 4300, 3600, 900]
    )[::-1]

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["oxygen_green"] = sk.constituent.MonochromaticVolumeEmissionRate(
        altitude, VER, 557.7
    )

    scen = []

    scen.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geo,
            "atmosphere": atmosphere,
        }
    )

    return scen


def test_ver_wf_native_grid():
    """
    Tests that the VMRAltitudeAbsorber class calculates derivatives with respect to VMR correctly when
    VMR is specified on the native grid
    """

    scens = _test_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["oxygen_green"].ver,
            0.0001,
            engine,
            atmosphere,
            "wf_oxygen_green_ver",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_oxygen_green_ver"],
            radiance["wf_oxygen_green_ver_numeric"],
            wf_dim="oxygen_green_altitude",
            decimal=5,
        )
