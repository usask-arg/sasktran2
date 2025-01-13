from __future__ import annotations

import numpy as np
import pandas as pd
import sasktran2 as sk


def test_limb_viewing_ray_conversion():
    time = pd.Timestamp.now()

    solar = sk.solar.SolarGeometryHandlerForced(45, 180)

    tan_alts = np.array([10, 20, 30, 40]) * 1000

    obs_geo = sk.WGS84()
    obs_geo.from_lat_lon_alt(20, 30, 700000)

    tp_geo = sk.WGS84()

    for tan_alt in tan_alts:
        lv = tp_geo.from_tangent_altitude(tan_alt, obs_geo.location, [0, 0, 1])

        ray = sk.viewinggeo.ecef.ecef_to_sasktran2_ray(
            obs_geo.location, lv, time, sk.WGS84(), solar
        )

        assert isinstance(ray, sk.TangentAltitudeSolar)

        string = str(ray).split(": ", 1)[1]

        # Split into individual parameter pairs
        pairs = string.split(", ")

        # Convert to dictionary
        param_dict = {}
        for pair in pairs:
            key, value = pair.split(": ")
            param_dict[key.strip()] = float(value)

        np.testing.assert_allclose(param_dict["cos_sza"], np.cos(np.deg2rad(45)))
        np.testing.assert_allclose(
            param_dict["relative_azimuth_angle"], np.deg2rad(180)
        )
        np.testing.assert_allclose(
            param_dict["tangentaltitude"], tan_alt, atol=0.1
        )  # Check to accurate within 0.1 m
        np.testing.assert_allclose(param_dict["observeraltitude"], 700000)
