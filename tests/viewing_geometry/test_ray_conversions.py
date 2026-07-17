from __future__ import annotations

import numpy as np
import pandas as pd
import sasktran2 as sk


def ray_parameters(ray) -> dict[str, float]:
    string = str(ray).split(": ", 1)[1]
    return {
        key.strip(): float(value)
        for pair in string.split(", ")
        for key, value in [pair.split(": ")]
    }


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

        param_dict = ray_parameters(ray)

        np.testing.assert_allclose(param_dict["cos_sza"], np.cos(np.deg2rad(45)))
        np.testing.assert_allclose(
            param_dict["relative_azimuth_angle"], np.deg2rad(180)
        )
        np.testing.assert_allclose(
            param_dict["tangentaltitude"], tan_alt, atol=0.1
        )  # Check to accurate within 0.1 m
        np.testing.assert_allclose(param_dict["observeraltitude"], 700000)


def test_ground_viewing_ray_conversion():
    time = pd.Timestamp("2020-01-01")
    solar_zenith = 48.0
    solar_azimuth = 137.0
    solar = sk.solar.SolarGeometryHandlerForced(solar_zenith, solar_azimuth)

    observer_geo = sk.WGS84()
    observer_geo.from_lat_lon_alt(20.0, 30.0, 700000.0)

    target_geo = sk.WGS84()
    target_geo.from_lat_lon_alt(18.0, 34.0, 0.0)

    look_vector = target_geo.location - observer_geo.location
    look_vector /= np.linalg.norm(look_vector)

    # A scaled look vector verifies that the conversion normalizes its input.
    ray = sk.viewinggeo.ecef.ecef_to_sasktran2_ray(
        observer_geo.location,
        3.0 * look_vector,
        time,
        sk.WGS84(),
        solar,
    )

    assert isinstance(ray, sk.GroundViewingSolar)
    param_dict = ray_parameters(ray)

    intercept_geo = sk.WGS84()
    intercept = intercept_geo.altitude_intercepts(
        0.0, observer_geo.location, look_vector
    )[0]
    intercept_geo.from_xyz(intercept)

    expected_viewing_azimuth = -np.rad2deg(
        np.arctan2(
            np.dot(look_vector, intercept_geo.local_west),
            -np.dot(look_vector, intercept_geo.local_south),
        )
    )

    np.testing.assert_allclose(param_dict["cos_sza"], np.cos(np.deg2rad(solar_zenith)))
    np.testing.assert_allclose(
        param_dict["relative_azimuth_angle"],
        np.deg2rad(solar_azimuth - expected_viewing_azimuth),
    )
    np.testing.assert_allclose(
        param_dict["cos_viewing_zenith"],
        -np.dot(look_vector, intercept_geo.local_up),
    )
    assert param_dict["cos_viewing_zenith"] > 0
    np.testing.assert_allclose(param_dict["observer_altitude_m"], 700000.0)


def test_ray_conversion_uses_default_solar_handler():
    time = pd.Timestamp("2020-01-01")

    observer_geo = sk.WGS84()
    observer_geo.from_lat_lon_alt(20.0, 30.0, 700000.0)

    target_geo = sk.WGS84()
    target_geo.from_lat_lon_alt(20.0, 30.0, 0.0)

    look_vector = target_geo.location - observer_geo.location
    ray = sk.viewinggeo.ecef.ecef_to_sasktran2_ray(
        observer_geo.location, look_vector, time
    )

    assert isinstance(ray, sk.GroundViewingSolar)
    param_dict = ray_parameters(ray)
    np.testing.assert_allclose(param_dict["cos_sza"], 1.0)
    np.testing.assert_allclose(param_dict["relative_azimuth_angle"], 0.0)
