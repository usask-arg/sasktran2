from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def test_viewing_container():
    """
    Tests that we can create each different type of viewing ray and add them to the container without issue
    """
    viewing_geo = sk.ViewingGeometry()

    ray = sk.TangentAltitudeSolar(10000, 0, 600000, 0.6)

    viewing_geo.add_ray(ray)

    ray = sk.GroundViewingSolar(0.6, 0, -0.8, 200000)

    viewing_geo.add_ray(ray)

    assert len(viewing_geo.observer_rays) == 2


def test_geometry_relative_tangent_ray_is_public_and_introspectable():
    ray = sk.TangentAltitude(
        tangent_altitude_m=20_000.0,
        observer_altitude_m=200_000.0,
        horizontal_angle_radians=0.25,
        viewing_azimuth_radians=np.pi / 3.0,
    )
    viewing_geo = sk.ViewingGeometry()
    viewing_geo.add_ray(ray)

    assert viewing_geo.observer_rays == [ray]
    assert ray.tangent_altitude_m == 20_000.0
    assert ray.observer_altitude_m == 200_000.0
    assert ray.horizontal_angle_radians == 0.25
    assert ray.viewing_azimuth_radians == np.pi / 3.0
    assert "horizontal_angle_radians: 0.25" in repr(ray)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tangent_altitude_m": -1.0},
        {"tangent_altitude_m": np.nan},
        {"observer_altitude_m": 10_000.0},
        {"observer_altitude_m": np.inf},
        {"horizontal_angle_radians": np.nan},
        {"viewing_azimuth_radians": np.inf},
    ],
)
def test_geometry_relative_tangent_ray_rejects_invalid_inputs(kwargs):
    parameters = {
        "tangent_altitude_m": 20_000.0,
        "observer_altitude_m": 200_000.0,
        "horizontal_angle_radians": 0.0,
        "viewing_azimuth_radians": 0.0,
    }
    parameters.update(kwargs)

    with pytest.raises(ValueError, match="must be"):
        sk.TangentAltitude(**parameters)
