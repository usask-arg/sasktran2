from __future__ import annotations

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
