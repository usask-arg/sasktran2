
import numpy as np
import sasktran2 as sk


def test_viewing_container():
    viewing_geo = sk.ViewingGeometry()

    ray = sk.TangentAltitudeSolar(10000, 0, 600000, 0.6)

    viewing_geo.add_ray(ray)
