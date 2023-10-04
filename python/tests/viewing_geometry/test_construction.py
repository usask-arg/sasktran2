
import sasktran2 as sk
import numpy as np


def test_viewing_container():
    viewing_geo = sk.ViewingGeometry()

    ray = sk.TangentAltitudeSolar(10000, 0, 600000, 0.6)

    viewing_geo.add_ray(ray)
    