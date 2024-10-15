import numpy as np
from sasktran2.cli.fer import FERGenerator


def test_full():
    acquisitions = {
        "acq_1": {
            "type": "limb_vertical_image",
            "time": {
                "year": 2024,
                "month": 1,
                "day": 1,
            },
            "tangent_altitudes": np.arange(5.0, 60.0, 1.0) * 1000,
            "tangent_latitude": 20,
            "tangent_longitude": 15,
            "inst_loc": {
                "type": "inferred",
                "altitude": 700000,
            },
            "solar_geometry": {
                "type": "forced",
                "solar_zenith": 60,
                "solar_azimuth": 0,
            },
        }
    }

    atmosphere = {
        "atmospheric_state": {
            "type": "us76",
            "include_rayleigh_scattering": True,
            "include_o2o2_absorption": False,
        },
        "constituents": {"o3": {"type": "mipas_climatology"}},
    }

    gen = FERGenerator()

    gen.generate(acquistions=acquisitions, atmosphere=atmosphere)
