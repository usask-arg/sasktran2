from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_fascode():
    fascode_valid = [
        "H2O",
        "CO2",
        "O3",
        "N2O",
        "CO",
        "CH4",
        "O2",
        "NO",
        "SO2",
        "NO2",
        "NH3",
        "HNO3",
        "OH",
        "HF",
        "HCl",
        "HBr",
        "HI",
        "ClO",
        "OCS",
        "H2CO",
        "HOCl",
        "N2",
        "HCN",
        "CH3Cl",
        "H2O2",
        "C2H2",
        "C2H6",
        "PH3",
        "COF2",
        "SF6",
        "H2S",
        "CFCl3",
        "CF2Cl2",
        "CClF3",
        "CF4",
        "CHCl2F",
        "CHClF2",
        "C2Cl3F3",
        "C2Cl2F4",
        "C2ClF5",
        "CCl4",
        "ClONO2",
        "N2O5",
        "HNO4",
    ]

    fascode_states = ["tro", "mls", "mlw", "sas", "saw", "std"]

    for fascode in fascode_valid:
        for fascode_state in fascode_states:
            sk.climatology.mipas.constituent(
                fascode, None, dataset="fascode", climatology=fascode_state
            )


def test_mipas1998():
    mipas_valid = [
        "N2",
        "O2",
        "O3P",
        "CO2",
        "O3",
        "H2O",
        "CH4",
        "N2O",
        "HNO3",
        "CO",
        "NO2",
        "N2O5",
        "ClO",
        "HOCl",
        "ClONO2",
        "NO",
        "HCN",
        "H2O2",
        "F12",
        "F14",
        "F22",
        "COF2",
        "OCS",
        "NH3",
        "SO2",
        "CFCl3",
        "C2H2",
        "C2H6",
        "CCl4",
        "SF6",
        "HNO4",
        "CH3Cl",
        "CClF3",
        "CHCl2F",
        "C2Cl3F3",
        "C2Cl2F4",
    ]

    mipas_states = ["day_imk", "ngt_imk", "win_imk", "sum_imk"]

    for mipas in mipas_valid:
        for mipas_state in mipas_states:
            sk.climatology.mipas.constituent(
                mipas, None, dataset="mipas_1998", climatology=mipas_state
            )


def test_mipas2001():
    mipas_valid = [
        "N2",
        "O2",
        "CO2",
        "O3",
        "H2O",
        "CH4",
        "N2O",
        "HNO3",
        "CO",
        "NO2",
        "N2O5",
        "ClO",
        "HOCl",
        "ClONO2",
        "NO",
        "HNO4",
        "HCN",
        "NH3",
        "F11",
        "F12",
        "F14",
        "F22",
        "CCl4",
        "COF2",
        "H2O2",
        "C2H2",
        "C2H6",
        "OCS",
        "SO2",
        "SF6",
        "CClF3",
        "CHCl2F",
        "C2Cl3F3",
        "C2Cl2F4",
        "C2ClF5",
        "CH3Cl",
        "H2S",
    ]

    mipas_states = ["day", "ngt", "win", "sum", "equ"]

    for mipas in mipas_valid:
        for mipas_state in mipas_states:
            sk.climatology.mipas.constituent(
                mipas, None, dataset="mipas_2001", climatology=mipas_state
            )


def test_add_to_atmosphere():
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 2

    model_geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.arange(0, 65001, 1000),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    for alt in [10000, 20000, 30000, 40000]:
        ray = sk.TangentAltitudeSolar(
            tangent_altitude_m=alt,
            relative_azimuth=0,
            observer_altitude_m=200000,
            cos_sza=0.6,
        )
        viewing_geo.add_ray(ray)

    wavel = np.arange(280.0, 800.0, 10)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    sk.climatology.mipas.add_to_atmosphere(
        atmosphere, {"O3": sk.optical.O3DBM(), "NO2": sk.optical.NO2Vandaele()}
    )

    engine = sk.Engine(config, model_geometry, viewing_geo)

    _ = engine.calculate_radiance(atmosphere)
