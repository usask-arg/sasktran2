from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def _default_settings():
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

    atmosphere.storage.ssa[:] = 0.9
    atmosphere.storage.total_extinction[:] = 0.1

    atmosphere.leg_coeff.a1[0, :, :] = 1

    return config, model_geometry, viewing_geo, atmosphere


def test_nan_extinction():
    config, model_geometry, viewing_geo, atmosphere = _default_settings()

    atmosphere.storage.total_extinction[10, 10] = np.nan

    engine = sk.Engine(config, model_geometry, viewing_geo)

    try:
        _ = engine.calculate_radiance(atmosphere)
        pytest.fail("Expected an exception")
    except Exception:
        pass


def test_negative_extinction():
    config, model_geometry, viewing_geo, atmosphere = _default_settings()

    atmosphere.storage.total_extinction[0, 0] = -0.1

    engine = sk.Engine(config, model_geometry, viewing_geo)

    try:
        _ = engine.calculate_radiance(atmosphere)
        pytest.fail("Expected an exception")
    except Exception:
        pass


def test_nan_ssa():
    config, model_geometry, viewing_geo, atmosphere = _default_settings()

    atmosphere.storage.ssa[10, 10] = np.nan

    engine = sk.Engine(config, model_geometry, viewing_geo)

    try:
        _ = engine.calculate_radiance(atmosphere)
        pytest.fail("Expected an exception")
    except Exception:
        pass


def test_negative_ssa():
    config, model_geometry, viewing_geo, atmosphere = _default_settings()

    atmosphere.storage.ssa[0, 0] = -0.1

    engine = sk.Engine(config, model_geometry, viewing_geo)

    try:
        _ = engine.calculate_radiance(atmosphere)
        pytest.fail("Expected an exception")
    except Exception:
        pass
