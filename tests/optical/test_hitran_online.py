from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def test_hitran_online():
    pytest.importorskip("hapi")

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

    wavel = np.arange(1360, 1370, 0.01)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geometry.altitudes(),
        np.ones_like(model_geometry.altitudes()) * 1e-6,
    )

    atmosphere["h2o"] = sk.climatology.mipas.constituent(
        "H2O", sk.optical.HITRANAbsorber("H2O")
    )

    engine = sk.Engine(config, model_geometry, viewing_geo)

    _ = engine.calculate_radiance(atmosphere)
