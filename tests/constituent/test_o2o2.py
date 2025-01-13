from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _test_scenarios():
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 2

    altitude_grid = np.arange(0, 65001, 1000.0)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    for tan_alt in np.arange(10000, 60000, 2000):
        viewing_geo.add_ray(sk.TangentAltitudeSolar(tan_alt, 0, 600000, 0.6))

    wavel = np.array([310, 330, 350, 600])

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    scen = []

    scen.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geo,
            "atmosphere": atmosphere,
        }
    )

    return scen


def test_o2o2_construction():
    """
    Test that the CollisionInducedAbsorber class can be constructed.
    """
    sk.constituent.CollisionInducedAbsorber(
        sk.optical.HITRANCollision("O2O2"), name="O2O2"
    )


def test_o2o2_wf_temperature():
    """
    Tests that the CollisionInducedAbsorber class calculates derivatives with respect to temperature correctly.
    """

    scens = _test_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

        atmosphere["o2o2"] = sk.constituent.CollisionInducedAbsorber(
            sk.optical.HITRANCollision("O2O2"), name="O2O2"
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere.temperature_k, 0.00001, engine, atmosphere, "wf_temperature_k"
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_temperature_k"],
            radiance["wf_temperature_k_numeric"],
            wf_dim="altitude",
            decimal=4,
        )


def test_o2o2_wf_pressure():
    """
    Tests that the CollisionInducedAbsorber class calculates derivatives with respect to pressure correctly.
    """

    scens = _test_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

        atmosphere["o2o2"] = sk.constituent.CollisionInducedAbsorber(
            sk.optical.HITRANCollision("O2O2"), name="O2O2"
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere.pressure_pa, 0.01, engine, atmosphere, "wf_pressure_pa"
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_pressure_pa"],
            radiance["wf_pressure_pa_numeric"],
            wf_dim="altitude",
            decimal=4,
        )
