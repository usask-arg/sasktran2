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

    atmos = []

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmos.append(atmosphere)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    vmr_altitude_grid = np.array([10000.0, 30000.0, 60000.0])
    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(), vmr_altitude_grid, np.ones_like(vmr_altitude_grid) * 1e-6
    )

    atmos.append(atmosphere)

    scen = []

    for atmo in atmos:
        scen.append(
            {
                "config": config,
                "geometry": geometry,
                "viewing_geo": viewing_geo,
                "atmosphere": atmo,
            }
        )

    return scen


def test_pressure_wf():
    """
    Checks that the pressure weighting function is working correctly
    """

    scens = _test_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

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


def test_temperature_wf():
    """
    Checks that the temperature weighting function is working correctly
    """

    scens = _test_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere.temperature_k, 0.0001, engine, atmosphere, "wf_temperature_k"
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_temperature_k"],
            radiance["wf_temperature_k_numeric"],
            wf_dim="altitude",
            decimal=4,
        )


def test_specific_humidity_wf():
    """
    Checks that the specific_humidity weighting function is working correctly
    """

    scens = _test_scenarios()

    # Only the second scenario has specific humidity derivatives
    for scen in scens[1:]:
        atmosphere = scen["atmosphere"]

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])
        atmosphere.specific_humidity = np.zeros_like(atmosphere.temperature_k) + 0.01

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere.specific_humidity,
            0.1,
            engine,
            atmosphere,
            "wf_specific_humidity",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_specific_humidity"],
            radiance["wf_specific_humidity_numeric"],
            wf_dim="altitude",
            decimal=4,
        )
