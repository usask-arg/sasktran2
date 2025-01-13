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


def test_vmr_altitude_construction():
    """
    Test that the VMRAltitude class can be constructed
    """
    alts = np.arange(0, 100001, 1000)

    sk.constituent.VMRAltitudeAbsorber(sk.optical.O3DBM(), alts, np.ones_like(alts))


def test_vmr_altitude_wf_native_grid():
    """
    Tests that the VMRAltitudeAbsorber class calculates derivatives with respect to VMR correctly when
    VMR is specified on the native grid
    """

    scens = _test_scenarios()

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["ozone"].vmr, 0.0001, engine, atmosphere, "wf_ozone_vmr"
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_ozone_vmr"],
            radiance["wf_ozone_vmr_numeric"],
            wf_dim="ozone_altitude",
            decimal=5,
        )


def test_vmr_altitude_wf_interpolated_grid():
    """
    Tests that the VMRAltitudeAbsorber class calculates derivatives with respect to VMR correctly when
    VMR is specified on an interpolated grid
    """

    scens = _test_scenarios()

    for scen in scens:
        altitude_grid = np.array([10000.0, 30000.0, 60000.0])
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["ozone"].vmr, 0.001, engine, atmosphere, "wf_ozone_vmr"
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_ozone_vmr"],
            radiance["wf_ozone_vmr_numeric"],
            wf_dim="ozone_altitude",
            decimal=5,
        )
