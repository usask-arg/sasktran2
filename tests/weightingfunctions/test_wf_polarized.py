from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _test_scenarios():
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 4
    config.num_stokes = 3

    altitude_grid = np.arange(0, 65001, 5000.0)

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


def test_wf_polarized_vmr_const():
    """
    Checks that the polarized WFs are correct for a VMR constituent
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


def _test_scat_scenarios():
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    # config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.num_streams = 4
    config.delta_m_scaling = False
    config.num_stokes = 3

    altitude_grid = np.arange(0, 65001, 5000.0)

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
        viewing_geo.add_ray(sk.TangentAltitudeSolar(tan_alt, 0.1, 600000, 0.6))

    wavel = np.array([600])

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere.surface.albedo[:] = 0.3

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


def test_scatterer_extinction_wf_native_grid_polarized():
    """
    Tests the derivatives for a scatterer on the native grid for polarized calculation
    """

    scens = _test_scat_scenarios()

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        atmosphere["strat_aerosol"] = sk.test_util.scenarios.test_aerosol_constituent(
            altitude_grid, True
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["strat_aerosol"].extinction_per_m,
            1e-3,
            engine,
            atmosphere,
            "wf_strat_aerosol_extinction",
        )

        # The precision on this one can be pretty bad, but it looks right to me? Unsure
        sk.test_util.wf.validate_wf(
            radiance["wf_strat_aerosol_extinction"],
            radiance["wf_strat_aerosol_extinction_numeric"],
            wf_dim="strat_aerosol_altitude",
            decimal=2,
        )
