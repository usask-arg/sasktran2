from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _test_scenarios():
    config = sk.Config()
    config.emission_source = sk.EmissionSource.Standard
    config.single_scatter_source = sk.SingleScatterSource.NoSource

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

    wavelengths = np.arange(7370, 7380, 0.01)
    hitran_db = sk.optical.database.OpticalDatabaseGenericAbsorber(
        sk.database.StandardDatabase().path(
            "hitran/CH4/sasktran2/6f006bcb051f81fc57d1bd09315589bfe77b4348.nc"
        )
    )

    hitran_db = sk.database.HITRANDatabase(
        molecule="CH4",
        start_wavenumber=1355,
        end_wavenumber=1357,
        wavenumber_resolution=0.01,
        reduction_factor=1,
        backend="sasktran2",
        profile="voigt",
    )

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["ch4"] = sk.climatology.mipas.constituent("CH4", hitran_db)
    atmosphere["emission"] = sk.constituent.ThermalEmission()

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


def test_wf_extinction_ssa_with_emission():
    """
    Checks that the WFs are correct for a VMR constituent When emissions are present
    """

    scens = _test_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["ch4"].vmr, 0.0001, engine, atmosphere, "wf_ch4_vmr"
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_ch4_vmr"],
            radiance["wf_ch4_vmr_numeric"],
            wf_dim="ch4_altitude",
            decimal=5,
        )
