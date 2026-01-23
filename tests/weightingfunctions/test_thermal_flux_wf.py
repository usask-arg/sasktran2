from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _raw_scenarios() -> list:
    """ """
    configs = []

    configs.append(sk.Config())
    configs[-1].multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    configs[-1].num_streams = 4
    configs[-1].single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    configs[-1].emission_source = sk.EmissionSource.DiscreteOrdinates
    configs[-1].num_stokes = 1

    geometrys = []
    altitude_grid = np.arange(0, 65001, 5000.0)

    geometrys.append(
        sk.Geometry1D(
            0.6,
            0,
            6327000,
            altitude_grid,
            sk.InterpolationMethod.LinearInterpolation,
            sk.GeometryType.PlaneParallel,
        )
    )

    viewing_geos = []
    viewing_geos.append(sk.ViewingGeometry())

    # for alt in np.arange(10000, 60000, 2000):
    #    viewing_geos[-1].add_flux_observer(sk.FluxObserverSolar(0.6, alt))
    viewing_geos[-1].add_flux_observer(sk.FluxObserverSolar(0.6, 0))
    viewing_geos[-1].add_flux_observer(sk.FluxObserverSolar(0.6, 65000))
    scen = []

    wavelengths = np.arange(7370, 7380, 0.01)
    wavelengths = np.array([7340.0])
    hitran_db = sk.optical.database.OpticalDatabaseGenericAbsorber(
        sk.database.StandardDatabase().path(
            "hitran/CH4/sasktran2/60fc6c547a1b2e1181f3296dead288d84fa7c178.nc"
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

    atmosphere = sk.Atmosphere(geometrys[-1], configs[-1], wavelengths_nm=wavelengths)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["ch4"] = sk.climatology.mipas.constituent("CH4", hitran_db)
    atmosphere["emission"] = sk.constituent.ThermalEmission()
    atmosphere["surface_emission"] = sk.constituent.SurfaceThermalEmission(300, 0.9)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["solar_irradiance"] = sk.constituent.SolarIrradiance()

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    for config in configs:
        for geometry in geometrys:
            for viewing_geo in viewing_geos:
                scen.append(
                    {
                        "config": config,
                        "geometry": geometry,
                        "viewing_geo": viewing_geo,
                        "atmosphere": atmosphere,
                    }
                )

    return scen


def test_thermal_flux_wf_temperature_with_emission():
    """
    Checks that the WFs are correct for a VMR constituent When emissions are present
    """

    scens = _raw_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere.temperature_k,
            0.01,
            engine,
            atmosphere,
            "wf_temperature_k",
            calc_vars=["radiance", "upwelling_flux", "downwelling_flux"],
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_temperature_k"],
            radiance["wf_temperature_k_numeric"],
            wf_dim="altitude",
            decimal=5,
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_temperature_k_upwelling_flux"].isel(flux_location=-1),
            radiance["wf_temperature_k_upwelling_flux_numeric"].isel(flux_location=-1),
            wf_dim="altitude",
            decimal=3,
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_temperature_k_downwelling_flux"].isel(flux_location=0),
            radiance["wf_temperature_k_downwelling_flux_numeric"].isel(flux_location=0),
            wf_dim="altitude",
            decimal=3,
        )
