from __future__ import annotations

import numpy as np
import sasktran2 as sk
import xarray as xr
from scipy.special import roots_legendre


def test_upwelling_flux():
    alts = [64750, 30000, 65000, 999.99]

    for alt in alts:
        config = sk.Config()
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        config.num_streams = 10
        config.num_forced_azimuth = 1

        model_geometry = sk.Geometry1D(
            cos_sza=0.6,
            solar_azimuth=0,
            earth_radius_m=6372000,
            altitude_grid_m=np.arange(0, 65001, 1000.0),
            interpolation_method=sk.InterpolationMethod.LinearInterpolation,
            geometry_type=sk.GeometryType.PlaneParallel,
        )

        viewing_geo = sk.ViewingGeometry()

        viewing_geo.add_flux_observer(sk.FluxObserverSolar(0.6, alt))

        nodes, weights = roots_legendre(config.num_streams / 2)
        nodes = 0.5 * nodes + 0.5
        weights *= 0.5

        for n in nodes:
            viewing_geo.add_ray(sk.SolarAnglesObserverLocation(0.6, 0.0, -n, alt))

        integrand = xr.DataArray(data=weights * nodes * 2.0 * np.pi, dims=["los"])

        wavel = np.arange(300.0, 800.0, 10)
        atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

        sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

        atmosphere["rayleigh"] = sk.constituent.Rayleigh()

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(),
            model_geometry.altitudes(),
            np.ones_like(model_geometry.altitudes()) * 1e-6,
        )

        atmosphere["surface"] = sk.constituent.LambertianSurface(0.3)

        engine = sk.Engine(config, model_geometry, viewing_geo)

        rad = engine.calculate_radiance(atmosphere)

        rad["upwelling_flux_numeric"] = rad["radiance"].isel(stokes=0) @ integrand

        np.testing.assert_allclose(
            rad["upwelling_flux"].isel(flux_location=0).values,
            rad["upwelling_flux_numeric"].values,
            rtol=2e-8,
            atol=0.0,
        )

        assert "upwelling_flux" in rad
        assert "downwelling_flux" in rad


def test_upwelling_flux_thermal():
    alts = [64750, 30000, 65000, 999.99]

    for alt in alts:
        config = sk.Config()
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        config.emission_source = sk.EmissionSource.DiscreteOrdinates
        config.num_streams = 10
        config.num_forced_azimuth = 1

        model_geometry = sk.Geometry1D(
            cos_sza=0.6,
            solar_azimuth=0,
            earth_radius_m=6372000,
            altitude_grid_m=np.arange(0, 65001, 1000.0),
            interpolation_method=sk.InterpolationMethod.LinearInterpolation,
            geometry_type=sk.GeometryType.PlaneParallel,
        )

        viewing_geo = sk.ViewingGeometry()

        viewing_geo.add_flux_observer(sk.FluxObserverSolar(0.6, alt))

        nodes, weights = roots_legendre(config.num_streams / 2)
        nodes = 0.5 * nodes + 0.5
        weights *= 0.5

        for n in nodes:
            viewing_geo.add_ray(sk.SolarAnglesObserverLocation(0.6, 0.0, -n, alt))

        integrand = xr.DataArray(data=weights * nodes * 2.0 * np.pi, dims=["los"])

        wavel = np.arange(7370, 7380, 0.1)
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

        atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)
        atmosphere["ch4"] = sk.climatology.mipas.constituent("CH4", hitran_db)

        sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

        atmosphere["rayleigh"] = sk.constituent.Rayleigh()

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(),
            model_geometry.altitudes(),
            np.ones_like(model_geometry.altitudes()) * 1e-6,
        )

        atmosphere["surface"] = sk.constituent.LambertianSurface(0.3)
        atmosphere["thermal"] = sk.constituent.ThermalEmission()
        atmosphere["surface_thermal"] = sk.constituent.SurfaceThermalEmission(290, 0.9)
        atmosphere["irradiance"] = sk.constituent.SolarIrradiance()

        engine = sk.Engine(config, model_geometry, viewing_geo)

        rad = engine.calculate_radiance(atmosphere)

        rad["upwelling_flux_numeric"] = rad["radiance"].isel(stokes=0) @ integrand

        np.testing.assert_allclose(
            rad["upwelling_flux"].isel(flux_location=0).values,
            rad["upwelling_flux_numeric"].values,
            rtol=2e-8,
            atol=0.0,
        )

        assert "upwelling_flux" in rad
        assert "downwelling_flux" in rad


def test_actinic_flux_output():
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.flux_types = [sk.FluxType.Actinic]
    config.num_streams = 4
    config.num_forced_azimuth = 1
    config.wavelength_batch_size = 4

    model_geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.arange(0, 30001, 10000.0),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.PlaneParallel,
    )

    viewing_geo = sk.ViewingGeometry()
    viewing_geo.add_flux_observer(sk.FluxObserverSolar(0.6, 20000.0))

    atmosphere = sk.Atmosphere(
        model_geometry,
        config,
        wavelengths_nm=np.array([310.0, 330.0, 350.0]),
    )
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["solar"] = sk.constituent.SolarIrradiance()

    engine = sk.Engine(config, model_geometry, viewing_geo)
    rad = engine.calculate_radiance(atmosphere)

    assert "actinic_flux" in rad
    assert "upwelling_flux" not in rad
    assert "downwelling_flux" not in rad
    assert rad["actinic_flux"].shape == (3, 1)
    assert np.all(np.isfinite(rad["actinic_flux"].values))
    assert np.any(rad["actinic_flux"].values > 0.0)


def test_direct_beam_flux_is_attenuated_to_the_observer():
    # A pure absorber with constant extinction, so the optical depth above
    # any altitude is exactly extinction * (top - altitude). The flux is then
    # only the direct beam, at grid levels and between them.
    extinction_per_m = 2.0e-5
    altitude_grid = np.arange(0.0, 50001.0, 1000.0)
    observers = np.array([0.0, 2500.0, 10000.0, 17300.0, 30000.0, 49500.0])

    for cos_sza in (1.0, 0.5):
        config = sk.Config()
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        config.flux_types = [sk.FluxType.Actinic, sk.FluxType.Downwelling]
        config.num_streams = 4
        config.num_forced_azimuth = 1

        geometry = sk.Geometry1D(
            cos_sza,
            0.0,
            6371000.0,
            altitude_grid,
            sk.InterpolationMethod.LinearInterpolation,
            sk.GeometryType.PlaneParallel,
        )
        viewing = sk.ViewingGeometry()
        for altitude in observers:
            viewing.add_flux_observer(sk.FluxObserverSolar(cos_sza, altitude))

        atmosphere = sk.Atmosphere(
            geometry, config, np.array([500.0]), calculate_derivatives=False
        )
        atmosphere.storage.total_extinction[:] = extinction_per_m
        atmosphere.storage.ssa[:] = 0.0
        atmosphere.storage.solar_irradiance[:] = 1.0
        atmosphere.leg_coeff.a1[0] = 1.0

        rad = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

        transmission = np.exp(
            -extinction_per_m * (altitude_grid[-1] - observers) / cos_sza
        )
        np.testing.assert_allclose(
            rad["actinic_flux"].to_numpy()[0], transmission, rtol=1e-10
        )
        np.testing.assert_allclose(
            rad["downwelling_flux"].to_numpy()[0],
            cos_sza * transmission,
            rtol=1e-10,
        )
