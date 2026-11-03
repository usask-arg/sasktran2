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


def _plume_flux_scenario(
    single_scatter_source, flux_alts, dz=1000.0, cos_sza=0.8, ssa=1.0
):
    """
    Rayleigh atmosphere with a Gaussian Henyey-Greenstein aerosol plume (AOD 1 at 4 km),
    plane parallel, with flux observers at flux_alts.
    """
    config = sk.Config()
    config.single_scatter_source = single_scatter_source
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.flux_types = [
        sk.FluxType.Upwelling,
        sk.FluxType.Downwelling,
        sk.FluxType.Actinic,
    ]
    config.num_streams = 8
    config.num_forced_azimuth = 1

    alt_grid = np.arange(0.0, 30001.0, dz)
    model_geometry = sk.Geometry1D(
        cos_sza=cos_sza,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=alt_grid,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.PlaneParallel,
    )

    viewing_geo = sk.ViewingGeometry()
    for alt in flux_alts:
        viewing_geo.add_flux_observer(sk.FluxObserverSolar(cos_sza, alt))

    wavel = np.array([500.0, 550.0])
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    hg = sk.optical.HenyeyGreenstein.from_parameters(
        wavel, np.ones(2), np.ones(2) * ssa, np.ones(2) * 0.7
    )
    profile = np.exp(-0.5 * ((alt_grid - 4000.0) / 1000.0) ** 2)
    profile /= np.trapezoid(profile, alt_grid)
    atmosphere["plume"] = sk.constituent.ExtinctionScatterer(
        hg, alt_grid, profile, 550.0
    )
    atmosphere["surface"] = sk.constituent.LambertianSurface(0.1)

    engine = sk.Engine(config, model_geometry, viewing_geo)
    rad = engine.calculate_radiance(atmosphere)

    # Copy, the storage arrays are views into the atmosphere
    extinction = np.array(atmosphere.storage.total_extinction)

    return rad, alt_grid, extinction


def test_direct_flux_inside_layer():
    """
    The direct solar beam included in the downwelling and actinic fluxes must be attenuated
    to the flux observer altitude, both for observers on and between grid levels.
    """
    cos_sza = 0.8
    flux_alts = np.array([0.0, 999.99, 2500.0, 3000.0, 3700.0, 4000.0, 5250.0, 20000.0])

    rad_do, alt_grid, extinction = _plume_flux_scenario(
        sk.SingleScatterSource.DiscreteOrdinates, flux_alts, cos_sza=cos_sza
    )
    rad_diffuse, _, _ = _plume_flux_scenario(
        sk.SingleScatterSource.NoSource, flux_alts, cos_sza=cos_sza
    )

    # Optical depth from TOA at the grid levels, integrating the linearly interpolated
    # extinction exactly.  The discrete ordinates layers are homogeneous, so inside a
    # layer the optical depth is linear in altitude.
    layer_od = 0.5 * (extinction[1:] + extinction[:-1]) * np.diff(alt_grid)[:, None]
    od_from_top = np.concatenate(
        [np.cumsum(layer_od[::-1], axis=0)[::-1], np.zeros((1, extinction.shape[1]))]
    )
    od_obs = np.stack(
        [np.interp(flux_alts, alt_grid, od) for od in od_from_top.T], axis=1
    )

    # (wavelength, flux_location)
    direct_actinic = np.exp(-od_obs / cos_sza).T

    np.testing.assert_allclose(
        rad_do["downwelling_flux"].to_numpy()
        - rad_diffuse["downwelling_flux"].to_numpy(),
        direct_actinic * cos_sza,
        rtol=1e-10,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        rad_do["actinic_flux"].to_numpy() - rad_diffuse["actinic_flux"].to_numpy(),
        direct_actinic,
        rtol=1e-10,
        atol=1e-12,
    )
    # The single scatter source only changes the direct beam, not the diffuse upwelling
    np.testing.assert_allclose(
        rad_do["upwelling_flux"].to_numpy(),
        rad_diffuse["upwelling_flux"].to_numpy(),
        rtol=1e-12,
    )


def test_flux_energy_conservation():
    """
    In a conservatively scattering atmosphere the net flux must be constant with altitude,
    i.e., there should be no spurious heating/cooling around a scattering layer.
    """
    flux_alts = np.arange(0.0, 15001.0, 500.0)

    rad, _, _ = _plume_flux_scenario(
        sk.SingleScatterSource.DiscreteOrdinates, flux_alts, dz=1000.0, ssa=1.0
    )

    net = rad["downwelling_flux"].to_numpy() - rad["upwelling_flux"].to_numpy()

    np.testing.assert_allclose(
        net, np.broadcast_to(net[:, -1:], net.shape), rtol=1e-8, atol=0.0
    )
