from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_mie_simple_database():
    _ = sk.database.MieDatabase(
        sk.mie.distribution.LogNormalDistribution(),
        sk.mie.refractive.H2SO4(),
        np.array([532, 1020]),
        median_radius=np.array([100, 200]),
        mode_width=np.array([1.5, 1.6, 1.7]),
    )


def test_mie_db_in_engine():
    config = sk.Config()

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

    wavel = np.arange(280.0, 800.0, 0.1)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["ozone"] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
    atmosphere["no2"] = sk.climatology.mipas.constituent(
        "NO2", sk.optical.NO2Vandaele()
    )

    mie_db = sk.database.MieDatabase(
        sk.mie.distribution.LogNormalDistribution().freeze(
            median_radius=80, mode_width=1.6
        ),
        sk.mie.refractive.H2SO4(),
        wavelengths_nm=np.arange(270, 1000, 50.0),
    )

    aero_ext = np.zeros(len(model_geometry.altitudes()))
    aero_ext[0:30] = 1e-7
    atmosphere["aerosol"] = sk.constituent.ExtinctionScatterer(
        mie_db,
        altitudes_m=model_geometry.altitudes(),
        extinction_per_m=aero_ext,
        extinction_wavelength_nm=745,
    )

    engine = sk.Engine(config, model_geometry, viewing_geo)
    _ = engine.calculate_radiance(atmosphere)


def test_mie_database_against_online():
    config = sk.Config()
    config.num_stokes = 3
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates

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

    wavel = np.array([750])
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["ozone"] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
    atmosphere["no2"] = sk.climatology.mipas.constituent(
        "NO2", sk.optical.NO2Vandaele()
    )

    mie_db = sk.database.MieDatabase(
        sk.mie.distribution.LogNormalDistribution().freeze(
            median_radius=80, mode_width=1.6
        ),
        sk.mie.refractive.H2SO4(),
        wavelengths_nm=np.array([700, 750, 800]),
    )

    aero_ext = np.zeros(len(model_geometry.altitudes()))
    aero_ext[0:30] = 1e-7
    atmosphere["aerosol"] = sk.constituent.ExtinctionScatterer(
        mie_db,
        altitudes_m=model_geometry.altitudes(),
        extinction_per_m=aero_ext,
        extinction_wavelength_nm=750,
    )

    engine = sk.Engine(config, model_geometry, viewing_geo)
    rad_with_db = engine.calculate_radiance(atmosphere)

    atmosphere["aerosol"] = sk.constituent.ExtinctionScatterer(
        sk.optical.Mie(
            sk.mie.distribution.LogNormalDistribution().freeze(
                median_radius=80, mode_width=1.6
            ),
            sk.mie.refractive.H2SO4(),
        ),
        altitudes_m=model_geometry.altitudes(),
        extinction_per_m=aero_ext,
        extinction_wavelength_nm=750,
    )

    rad_without_db = engine.calculate_radiance(atmosphere)

    p_diff = (
        np.abs(
            rad_with_db["radiance"].to_numpy() - rad_without_db["radiance"].to_numpy()
        )
        / rad_with_db["radiance"].to_numpy()
        * 100
    )

    np.testing.assert_array_almost_equal(p_diff, 0.0, decimal=3)
