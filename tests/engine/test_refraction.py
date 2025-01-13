from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_los_refraction_refractive_one():
    # Tests that the model gives the same results with LOS refraction enabled if the refractive
    # index is forced to 1
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

    config = sk.Config()
    config.los_refraction = True

    wavel = np.arange(280.0, 800.0, 10)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geometry.altitudes(),
        np.ones_like(model_geometry.altitudes()) * 1e-6,
    )

    engine_refraction = sk.Engine(config, model_geometry, viewing_geo)

    radiance_refracted = engine_refraction.calculate_radiance(atmosphere)

    config = sk.Config()
    config.los_refraction = False

    wavel = np.arange(280.0, 800.0, 10)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geometry.altitudes(),
        np.ones_like(model_geometry.altitudes()) * 1e-6,
    )

    engine = sk.Engine(config, model_geometry, viewing_geo)

    radiance = engine.calculate_radiance(atmosphere)

    np.testing.assert_allclose(
        radiance_refracted["radiance"].to_numpy(),
        radiance["radiance"].to_numpy(),
        rtol=1e-4,
    )


def test_multiple_scatter_refraction_refractive_one():
    # Tests that the model gives the same results with multiple scatter refraction enabled if the refractive
    # index is forced to 1
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

    config = sk.Config()
    config.multiple_scatter_refraction = True
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders

    wavel = np.array([500.0])
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geometry.altitudes(),
        np.ones_like(model_geometry.altitudes()) * 1e-6,
    )

    engine_refraction = sk.Engine(config, model_geometry, viewing_geo)

    radiance_refracted = engine_refraction.calculate_radiance(atmosphere)

    config = sk.Config()
    config.multiple_scatter_refraction = False
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders

    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geometry.altitudes(),
        np.ones_like(model_geometry.altitudes()) * 1e-6,
    )

    engine = sk.Engine(config, model_geometry, viewing_geo)

    radiance = engine.calculate_radiance(atmosphere)

    np.testing.assert_allclose(
        radiance_refracted["radiance"].to_numpy(),
        radiance["radiance"].to_numpy(),
        rtol=1e-4,
    )


def test_solar_refraction_refractive_one_discrete_ordinates():
    # Tests that the model gives the same results with solar refraction enabled if the refractive
    # index is forced to 1 and we are using the discrete ordinates source
    csz = 0.1
    model_geometry = sk.Geometry1D(
        cos_sza=csz,
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
            cos_sza=csz,
        )
        viewing_geo.add_ray(ray)

    config = sk.Config()
    config.solar_refraction = True
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 2

    wavel = np.arange(280.0, 800.0, 10)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    # model_geometry.refractive_index = sk.optical.refraction.ciddor_index_of_refraction(
    #    atmosphere.temperature_k, atmosphere.pressure_pa, 0.0, 450, 600
    # )

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geometry.altitudes(),
        np.ones_like(model_geometry.altitudes()) * 1e-6,
    )

    engine_refraction = sk.Engine(config, model_geometry, viewing_geo)

    radiance_refracted = engine_refraction.calculate_radiance(atmosphere)

    config = sk.Config()
    config.solar_refraction = False
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 2

    wavel = np.arange(280.0, 800.0, 10)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geometry.altitudes(),
        np.ones_like(model_geometry.altitudes()) * 1e-6,
    )

    engine = sk.Engine(config, model_geometry, viewing_geo)

    radiance = engine.calculate_radiance(atmosphere)

    np.testing.assert_allclose(
        radiance_refracted["radiance"].to_numpy(),
        radiance["radiance"].to_numpy(),
        rtol=1e-4,
    )
