from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def test_los_refraction_refractive_one():
    # Tests that the model gives the same results with LOS refraction enabled if the refractive
    # index is forced to 1
    model_geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.arange(0, 65001, 1000.0),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    for alt in [10000, 20000, 30000, 40000]:
        ray = sk.TangentAltitudeSolar(
            tangent_altitude_m=alt,
            relative_azimuth=0.1,
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
        altitude_grid_m=np.arange(0, 65001, 1000.0),
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
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersLegacy

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
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersLegacy

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
        altitude_grid_m=np.arange(0, 65001, 1000.0),
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


def test_refraction_enabling():
    # Tests that the model gives different radiances when refraction is enabled
    csz = 0.1
    model_geometry = sk.Geometry1D(
        cos_sza=csz,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.arange(0, 65001, 1000.0),
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
    config.los_refraction = True
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 2

    wavel = np.arange(280.0, 800.0, 10)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    model_geometry.refractive_index = sk.optical.refraction.ciddor_index_of_refraction(
        atmosphere.temperature_k, atmosphere.pressure_pa, 0.0, 450, 600
    )

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

    # Verify that this should fail
    try:
        np.testing.assert_allclose(
            radiance_refracted["radiance"].to_numpy(),
            radiance["radiance"].to_numpy(),
            rtol=1e-4,
        )
        pytest.fail(
            "Refraction enabled and disabled should give different results, but they are the same."
        )
    except AssertionError:
        # This is expected
        pass


def _thin_shell_single_scatter(los_refraction, num_stokes, single_scatter_source):
    altitude_grid_m = np.arange(0.0, 100001.0, 1000.0)
    cos_sza = 0.3

    model_geometry = sk.Geometry1D(
        cos_sza=cos_sza,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=altitude_grid_m,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    # Grazing rays that hit the ground, the scattering angle changes by a few tenths of a degree
    # along the refracted path, almost all of it below the scattering shell
    viewing_geo = sk.ViewingGeometry()
    for relative_azimuth in np.deg2rad([0.0, 60.0, 180.0]):
        viewing_geo.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m=-5000,
                relative_azimuth=relative_azimuth,
                observer_altitude_m=200000,
                cos_sza=cos_sza,
            )
        )

    config = sk.Config()
    config.num_stokes = num_stokes
    config.single_scatter_source = single_scatter_source
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.num_singlescatter_moments = 64
    config.los_refraction = los_refraction

    state_atmosphere = sk.Atmosphere(
        model_geometry, sk.Config(), wavelengths_nm=np.array([500.0])
    )
    sk.climatology.us76.add_us76_standard_atmosphere(state_atmosphere)
    model_geometry.refractive_index = sk.optical.refraction.ciddor_index_of_refraction(
        state_atmosphere.temperature_k, state_atmosphere.pressure_pa, 0.0, 450, 500
    )

    # Scattering only takes place in a thin shell at 40 km above a black surface
    atmosphere = sk.Atmosphere(model_geometry, config, numwavel=1)
    atmosphere.storage.total_extinction[:] = np.where(
        altitude_grid_m == 40000.0, 1.0e-6, 0.0
    )[:, np.newaxis]
    atmosphere.storage.ssa[:] = 1.0
    atmosphere.leg_coeff.a1[:] = 0.0
    if num_stokes == 1:
        # Strongly forward peaked Henyey-Greenstein phase function
        g = 0.75
        order = np.arange(atmosphere.leg_coeff.a1.shape[0])
        atmosphere.leg_coeff.a1[:] = ((2 * order + 1) * g**order)[:, None, None]
    else:
        # Rayleigh phase matrix
        atmosphere.leg_coeff.a1[0] = 1.0
        atmosphere.leg_coeff.a1[2] = 0.5
        atmosphere.leg_coeff.a2[2] = 3.0
        atmosphere.leg_coeff.b1[2] = -np.sqrt(6.0) / 2.0
    atmosphere.surface.albedo[:] = 0.0

    engine = sk.Engine(config, model_geometry, viewing_geo)

    return engine.calculate_radiance(atmosphere)["radiance"].to_numpy()


@pytest.mark.parametrize(
    "single_scatter_source",
    [sk.SingleScatterSource.Exact, sk.SingleScatterSource.Table],
)
@pytest.mark.parametrize("num_stokes", [1, 3])
def test_los_refraction_single_scatter_uses_layer_direction(
    num_stokes, single_scatter_source
):
    # Tests that each layer of a refracted line of sight uses its own direction for the single
    # scatter phase function and Stokes rotation.  Since there is nothing to scatter below the
    # shell, the refracted radiance must match the straight ray radiance up to the small amount
    # of bending above the shell.  Using the direction of the far end of the ray for every layer
    # causes errors of ~2% for the scalar case and ~0.1-0.3% for I, Q, U in the vector case.
    radiance_refracted = _thin_shell_single_scatter(
        True, num_stokes, single_scatter_source
    )
    radiance = _thin_shell_single_scatter(False, num_stokes, single_scatter_source)

    np.testing.assert_allclose(
        radiance_refracted[..., 0], radiance[..., 0], rtol=5e-4, atol=0
    )
    if num_stokes == 3:
        np.testing.assert_allclose(
            radiance_refracted[..., 1:] / radiance[..., :1],
            radiance[..., 1:] / radiance[..., :1],
            rtol=0,
            atol=2e-4,
        )
