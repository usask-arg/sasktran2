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


def _log_linear_bending(altitudes, refractive_index, earth_radius, impact, legs):
    """Independent bending-angle reference for a log-linear refractive profile.

    Each leg is a (lower, upper) radius pair traversed monotonically. The bending
    increment is -b dln(n)/dr / sqrt(n^2 r^2 - b^2) dr, integrated with the
    substitution r = r0 + x^2 so the tangent-point singularity is removed.
    """
    from scipy.integrate import quad

    log_n = np.log(refractive_index)
    radii = earth_radius + altitudes
    slopes = np.diff(log_n) / np.diff(radii)

    def n_times_r(r):
        return np.exp(np.interp(r - earth_radius, altitudes, log_n)) * r

    total = 0.0
    for r_low, r_high, r0 in legs:

        def integrand(x, r0=r0):
            r = r0 + x * x
            i = np.clip(np.searchsorted(radii, r) - 1, 0, len(slopes) - 1)
            nr = n_times_r(r)
            return -impact * slopes[i] * 2 * x / np.sqrt((nr - impact) * (nr + impact))

        x_low, x_high = np.sqrt(r_low - r0), np.sqrt(r_high - r0)
        nodes = np.sqrt(radii[(radii > r_low) & (radii < r_high)] - r0)
        total += quad(integrand, x_low, x_high, points=nodes, limit=500)[0]
    return total


def _deflection_scene():
    earth_radius = 6372000.0
    altitudes = np.arange(0, 100001, 1000.0)
    model_geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0,
        earth_radius_m=earth_radius,
        altitude_grid_m=altitudes,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    # Exponential refractivity that reaches exactly 1 at the top of the grid
    refractive_index = 1.0 + 2.8e-4 * (
        np.exp(-altitudes / 7500.0) - np.exp(-altitudes[-1] / 7500.0)
    )
    model_geometry.refractive_index = refractive_index
    return model_geometry, altitudes, refractive_index, earth_radius


def _calculate_deflection(model_geometry, viewing_geo, los_refraction=True):
    config = sk.Config()
    config.los_refraction = los_refraction
    config.single_scatter_source = sk.SingleScatterSource.NoSource

    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=np.array([600.0]))
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    return sk.Engine(config, model_geometry, viewing_geo).calculate_radiance(atmosphere)


def test_los_refraction_deflection_angle_matches_bending_integral():
    from scipy.optimize import brentq

    model_geometry, altitudes, refractive_index, earth_radius = _deflection_scene()
    r_top = earth_radius + altitudes[-1]
    log_n = np.log(refractive_index)

    def n_at(r):
        return np.exp(np.interp(r - earth_radius, altitudes, log_n))

    viewing_geo = sk.ViewingGeometry()
    expected = []

    # Observer outside the atmosphere, limb viewing: bending accumulates on both
    # sides of the refracted tangent point
    for tangent_altitude in [2000.0, 5000.0, 10000.0, 20000.0, 35000.0]:
        viewing_geo.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m=tangent_altitude,
                relative_azimuth=0.0,
                observer_altitude_m=200000,
                cos_sza=0.6,
            )
        )
        impact = earth_radius + tangent_altitude
        rt = brentq(lambda r, b=impact: n_at(r) * r - b, earth_radius, impact)
        expected.append(
            _log_linear_bending(
                altitudes,
                refractive_index,
                earth_radius,
                impact,
                [(rt, r_top, rt), (rt, r_top, rt)],
            )
        )

    # Refracted tangent points just below a grid node, where the tangent layer
    # is only a few centimetres thick
    for node, offset in [(5000.0, 0.05), (20000.0, 0.02)]:
        rt = earth_radius + node - offset
        impact = n_at(rt) * rt
        viewing_geo.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m=impact - earth_radius,
                relative_azimuth=0.0,
                observer_altitude_m=200000,
                cos_sza=0.6,
            )
        )
        expected.append(
            _log_linear_bending(
                altitudes,
                refractive_index,
                earth_radius,
                impact,
                [(rt, r_top, rt), (rt, r_top, rt)],
            )
        )

    # Observer outside the atmosphere looking at the ground
    observer_radius = earth_radius + 200000.0
    cos_viewing = -0.3
    viewing_geo.add_ray(sk.SolarAnglesObserverLocation(0.6, 0.0, cos_viewing, 200000.0))
    impact = observer_radius * np.sqrt(1 - cos_viewing**2)
    expected.append(
        _log_linear_bending(
            altitudes,
            refractive_index,
            earth_radius,
            impact,
            [(earth_radius, r_top, impact / n_at(earth_radius) - 1.0)],
        )
    )

    # Observer inside the atmosphere looking up at 10 degrees elevation, where
    # the ray invariant includes the refractive index at the observer
    observer_altitude = 1000.0
    observer_radius = earth_radius + observer_altitude
    cos_viewing = np.sin(np.deg2rad(10.0))
    viewing_geo.add_ray(
        sk.SolarAnglesObserverLocation(0.6, 0.0, cos_viewing, observer_altitude)
    )
    impact = n_at(observer_radius) * observer_radius * np.sqrt(1 - cos_viewing**2)
    expected.append(
        _log_linear_bending(
            altitudes,
            refractive_index,
            earth_radius,
            impact,
            [(observer_radius, r_top, impact / n_at(observer_radius) - 1.0)],
        )
    )

    radiance = _calculate_deflection(model_geometry, viewing_geo)

    deflection = radiance["los_refraction_deflection_angle"]
    assert deflection.dims == ("los",)
    assert deflection.attrs["units"] == "rad"
    assert np.all(deflection.to_numpy() > 0)
    np.testing.assert_allclose(deflection.to_numpy(), expected, rtol=1e-5, atol=2e-8)


def test_los_refraction_deflection_angle_unity_index_and_disabled():
    model_geometry, altitudes, _, _ = _deflection_scene()

    viewing_geo = sk.ViewingGeometry()
    for tangent_altitude in [5000.0, 20000.0]:
        viewing_geo.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m=tangent_altitude,
                relative_azimuth=0.0,
                observer_altitude_m=200000,
                cos_sza=0.6,
            )
        )

    radiance = _calculate_deflection(model_geometry, viewing_geo, los_refraction=False)
    assert "los_refraction_deflection_angle" not in radiance

    model_geometry.refractive_index = np.ones_like(altitudes)
    radiance = _calculate_deflection(model_geometry, viewing_geo)
    np.testing.assert_allclose(
        radiance["los_refraction_deflection_angle"].to_numpy(), 0.0, atol=1e-9
    )
