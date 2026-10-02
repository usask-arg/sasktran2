from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.arange(0.0, 60_001.0, 1000.0)
TOP_RADIUS_M = EARTH_RADIUS_M + ALTITUDES_M[-1]
# Logarithm of an approximately visible wavelength dry air refractive index
LOG_REFRACTIVE_INDEX = np.log(1.0 + 2.8e-4 * np.exp(-ALTITUDES_M / 8000.0))
EXTINCTION_PER_M = 6.0e-6 * np.exp(-ALTITUDES_M / 7000.0)

SINGLE_SCATTER_SOURCES = [sk.SingleScatterSource.Exact, sk.SingleScatterSource.Table]


def geometry1d(cos_sza: float, refracting: bool = True) -> sk.Geometry1D:
    geometry = sk.Geometry1D(
        cos_sza=cos_sza,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    if refracting:
        geometry.refractive_index = np.exp(LOG_REFRACTIVE_INDEX)
    return geometry


def config(
    single_scatter_source: sk.SingleScatterSource,
    solar_refraction: bool,
    num_stokes: int = 1,
    multiple_scatter_source: sk.MultipleScatterSource = (
        sk.MultipleScatterSource.NoSource
    ),
) -> sk.Config:
    result = sk.Config()
    result.num_stokes = num_stokes
    result.single_scatter_source = single_scatter_source
    result.multiple_scatter_source = multiple_scatter_source
    result.num_singlescatter_moments = 6
    result.solar_refraction = solar_refraction
    return result


def limb_viewing(
    cos_sza: float,
    tangent_altitudes_m=(5000.0, 12000.0, 25000.0),
    relative_azimuths=(0.0, 0.4, 1.2),
):
    viewing = sk.ViewingGeometry()
    for tangent_altitude, relative_azimuth in zip(
        tangent_altitudes_m, relative_azimuths, strict=True
    ):
        viewing.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m=tangent_altitude,
                relative_azimuth=relative_azimuth,
                observer_altitude_m=200_000.0,
                cos_sza=cos_sza,
            )
        )
    return viewing


def nadir_viewing(cos_sza: float) -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=cos_sza,
            relative_azimuth=0.0,
            cos_viewing_zenith=1.0,
            observer_altitude_m=200_000.0,
        )
    )
    return viewing


def rayleigh_atmosphere(
    geometry: sk.Geometry1D, model_config: sk.Config, calculate_derivatives=False
) -> sk.Atmosphere:
    atmosphere = sk.Atmosphere(
        geometry,
        model_config,
        wavelengths_nm=np.array([350.0, 600.0]),
        calculate_derivatives=calculate_derivatives,
    )
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    return atmosphere


def radiance(geometry, viewing, model_config) -> np.ndarray:
    atmosphere = rayleigh_atmosphere(geometry, model_config)
    return (
        sk.Engine(model_config, geometry, viewing)
        .calculate_radiance(atmosphere)["radiance"]
        .to_numpy()
    )


def single_level_radiance(
    sza_deg: float,
    level_altitude_m: float,
    solar_refraction: bool,
    single_scatter_source: sk.SingleScatterSource,
) -> float:
    """Nadir radiance with an isotropic scatterer at a single altitude level.

    The radiance is then proportional to the solar transmission at one point.
    """
    cos_sza = np.cos(np.radians(sza_deg))
    geometry = geometry1d(cos_sza)
    model_config = config(single_scatter_source, solar_refraction)
    model_config.num_singlescatter_moments = 2
    atmosphere = sk.Atmosphere(
        geometry,
        model_config,
        wavelengths_nm=np.array([600.0]),
        calculate_derivatives=False,
    )
    ssa = np.zeros_like(ALTITUDES_M)
    ssa[np.argmin(np.abs(ALTITUDES_M - level_altitude_m))] = 1.0
    atmosphere.storage.total_extinction[:] = EXTINCTION_PER_M[:, None]
    atmosphere.storage.ssa[:] = ssa[:, None]
    atmosphere.leg_coeff.a1[0] = 1.0
    result = sk.Engine(model_config, geometry, nadir_viewing(cos_sza))
    return float(
        result.calculate_radiance(atmosphere)["radiance"].to_numpy().ravel()[0]
    )


def reference_ray(radius: float, zenith: float, refract: bool = True):
    """Integrates the ray equation with an ODE solver, independent of the model.

    Returns the angle between the start position and the direction of the ray
    after it leaves the atmosphere and the optical depth along the ray, or None
    if the ray intersects the surface.
    """
    spacing = ALTITUDES_M[1] - ALTITUDES_M[0]

    def log_index_gradient(altitude):
        if not refract or not (0.0 <= altitude < ALTITUDES_M[-1]):
            return 0.0
        i = min(int(altitude // spacing), len(ALTITUDES_M) - 2)
        return (LOG_REFRACTIVE_INDEX[i + 1] - LOG_REFRACTIVE_INDEX[i]) / spacing

    def derivative(_, y):
        r, _, zenith_angle, _ = y
        altitude = r - EARTH_RADIUS_M
        return [
            np.cos(zenith_angle),
            np.sin(zenith_angle) / r,
            -np.sin(zenith_angle) * (1.0 / r + log_index_gradient(altitude)),
            np.interp(altitude, ALTITUDES_M, EXTINCTION_PER_M),
        ]

    def leaves_atmosphere(_, y):
        return y[0] - TOP_RADIUS_M

    def hits_surface(_, y):
        return y[0] - EARTH_RADIUS_M

    leaves_atmosphere.terminal = True
    leaves_atmosphere.direction = 1
    hits_surface.terminal = True
    hits_surface.direction = -1
    solution = solve_ivp(
        derivative,
        [0.0, 5e6],
        [radius, 0.0, zenith, 0.0],
        events=[leaves_atmosphere, hits_surface],
        rtol=1e-11,
        atol=1e-9,
        max_step=500.0,
    )
    if solution.t_events[1].size:
        return None
    _, swept, exit_zenith, optical_depth = solution.y_events[0][0]
    # The refractive index is one at the top of the atmosphere to ~1e-11
    return swept + exit_zenith, optical_depth


def reference_transmission(sza_deg: float, level_altitude_m: float, refract=True):
    radius = EARTH_RADIUS_M + level_altitude_m
    geometric_zenith = np.radians(sza_deg)
    if not refract:
        return np.exp(-reference_ray(radius, geometric_zenith, refract=False)[1])

    def residual(zenith):
        return reference_ray(radius, zenith)[0] - geometric_zenith

    apparent_zenith = brentq(
        residual, geometric_zenith - 0.03, geometric_zenith, xtol=1e-12
    )
    return np.exp(-reference_ray(radius, apparent_zenith)[1])


@pytest.mark.parametrize("single_scatter_source", SINGLE_SCATTER_SOURCES)
@pytest.mark.parametrize("num_stokes", [1, 3])
def test_1d_unity_solar_refraction_matches_straight(single_scatter_source, num_stokes):
    cos_sza = 0.05
    geometry = geometry1d(cos_sza, refracting=False)
    viewing = limb_viewing(cos_sza)

    straight = radiance(
        geometry, viewing, config(single_scatter_source, False, num_stokes)
    )
    refracted = radiance(
        geometry, viewing, config(single_scatter_source, True, num_stokes)
    )

    # Refracted rays rebuild their geometry from integrated deflection angles,
    # which agrees with the straight line construction to roundoff
    np.testing.assert_allclose(refracted, straight, rtol=1e-6, atol=1e-10)


@pytest.mark.parametrize("single_scatter_source", SINGLE_SCATTER_SOURCES)
def test_1d_solar_refraction_brightens_low_sun_limb(single_scatter_source):
    cos_sza = 0.05
    geometry = geometry1d(cos_sza)
    viewing = limb_viewing(cos_sza)

    straight = radiance(geometry, viewing, config(single_scatter_source, False))
    refracted = radiance(geometry, viewing, config(single_scatter_source, True))

    assert np.all(np.isfinite(refracted))
    # The sun appears higher, shortening the solar paths at low altitudes
    assert np.all(refracted[:, 0, 0] > straight[:, 0, 0] * 1.005)


@pytest.mark.parametrize("single_scatter_source", SINGLE_SCATTER_SOURCES)
@pytest.mark.parametrize(
    ("sza_deg", "level_altitude_m"), [(88.0, 10_000.0), (91.5, 10_000.0)]
)
def test_1d_refracted_solar_transmission_matches_ray_equation(
    single_scatter_source, sza_deg, level_altitude_m
):
    # Straight solar rays calibrate the factor relating the radiance to the
    # solar transmission at the scattering point
    straight = single_level_radiance(
        sza_deg, level_altitude_m, False, single_scatter_source
    )
    refracted = single_level_radiance(
        sza_deg, level_altitude_m, True, single_scatter_source
    )
    straight_transmission = reference_transmission(
        sza_deg, level_altitude_m, refract=False
    )
    model_transmission = refracted / straight * straight_transmission

    # The remaining difference is the straight chord approximation of the ray
    # tracer within refracted layers
    np.testing.assert_allclose(
        model_transmission,
        reference_transmission(sza_deg, level_altitude_m),
        rtol=2e-3,
    )


@pytest.mark.parametrize("single_scatter_source", SINGLE_SCATTER_SOURCES)
def test_1d_solar_refraction_illuminates_beyond_geometric_shadow(
    single_scatter_source,
):
    # At 4 km the straight shadow begins at 92.03 deg and the refracted one at
    # 92.87 deg
    def calculate(sza_deg, solar_refraction):
        return single_level_radiance(
            sza_deg, 4000.0, solar_refraction, single_scatter_source
        )

    assert calculate(92.5, False) == 0.0
    assert calculate(92.5, True) > 0.0
    assert calculate(93.0, True) == 0.0


@pytest.mark.parametrize("solar_refraction", [False, True])
@pytest.mark.parametrize("cos_sza", [0.05, -0.05])
def test_1d_table_single_scatter_matches_exact_in_twilight(cos_sza, solar_refraction):
    geometry = geometry1d(cos_sza)
    viewing = limb_viewing(cos_sza)

    exact = radiance(
        geometry, viewing, config(sk.SingleScatterSource.Exact, solar_refraction)
    )
    table = radiance(
        geometry, viewing, config(sk.SingleScatterSource.Table, solar_refraction)
    )

    assert np.all(exact > 0.0)
    np.testing.assert_allclose(table, exact, rtol=5e-3)


@pytest.mark.parametrize("single_scatter_source", SINGLE_SCATTER_SOURCES)
def test_1d_empty_ray_does_not_shift_following_rays(single_scatter_source):
    cos_sza = 0.3
    geometry = geometry1d(cos_sza, refracting=False)
    model_config = config(single_scatter_source, False)

    alone = radiance(geometry, limb_viewing(cos_sza, [20_000.0], [0.4]), model_config)
    # The 70 km tangent altitude is above the top of the atmosphere
    with_empty = radiance(
        geometry,
        limb_viewing(cos_sza, [70_000.0, 20_000.0], [0.4, 0.4]),
        model_config,
    )

    assert np.all(with_empty[:, 0] == 0.0)
    np.testing.assert_allclose(with_empty[:, 1], alone[:, 0], rtol=1e-12)


@pytest.mark.parametrize("solar_refraction", [False, True])
def test_1d_table_single_scatter_derivatives_match_exact(solar_refraction):
    cos_sza = 0.3
    geometry = geometry1d(cos_sza)
    viewing = limb_viewing(cos_sza)

    def calculate(single_scatter_source):
        model_config = config(single_scatter_source, solar_refraction)
        atmosphere = rayleigh_atmosphere(
            geometry, model_config, calculate_derivatives=True
        )
        return sk.Engine(model_config, geometry, viewing).calculate_radiance(atmosphere)

    exact = calculate(sk.SingleScatterSource.Exact)
    table = calculate(sk.SingleScatterSource.Table)

    for name in ["wf_pressure_pa", "wf_temperature_k"]:
        exact_wf = exact[name].to_numpy()
        table_wf = table[name].to_numpy()
        assert np.all(np.isfinite(table_wf))
        np.testing.assert_allclose(
            table_wf, exact_wf, rtol=0.0, atol=2e-3 * np.abs(exact_wf).max()
        )


@pytest.mark.parametrize(
    "multiple_scatter_source",
    [
        sk.MultipleScatterSource.SuccessiveOrders,
        sk.MultipleScatterSource.DiscreteOrdinates,
    ],
)
def test_1d_multiple_scatter_unity_solar_refraction_matches_straight(
    multiple_scatter_source,
):
    cos_sza = 0.1
    geometry = geometry1d(cos_sza, refracting=False)
    viewing = limb_viewing(cos_sza)

    def model_config(solar_refraction):
        result = config(
            sk.SingleScatterSource.Exact,
            solar_refraction,
            multiple_scatter_source=multiple_scatter_source,
        )
        result.num_streams = 4
        return result

    straight = radiance(geometry, viewing, model_config(False))
    refracted = radiance(geometry, viewing, model_config(True))

    np.testing.assert_allclose(refracted, straight, rtol=1e-8)


@pytest.mark.parametrize(
    "multiple_scatter_source",
    [
        sk.MultipleScatterSource.SuccessiveOrders,
        sk.MultipleScatterSource.DiscreteOrdinates,
    ],
)
def test_1d_multiple_scatter_solar_refraction_brightens_low_sun(
    multiple_scatter_source,
):
    cos_sza = 0.05
    geometry = geometry1d(cos_sza)
    viewing = limb_viewing(cos_sza)

    def model_config(solar_refraction):
        result = config(
            sk.SingleScatterSource.NoSource,
            solar_refraction,
            multiple_scatter_source=multiple_scatter_source,
        )
        result.num_streams = 4
        return result

    straight = radiance(geometry, viewing, model_config(False))
    refracted = radiance(geometry, viewing, model_config(True))

    assert np.all(np.isfinite(refracted))
    assert np.all(refracted[:, 0, 0] > straight[:, 0, 0])
