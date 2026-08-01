from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.arange(0.0, 30_001.0, 5_000.0)
WAVELENGTHS_NM = np.array([500.0, 750.0])


def successive_orders_config(num_stokes: int) -> sk.Config:
    config = sk.Config()
    config.num_stokes = num_stokes
    config.num_streams = 8
    config.num_singlescatter_moments = 16
    config.num_sza = 5
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersRust
    config.num_successive_orders_incoming = 26
    config.num_successive_orders_outgoing = 26
    config.successive_orders_max_iterations = 3
    config.successive_orders_relative_tolerance = 0.0
    config.successive_orders_absolute_tolerance = 0.0
    config.successive_orders_anderson_depth = 0
    return config


def geometry1d() -> sk.Geometry1D:
    return sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.2,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )


def geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.2,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=np.linspace(-0.2, 0.2, 9),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
    )


def viewing_geometry() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    # Cancelling the geometry's solar azimuth keeps the LOS in the modeled 2D
    # plane while retaining an oblique solar plane and therefore nonzero U.
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.6,
            relative_azimuth=-0.2,
            cos_viewing_zenith=0.7,
            observer_altitude_m=100_000.0,
        )
    )
    return viewing


def homogeneous_atmosphere(
    geometry: sk.Geometry1D | sk.Geometry2D,
    config: sk.Config,
    *,
    calculate_derivatives: bool = False,
) -> sk.Atmosphere:
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=calculate_derivatives,
    )
    altitude_factor = np.exp(-ALTITUDES_M / 7_500.0)[:, np.newaxis]
    spectral_factor = np.array([[0.8, 1.2]])
    extinction = (2.4e-5 * altitude_factor + 1.0e-9) * spectral_factor
    ssa = (
        0.91
        + 0.025 * np.exp(-ALTITUDES_M / 18_000.0)[:, np.newaxis]
        - np.array([[0.0, 0.01]])
    )
    if isinstance(geometry, sk.Geometry2D):
        extinction = np.tile(extinction, (geometry.shape[0], 1))
        ssa = np.tile(ssa, (geometry.shape[0], 1))

    atmosphere.storage.total_extinction[:] = extinction
    atmosphere.storage.ssa[:] = ssa
    atmosphere.leg_coeff.a1[0] = 1.0
    atmosphere.leg_coeff.a1[1] = 0.08
    atmosphere.leg_coeff.a1[2] = 0.5
    if atmosphere.nstokes == 3:
        atmosphere.leg_coeff.a2[2] = 3.0
        atmosphere.leg_coeff.b1[2] = -np.sqrt(6.0) / 2.0
    # Isolate atmospheric vector transport from the differently discretized
    # 1D SZA and 2D horizontal surface-source grids.
    atmosphere.surface.albedo[:] = 0.0
    return atmosphere


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_homogeneous_2d_matches_1d(num_stokes: int):
    config = successive_orders_config(num_stokes)
    geometry_1d = geometry1d()
    geometry_2d = geometry2d()
    viewing = viewing_geometry()

    result_1d = sk.Engine(config, geometry_1d, viewing).calculate_radiance(
        homogeneous_atmosphere(geometry_1d, config)
    )
    result_2d = sk.Engine(config, geometry_2d, viewing).calculate_radiance(
        homogeneous_atmosphere(geometry_2d, config)
    )

    if num_stokes == 3:
        assert np.any(np.abs(result_1d.radiance.sel(stokes="Q")) > 1.0e-10)
        assert np.any(np.abs(result_1d.radiance.sel(stokes="U")) > 1.0e-10)
    np.testing.assert_allclose(
        result_2d.radiance, result_1d.radiance, rtol=1.1e-2, atol=2.0e-13
    )


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_twostream_initialization_preserves_converged_2d_solution(num_stokes: int):
    baseline_config = successive_orders_config(num_stokes)
    warm_config = successive_orders_config(num_stokes)
    for config in (baseline_config, warm_config):
        config.successive_orders_max_iterations = 50
        config.successive_orders_relative_tolerance = 1.0e-9
        config.successive_orders_absolute_tolerance = 1.0e-12
        config.successive_orders_anderson_depth = 3
    warm_config.successive_orders_initialization = (
        sk.SuccessiveOrdersInitialization.TwoStream
    )

    geometry = geometry2d()
    viewing = viewing_geometry()
    baseline = sk.Engine(baseline_config, geometry, viewing).calculate_radiance(
        homogeneous_atmosphere(geometry, baseline_config)
    )
    warm = sk.Engine(warm_config, geometry, viewing).calculate_radiance(
        homogeneous_atmosphere(geometry, warm_config)
    )

    np.testing.assert_allclose(
        warm.radiance, baseline.radiance, rtol=5.0e-7, atol=2.0e-12
    )


@pytest.mark.parametrize(
    "source_altitudes",
    [
        np.array([2_500.0, 9_000.0, 18_000.0, 27_500.0]),
        np.array([15_000.0]),
    ],
)
def test_explicit_successive_orders_altitude_grid_m_supports_1d_and_2d(
    source_altitudes: np.ndarray,
):
    config = successive_orders_config(1)
    config.successive_orders_altitude_grid_m = source_altitudes
    geometry_1d = geometry1d()
    geometry_2d = geometry2d()
    viewing = viewing_geometry()

    result_1d = sk.Engine(config, geometry_1d, viewing).calculate_radiance(
        homogeneous_atmosphere(geometry_1d, config)
    )
    result_2d = sk.Engine(config, geometry_2d, viewing).calculate_radiance(
        homogeneous_atmosphere(geometry_2d, config)
    )

    assert np.all(np.isfinite(result_1d.radiance))
    assert np.all(np.isfinite(result_2d.radiance))
    np.testing.assert_allclose(
        result_2d.radiance, result_1d.radiance, rtol=1.1e-2, atol=2.0e-13
    )


def test_explicit_midpoint_grid_preserves_default_solution():
    geometry = geometry2d()
    viewing = viewing_geometry()
    default_config = successive_orders_config(1)
    explicit_config = successive_orders_config(1)
    explicit_config.successive_orders_altitude_grid_m = (
        ALTITUDES_M[:-1] + ALTITUDES_M[1:]
    ) / 2

    default_result = sk.Engine(default_config, geometry, viewing).calculate_radiance(
        homogeneous_atmosphere(geometry, default_config)
    )
    explicit_result = sk.Engine(explicit_config, geometry, viewing).calculate_radiance(
        homogeneous_atmosphere(geometry, explicit_config)
    )

    np.testing.assert_allclose(
        explicit_result.radiance, default_result.radiance, rtol=0.0, atol=0.0
    )


def test_explicit_successive_orders_altitudes_must_be_inside_atmosphere():
    config = successive_orders_config(1)
    config.successive_orders_altitude_grid_m = np.array([-1.0, 10_000.0])

    with pytest.raises(
        ValueError, match="successive_orders_altitude_grid_m.*atmospheric"
    ):
        sk.Engine(config, geometry2d(), viewing_geometry())


def test_native_products_resolve_horizontal_atmosphere_nodes():
    config = successive_orders_config(1)
    config.successive_orders_max_iterations = 50
    config.successive_orders_relative_tolerance = 1.0e-10
    config.successive_orders_absolute_tolerance = 1.0e-12
    config.successive_orders_anderson_depth = 3
    config.successive_orders_altitude_grid_m = np.array(
        [2_500.0, 9_000.0, 18_000.0, 27_500.0]
    )
    geometry = geometry2d()
    atmosphere = homogeneous_atmosphere(geometry, config, calculate_derivatives=True)
    horizontal_factor = np.linspace(0.9, 1.1, geometry.shape[0])[:, np.newaxis]
    location_factor = np.repeat(horizontal_factor, geometry.shape[1], axis=1).reshape(
        -1, 1
    )
    atmosphere.storage.total_extinction[:] *= location_factor

    engine = sk.Engine(config, geometry, viewing_geometry())
    linearization = engine.linearize(atmosphere)
    assert linearization.backends == {
        "jvp": sk.LinearizationBackend.Native,
        "vjp": sk.LinearizationBackend.Native,
    }

    tangent = linearization.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = np.linspace(
        -2.0e-7, 3.0e-7, tangent["extinction"].size
    ).reshape(tangent["extinction"].shape)
    tangent["ssa"].data[:] = np.linspace(-0.02, 0.03, tangent["ssa"].size).reshape(
        tangent["ssa"].shape
    )
    jvp = linearization.jvp(tangent)

    cotangent = xr.ones_like(linearization.value)
    cotangent.data[:] = np.linspace(0.4, 1.0, cotangent.size).reshape(cotangent.shape)
    gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))
    lhs = float((jvp * cotangent).sum())
    rhs = float(
        (tangent["extinction"] * gradient["extinction"]).sum()
        + (tangent["ssa"] * gradient["ssa"]).sum()
    )
    np.testing.assert_allclose(lhs, rhs, rtol=3.0e-9, atol=2.0e-11)

    jacobian = linearization.jacobian
    assert jacobian["extinction"].dims == (
        "horizontal_angle",
        "altitude",
        "wavelength",
        "los",
        "stokes",
    )
    jacobian_jvp = (
        jacobian["extinction"] * tangent["extinction"]
        + jacobian["ssa"] * tangent["ssa"]
    ).sum(("horizontal_angle", "altitude"))
    xr.testing.assert_allclose(
        jacobian_jvp.transpose(*linearization.value.dims),
        jvp,
        rtol=3.0e-9,
        atol=2.0e-11,
    )
    for name in ("extinction", "ssa"):
        xr.testing.assert_allclose(
            (jacobian[name] * cotangent)
            .sum(linearization.value.dims)
            .transpose("horizontal_angle", "altitude"),
            gradient[name],
            rtol=1.0e-6,
            atol=2.0e-9,
        )

    epsilon = 2.0e-5
    perturbed_radiances = []
    for sign in (-1.0, 1.0):
        perturbed = homogeneous_atmosphere(geometry, config)
        perturbed.storage.total_extinction[:] *= location_factor
        perturbed.storage.total_extinction[:] += (
            sign * epsilon * tangent["extinction"].values.reshape(-1, 1)
        )
        perturbed.storage.ssa[:] += (
            sign * epsilon * tangent["ssa"].values.reshape(-1, 1)
        )
        perturbed_radiances.append(engine.calculate_radiance(perturbed).radiance.values)
    finite_difference = (perturbed_radiances[1] - perturbed_radiances[0]) / (
        2.0 * epsilon
    )
    np.testing.assert_allclose(jvp.values, finite_difference, rtol=5.0e-5, atol=2.0e-9)

    assert np.any(np.abs(gradient["extinction"].isel(horizontal_angle=0)) > 0.0)
    assert not np.allclose(
        gradient["extinction"].isel(horizontal_angle=0),
        gradient["extinction"].isel(horizontal_angle=-1),
        rtol=1.0e-4,
        atol=0.0,
    )


def test_rayon_native_products_match_serial_in_2d():
    serial_config = successive_orders_config(1)
    rayon_config = successive_orders_config(1)
    for config in (serial_config, rayon_config):
        config.successive_orders_max_iterations = 50
        config.successive_orders_relative_tolerance = 1.0e-10
        config.successive_orders_absolute_tolerance = 1.0e-12
        config.successive_orders_anderson_depth = 3
        config.successive_orders_initialization = (
            sk.SuccessiveOrdersInitialization.TwoStream
        )
    rayon_config.num_threads = 2
    rayon_config.threading_lib = sk.ThreadingLib.Rayon
    rayon_config.threading_model = sk.ThreadingModel.Wavelength

    serial_geometry = geometry2d()
    rayon_geometry = geometry2d()
    serial_atmosphere = homogeneous_atmosphere(
        serial_geometry, serial_config, calculate_derivatives=True
    )
    rayon_atmosphere = homogeneous_atmosphere(
        rayon_geometry, rayon_config, calculate_derivatives=True
    )
    serial = sk.Engine(serial_config, serial_geometry, viewing_geometry()).linearize(
        serial_atmosphere
    )
    rayon = sk.Engine(rayon_config, rayon_geometry, viewing_geometry()).linearize(
        rayon_atmosphere
    )
    xr.testing.assert_identical(rayon.value, serial.value)

    tangent = serial.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = np.linspace(
        -2.0e-7, 3.0e-7, tangent["extinction"].size
    ).reshape(tangent["extinction"].shape)
    tangent["ssa"].data[:] = np.linspace(-0.02, 0.03, tangent["ssa"].size).reshape(
        tangent["ssa"].shape
    )
    xr.testing.assert_allclose(rayon.jvp(tangent), serial.jvp(tangent))

    cotangent = xr.ones_like(serial.value)
    cotangent.data[:] = np.linspace(0.4, 1.0, cotangent.size).reshape(cotangent.shape)
    serial_gradient = serial.vjp(cotangent, parameters=("extinction", "ssa"))
    rayon_gradient = rayon.vjp(cotangent, parameters=("extinction", "ssa"))
    for parameter in ("extinction", "ssa"):
        xr.testing.assert_allclose(
            rayon_gradient[parameter],
            serial_gradient[parameter],
            rtol=2.0e-12,
            atol=2.0e-12,
        )


def test_geometry2d_rejects_diffuse_refraction():
    config = successive_orders_config(1)
    config.multiple_scatter_refraction = True

    with pytest.raises(NotImplementedError, match="refracted diffuse rays"):
        sk.Engine(config, geometry2d(), viewing_geometry())
