from __future__ import annotations

import gc

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.arange(0.0, 30_001.0, 5_000.0)
WAVELENGTHS_NM = np.array([510.0, 690.0])


def _config(
    source: sk.MultipleScatterSource,
    *,
    num_stokes: int = 1,
    iterations: int = 4,
    relative_tolerance: float = 0.0,
    absolute_tolerance: float = 0.0,
    anderson_depth: int = 0,
) -> sk.Config:
    config = sk.Config()
    config.num_threads = 1
    config.num_stokes = num_stokes
    config.num_streams = 8
    config.num_singlescatter_moments = 8
    config.num_sza = 1
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = source
    config.num_successive_orders_incoming = 26
    config.num_successive_orders_outgoing = 26
    config.num_successive_orders_iterations = iterations
    config.successive_orders_relative_tolerance = relative_tolerance
    config.successive_orders_absolute_tolerance = absolute_tolerance
    config.successive_orders_anderson_depth = anderson_depth
    config.successive_orders_damping = 1.0
    config.delta_m_scaling = False
    return config


def _geometry() -> sk.Geometry1D:
    return sk.Geometry1D(
        cos_sza=0.55,
        solar_azimuth=0.15,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )


def _viewing_geometry() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.55,
            relative_azimuth=0.35,
            cos_viewing_zenith=0.72,
            observer_altitude_m=100_000.0,
        )
    )
    viewing.add_ray(
        sk.TangentAltitudeSolar(
            tangent_altitude_m=11_000.0,
            relative_azimuth=-0.4,
            observer_altitude_m=100_000.0,
            cos_sza=0.55,
        )
    )
    return viewing


def _atmosphere(
    geometry: sk.Geometry1D,
    config: sk.Config,
    *,
    calculate_derivatives: bool = False,
    extinction_shift: np.ndarray | None = None,
    ssa_shift: np.ndarray | None = None,
    surface_albedo: float = 0.12,
) -> sk.Atmosphere:
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=calculate_derivatives,
    )

    altitude_factor = np.exp(-ALTITUDES_M / 7_600.0)[:, np.newaxis]
    spectral_factor = np.array([0.9, 1.1])[np.newaxis, :]
    atmosphere.storage.total_extinction[:] = (
        1.4e-5 * altitude_factor + 2.0e-9
    ) * spectral_factor
    atmosphere.storage.ssa[:] = (
        0.86 + 0.04 * np.exp(-ALTITUDES_M / 14_000.0)[:, np.newaxis]
    ) - np.array([0.0, 0.015])[np.newaxis, :]

    if extinction_shift is not None:
        atmosphere.storage.total_extinction[:] += extinction_shift[:, np.newaxis]
    if ssa_shift is not None:
        atmosphere.storage.ssa[:] += ssa_shift[:, np.newaxis]

    atmosphere.leg_coeff.a1[0, :, :] = 1.0
    atmosphere.leg_coeff.a1[1, :, :] = 0.08
    atmosphere.leg_coeff.a1[2, :, :] = 0.5
    if config.num_stokes == 3:
        atmosphere.leg_coeff.a2[2, :, :] = 3.0
        atmosphere.leg_coeff.b1[2, :, :] = -np.sqrt(6.0) / 2.0
    atmosphere.surface.albedo[:] = surface_albedo
    return atmosphere


def _add_forward_delta_peak(
    atmosphere: sk.Atmosphere, num_stokes: int, fraction: float
) -> None:
    for degree in range(atmosphere.leg_coeff.a1.shape[0]):
        delta_coefficient = (2 * degree + 1) * fraction
        atmosphere.leg_coeff.a1[degree] = (1 - fraction) * atmosphere.leg_coeff.a1[
            degree
        ] + delta_coefficient
        if num_stokes == 3:
            atmosphere.leg_coeff.a2[degree] = (1 - fraction) * atmosphere.leg_coeff.a2[
                degree
            ] + delta_coefficient
            atmosphere.leg_coeff.a3[degree] = (1 - fraction) * atmosphere.leg_coeff.a3[
                degree
            ] + delta_coefficient
            atmosphere.leg_coeff.b1[degree] *= 1 - fraction


def _calculate(
    config: sk.Config,
    geometry: sk.Geometry1D | None = None,
    **atmosphere_kwargs,
) -> xr.Dataset:
    if geometry is None:
        geometry = _geometry()
    engine = sk.Engine(config, geometry, _viewing_geometry())
    return engine.calculate_radiance(_atmosphere(geometry, config, **atmosphere_kwargs))


def test_scalar_1d_primal_is_finite_and_nonnegative():
    config = _config(sk.MultipleScatterSource.SuccessiveOrders)
    result = _calculate(config)

    assert result.radiance.shape == (WAVELENGTHS_NM.size, 2, 1)
    assert np.all(np.isfinite(result.radiance))
    assert np.min(result.radiance.values) >= -2.0e-13
    assert np.max(result.radiance.values) > 1.0e-6


@pytest.mark.parametrize(
    ("num_stokes", "relative_tolerance", "absolute_tolerance"),
    [(1, 2.0e-3, 2.0e-7), (3, 2.0e-3, 2.0e-7)],
)
def test_fixed_iteration_solution_agrees_with_legacy_successive_orders(
    num_stokes: int,
    relative_tolerance: float,
    absolute_tolerance: float,
):
    cpp_config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        num_stokes=num_stokes,
        iterations=3,
    )
    legacy_config = _config(
        sk.MultipleScatterSource.SuccessiveOrdersLegacy,
        num_stokes=num_stokes,
        iterations=3,
    )

    cpp = _calculate(cpp_config, surface_albedo=0.0).radiance
    legacy = _calculate(legacy_config, surface_albedo=0.0).radiance

    xr.testing.assert_allclose(
        cpp,
        legacy,
        rtol=relative_tolerance,
        atol=absolute_tolerance,
    )
    if num_stokes == 3:
        assert np.any(np.abs(cpp.sel(stokes="Q")) > 1.0e-8)
        assert np.any(np.abs(cpp.sel(stokes="U")) > 1.0e-8)


def test_explicit_midpoint_altitudes_preserve_default_and_coarse_grid_is_used():
    default_config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=3,
    )
    midpoint_config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=3,
    )
    midpoint_config.successive_orders_altitude_grid_m = (
        ALTITUDES_M[:-1] + ALTITUDES_M[1:]
    ) / 2.0
    coarse_config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=3,
    )
    coarse_config.successive_orders_altitude_grid_m = np.array(
        [2_000.0, 9_000.0, 19_000.0, 28_000.0]
    )

    default = _calculate(default_config).radiance
    midpoint = _calculate(midpoint_config).radiance
    coarse = _calculate(coarse_config).radiance

    xr.testing.assert_allclose(midpoint, default, rtol=2.0e-12, atol=2.0e-14)
    assert np.all(np.isfinite(coarse))
    assert not np.allclose(coarse, default, rtol=1.0e-7, atol=1.0e-12)


def test_explicit_altitudes_must_be_inside_the_model_atmosphere():
    config = _config(sk.MultipleScatterSource.SuccessiveOrders)
    config.successive_orders_altitude_grid_m = np.array([-1.0, 10_000.0])

    with pytest.raises(RuntimeError, match="Failed to create Engine"):
        sk.Engine(config, _geometry(), _viewing_geometry())


def test_explicit_horizontal_source_angles_require_geometry2d():
    config = _config(sk.MultipleScatterSource.SuccessiveOrders)
    config.successive_orders_horizontal_angle_grid_radians = np.array([-0.1, 0.1])

    with pytest.raises(ValueError, match="supported only with Geometry2D"):
        sk.Engine(config, _geometry(), _viewing_geometry())


def test_zero_tolerances_use_the_configured_fixed_iteration_count():
    one_iteration = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=1,
    )
    two_iterations = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=2,
    )
    early_exit = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=8,
        absolute_tolerance=1.0e100,
    )

    one = _calculate(one_iteration).radiance
    two = _calculate(two_iterations).radiance
    early = _calculate(early_exit).radiance

    xr.testing.assert_allclose(early, one, rtol=2.0e-13, atol=2.0e-14)
    difference = np.linalg.norm((two - one).values)
    assert difference > 1.0e-7 * np.linalg.norm(one.values)


def test_fixed_iteration_result_does_not_depend_on_engine_history():
    config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=3,
    )
    geometry = _geometry()
    viewing = _viewing_geometry()
    reused_engine = sk.Engine(config, geometry, viewing)

    reused_engine.calculate_radiance(_atmosphere(geometry, config))
    shifted_atmosphere = _atmosphere(
        geometry,
        config,
        extinction_shift=np.linspace(0.0, 3.0e-7, ALTITUDES_M.size),
        ssa_shift=np.linspace(-0.02, 0.01, ALTITUDES_M.size),
    )
    reused = reused_engine.calculate_radiance(shifted_atmosphere).radiance
    fresh = (
        sk.Engine(config, geometry, viewing)
        .calculate_radiance(shifted_atmosphere)
        .radiance
    )

    xr.testing.assert_allclose(reused, fresh, rtol=2.0e-12, atol=2.0e-14)


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_native_jvp_vjp_match_finite_difference_jacobian_and_adjoint(num_stokes):
    config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        num_stokes=num_stokes,
        iterations=60,
        relative_tolerance=1.0e-11,
        absolute_tolerance=1.0e-13,
        anderson_depth=3,
    )
    geometry = _geometry()
    atmosphere = _atmosphere(geometry, config, calculate_derivatives=True)
    engine = sk.Engine(config, geometry, _viewing_geometry())
    linearization = engine.linearize(atmosphere)

    assert linearization.backends == {
        "jvp": sk.LinearizationBackend.Native,
        "vjp": sk.LinearizationBackend.Native,
    }

    tangent = linearization.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = 0.08 * np.mean(
        atmosphere.storage.total_extinction, axis=1
    )
    tangent["ssa"].data[:] = np.linspace(-0.015, 0.02, ALTITUDES_M.size)
    jvp = linearization.jvp(tangent)

    cotangent = xr.ones_like(linearization.value)
    cotangent.data[:] = np.linspace(0.35, 1.1, cotangent.size).reshape(cotangent.shape)
    gradient = linearization.vjp(
        cotangent,
        parameters=("extinction", "ssa"),
    )

    lhs = float((jvp * cotangent).sum())
    rhs = float(
        (tangent["extinction"] * gradient["extinction"]).sum()
        + (tangent["ssa"] * gradient["ssa"]).sum()
    )
    np.testing.assert_allclose(lhs, rhs, rtol=2.0e-8, atol=2.0e-11)

    jacobian = linearization.jacobian
    jacobian_jvp = (
        jacobian["extinction"] * tangent["extinction"]
        + jacobian["ssa"] * tangent["ssa"]
    ).sum("altitude")
    xr.testing.assert_allclose(
        jacobian_jvp.transpose(*jvp.dims),
        jvp,
        rtol=2.0e-7,
        atol=2.0e-11,
    )
    for parameter in ("extinction", "ssa"):
        expected_gradient = (jacobian[parameter] * cotangent).sum(
            linearization.value.dims
        )
        xr.testing.assert_allclose(
            gradient[parameter],
            expected_gradient,
            rtol=2.0e-7,
            atol=2.0e-10,
        )

    epsilon = 2.0e-3
    perturbed = []
    for sign in (-1.0, 1.0):
        shifted_atmosphere = _atmosphere(
            geometry,
            config,
            extinction_shift=(sign * epsilon * tangent["extinction"].values),
            ssa_shift=sign * epsilon * tangent["ssa"].values,
        )
        perturbed.append(engine.calculate_radiance(shifted_atmosphere).radiance.values)
    finite_difference = (perturbed[1] - perturbed[0]) / (2.0 * epsilon)
    np.testing.assert_allclose(
        jvp.values,
        finite_difference,
        rtol=2.0e-4,
        atol=2.0e-9,
    )


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_delta_m_native_products_match_finite_difference_and_are_adjoint(num_stokes):
    config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        num_stokes=num_stokes,
        iterations=80,
        relative_tolerance=1.0e-11,
        absolute_tolerance=1.0e-13,
        anderson_depth=3,
    )
    config.num_singlescatter_moments = 9
    config.delta_m_scaling = True
    geometry = _geometry()
    atmosphere = _atmosphere(geometry, config, calculate_derivatives=True)
    _add_forward_delta_peak(atmosphere, num_stokes, 0.1)
    engine = sk.Engine(config, geometry, _viewing_geometry())
    linearization = engine.linearize(atmosphere)

    coefficient_index = 8 if num_stokes == 1 else 4 * 8
    parameter = f"leg_coeff_{coefficient_index}"
    tangent = linearization.tangent_template[[parameter]]
    tangent[parameter].data[:] = np.linspace(-0.04, 0.06, ALTITUDES_M.size)
    jvp = linearization.jvp(tangent)

    cotangent = xr.ones_like(linearization.value)
    cotangent.data[:] = np.linspace(0.3, 1.2, cotangent.size).reshape(cotangent.shape)
    gradient = linearization.vjp(cotangent, parameters=(parameter,))
    lhs = float((jvp * cotangent).sum())
    rhs = float((tangent[parameter] * gradient[parameter]).sum())
    np.testing.assert_allclose(lhs, rhs, rtol=3.0e-8, atol=3.0e-11)

    epsilon = 2.0e-4
    perturbed = []
    for sign in (-1.0, 1.0):
        shifted_atmosphere = _atmosphere(geometry, config)
        _add_forward_delta_peak(shifted_atmosphere, num_stokes, 0.1)
        shifted_atmosphere.leg_coeff.a1[8] += (
            sign * epsilon * tangent[parameter].values[:, np.newaxis]
        )
        perturbed.append(engine.calculate_radiance(shifted_atmosphere).radiance.values)
    finite_difference = (perturbed[1] - perturbed[0]) / (2.0 * epsilon)
    np.testing.assert_allclose(
        jvp.values,
        finite_difference,
        rtol=3.0e-4,
        atol=3.0e-9,
    )


@pytest.mark.parametrize(
    "threading_model",
    [sk.ThreadingModel.Wavelength, sk.ThreadingModel.Source],
)
def test_scalar_threading_models_match_single_thread(threading_model):
    reference_config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=4,
    )
    threaded_config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=4,
    )
    threaded_config.num_threads = 2
    threaded_config.threading_model = threading_model

    reference = _calculate(reference_config).radiance
    threaded = _calculate(threaded_config).radiance
    xr.testing.assert_allclose(threaded, reference, rtol=2.0e-12, atol=2.0e-14)


@pytest.mark.parametrize(
    "threading_model",
    [sk.ThreadingModel.Wavelength, sk.ThreadingModel.Source],
)
@pytest.mark.parametrize("num_stokes", [1, 3])
def test_threaded_native_products_are_adjoint(threading_model, num_stokes):
    config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        num_stokes=num_stokes,
        iterations=60,
        relative_tolerance=1.0e-11,
        absolute_tolerance=1.0e-13,
        anderson_depth=3,
    )
    config.num_threads = 2
    config.threading_model = threading_model
    geometry = _geometry()
    atmosphere = _atmosphere(geometry, config, calculate_derivatives=True)
    linearization = sk.Engine(config, geometry, _viewing_geometry()).linearize(
        atmosphere
    )

    tangent = linearization.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = 0.03 * np.mean(
        atmosphere.storage.total_extinction, axis=1
    )
    tangent["ssa"].data[:] = np.linspace(-0.01, 0.012, ALTITUDES_M.size)
    jvp = linearization.jvp(tangent)
    cotangent = xr.ones_like(linearization.value)
    cotangent.data[:] = np.linspace(0.4, 1.2, cotangent.size).reshape(cotangent.shape)
    gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))

    lhs = float((jvp * cotangent).sum())
    rhs = float(
        (tangent["extinction"] * gradient["extinction"]).sum()
        + (tangent["ssa"] * gradient["ssa"]).sum()
    )
    assert np.all(np.isfinite(jvp))
    assert all(np.all(np.isfinite(value)) for value in gradient.data_vars.values())
    np.testing.assert_allclose(lhs, rhs, rtol=3.0e-8, atol=3.0e-11)


def test_native_vjp_returns_a_finite_result_without_convergence():
    config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=1,
        relative_tolerance=1.0e-15,
        absolute_tolerance=1.0e-15,
        anderson_depth=0,
    )
    geometry = _geometry()
    atmosphere = _atmosphere(geometry, config, calculate_derivatives=True)
    linearization = sk.Engine(config, geometry, _viewing_geometry()).linearize(
        atmosphere
    )
    gradient = linearization.vjp(
        xr.ones_like(linearization.value), parameters=("extinction", "ssa")
    )

    assert all(np.all(np.isfinite(value)) for value in gradient.data_vars.values())


def test_homogeneous_scalar_fast_paths_preserve_native_adjoint():
    config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=60,
        relative_tolerance=1.0e-11,
        absolute_tolerance=1.0e-13,
        anderson_depth=3,
    )
    geometry = _geometry()
    atmosphere = _atmosphere(geometry, config, calculate_derivatives=True)
    atmosphere.storage.ssa[:] = 0.91
    atmosphere.mark_changed()
    linearization = sk.Engine(config, geometry, _viewing_geometry()).linearize(
        atmosphere
    )

    tangent = linearization.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = 0.04 * np.mean(
        atmosphere.storage.total_extinction, axis=1
    )
    tangent["ssa"].data[:] = 0.012
    jvp = linearization.jvp(tangent)
    cotangent = xr.ones_like(linearization.value)
    gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))

    lhs = float((jvp * cotangent).sum())
    rhs = float(
        (tangent["extinction"] * gradient["extinction"]).sum()
        + (tangent["ssa"] * gradient["ssa"]).sum()
    )
    np.testing.assert_allclose(lhs, rhs, rtol=3.0e-8, atol=3.0e-11)


def test_scalar_repeated_vjp_survives_atmosphere_and_engine_lifetimes():
    """Exercise the compact LOS metadata across retrieval-like lifetimes.

    This intentionally performs enough primal/VJP pairs to cover the failure
    window of the original repeated-evaluation crash. It swaps atmosphere
    objects without rebuilding the engine, then destroys and recreates the
    engine, so source caches and derivative scratch must remain scoped to the
    active calculation objects.
    """
    config = _config(
        sk.MultipleScatterSource.SuccessiveOrders,
        iterations=40,
        relative_tolerance=1.0e-9,
        absolute_tolerance=1.0e-12,
        anderson_depth=3,
    )
    config.num_threads = 1
    config.threading_model = sk.ThreadingModel.Wavelength
    geometry = _geometry()
    viewing = _viewing_geometry()

    product_count = 0
    for engine_index in range(2):
        engine = sk.Engine(config, geometry, viewing)
        for atmosphere_index in range(2):
            atmosphere = _atmosphere(
                geometry,
                config,
                calculate_derivatives=True,
                surface_albedo=0.08 + 0.01 * atmosphere_index,
            )
            base_extinction = atmosphere.storage.total_extinction.copy()
            base_ssa = atmosphere.storage.ssa.copy()

            for evaluation in range(50):
                scale = 1.0 + 2.0e-3 * np.sin(0.37 * evaluation + 0.11 * engine_index)
                atmosphere.storage.total_extinction[:] = base_extinction * scale
                atmosphere.storage.ssa[:] = np.clip(
                    base_ssa
                    + 5.0e-4 * np.cos(0.23 * evaluation + 0.17 * atmosphere_index),
                    0.0,
                    1.0,
                )
                atmosphere.mark_changed()

                linearization = engine.linearize(atmosphere)
                cotangent = xr.ones_like(linearization.value)
                cotangent.data[:] = np.linspace(
                    0.4 + 1.0e-4 * evaluation,
                    1.1 + 1.0e-4 * evaluation,
                    cotangent.size,
                ).reshape(cotangent.shape)
                gradient = linearization.vjp(
                    cotangent, parameters=("extinction", "ssa")
                )

                assert np.all(np.isfinite(linearization.value))
                assert np.all(np.isfinite(gradient["extinction"]))
                assert np.all(np.isfinite(gradient["ssa"]))
                product_count += 1

            del atmosphere, linearization, gradient

        del engine
        gc.collect()

    assert product_count == 200
