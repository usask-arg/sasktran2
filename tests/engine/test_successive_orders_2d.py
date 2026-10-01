from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.array([0.0, 10_000.0, 30_000.0])
HORIZONTAL_ANGLES_RAD = np.array([-0.6, 0.0, 0.6])


def successive_orders_config(
    *,
    num_stokes: int = 1,
    single_scatter_source: sk.SingleScatterSource = sk.SingleScatterSource.Exact,
    multiple_scatter_source: sk.MultipleScatterSource = (
        sk.MultipleScatterSource.SuccessiveOrders
    ),
) -> sk.Config:
    config = sk.Config()
    config.num_threads = 1
    config.num_stokes = num_stokes
    config.num_streams = 4
    config.num_singlescatter_moments = 4
    config.single_scatter_source = single_scatter_source
    config.multiple_scatter_source = multiple_scatter_source
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    config.num_successive_orders_incoming = 6
    config.num_successive_orders_outgoing = 6
    config.num_successive_orders_iterations = 2
    config.successive_orders_relative_tolerance = 0.0
    config.successive_orders_absolute_tolerance = 0.0
    config.delta_m_scaling = False
    return config


def geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=HORIZONTAL_ANGLES_RAD,
    )


def viewing_geometry() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.6,
            relative_azimuth=0.2,
            cos_viewing_zenith=0.7,
            observer_altitude_m=100_000.0,
        )
    )
    return viewing


def atmosphere(
    geometry: sk.Geometry2D,
    config: sk.Config,
    *,
    horizontal_slope: float = 0.0,
    calculate_derivatives: bool = False,
) -> sk.Atmosphere:
    result = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        calculate_derivatives=calculate_derivatives,
    )
    horizontal, altitude = np.meshgrid(
        HORIZONTAL_ANGLES_RAD, ALTITUDES_M, indexing="ij"
    )
    result.storage.total_extinction[:, 0] = (
        1.5e-5 * np.exp(-altitude / 15_000.0) * (1.0 + horizontal_slope * horizontal)
    ).ravel()
    result.storage.ssa[:, 0] = (0.88 - 0.08 * altitude / ALTITUDES_M[-1]).ravel()
    result.leg_coeff.a1[0] = 1.0
    result.leg_coeff.a1[2] = 0.3
    if config.num_stokes == 3:
        result.leg_coeff.a2[2] = 3.0
        result.leg_coeff.b1[2] = -np.sqrt(6.0) / 2.0
    result.surface.albedo[:] = 0.1
    return result


def calculate(
    geometry: sk.Geometry2D,
    *,
    num_stokes: int = 1,
    horizontal_slope: float = 0.0,
    single_source: sk.SingleScatterSource = sk.SingleScatterSource.Exact,
    source: sk.MultipleScatterSource = sk.MultipleScatterSource.SuccessiveOrders,
) -> np.ndarray:
    config = successive_orders_config(
        num_stokes=num_stokes,
        single_scatter_source=single_source,
        multiple_scatter_source=source,
    )
    result = sk.Engine(config, geometry, viewing_geometry()).calculate_radiance(
        atmosphere(geometry, config, horizontal_slope=horizontal_slope)
    )
    return result.radiance.values


def _spectral_cache_atmosphere(
    geometry: sk.Geometry2D,
    config: sk.Config,
    spectral_indices: np.ndarray,
    *,
    calculate_derivatives: bool,
    uniform_phase: bool = False,
    legendre_derivative: bool = False,
) -> sk.Atmosphere:
    """Use distinct optics so reusing another wavelength's cache is visible."""
    wavelengths_nm = np.array([410.0, 530.0, 690.0])[spectral_indices]
    result = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=wavelengths_nm,
        calculate_derivatives=calculate_derivatives,
        legendre_derivative=legendre_derivative,
    )
    horizontal, altitude = np.meshgrid(
        HORIZONTAL_ANGLES_RAD, ALTITUDES_M, indexing="ij"
    )
    result.storage.total_extinction[:] = (
        1.5e-5 * np.exp(-altitude / 15_000.0) * (1.0 + 0.3 * horizontal)
    ).reshape(-1, 1) * np.array([0.8, 1.1, 1.5])[spectral_indices]
    result.storage.ssa[:] = (
        0.89 - 0.08 * altitude / ALTITUDES_M[-1] + 0.02 * horizontal
    ).reshape(-1, 1) - np.array([0.0, 0.03, 0.07])[spectral_indices]
    result.leg_coeff.a1[0] = 1.0
    phase_profile = 0.3 if uniform_phase else (0.3 + 0.04 * horizontal).reshape(-1, 1)
    result.leg_coeff.a1[2] = phase_profile + 0.05 * spectral_indices
    result.surface.albedo[:] = np.array([0.06, 0.12, 0.23])[spectral_indices]
    result.mark_changed()
    return result


@pytest.mark.parametrize("num_threads", [1, 2])
def test_2d_mixed_uniform_phase_history_preserves_native_products(num_threads):
    config = successive_orders_config()
    config.num_threads = num_threads
    config.threading_model = sk.ThreadingModel.Wavelength
    config.num_sza = 3
    geometry = geometry2d()
    viewing = viewing_geometry()
    indices = np.arange(2)
    scene = _spectral_cache_atmosphere(
        geometry,
        config,
        indices,
        calculate_derivatives=True,
        uniform_phase=True,
        legendre_derivative=True,
    )
    engine = sk.Engine(config, geometry, viewing)
    initial_extinction = scene.storage.total_extinction.copy()
    initial_ssa = scene.storage.ssa.copy()
    initial_albedo = scene.surface.albedo.copy()
    uniform_phase = scene.leg_coeff.a1[2].copy()
    variation = np.linspace(-0.04, 0.04, initial_extinction.shape[0]).reshape(-1, 1)
    scene.leg_coeff.a1[2] = uniform_phase + variation * [0.0, 1.0]
    initial_phase = scene.leg_coeff.a1[2].copy()
    scene.mark_changed()
    parameters = ("extinction", "ssa", "leg_coeff_3", "albedo")
    initial_products = None

    def assert_gradient_equal(actual, expected):
        try:
            if num_threads == 1:
                xr.testing.assert_identical(actual, expected)
            else:
                # Parallel wavelength accumulation already permits final-bit
                # rounding differences when the worker completion order changes.
                xr.testing.assert_allclose(
                    actual, expected, rtol=8 * np.finfo(float).eps, atol=0
                )
        except AssertionError as error:
            diagnostics = [
                f"Native gradient mismatch: update={update!r}, "
                f"num_threads={num_threads}",
                str(error),
            ]
            for name in actual.data_vars:
                if name not in expected:
                    diagnostics.append(f"{name}: missing from expected gradient")
                    continue
                actual_values = actual[name].values
                expected_values = expected[name].values
                if actual_values.shape == expected_values.shape:
                    matches = (
                        np.array_equal(actual_values, expected_values, equal_nan=True)
                        if num_threads == 1
                        else np.allclose(
                            actual_values,
                            expected_values,
                            rtol=8 * np.finfo(float).eps,
                            atol=0,
                            equal_nan=True,
                        )
                    )
                    if matches:
                        continue
                diagnostics.append(
                    f"{name}: actual shape={actual_values.shape}, "
                    f"expected shape={expected_values.shape}"
                )
                for label, values in (
                    ("actual", actual_values),
                    ("expected", expected_values),
                ):
                    diagnostics.append(
                        f"{label}=\n"
                        + np.array2string(
                            values,
                            precision=17,
                            threshold=values.size,
                            floatmode="maxprec_equal",
                        )
                    )
                if (
                    actual_values.shape != expected_values.shape
                    or actual_values.size == 0
                ):
                    continue
                absolute_error = np.abs(actual_values - expected_values)
                ranked_error = np.nan_to_num(
                    absolute_error, nan=np.inf, posinf=np.inf, neginf=np.inf
                )
                worst_index = np.unravel_index(
                    ranked_error.argmax(), actual_values.shape
                )
                nonzero = expected_values != 0
                relative_error = np.full(actual_values.shape, np.nan)
                np.divide(
                    absolute_error,
                    np.abs(expected_values),
                    out=relative_error,
                    where=nonzero,
                )
                max_relative = (
                    relative_error[nonzero].max() if nonzero.any() else np.nan
                )
                diagnostics.append(
                    f"max_absolute_error={absolute_error.max():.17g}, "
                    f"worst_index={worst_index}, "
                    f"actual_at_worst={actual_values[worst_index]:.17g}, "
                    f"expected_at_worst={expected_values[worst_index]:.17g}, "
                    f"relative_error_at_worst={relative_error[worst_index]:.17g}, "
                    f"max_relative_error_where_expected_nonzero={max_relative:.17g}"
                )
            raise AssertionError("\n".join(diagnostics)) from error

    for update in ("initial", "swap", "surface", "both_varying", "restore"):
        if update == "swap":
            scene.leg_coeff.a1[2] = uniform_phase + variation * [1.0, 0.0]
            scene.storage.total_extinction[:] *= 1.12
            scene.storage.ssa[:] -= 0.015
        elif update == "surface":
            scene.surface.albedo[:] *= 1.3
        elif update == "both_varying":
            scene.leg_coeff.a1[2] = uniform_phase + variation * [1.0, 0.7]
        elif update == "restore":
            scene.storage.total_extinction[:] = initial_extinction
            scene.storage.ssa[:] = initial_ssa
            scene.surface.albedo[:] = initial_albedo
            scene.leg_coeff.a1[2] = initial_phase
        if update != "initial":
            # The raw Python API uses a full revision even when only albedo
            # changes; the native source fixture covers surface-only revisions.
            scene.mark_changed()

        linearization = engine.linearize(scene)
        assert linearization.backends == {
            "jvp": sk.LinearizationBackend.Native,
            "vjp": sk.LinearizationBackend.Native,
        }
        tangent = linearization.tangent_template[list(parameters)]
        for name, scale in zip(parameters, (2.0e-7, 0.015, 0.03, 0.04), strict=True):
            tangent[name].values[...] = scale * np.linspace(
                -0.7, 1.1, tangent[name].size
            ).reshape(tangent[name].shape)
        # Degree three is absent from the primal phase but present in its JVP.
        cotangent = xr.ones_like(linearization.value)
        cotangent.values[:] = np.linspace(0.35, 1.1, cotangent.size).reshape(
            cotangent.shape
        )
        gradient = linearization.vjp(cotangent, parameters=parameters)
        jvp = linearization.jvp(tangent)
        assert_gradient_equal(
            linearization.vjp(cotangent, parameters=parameters), gradient
        )
        xr.testing.assert_identical(linearization.jvp(tangent), jvp)

        reference_scene = _spectral_cache_atmosphere(
            geometry,
            config,
            indices,
            calculate_derivatives=True,
            legendre_derivative=True,
        )
        reference_scene.storage.total_extinction[:] = scene.storage.total_extinction
        reference_scene.storage.ssa[:] = scene.storage.ssa
        reference_scene.leg_coeff.a1[:] = scene.leg_coeff.a1
        reference_scene.surface.albedo[:] = scene.surface.albedo
        reference_scene.mark_changed()
        reference = sk.Engine(config, geometry, viewing).linearize(reference_scene)
        xr.testing.assert_identical(linearization.value, reference.value)
        xr.testing.assert_identical(jvp, reference.jvp(tangent))
        assert_gradient_equal(gradient, reference.vjp(cotangent, parameters=parameters))
        products = tuple(
            product.copy(deep=True) for product in (linearization.value, jvp, gradient)
        )
        if update == "initial":
            initial_products = products
        elif update == "restore":
            for actual, initial in zip(products[:2], initial_products[:2], strict=True):
                xr.testing.assert_identical(actual, initial)
            assert_gradient_equal(products[2], initial_products[2])


def test_2d_native_phase_products_survive_derivative_group_transitions():
    config = successive_orders_config()
    config.num_sza = 3
    geometry = geometry2d()
    viewing = viewing_geometry()
    engine = sk.Engine(config, geometry, viewing)
    initial_products = None
    # All phases vary spatially. The same engine must allocate phase scratch
    # for phase mappings and release it when those mappings disappear.
    for evaluation, phase_derivatives in enumerate((False, True, False)):
        scene = _spectral_cache_atmosphere(
            geometry,
            config,
            np.arange(2),
            calculate_derivatives=True,
            legendre_derivative=phase_derivatives,
        )
        linearization = engine.linearize(scene)
        parameters = ("extinction", "ssa", "albedo")
        if phase_derivatives:
            parameters += ("leg_coeff_3",)
        tangent = linearization.tangent_template[list(parameters)]
        for name, scale in zip(parameters, (2.0e-7, 0.015, 0.04, 0.03), strict=False):
            tangent[name].values[...] = scale * np.linspace(
                -0.7, 1.1, tangent[name].size
            ).reshape(tangent[name].shape)
        cotangent = xr.ones_like(linearization.value)
        cotangent.values[:] = np.linspace(0.35, 1.1, cotangent.size).reshape(
            cotangent.shape
        )
        jvp = linearization.jvp(tangent)
        gradient = linearization.vjp(cotangent, parameters=parameters)
        reference = sk.Engine(config, geometry, viewing).linearize(scene)
        xr.testing.assert_identical(linearization.value, reference.value)
        xr.testing.assert_identical(jvp, reference.jvp(tangent))
        xr.testing.assert_identical(
            gradient, reference.vjp(cotangent, parameters=parameters)
        )
        products = tuple(
            product.copy(deep=True) for product in (linearization.value, jvp, gradient)
        )
        if evaluation == 0:
            initial_products = products
        elif evaluation == 2:
            for actual, initial in zip(products, initial_products, strict=True):
                xr.testing.assert_identical(actual, initial)


@pytest.mark.parametrize("uniform_phase", [False, True])
@pytest.mark.parametrize(
    ("num_threads", "threading_model"),
    [
        (1, sk.ThreadingModel.Wavelength),
        (2, sk.ThreadingModel.Wavelength),
        (2, sk.ThreadingModel.Source),
    ],
)
def test_2d_spectral_worker_cache_preserves_native_products_after_updates(
    uniform_phase,
    num_threads,
    threading_model,
):
    config = successive_orders_config()
    config.num_threads = num_threads
    config.threading_model = threading_model
    config.num_sza = 3
    config.num_successive_orders_iterations = 80
    config.successive_orders_relative_tolerance = 1.0e-12
    config.successive_orders_absolute_tolerance = 1.0e-14
    config.successive_orders_anderson_depth = 3
    geometry = geometry2d()
    viewing = viewing_geometry()
    engine = sk.Engine(config, geometry, viewing)
    scene = _spectral_cache_atmosphere(
        geometry,
        config,
        np.arange(3),
        calculate_derivatives=True,
        uniform_phase=uniform_phase,
    )
    # Three wavelengths exceed the worker count, forcing cache eviction.
    # Independent engines retain their single wavelength throughout.
    reference_engines = [sk.Engine(config, geometry, viewing) for _ in range(3)]
    reference_scenes = [
        _spectral_cache_atmosphere(
            geometry,
            config,
            np.array([index]),
            calculate_derivatives=True,
            uniform_phase=uniform_phase,
        )
        for index in range(3)
    ]

    for evaluation in range(3):
        if evaluation:
            for current in [scene, *reference_scenes]:
                if evaluation == 1:
                    current.surface.albedo[:] *= 1.3
                else:
                    current.storage.total_extinction[:] *= 1.12
                    current.storage.ssa[:] -= 0.015
                current.mark_changed()

        linearization = engine.linearize(scene)
        tangent = linearization.tangent_template[["extinction", "ssa"]]
        tangent.extinction.values[:] = np.linspace(
            -2.0e-7, 3.0e-7, tangent.extinction.size
        ).reshape(tangent.extinction.shape)
        tangent.ssa.values[:] = np.linspace(-0.015, 0.02, tangent.ssa.size).reshape(
            tangent.ssa.shape
        )
        cotangent = xr.ones_like(linearization.value)
        cotangent.values[:] = np.linspace(0.35, 1.1, cotangent.size).reshape(
            cotangent.shape
        )
        gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))
        jvp = linearization.jvp(tangent)
        xr.testing.assert_allclose(
            linearization.vjp(cotangent, parameters=("extinction", "ssa")),
            gradient,
            rtol=2.0e-12,
            atol=2.0e-13,
        )
        xr.testing.assert_allclose(
            linearization.jvp(tangent), jvp, rtol=2.0e-12, atol=2.0e-14
        )

        reference_values = []
        reference_jvps = []
        reference_gradient = xr.zeros_like(gradient)
        for index, (reference_engine, reference_scene) in enumerate(
            zip(reference_engines, reference_scenes, strict=True)
        ):
            reference = reference_engine.linearize(reference_scene)
            reference_values.append(reference.value)
            reference_gradient += reference.vjp(
                cotangent.isel(wavelength=[index]), parameters=("extinction", "ssa")
            )
            reference_jvps.append(reference.jvp(tangent))
        xr.testing.assert_allclose(
            linearization.value,
            xr.concat(reference_values, dim="wavelength"),
            rtol=2.0e-12,
            atol=2.0e-14,
        )
        xr.testing.assert_allclose(
            jvp,
            xr.concat(reference_jvps, dim="wavelength"),
            rtol=2.0e-12,
            atol=2.0e-14,
        )
        xr.testing.assert_allclose(
            gradient, reference_gradient, rtol=2.0e-12, atol=2.0e-13
        )
        np.testing.assert_allclose(
            float((jvp * cotangent).sum()),
            float((tangent * gradient).to_array().sum()),
            rtol=3.0e-8,
            atol=3.0e-11,
        )


@pytest.mark.parametrize("iterations", [0, 2])
def test_2d_spectral_worker_cache_preserves_fixed_iteration_history(iterations):
    config = successive_orders_config()
    config.num_sza = 3
    config.num_successive_orders_iterations = iterations
    geometry = geometry2d()
    viewing = viewing_geometry()
    engine = sk.Engine(config, geometry, viewing)
    scene = _spectral_cache_atmosphere(
        geometry, config, np.arange(3), calculate_derivatives=True
    )

    for evaluation in range(2):
        if evaluation:
            scene.storage.total_extinction[:] *= 1.2
            scene.storage.ssa[:] -= 0.025
            scene.surface.albedo[:] *= 1.3
            scene.mark_changed()
        linearization = engine.linearize(scene)
        result = linearization.value
        xr.testing.assert_identical(engine.linearize(scene).value, result)
        tangent = linearization.tangent_template[["extinction", "ssa"]]
        tangent.extinction.values[:] = np.linspace(
            -2.0e-7, 3.0e-7, tangent.extinction.size
        ).reshape(tangent.extinction.shape)
        tangent.ssa.values[:] = 0.012
        cotangent = xr.ones_like(result)
        cotangent.values[:] = np.linspace(0.35, 1.1, cotangent.size).reshape(
            cotangent.shape
        )
        gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))
        jvp = linearization.jvp(tangent)
        # Every native pass revisits wavelengths evicted by the previous pass.
        # These comparisons exercise the saved forcing even at zero orders,
        # where radiance alone does not constrain the implicit native product.
        xr.testing.assert_identical(
            linearization.vjp(cotangent, parameters=("extinction", "ssa")), gradient
        )
        xr.testing.assert_identical(linearization.jvp(tangent), jvp)
        reference_values = []
        reference_jvps = []
        reference_gradient = xr.zeros_like(gradient)
        for index in range(3):
            reference_scene = _spectral_cache_atmosphere(
                geometry, config, np.array([index]), calculate_derivatives=True
            )
            reference_scene.storage.total_extinction[:] = (
                scene.storage.total_extinction[:, [index]]
            )
            reference_scene.storage.ssa[:] = scene.storage.ssa[:, [index]]
            reference_scene.surface.albedo[:] = scene.surface.albedo[:, [index]]
            reference_scene.mark_changed()
            reference = sk.Engine(config, geometry, viewing).linearize(reference_scene)
            reference_values.append(reference.value)
            reference_jvps.append(reference.jvp(tangent))
            reference_gradient += reference.vjp(
                cotangent.isel(wavelength=[index]), parameters=("extinction", "ssa")
            )
        xr.testing.assert_allclose(
            result,
            xr.concat(reference_values, dim="wavelength"),
            rtol=2.0e-12,
            atol=2.0e-14,
        )

        xr.testing.assert_allclose(
            jvp,
            xr.concat(reference_jvps, dim="wavelength"),
            rtol=2.0e-12,
            atol=2.0e-14,
        )
        xr.testing.assert_allclose(
            gradient, reference_gradient, rtol=2.0e-12, atol=2.0e-13
        )


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_2d_successive_orders_is_finite_and_adds_multiple_scatter(num_stokes: int):
    geometry = geometry2d()
    multiple = calculate(
        geometry,
        num_stokes=num_stokes,
        single_source=sk.SingleScatterSource.NoSource,
    )

    assert multiple.shape == (1, 1, num_stokes)
    assert np.all(np.isfinite(multiple))
    assert multiple[0, 0, 0] > 0.0
    if num_stokes == 1:
        combined = calculate(geometry, num_stokes=num_stokes)
        assert np.all(np.isfinite(combined))
        assert combined[0, 0, 0] > multiple[0, 0, 0]


def test_2d_successive_orders_uses_horizontal_atmospheric_structure():
    geometry = geometry2d()
    config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.NoSource
    )
    engine = sk.Engine(config, geometry, viewing_geometry())
    uniform_multiple = engine.calculate_radiance(
        atmosphere(geometry, config)
    ).radiance.values
    varying_multiple = engine.calculate_radiance(
        atmosphere(geometry, config, horizontal_slope=0.8)
    ).radiance.values

    assert not np.isclose(
        varying_multiple.item(), uniform_multiple.item(), rtol=1.0e-4, atol=0.0
    )


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_2d_reduced_horizon_supports_arbitrary_incoming_count(num_stokes: int):
    geometry = geometry2d()
    config = successive_orders_config(
        num_stokes=num_stokes,
        single_scatter_source=sk.SingleScatterSource.NoSource,
    )
    config.num_sza = 3
    config.num_successive_orders_incoming = 37
    config.num_successive_orders_outgoing = 14
    config.successive_orders_reduced_horizon_quadrature = True
    engine = sk.Engine(config, geometry, viewing_geometry())
    uniform = engine.calculate_radiance(atmosphere(geometry, config)).radiance.values
    varying = engine.calculate_radiance(
        atmosphere(geometry, config, horizontal_slope=0.8)
    ).radiance.values

    assert uniform.shape == (1, 1, num_stokes)
    assert np.all(np.isfinite(uniform))
    assert np.all(np.isfinite(varying))
    assert uniform[0, 0, 0] > 0.0
    assert not np.isclose(varying[0, 0, 0], uniform[0, 0, 0], rtol=1.0e-4)


def test_2d_successive_orders_accepts_explicit_horizontal_source_angles():
    geometry = geometry2d()
    config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.NoSource
    )
    config.num_sza = 99
    config.successive_orders_horizontal_angle_grid_radians = np.array(
        [-0.55, -0.1, 0.2, 0.5]
    )

    result = sk.Engine(config, geometry, viewing_geometry()).calculate_radiance(
        atmosphere(geometry, config)
    )

    assert np.all(np.isfinite(result.radiance.values))
    assert result.radiance.values.item() > 0.0


def test_2d_successive_orders_rejects_source_angles_outside_geometry():
    config = successive_orders_config()
    config.successive_orders_horizontal_angle_grid_radians = np.array([-0.7, 0.0])

    with pytest.raises(ValueError, match="must lie inside the Geometry2D"):
        sk.Engine(config, geometry2d(), viewing_geometry())


@pytest.mark.parametrize(
    ("num_stokes", "reduced_horizon"), [(1, False), (1, True), (3, True)]
)
def test_2d_successive_orders_native_products_are_adjoint(
    num_stokes: int, reduced_horizon: bool
):
    geometry = geometry2d()
    config = successive_orders_config(
        num_stokes=num_stokes, single_scatter_source=sk.SingleScatterSource.NoSource
    )
    if reduced_horizon:
        config.num_sza = 3
        config.num_successive_orders_incoming = 37
        config.num_successive_orders_outgoing = 14
        config.successive_orders_reduced_horizon_quadrature = True
    config.num_successive_orders_iterations = 60
    config.successive_orders_relative_tolerance = 1.0e-11
    config.successive_orders_absolute_tolerance = 1.0e-13
    config.successive_orders_anderson_depth = 3
    linearization = sk.Engine(config, geometry, viewing_geometry()).linearize(
        atmosphere(geometry, config, horizontal_slope=0.4, calculate_derivatives=True)
    )

    assert linearization.backends == {
        "jvp": sk.LinearizationBackend.Native,
        "vjp": sk.LinearizationBackend.Native,
    }
    tangent = linearization.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = np.linspace(
        -2.0e-7, 3.0e-7, tangent["extinction"].size
    ).reshape(tangent["extinction"].shape)
    tangent["ssa"].data[:] = np.linspace(-0.015, 0.02, tangent["ssa"].size).reshape(
        tangent["ssa"].shape
    )
    jvp = linearization.jvp(tangent)

    cotangent = xr.ones_like(linearization.value)
    gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))
    lhs = float((jvp * cotangent).sum())
    rhs = float(
        (tangent["extinction"] * gradient["extinction"]).sum()
        + (tangent["ssa"] * gradient["ssa"]).sum()
    )

    assert np.all(np.isfinite(jvp))
    assert np.all(np.isfinite(gradient["extinction"]))
    assert np.all(np.isfinite(gradient["ssa"]))
    np.testing.assert_allclose(lhs, rhs, rtol=3.0e-8, atol=3.0e-11)


def test_2d_successive_orders_rejects_diffuse_refraction():
    config = successive_orders_config()
    config.multiple_scatter_refraction = True

    with pytest.raises(NotImplementedError, match="diffuse-ray refraction"):
        sk.Engine(config, geometry2d(), viewing_geometry())


def test_2d_successive_orders_supports_solar_refraction():
    geometry = geometry2d()
    geometry.refractive_index = np.array([1.001, 1.0004, 1.0])

    straight_config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.NoSource
    )
    refracted_config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.NoSource
    )
    refracted_config.solar_refraction = True

    straight = (
        sk.Engine(straight_config, geometry, viewing_geometry())
        .calculate_radiance(atmosphere(geometry, straight_config))
        .radiance.values
    )
    refracted = (
        sk.Engine(refracted_config, geometry, viewing_geometry())
        .calculate_radiance(atmosphere(geometry, refracted_config))
        .radiance.values
    )

    assert np.all(np.isfinite(refracted))
    assert not np.allclose(refracted, straight, rtol=1.0e-8, atol=0.0)

    refracted_config.num_successive_orders_iterations = 60
    refracted_config.successive_orders_relative_tolerance = 1.0e-11
    refracted_config.successive_orders_absolute_tolerance = 1.0e-13
    refracted_config.successive_orders_anderson_depth = 3
    linearization = sk.Engine(refracted_config, geometry, viewing_geometry()).linearize(
        atmosphere(
            geometry,
            refracted_config,
            horizontal_slope=0.3,
            calculate_derivatives=True,
        )
    )
    tangent = linearization.tangent_template[["extinction"]]
    tangent["extinction"].data[:] = np.linspace(
        -2.0e-7, 3.0e-7, tangent["extinction"].size
    ).reshape(tangent["extinction"].shape)
    cotangent = xr.ones_like(linearization.value)
    jvp = linearization.jvp(tangent)
    gradient = linearization.vjp(cotangent, parameters=("extinction",))
    np.testing.assert_allclose(
        float((jvp * cotangent).sum()),
        float((tangent["extinction"] * gradient["extinction"]).sum()),
        rtol=3.0e-8,
        atol=3.0e-11,
    )


def test_2d_successive_orders_unity_solar_refraction_matches_straight_paths():
    geometry = geometry2d()
    straight_config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.NoSource
    )
    refracted_config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.NoSource
    )
    refracted_config.solar_refraction = True

    straight = (
        sk.Engine(straight_config, geometry, viewing_geometry())
        .calculate_radiance(atmosphere(geometry, straight_config))
        .radiance.values
    )
    refracted = (
        sk.Engine(refracted_config, geometry, viewing_geometry())
        .calculate_radiance(atmosphere(geometry, refracted_config))
        .radiance.values
    )

    np.testing.assert_allclose(refracted, straight, rtol=2.0e-12, atol=1.0e-14)


@pytest.mark.parametrize(
    "single_source",
    [sk.SingleScatterSource.Exact, sk.SingleScatterSource.Table],
)
@pytest.mark.parametrize("num_stokes", [1, 3])
def test_2d_single_scatter_supports_solar_refraction(single_source, num_stokes):
    geometry = geometry2d()
    geometry.refractive_index = np.array([1.001, 1.0004, 1.0])
    straight_config = successive_orders_config(
        num_stokes=num_stokes,
        single_scatter_source=single_source,
        multiple_scatter_source=sk.MultipleScatterSource.NoSource,
    )
    refracted_config = successive_orders_config(
        num_stokes=num_stokes,
        single_scatter_source=single_source,
        multiple_scatter_source=sk.MultipleScatterSource.NoSource,
    )
    refracted_config.solar_refraction = True

    straight = (
        sk.Engine(straight_config, geometry, viewing_geometry())
        .calculate_radiance(atmosphere(geometry, straight_config))
        .radiance.values
    )
    refracted = (
        sk.Engine(refracted_config, geometry, viewing_geometry())
        .calculate_radiance(atmosphere(geometry, refracted_config))
        .radiance.values
    )

    assert np.all(np.isfinite(refracted))
    assert refracted[0, 0, 0] > 0.0
    assert not np.allclose(refracted, straight, rtol=1.0e-8, atol=0.0)


@pytest.mark.parametrize(
    "single_source",
    [sk.SingleScatterSource.Exact, sk.SingleScatterSource.Table],
)
@pytest.mark.parametrize("num_stokes", [1, 3])
def test_2d_single_scatter_unity_solar_refraction_matches_straight(
    single_source, num_stokes
):
    geometry = geometry2d()
    straight_config = successive_orders_config(
        num_stokes=num_stokes,
        single_scatter_source=single_source,
        multiple_scatter_source=sk.MultipleScatterSource.NoSource,
    )
    refracted_config = successive_orders_config(
        num_stokes=num_stokes,
        single_scatter_source=single_source,
        multiple_scatter_source=sk.MultipleScatterSource.NoSource,
    )
    refracted_config.solar_refraction = True

    straight = (
        sk.Engine(straight_config, geometry, viewing_geometry())
        .calculate_radiance(atmosphere(geometry, straight_config))
        .radiance.values
    )
    refracted = (
        sk.Engine(refracted_config, geometry, viewing_geometry())
        .calculate_radiance(atmosphere(geometry, refracted_config))
        .radiance.values
    )

    np.testing.assert_allclose(refracted, straight, rtol=2.0e-12, atol=1.0e-14)


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_2d_table_single_scatter_native_products_are_adjoint(num_stokes):
    geometry = geometry2d()
    geometry.refractive_index = np.array([1.001, 1.0004, 1.0])
    config = successive_orders_config(
        num_stokes=num_stokes,
        single_scatter_source=sk.SingleScatterSource.Table,
        multiple_scatter_source=sk.MultipleScatterSource.NoSource,
    )
    config.solar_refraction = True
    linearization = sk.Engine(config, geometry, viewing_geometry()).linearize(
        atmosphere(
            geometry,
            config,
            horizontal_slope=0.35,
            calculate_derivatives=True,
        )
    )

    assert linearization.backends == {
        "jvp": sk.LinearizationBackend.Native,
        "vjp": sk.LinearizationBackend.Native,
    }
    with pytest.raises(NotImplementedError, match="cannot materialize"):
        _ = linearization.jacobian
    tangent = linearization.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = np.linspace(
        -2.0e-7, 3.0e-7, tangent["extinction"].size
    ).reshape(tangent["extinction"].shape)
    tangent["ssa"].data[:] = np.linspace(-0.015, 0.02, tangent["ssa"].size).reshape(
        tangent["ssa"].shape
    )
    cotangent = xr.ones_like(linearization.value)
    jvp = linearization.jvp(tangent)

    finite_difference_engine = sk.Engine(config, geometry, viewing_geometry())
    epsilon = 1.0e-3

    def perturbed_radiance(sign: float) -> xr.DataArray:
        perturbed = atmosphere(
            geometry,
            config,
            horizontal_slope=0.35,
            calculate_derivatives=False,
        )
        perturbed.storage.total_extinction[:] += (
            sign
            * epsilon
            * tangent["extinction"].values.reshape(
                perturbed.storage.total_extinction.shape
            )
        )
        perturbed.storage.ssa[:] += (
            sign * epsilon * tangent["ssa"].values.reshape(perturbed.storage.ssa.shape)
        )
        return finite_difference_engine.calculate_radiance(perturbed).radiance

    finite_difference = (perturbed_radiance(1.0) - perturbed_radiance(-1.0)) / (
        2.0 * epsilon
    )
    xr.testing.assert_allclose(jvp, finite_difference, rtol=2.0e-7, atol=1.0e-12)

    gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))

    np.testing.assert_allclose(
        float((jvp * cotangent).sum()),
        float(
            (tangent["extinction"] * gradient["extinction"]).sum()
            + (tangent["ssa"] * gradient["ssa"]).sum()
        ),
        rtol=3.0e-8,
        atol=3.0e-11,
    )


def test_2d_table_single_scatter_can_share_solar_table_with_successive_orders():
    geometry = geometry2d()
    geometry.refractive_index = np.array([1.001, 1.0004, 1.0])
    config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.Table
    )
    config.solar_refraction = True
    result = sk.Engine(config, geometry, viewing_geometry()).calculate_radiance(
        atmosphere(geometry, config, horizontal_slope=0.2)
    )

    assert np.all(np.isfinite(result.radiance.values))
    assert result.radiance.values.item() > 0.0


@pytest.mark.parametrize(
    ("single_source", "multiple_source"),
    [
        (sk.SingleScatterSource.Table, sk.MultipleScatterSource.NoSource),
        (
            sk.SingleScatterSource.NoSource,
            sk.MultipleScatterSource.SuccessiveOrders,
        ),
    ],
)
def test_2d_solar_table_has_stable_tangent_topology(single_source, multiple_source):
    altitude_grid_m = np.linspace(0.0, 80_000.0, 80)
    horizontal_grid = np.linspace(-0.4, 0.4, 40)
    geometry = sk.Geometry2D(
        cos_sza=0.15,
        solar_azimuth=0.25,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=altitude_grid_m,
        horizontal_angle_grid_radians=horizontal_grid,
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.TangentAltitude(
            tangent_altitude_m=20_000.0,
            observer_altitude_m=150_000.0,
            horizontal_angle_radians=-0.3,
            viewing_azimuth_radians=0.0,
        )
    )
    config = successive_orders_config(
        single_scatter_source=single_source,
        multiple_scatter_source=multiple_source,
    )
    config.num_sza = 5
    config.successive_orders_altitude_grid_m = np.linspace(1_000.0, 79_000.0, 25)

    sk.Engine(config, geometry, viewing)


def test_2d_rejects_legacy_successive_orders_source():
    config = successive_orders_config(
        multiple_scatter_source=sk.MultipleScatterSource.SuccessiveOrdersLegacy
    )

    with pytest.raises(NotImplementedError, match="successive-orders"):
        sk.Engine(config, geometry2d(), viewing_geometry())
