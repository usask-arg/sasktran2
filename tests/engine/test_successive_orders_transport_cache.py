from __future__ import annotations

import json

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr


def _config(cache_count: int, *, threads: int = 1, converged: bool = False):
    config = sk.Config()
    config.num_threads = threads
    config.num_streams = 4
    config.num_singlescatter_moments = 4
    config.num_sza = 3
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    config.num_successive_orders_incoming = 6
    config.num_successive_orders_outgoing = 6
    config.num_successive_orders_iterations = 60 if converged else 2
    config.successive_orders_relative_tolerance = 1.0e-12 if converged else 0.0
    config.successive_orders_absolute_tolerance = 1.0e-14 if converged else 0.0
    config.delta_m_scaling = False
    config.successive_orders_transport_cache_wavelengths = cache_count
    return config


def _geometry():
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=np.array([0.0, 10_000.0, 30_000.0]),
        horizontal_angle_grid_radians=np.array([-0.6, 0.0, 0.6]),
    )


def _viewing():
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0.2, 0.7, 100_000.0))
    viewing.add_ray(sk.TangentAltitudeSolar(10_000.0, -0.3, 100_000.0, 0.6))
    return viewing


def _scene(geometry, config, *, uniform_phase: bool = False):
    scene = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([410.0, 530.0, 690.0]),
        calculate_derivatives=True,
        legendre_derivative=False,
    )
    horizontal, altitude = np.meshgrid(
        [-0.6, 0.0, 0.6], [0.0, 10_000.0, 30_000.0], indexing="ij"
    )
    scene.storage.total_extinction[:] = (
        1.5e-5 * np.exp(-altitude / 15_000.0) * (1.0 + 0.3 * horizontal)
    ).reshape(-1, 1) * [0.8, 1.1, 1.5]
    scene.storage.ssa[:] = (
        0.89 - 0.08 * altitude / 30_000.0 + 0.02 * horizontal
    ).reshape(-1, 1) - [0.0, 0.03, 0.07]
    scene.leg_coeff.a1[0] = 1.0
    phase = 0.3 if uniform_phase else (0.3 + 0.04 * horizontal).reshape(-1, 1)
    scene.leg_coeff.a1[2] = phase + np.array([0.0, 0.05, 0.10])
    if config.num_stokes == 3:
        scene.leg_coeff.a2[2] = 3.0
        scene.leg_coeff.b1[2] = -np.sqrt(6.0) / 2.0
    scene.surface.albedo[:] = [0.06, 0.12, 0.23]
    scene.mark_changed()
    return scene


def _assert_bits(actual, expected):
    if isinstance(actual, xr.Dataset):
        assert set(actual.data_vars) == set(expected.data_vars)
        for name in actual.data_vars:
            _assert_bits(actual[name], expected[name])
        return
    assert actual.dims == expected.dims
    assert actual.shape == expected.shape
    np.testing.assert_array_equal(
        actual.values.view(np.uint64), expected.values.view(np.uint64)
    )


def _products(
    engine, scene, *, exact_repeats: bool = True, parallel_gradient: bool = False
):
    linearization = engine.linearize(scene)
    tangent = linearization.tangent_template[["extinction", "ssa"]]
    tangent.extinction.values[:] = np.linspace(
        -2.0e-7, 3.0e-7, tangent.extinction.size
    ).reshape(tangent.extinction.shape)
    tangent.ssa.values[:] = 0.012
    cotangent = xr.ones_like(linearization.value)
    cotangent.values[:] = np.linspace(0.35, 1.1, cotangent.size).reshape(
        cotangent.shape
    )
    gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))
    jvp = linearization.jvp(tangent)
    repeated_gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))
    repeated_jvp = linearization.jvp(tangent)
    if exact_repeats:
        if parallel_gradient:
            # Rayon can assign wavelengths to different workers between calls.
            # The existing mapped-gradient reduction then groups sums differently.
            xr.testing.assert_allclose(
                repeated_gradient, gradient, rtol=2.0e-12, atol=2.0e-13
            )
        else:
            _assert_bits(repeated_gradient, gradient)
        _assert_bits(repeated_jvp, jvp)
    else:
        xr.testing.assert_allclose(
            repeated_gradient, gradient, rtol=2.0e-12, atol=2.0e-13
        )
        xr.testing.assert_allclose(repeated_jvp, jvp, rtol=2.0e-12, atol=2.0e-13)
    return linearization.value.copy(), jvp, gradient


@pytest.mark.parametrize("cache_count", [1, 3, 10])
@pytest.mark.parametrize("threads", [1, 2])
@pytest.mark.parametrize("converged", [False, True])
def test_transport_cache_preserves_products_and_update_history(
    cache_count, threads, converged
):
    geometry = _geometry()
    cached_config = _config(cache_count, threads=threads, converged=converged)
    reference_config = _config(0, threads=threads, converged=converged)
    cached_engine = sk.Engine(cached_config, geometry, _viewing())
    reference_engine = sk.Engine(reference_config, geometry, _viewing())
    scenes = [_scene(geometry, cached_config), _scene(geometry, reference_config)]
    saved = [
        (
            scene.storage.total_extinction.copy(),
            scene.storage.ssa.copy(),
            scene.surface.albedo.copy(),
        )
        for scene in scenes
    ]
    for evaluation in range(4):
        for scene, initial in zip(scenes, saved, strict=True):
            if evaluation == 1:
                scene.surface.albedo[:] += 0.11
                scene.mark_changed()
            elif evaluation == 2:
                scene.storage.total_extinction[:] *= 1.25
                scene.storage.ssa[:] -= 0.04
                scene.mark_changed()
            elif evaluation == 3:
                scene.storage.total_extinction[:] = initial[0]
                scene.storage.ssa[:] = initial[1]
                scene.surface.albedo[:] = initial[2]
                scene.mark_changed()
        actual = _products(cached_engine, scenes[0], parallel_gradient=threads > 1)
        expected = _products(reference_engine, scenes[1], parallel_gradient=threads > 1)
        for index, (got, wanted) in enumerate(zip(actual, expected, strict=True)):
            if threads > 1 and index == 2:
                xr.testing.assert_allclose(got, wanted, rtol=2.0e-12, atol=2.0e-13)
            else:
                _assert_bits(got, wanted)


@pytest.mark.parametrize(("cache_count", "maximum_vectors"), [(1, 2), (3, 3), (10, 3)])
def test_transport_cache_hits_have_exclusive_bounded_buffers(
    cache_count, maximum_vectors, monkeypatch, capfd
):
    monkeypatch.setenv("SASKTRAN2_PROFILE_MEMORY", "1")
    config = _config(cache_count)
    geometry = _geometry()
    scene = _scene(geometry, config, uniform_phase=True)
    _products(sk.Engine(config, geometry, _viewing()), scene)
    diagnostic = capfd.readouterr().err
    events = [
        json.loads(line.removeprefix("SASKTRAN2_MEMORY "))
        for line in diagnostic.splitlines()
        if line.startswith("SASKTRAN2_MEMORY ")
    ]
    events = [event for event in events if event["kind"] == "scalar_transport_cache"]
    assert events
    assert any(event["hit"] for event in events)
    assert {event["cache_wavelengths"] for event in events} == {min(cache_count, 3)}
    assert max(event["retained_vectors"] for event in events) == maximum_vectors
    assert all(
        event["retained_bytes"] == event["inactive_bytes"] + event["active_bytes"]
        for event in events
    )


def test_transport_cache_defaults_to_disabled_and_rejects_negative_count():
    config = sk.Config()
    assert config.successive_orders_transport_cache_wavelengths == 0
    with pytest.raises(ValueError, match="non-negative"):
        config.successive_orders_transport_cache_wavelengths = -1


def test_transport_cache_zero_order_products_match_disabled_cache():
    geometry = _geometry()
    results = []
    for count in (0, 3):
        config = _config(count)
        config.num_successive_orders_iterations = 0
        results.append(
            _products(sk.Engine(config, geometry, _viewing()), _scene(geometry, config))
        )
    for actual, expected in zip(*results, strict=True):
        _assert_bits(actual, expected)


def test_transport_cache_preserves_source_thread_products():
    geometry = _geometry()
    results = []
    for count in (0, 3):
        config = _config(count, threads=2)
        config.threading_model = sk.ThreadingModel.Source
        results.append(
            _products(
                sk.Engine(config, geometry, _viewing()),
                _scene(geometry, config),
                exact_repeats=False,
            )
        )
    for actual, expected in zip(*results, strict=True):
        # Source-thread reductions may arrive in a different order between
        # separate engines. The cache cannot change their numerical result.
        xr.testing.assert_allclose(actual, expected, rtol=2.0e-12, atol=2.0e-13)


def test_transport_cache_positive_setting_keeps_polarized_fallback():
    geometry = _geometry()
    results = []
    for count in (0, 3):
        config = _config(count)
        config.num_stokes = 3
        result = sk.Engine(config, geometry, _viewing()).calculate_radiance(
            _scene(geometry, config)
        )
        results.append(result.radiance)
    _assert_bits(*results)
