from __future__ import annotations

import numpy as np
import pytest
from sasktran2.test_util.source_integration import (
    accuracy_metrics,
    benchmark_callable,
    minimum_layers_for_tolerance,
    normalized_error,
    observed_orders,
)


def test_normalized_error_uses_reference_floor_near_zero():
    reference = np.array([10.0, 1.0, 0.0])
    candidate = np.array([11.0, 0.9, 1e-11])

    error, floor = normalized_error(candidate, reference, floor_fraction=1e-12)

    assert floor == pytest.approx(1e-11)
    np.testing.assert_allclose(error, [0.1, 0.1, 1.0])


def test_accuracy_metrics_summarize_normalized_and_absolute_error():
    metrics = accuracy_metrics(
        np.array([1.1, 1.8]), np.array([1.0, 2.0]), floor_fraction=1e-12
    )

    assert metrics.max_absolute == pytest.approx(0.2)
    assert metrics.max_normalized == pytest.approx(0.1)
    assert metrics.rms_normalized == pytest.approx(0.1)
    assert metrics.p95_normalized == pytest.approx(0.1)


def test_normalized_error_rejects_incompatible_shapes():
    with pytest.raises(ValueError, match="same shape"):
        normalized_error(np.ones(2), np.ones(3))


def test_observed_orders_support_nonuniform_refinement():
    errors = [0.16, 0.04, 0.01]

    np.testing.assert_allclose(observed_orders(errors), [2.0, 2.0])
    np.testing.assert_allclose(
        observed_orders([0.09, 0.01], refinement_ratios=[3.0]), [2.0]
    )


def test_minimum_layers_for_tolerance_is_order_independent():
    assert minimum_layers_for_tolerance([80, 20, 40], [0.01, 0.2, 0.04], 0.05) == 40
    assert minimum_layers_for_tolerance([20, 40], [0.2, 0.1], 0.05) is None


def test_benchmark_callable_reports_per_call_samples():
    call_count = 0

    def increment() -> None:
        nonlocal call_count
        call_count += 1

    timing = benchmark_callable(
        increment,
        warmups=1,
        samples=3,
        minimum_sample_seconds=0,
    )

    assert timing.samples == 3
    assert timing.calls_per_sample == 1
    assert timing.minimum_seconds >= 0
    assert timing.median_seconds >= timing.minimum_seconds
    assert call_count == 5
