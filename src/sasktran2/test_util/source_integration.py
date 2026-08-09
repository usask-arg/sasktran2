from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass
from math import ceil
from time import perf_counter
from typing import Any

import numpy as np


@dataclass(frozen=True)
class AccuracyMetrics:
    """Error statistics normalized by the magnitude of a dense reference."""

    max_absolute: float
    max_normalized: float
    rms_normalized: float
    p95_normalized: float
    reference_floor: float

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class TimingMetrics:
    """Robust per-call timing statistics."""

    median_seconds: float
    mad_seconds: float
    minimum_seconds: float
    samples: int
    calls_per_sample: int

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


def normalized_error(
    candidate: np.ndarray,
    reference: np.ndarray,
    *,
    floor_fraction: float = 1e-12,
) -> tuple[np.ndarray, float]:
    """Return pointwise error with a stable relative-error denominator.

    The denominator is ``max(abs(reference), floor_fraction * max(abs(reference)))``.
    This retains sensitivity away from zero without allowing insignificant reference
    values to dominate the reported error.
    """

    candidate = np.asarray(candidate, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if candidate.shape != reference.shape:
        msg = (
            "candidate and reference must have the same shape, got "
            f"{candidate.shape} and {reference.shape}"
        )
        raise ValueError(msg)
    if floor_fraction <= 0:
        msg = "floor_fraction must be positive"
        raise ValueError(msg)
    if not np.all(np.isfinite(candidate)) or not np.all(np.isfinite(reference)):
        msg = "candidate and reference must contain only finite values"
        raise ValueError(msg)

    reference_scale = float(np.max(np.abs(reference), initial=0.0))
    reference_floor = max(reference_scale * floor_fraction, np.finfo(float).tiny)
    denominator = np.maximum(np.abs(reference), reference_floor)
    return np.abs(candidate - reference) / denominator, reference_floor


def accuracy_metrics(
    candidate: np.ndarray,
    reference: np.ndarray,
    *,
    floor_fraction: float = 1e-12,
) -> AccuracyMetrics:
    """Summarize absolute and normalized error against a dense reference."""

    candidate = np.asarray(candidate, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    error, reference_floor = normalized_error(
        candidate, reference, floor_fraction=floor_fraction
    )
    absolute_error = np.abs(candidate - reference)
    return AccuracyMetrics(
        max_absolute=float(np.max(absolute_error, initial=0.0)),
        max_normalized=float(np.max(error, initial=0.0)),
        rms_normalized=float(np.sqrt(np.mean(np.square(error)))),
        p95_normalized=float(np.percentile(error, 95)),
        reference_floor=reference_floor,
    )


def observed_orders(
    errors: Sequence[float], refinement_ratios: float | Sequence[float] = 2.0
) -> np.ndarray:
    """Estimate convergence order between successive grid refinements."""

    error_values = np.asarray(errors, dtype=np.float64)
    if error_values.ndim != 1 or error_values.size < 2:
        msg = "errors must be a one-dimensional sequence with at least two values"
        raise ValueError(msg)

    ratios = np.asarray(refinement_ratios, dtype=np.float64)
    if ratios.ndim == 0:
        ratios = np.full(error_values.size - 1, ratios)
    if ratios.shape != (error_values.size - 1,):
        msg = "refinement_ratios must be scalar or have len(errors) - 1 values"
        raise ValueError(msg)
    if np.any(ratios <= 1):
        msg = "refinement ratios must be greater than one"
        raise ValueError(msg)

    orders = np.full(error_values.size - 1, np.nan)
    valid = (error_values[:-1] > 0) & (error_values[1:] > 0)
    orders[valid] = np.log(error_values[:-1][valid] / error_values[1:][valid]) / np.log(
        ratios[valid]
    )
    return orders


def minimum_layers_for_tolerance(
    layer_counts: Sequence[int], errors: Sequence[float], tolerance: float
) -> int | None:
    """Return the smallest layer count whose error meets ``tolerance``."""

    layers = np.asarray(layer_counts, dtype=np.int64)
    error_values = np.asarray(errors, dtype=np.float64)
    if layers.ndim != 1 or layers.shape != error_values.shape:
        msg = "layer_counts and errors must be one-dimensional and have equal length"
        raise ValueError(msg)
    if tolerance < 0:
        msg = "tolerance must be non-negative"
        raise ValueError(msg)

    passing = layers[error_values <= tolerance]
    return None if passing.size == 0 else int(np.min(passing))


def benchmark_callable(
    function: Callable[[], Any],
    *,
    warmups: int = 3,
    samples: int = 15,
    minimum_sample_seconds: float = 0.05,
    maximum_calls_per_sample: int = 10_000,
) -> TimingMetrics:
    """Benchmark a callable using warmups and median/MAD statistics."""

    if warmups < 0 or samples < 1:
        msg = "warmups must be non-negative and samples must be positive"
        raise ValueError(msg)
    if minimum_sample_seconds < 0 or maximum_calls_per_sample < 1:
        msg = "minimum_sample_seconds must be non-negative and call limit positive"
        raise ValueError(msg)

    for _ in range(warmups):
        function()

    start = perf_counter()
    function()
    probe_seconds = max(perf_counter() - start, np.finfo(float).eps)
    calls_per_sample = min(
        maximum_calls_per_sample,
        max(1, ceil(minimum_sample_seconds / probe_seconds)),
    )

    per_call_seconds = np.empty(samples)
    for sample_index in range(samples):
        start = perf_counter()
        for _ in range(calls_per_sample):
            function()
        per_call_seconds[sample_index] = (perf_counter() - start) / calls_per_sample

    median = float(np.median(per_call_seconds))
    return TimingMetrics(
        median_seconds=median,
        mad_seconds=float(np.median(np.abs(per_call_seconds - median))),
        minimum_seconds=float(np.min(per_call_seconds)),
        samples=samples,
        calls_per_sample=calls_per_sample,
    )
