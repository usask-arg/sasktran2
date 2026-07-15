"""Build the SASKTRAN2 Baum V3.6 ice-crystal scattering database.

The Baum files store P11 as an absolute normalized phase function and the other
five matrix elements as ratios to P11. This converter restores the absolute
matrix, computes the six Greek-coefficient families, and streams a rectangular,
zero-padded database to NetCDF4. When ``--default-output`` is supplied, the first
256 moments (configurable with ``--default-moments``) are then copied exactly into
a second lightweight runtime database.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr
from netCDF4 import Dataset
from numpy.polynomial.legendre import leggauss
from sasktran2._core_rust import WignerD
from scipy.sparse import csr_matrix

LOGGER = logging.getLogger(__name__)

DEFAULT_DATABASE_MOMENTS = 256
DEFAULT_MAX_NORMALIZED_TOLERANCE = 1e-3
DEFAULT_STANDARD_DATABASE_KEY = "cross_sections/ice/baum_ice_crystals_v3_6.nc"
FULL_STANDARD_DATABASE_KEY = "cross_sections/ice/baum_ice_crystals_v3_6_full.nc"

MODEL_FILES = {
    "general_habit_mixture": (
        "GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc"
    ),
    "aggregate_solid_columns": (
        "AggregateSolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc"
    ),
    "solid_columns": "SolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc",
}

PHASE_MATRIX_ELEMENTS = ("p11", "p12", "p22", "p33", "p34", "p44")
COEFFICIENT_NAMES = ("lm_a1", "lm_a2", "lm_a3", "lm_a4", "lm_b1", "lm_b2")
_SOURCE_PHASE_VARIABLES = (
    "p11_phase_function",
    "p21_phase_function",
    "p22_phase_function",
    "p33_phase_function",
    "p43_phase_function",
    "p44_phase_function",
)
_WIGNER_KEYS = ("00", "22", "2m2", "02")


class ConvergenceError(RuntimeError):
    """Raised when a phase matrix does not converge before the safety cap."""

    def __init__(
        self,
        cell_index: int,
        element: str,
        rms_error: float,
        max_error: float,
        safety_cap: int,
    ) -> None:
        self.cell_index = cell_index
        self.element = element
        self.rms_error = rms_error
        self.max_error = max_error
        self.safety_cap = safety_cap
        super().__init__(
            f"cell {cell_index} element {element} did not converge by "
            f"{safety_cap} moments (relative RMS={rms_error:.6g}, "
            f"max normalized error={max_error:.6g})"
        )


@dataclass
class TransformOperator:
    """Precomputed quadrature transforms for one native phase-angle grid."""

    phase_angles_deg: np.ndarray
    cos_angles: np.ndarray
    solid_angle_weights: np.ndarray
    kernels: dict[str, np.ndarray]
    native_basis: dict[str, np.ndarray]

    @property
    def max_moments(self) -> int:
        return self.kernels["00"].shape[0]


@dataclass
class TransformResult:
    coefficients: np.ndarray
    required_moments: np.ndarray
    converged: np.ndarray
    rms_error: np.ndarray
    max_normalized_error: np.ndarray
    a1_normalization: np.ndarray


def reconstruct_absolute_phase_matrix(
    p11: np.ndarray,
    p21_over_p11: np.ndarray,
    p22_over_p11: np.ndarray,
    p33_over_p11: np.ndarray,
    p43_over_p11: np.ndarray,
    p44_over_p11: np.ndarray,
) -> np.ndarray:
    """Restore absolute SASKTRAN2 phase elements from the Baum representation."""
    arrays = np.broadcast_arrays(
        np.asarray(p11, dtype=np.float64),
        np.asarray(p21_over_p11, dtype=np.float64),
        np.asarray(p22_over_p11, dtype=np.float64),
        np.asarray(p33_over_p11, dtype=np.float64),
        np.asarray(p43_over_p11, dtype=np.float64),
        np.asarray(p44_over_p11, dtype=np.float64),
    )
    p11_abs, p21, p22, p33, p43, p44 = arrays
    return np.stack(
        (
            p11_abs,
            p11_abs * p21,
            p11_abs * p22,
            p11_abs * p33,
            -p11_abs * p43,
            p11_abs * p44,
        ),
        axis=-2,
    )


def _validate_angle_grid(phase_angles_deg: np.ndarray) -> np.ndarray:
    angles = np.asarray(phase_angles_deg, dtype=np.float64)
    if angles.ndim != 1 or len(angles) < 2:
        msg = "phase_angles must be a one-dimensional array with at least two values"
        raise ValueError(msg)
    if not np.all(np.isfinite(angles)) or np.any(np.diff(angles) <= 0):
        msg = "phase_angles must be finite and strictly increasing"
        raise ValueError(msg)
    if not np.isclose(angles[0], 0.0) or not np.isclose(angles[-1], 180.0):
        msg = "phase_angles must span 0 to 180 degrees"
        raise ValueError(msg)
    return angles


def _quadrature_projection(
    phase_angles_deg: np.ndarray,
    max_moments: int,
    oversampling: float,
    minimum_order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, csr_matrix]:
    """Construct composite Gauss-Legendre nodes and a linear interpolation map."""
    angles = _validate_angle_grid(phase_angles_deg)
    if max_moments <= 0 or oversampling <= 0 or minimum_order <= 0:
        msg = "quadrature settings must be positive"
        raise ValueError(msg)

    # Work in ascending cos(theta), matching phase arrays reversed from the Baum
    # angle dimension. Each quadrature interval is exactly one native interval.
    cos_native = np.cos(np.deg2rad(angles))[::-1]
    theta_native = np.arccos(np.clip(cos_native, -1.0, 1.0))
    node_parts: list[np.ndarray] = []
    weight_parts: list[np.ndarray] = []
    left_parts: list[np.ndarray] = []
    fraction_parts: list[np.ndarray] = []
    rule_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    for left in range(len(cos_native) - 1):
        x0 = cos_native[left]
        x1 = cos_native[left + 1]
        delta_theta = abs(theta_native[left + 1] - theta_native[left])
        order = max(
            minimum_order,
            int(np.ceil(oversampling * max_moments * delta_theta / np.pi)),
        )
        if order not in rule_cache:
            rule_cache[order] = leggauss(order)
        canonical_nodes, canonical_weights = rule_cache[order]
        nodes = 0.5 * (canonical_nodes + 1.0) * (x1 - x0) + x0
        weights = 0.5 * canonical_weights * (x1 - x0)

        node_parts.append(nodes)
        weight_parts.append(weights)
        left_parts.append(np.full(order, left, dtype=np.int32))
        fraction_parts.append((nodes - x0) / (x1 - x0))

    nodes = np.concatenate(node_parts)
    weights = np.concatenate(weight_parts)
    left_indices = np.concatenate(left_parts)
    fractions = np.concatenate(fraction_parts)
    rows = np.arange(len(nodes), dtype=np.int64)
    projection = csr_matrix(
        (
            np.concatenate((1.0 - fractions, fractions)),
            (
                np.concatenate((rows, rows)),
                np.concatenate((left_indices, left_indices + 1)),
            ),
        ),
        shape=(len(nodes), len(cos_native)),
    )

    if not np.isclose(weights.sum(), 2.0, rtol=1e-12, atol=1e-12):
        msg = "composite quadrature weights do not integrate the full cosine domain"
        raise RuntimeError(msg)

    return cos_native, nodes, weights, projection


def _basis_block(
    wigner: WignerD, theta: np.ndarray, start: int, stop: int
) -> np.ndarray:
    return np.vstack(
        [
            np.asarray(wigner.d(theta, int(order)), dtype=np.float64)
            for order in range(start, stop)
        ]
    )


class _WignerBasisSequence:
    """Generate consecutive Wigner-d orders with one recurrence traversal."""

    def __init__(self, m: int, n: int, theta: np.ndarray) -> None:
        self._m = m
        self._n = n
        self._lmin = max(abs(m), abs(n))
        self._next_order = 0
        self._x = np.cos(np.asarray(theta, dtype=np.float64))

        zeta = 1 if n >= m or (m - n) % 2 == 0 else -1
        factorial_ratio = math.factorial(2 * self._lmin) / (
            math.factorial(abs(m - n)) * math.factorial(abs(m + n))
        )
        start_factor = zeta * 2.0 ** (-self._lmin) * math.sqrt(factorial_ratio)
        self._value = (
            start_factor
            * np.power(1.0 - self._x, abs(m - n) / 2)
            * np.power(1.0 + self._x, abs(m + n) / 2)
        )
        self._previous = np.zeros_like(self._value)

    def block(self, start: int, stop: int) -> np.ndarray:
        if start != self._next_order or stop <= start:
            msg = "Wigner basis blocks must be nonempty, consecutive, and ordered"
            raise ValueError(msg)

        result = np.zeros((stop - start, len(self._x)), dtype=np.float64)
        for order in range(start, stop):
            if order < self._lmin:
                continue
            if order > self._lmin:
                order_float = float(order)
                mm = float(order * order - self._m * self._m)
                if self._n == 0:
                    multiplier = 1.0 / (math.sqrt(mm) * order_float)
                    current_factor = (2 * order - 1) * order_float * self._x
                    prior_factor = order_float * math.sqrt(
                        (order - 1) ** 2 - self._m * self._m
                    )
                else:
                    nn = float(order * order - self._n * self._n)
                    multiplier = 1.0 / ((order - 1) * math.sqrt(mm) * math.sqrt(nn))
                    current_factor = (2 * order - 1) * (
                        order_float * (order_float - 1.0) * self._x - self._n * self._m
                    )
                    prior_factor = (
                        order_float
                        * math.sqrt((order - 1) ** 2 - self._m * self._m)
                        * math.sqrt((order - 1) ** 2 - self._n * self._n)
                    )

                previous_value = self._value
                self._value = multiplier * (
                    current_factor * self._value - prior_factor * self._previous
                )
                self._previous = previous_value

            result[order - start] = self._value

        self._next_order = stop
        return result


def build_transform_operator(
    phase_angles_deg: np.ndarray,
    max_moments: int = 16384,
    block_size: int = 128,
    quadrature_oversampling: float = 2.0,
    minimum_quadrature_order: int = 4,
) -> TransformOperator:
    """Precompute float64 transforms shared by every Baum phase matrix."""
    if block_size <= 0:
        msg = "block_size must be positive"
        raise ValueError(msg)
    cos_native, nodes, weights, projection = _quadrature_projection(
        phase_angles_deg,
        max_moments,
        quadrature_oversampling,
        minimum_quadrature_order,
    )
    theta_nodes = np.arccos(np.clip(nodes, -1.0, 1.0))
    theta_native = np.arccos(np.clip(cos_native, -1.0, 1.0))
    wigner_indices = {"00": (0, 0), "22": (2, 2), "2m2": (2, -2), "02": (0, 2)}
    node_sequences = {
        key: _WignerBasisSequence(m, n, theta_nodes)
        for key, (m, n) in wigner_indices.items()
    }
    native_sequences = {
        key: _WignerBasisSequence(m, n, theta_native)
        for key, (m, n) in wigner_indices.items()
    }
    kernels = {
        key: np.empty((max_moments, len(cos_native)), dtype=np.float64)
        for key in _WIGNER_KEYS
    }
    native_basis = {
        key: np.empty((max_moments, len(cos_native)), dtype=np.float64)
        for key in _WIGNER_KEYS
    }

    LOGGER.info(
        "Building transforms for %d moments using %d composite quadrature nodes",
        max_moments,
        len(nodes),
    )
    for start in range(0, max_moments, block_size):
        stop = min(start + block_size, max_moments)
        orders = np.arange(start, stop, dtype=np.float64)
        coefficient_weight = (2.0 * orders + 1.0) / 2.0
        for key in _WIGNER_KEYS:
            basis_nodes = node_sequences[key].block(start, stop)
            basis_nodes *= coefficient_weight[:, np.newaxis] * weights
            kernels[key][start:stop] = projection.T.dot(basis_nodes.T).T
            native_basis[key][start:stop] = native_sequences[key].block(start, stop)
        LOGGER.info("Built moment transform [%d, %d)", start, stop)

    solid_angle_weights = np.empty_like(cos_native)
    delta = np.diff(cos_native)
    solid_angle_weights[0] = delta[0] / 2.0
    solid_angle_weights[-1] = delta[-1] / 2.0
    solid_angle_weights[1:-1] = (cos_native[2:] - cos_native[:-2]) / 2.0

    return TransformOperator(
        phase_angles_deg=np.asarray(phase_angles_deg, dtype=np.float64),
        cos_angles=cos_native,
        solid_angle_weights=solid_angle_weights,
        kernels=kernels,
        native_basis=native_basis,
    )


def calculate_greek_coefficients(
    phase_matrix: np.ndarray,
    cos_angles: np.ndarray,
    angle_weights: np.ndarray,
    num_moments: int,
    block_size: int = 128,
    normalize: bool = True,
) -> np.ndarray:
    """Calculate Greek coefficients directly on a supplied quadrature grid.

    This direct path is primarily an independent convention oracle for the Rust
    Mie integrator. The production converter uses ``TransformOperator`` to reuse
    its composite-quadrature kernels across all database cells.
    """
    phase = np.asarray(phase_matrix, dtype=np.float64)
    cos_angles = np.asarray(cos_angles, dtype=np.float64)
    angle_weights = np.asarray(angle_weights, dtype=np.float64)
    if num_moments <= 0 or block_size <= 0:
        msg = "num_moments and block_size must be positive"
        raise ValueError(msg)
    if phase.shape[-2:] != (6, len(cos_angles)):
        msg = "phase_matrix must have trailing dimensions (6, num_angles)"
        raise ValueError(msg)
    if angle_weights.shape != cos_angles.shape:
        msg = "angle_weights must match cos_angles"
        raise ValueError(msg)

    flat_phase = phase.reshape((-1, 6, len(cos_angles)))
    coefficients = np.zeros((len(flat_phase), 6, num_moments), dtype=np.float64)
    theta = np.arccos(np.clip(cos_angles, -1.0, 1.0))
    wigners = {
        "00": WignerD(0, 0),
        "22": WignerD(2, 2),
        "2m2": WignerD(2, -2),
        "02": WignerD(0, 2),
    }

    for start in range(0, num_moments, block_size):
        stop = min(start + block_size, num_moments)
        orders = np.arange(start, stop, dtype=np.float64)
        factor = (2.0 * orders + 1.0)[:, np.newaxis] / 2.0
        transforms = {
            key: _basis_block(wigner, theta, start, stop) * (factor * angle_weights)
            for key, wigner in wigners.items()
        }
        coefficients[:, 0, start:stop] = flat_phase[:, 0] @ transforms["00"].T
        coefficients[:, 3, start:stop] = flat_phase[:, 5] @ transforms["00"].T
        coefficients[:, 4, start:stop] = flat_phase[:, 1] @ transforms["02"].T
        coefficients[:, 5, start:stop] = -(flat_phase[:, 4] @ transforms["02"].T)
        plus = (flat_phase[:, 2] + flat_phase[:, 3]) @ transforms["22"].T
        minus = (flat_phase[:, 2] - flat_phase[:, 3]) @ transforms["2m2"].T
        coefficients[:, 1, start:stop] = (plus + minus) / 2.0
        coefficients[:, 2, start:stop] = (plus - minus) / 2.0

    if normalize:
        normalization = coefficients[:, 0, 0]
        if np.any(~np.isfinite(normalization)) or np.any(normalization == 0):
            msg = "A1[0] must be finite and nonzero"
            raise ValueError(msg)
        coefficients /= normalization[:, np.newaxis, np.newaxis]
    return coefficients.reshape((*phase.shape[:-2], 6, num_moments))


def reconstruct_phase_matrix(
    coefficients: np.ndarray, cos_angles: np.ndarray, block_size: int = 128
) -> np.ndarray:
    """Reconstruct the six absolute phase-matrix elements from Greek coefficients."""
    coeff = np.asarray(coefficients, dtype=np.float64)
    cos_angles = np.asarray(cos_angles, dtype=np.float64)
    if block_size <= 0:
        msg = "block_size must be positive"
        raise ValueError(msg)
    if coeff.shape[-2] != 6:
        msg = "coefficients must have a six-element coefficient-family dimension"
        raise ValueError(msg)
    num_moments = coeff.shape[-1]
    flat_coeff = coeff.reshape((-1, 6, num_moments))
    phase = np.zeros((len(flat_coeff), 6, len(cos_angles)), dtype=np.float64)
    theta = np.arccos(np.clip(cos_angles, -1.0, 1.0))
    wigners = {
        "00": WignerD(0, 0),
        "22": WignerD(2, 2),
        "2m2": WignerD(2, -2),
        "02": WignerD(0, 2),
    }

    for start in range(0, num_moments, block_size):
        stop = min(start + block_size, num_moments)
        basis = {
            key: _basis_block(wigner, theta, start, stop)
            for key, wigner in wigners.items()
        }
        local = flat_coeff[:, :, start:stop]
        phase[:, 0] += local[:, 0] @ basis["00"]
        phase[:, 5] += local[:, 3] @ basis["00"]
        phase[:, 1] += local[:, 4] @ basis["02"]
        phase[:, 4] -= local[:, 5] @ basis["02"]
        plus = (local[:, 1] + local[:, 2]) @ basis["22"]
        minus = (local[:, 1] - local[:, 2]) @ basis["2m2"]
        phase[:, 2] += (plus + minus) / 2.0
        phase[:, 3] += (plus - minus) / 2.0

    return phase.reshape((*coeff.shape[:-2], 6, len(cos_angles)))


def _phase_errors(
    reconstructed: np.ndarray,
    target: np.ndarray,
    solid_angle_weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    difference = reconstructed - target
    numerator = np.sum(difference**2 * solid_angle_weights, axis=-1)
    denominator = np.sum(target**2 * solid_angle_weights, axis=-1)
    rms = np.sqrt(
        np.divide(
            numerator,
            denominator,
            out=np.zeros_like(numerator),
            where=denominator > np.finfo(np.float64).tiny,
        )
    )
    rms = np.where(
        denominator > np.finfo(np.float64).tiny,
        rms,
        np.sqrt(numerator),
    )

    value_range = np.ptp(target, axis=-1)
    value_scale = np.max(np.abs(target), axis=-1)
    range_floor = np.finfo(np.float64).eps * np.maximum(value_scale, 1.0)
    max_denominator = np.where(value_range > range_floor, value_range, value_scale)
    max_difference = np.max(np.abs(difference), axis=-1)
    max_normalized = np.divide(
        max_difference,
        max_denominator,
        out=max_difference.copy(),
        where=max_denominator > np.finfo(np.float64).tiny,
    )
    return rms, max_normalized


def adaptive_transform_batch(
    phase_matrix: np.ndarray,
    operator: TransformOperator,
    block_size: int = 128,
    relative_rms_tolerance: float = 1e-4,
    max_normalized_tolerance: float = DEFAULT_MAX_NORMALIZED_TOLERANCE,
    consecutive_passes: int = 2,
    allow_unconverged_at_cap: bool = False,
) -> TransformResult:
    """Transform a batch adaptively, optionally retaining cells at the hard cap."""
    phase = np.asarray(phase_matrix, dtype=np.float64)
    if phase.ndim == 2:
        phase = phase[np.newaxis, ...]
    if phase.ndim != 3 or phase.shape[1:] != (6, len(operator.phase_angles_deg)):
        msg = "phase_matrix must have shape (cell, 6, phase_angle)"
        raise ValueError(msg)
    if block_size <= 0 or consecutive_passes <= 0:
        msg = "block_size and consecutive_passes must be positive"
        raise ValueError(msg)
    if relative_rms_tolerance <= 0 or max_normalized_tolerance <= 0:
        msg = "phase reconstruction tolerances must be positive"
        raise ValueError(msg)

    # The operator uses ascending cos(theta); Baum stores increasing theta.
    phase = phase[..., ::-1]
    normalization = phase[:, 0] @ operator.kernels["00"][0]
    if np.any(~np.isfinite(normalization)) or np.any(normalization == 0):
        msg = "Every phase matrix must have a finite, nonzero A1[0]"
        raise ValueError(msg)
    target = phase / normalization[:, np.newaxis, np.newaxis]

    num_cells = len(target)
    max_moments = operator.max_moments
    coefficients = np.zeros((num_cells, 6, max_moments), dtype=np.float64)
    reconstructed = np.zeros_like(target)
    required = np.zeros(num_cells, dtype=np.int32)
    consecutive = np.zeros(num_cells, dtype=np.int32)
    active = np.ones(num_cells, dtype=bool)
    rms_error = np.full((num_cells, 6), np.inf, dtype=np.float64)
    max_error = np.full((num_cells, 6), np.inf, dtype=np.float64)

    for start in range(0, max_moments, block_size):
        stop = min(start + block_size, max_moments)
        active_indices = np.flatnonzero(active)
        if len(active_indices) == 0:
            break
        local_phase = target[active_indices]
        k00 = operator.kernels["00"][start:stop]
        k22 = operator.kernels["22"][start:stop]
        k2m2 = operator.kernels["2m2"][start:stop]
        k02 = operator.kernels["02"][start:stop]

        local_coeff = np.empty((len(active_indices), 6, stop - start), dtype=np.float64)
        local_coeff[:, 0] = local_phase[:, 0] @ k00.T
        local_coeff[:, 3] = local_phase[:, 5] @ k00.T
        local_coeff[:, 4] = local_phase[:, 1] @ k02.T
        local_coeff[:, 5] = -(local_phase[:, 4] @ k02.T)
        plus = (local_phase[:, 2] + local_phase[:, 3]) @ k22.T
        minus = (local_phase[:, 2] - local_phase[:, 3]) @ k2m2.T
        local_coeff[:, 1] = (plus + minus) / 2.0
        local_coeff[:, 2] = (plus - minus) / 2.0
        coefficients[active_indices, :, start:stop] = local_coeff

        local_reconstruction = reconstructed[active_indices]
        b00 = operator.native_basis["00"][start:stop]
        b22 = operator.native_basis["22"][start:stop]
        b2m2 = operator.native_basis["2m2"][start:stop]
        b02 = operator.native_basis["02"][start:stop]
        local_reconstruction[:, 0] += local_coeff[:, 0] @ b00
        local_reconstruction[:, 5] += local_coeff[:, 3] @ b00
        local_reconstruction[:, 1] += local_coeff[:, 4] @ b02
        local_reconstruction[:, 4] -= local_coeff[:, 5] @ b02
        plus = (local_coeff[:, 1] + local_coeff[:, 2]) @ b22
        minus = (local_coeff[:, 1] - local_coeff[:, 2]) @ b2m2
        local_reconstruction[:, 2] += (plus + minus) / 2.0
        local_reconstruction[:, 3] += (plus - minus) / 2.0
        reconstructed[active_indices] = local_reconstruction

        local_rms, local_max = _phase_errors(
            local_reconstruction,
            local_phase,
            operator.solid_angle_weights,
        )
        rms_error[active_indices] = local_rms
        max_error[active_indices] = local_max
        passed = np.all(
            (local_rms <= relative_rms_tolerance)
            & (local_max <= max_normalized_tolerance),
            axis=1,
        )
        consecutive[active_indices] = np.where(
            passed, consecutive[active_indices] + 1, 0
        )
        completed = active_indices[consecutive[active_indices] >= consecutive_passes]
        required[completed] = stop
        active[completed] = False

    converged = ~active.copy()
    if np.any(active) and not allow_unconverged_at_cap:
        active_indices = np.flatnonzero(active)
        rms_ratio = rms_error[active_indices] / relative_rms_tolerance
        max_ratio = max_error[active_indices] / max_normalized_tolerance
        cell_scores = np.max(np.maximum(rms_ratio, max_ratio), axis=1)
        worst_local = int(np.argmax(cell_scores))
        worst_cell = int(active_indices[worst_local])
        element_scores = np.maximum(rms_ratio[worst_local], max_ratio[worst_local])
        worst_element = int(np.argmax(element_scores))
        raise ConvergenceError(
            worst_cell,
            PHASE_MATRIX_ELEMENTS[worst_element],
            float(rms_error[worst_cell, worst_element]),
            float(max_error[worst_cell, worst_element]),
            max_moments,
        )
    required[active] = max_moments

    used_moments = int(required.max())
    return TransformResult(
        coefficients=coefficients[..., :used_moments],
        required_moments=required,
        converged=converged,
        rms_error=rms_error,
        max_normalized_error=max_error,
        a1_normalization=normalization,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_metadata(
    input_dir: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, str]]:
    paths = {model: input_dir / filename for model, filename in MODEL_FILES.items()}
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        msg = f"Missing Baum source files: {missing}"
        raise FileNotFoundError(msg)

    phase_angles = None
    diameters = None
    wavelengths_nm = None
    hashes = {}
    for model, path in paths.items():
        hashes[model] = _sha256(path)
        with xr.open_dataset(path) as source:
            local_angles = np.asarray(source["phase_angles"], dtype=np.float64)
            local_diameters = np.asarray(source["effective_diameter"], dtype=np.float64)
            local_wavelengths = (
                np.asarray(source["wavelengths"], dtype=np.float64) * 1000.0
            )
            for variable in (
                *_SOURCE_PHASE_VARIABLES,
                "asymmetry_parameter",
                "single_scattering_albedo",
                "extinction_coefficient_over_iwc",
            ):
                if variable not in source:
                    msg = f"{path} is missing required variable {variable}"
                    raise ValueError(msg)

        if phase_angles is None:
            phase_angles = local_angles
            diameters = local_diameters
            wavelengths_nm = local_wavelengths
        else:
            if not np.array_equal(phase_angles, local_angles):
                msg = f"Phase-angle grid differs in {path}"
                raise ValueError(msg)
            if not np.array_equal(diameters, local_diameters):
                msg = f"Effective-diameter grid differs in {path}"
                raise ValueError(msg)
            if not np.array_equal(wavelengths_nm, local_wavelengths):
                msg = f"Wavelength grid differs in {path}"
                raise ValueError(msg)

    return phase_angles, diameters, wavelengths_nm, hashes


def _initialize_output(
    output: Dataset,
    diameters: np.ndarray,
    wavelengths_nm: np.ndarray,
    source_hashes: dict[str, str],
    settings: dict,
) -> dict[str, object]:
    num_models = len(MODEL_FILES)
    num_diameters = len(diameters)
    num_wavelengths = len(wavelengths_nm)
    output.createDimension("particle_model", num_models)
    output.createDimension("effective_diameter_um", num_diameters)
    output.createDimension("wavelength_nm", num_wavelengths)
    output.createDimension("phase_matrix_element", len(PHASE_MATRIX_ELEMENTS))
    output.createDimension("legendre", None)

    model_var = output.createVariable("particle_model", str, ("particle_model",))
    model_var[:] = np.asarray(tuple(MODEL_FILES), dtype=object)
    diameter_var = output.createVariable(
        "effective_diameter_um", "f8", ("effective_diameter_um",)
    )
    diameter_var[:] = diameters
    diameter_var.units = "micron"
    wavelength_var = output.createVariable("wavelength_nm", "f8", ("wavelength_nm",))
    wavelength_var[:] = wavelengths_nm
    wavelength_var.units = "nm"
    element_var = output.createVariable(
        "phase_matrix_element", str, ("phase_matrix_element",)
    )
    element_var[:] = np.asarray(PHASE_MATRIX_ELEMENTS, dtype=object)
    legendre_var = output.createVariable("legendre", "i4", ("legendre",))

    compression = {"zlib": True, "complevel": 4, "shuffle": True}
    scalar_chunks = (1, 1, min(445, num_wavelengths))
    coefficient_chunks = (*scalar_chunks, 64)
    scalar_dimensions = (
        "particle_model",
        "effective_diameter_um",
        "wavelength_nm",
    )
    coefficient_dimensions = (*scalar_dimensions, "legendre")
    error_dimensions = (*scalar_dimensions, "phase_matrix_element")
    variables: dict[str, object] = {"legendre": legendre_var}

    for name in ("xs_total", "xs_scattering"):
        variables[name] = output.createVariable(
            name,
            "f4",
            scalar_dimensions,
            chunksizes=scalar_chunks,
            **compression,
        )
        variables[name].units = "m2 g-1"
    for name in COEFFICIENT_NAMES:
        variables[name] = output.createVariable(
            name,
            "f4",
            coefficient_dimensions,
            chunksizes=coefficient_chunks,
            **compression,
        )
        variables[name].normalization = "all families divided by lm_a1[0]"

    for name in (
        "source_asymmetry_parameter",
        "coefficient_asymmetry_parameter",
        "asymmetry_parameter_mismatch",
        "a1_normalization",
    ):
        variables[name] = output.createVariable(
            name,
            "f4",
            scalar_dimensions,
            chunksizes=scalar_chunks,
            **compression,
        )
    variables["required_moments"] = output.createVariable(
        "required_moments",
        "i4",
        scalar_dimensions,
        chunksizes=scalar_chunks,
        **compression,
    )
    variables[
        "required_moments"
    ].description = (
        "moments at convergence, or the safety cap when convergence was not reached"
    )
    variables["phase_reconstruction_converged"] = output.createVariable(
        "phase_reconstruction_converged",
        "i1",
        scalar_dimensions,
        chunksizes=scalar_chunks,
        **compression,
    )
    variables[
        "phase_reconstruction_converged"
    ].description = (
        "1 when both phase reconstruction criteria passed for two consecutive blocks"
    )
    variables["phase_relative_rms_error"] = output.createVariable(
        "phase_relative_rms_error",
        "f4",
        error_dimensions,
        chunksizes=(*scalar_chunks, len(PHASE_MATRIX_ELEMENTS)),
        **compression,
    )
    variables["phase_max_normalized_error"] = output.createVariable(
        "phase_max_normalized_error",
        "f4",
        error_dimensions,
        chunksizes=(*scalar_chunks, len(PHASE_MATRIX_ELEMENTS)),
        **compression,
    )

    output.title = "Baum V3.6 severely rough ice-crystal optical properties"
    output.baum_version = "3.6"
    output.phase_matrix_source_convention = (
        "P11 is absolute; P21, P22, P33, P43, and P44 are stored divided by P11"
    )
    output.sasktran2_phase_matrix_convention = (
        "P12=P21 and P34=-P43; lm_b2 uses -integral(d02*P34)"
    )
    output.coefficient_normalization = "all six Greek families divided by A1[0]"
    output.phase_interpolation = "piecewise linear in cos(scattering angle)"
    output.source_files = json.dumps(MODEL_FILES, sort_keys=True)
    output.source_sha256 = json.dumps(source_hashes, sort_keys=True)
    output.generation_settings = json.dumps(settings, sort_keys=True)
    return variables


def generate_database(
    input_dir: Path,
    output_path: Path,
    *,
    relative_rms_tolerance: float = 1e-4,
    max_normalized_tolerance: float = DEFAULT_MAX_NORMALIZED_TOLERANCE,
    moment_block_size: int = 128,
    safety_cap: int = 16384,
    quadrature_oversampling: float = 2.0,
    minimum_quadrature_order: int = 4,
    wavelength_batch_size: int = 16,
    allow_unconverged_at_cap: bool = True,
    overwrite: bool = False,
) -> Path:
    """Generate the capped, chunked Baum ice-crystal database."""
    input_dir = Path(input_dir)
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        msg = f"Output already exists: {output_path}; pass overwrite=True to replace it"
        raise FileExistsError(msg)
    if wavelength_batch_size <= 0:
        msg = "wavelength_batch_size must be positive"
        raise ValueError(msg)
    if relative_rms_tolerance <= 0 or max_normalized_tolerance <= 0:
        msg = "phase reconstruction tolerances must be positive"
        raise ValueError(msg)

    phase_angles, diameters, wavelengths_nm, source_hashes = _source_metadata(input_dir)
    operator = build_transform_operator(
        phase_angles,
        max_moments=safety_cap,
        block_size=moment_block_size,
        quadrature_oversampling=quadrature_oversampling,
        minimum_quadrature_order=minimum_quadrature_order,
    )
    settings = {
        "relative_rms_tolerance": relative_rms_tolerance,
        "max_normalized_tolerance": max_normalized_tolerance,
        "moment_block_size": moment_block_size,
        "safety_cap": safety_cap,
        "quadrature_oversampling": quadrature_oversampling,
        "minimum_quadrature_order": minimum_quadrature_order,
        "wavelength_batch_size": wavelength_batch_size,
        "consecutive_passes": 2,
        "allow_unconverged_at_cap": allow_unconverged_at_cap,
        "coefficient_dtype": "float32",
        "calculation_dtype": "float64",
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent,
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)

    global_moments = 0
    unconverged_cells = 0
    maximum_rms_error = -np.inf
    maximum_normalized_error = -np.inf
    worst_rms_cell = ""
    worst_normalized_cell = ""
    try:
        with Dataset(temporary_path, "w", format="NETCDF4") as output:
            variables = _initialize_output(
                output, diameters, wavelengths_nm, source_hashes, settings
            )
            for model_index, (model, filename) in enumerate(MODEL_FILES.items()):
                LOGGER.info("Processing %s", model)
                with xr.open_dataset(input_dir / filename) as source:
                    xs_total = np.asarray(
                        source["extinction_coefficient_over_iwc"], dtype=np.float32
                    )
                    xs_scattering = xs_total * np.asarray(
                        source["single_scattering_albedo"], dtype=np.float32
                    )
                    variables["xs_total"][model_index] = xs_total
                    variables["xs_scattering"][model_index] = xs_scattering
                    variables["source_asymmetry_parameter"][model_index] = np.asarray(
                        source["asymmetry_parameter"], dtype=np.float32
                    )

                    for diameter_index in range(len(diameters)):
                        for wavelength_start in range(
                            0, len(wavelengths_nm), wavelength_batch_size
                        ):
                            wavelength_stop = min(
                                wavelength_start + wavelength_batch_size,
                                len(wavelengths_nm),
                            )
                            selection = {
                                "nDeff": diameter_index,
                                "nWaveLen": slice(wavelength_start, wavelength_stop),
                            }
                            source_phase = [
                                source[name]
                                .isel(selection)
                                .transpose("nWaveLen", "nScatAng")
                                .to_numpy()
                                for name in _SOURCE_PHASE_VARIABLES
                            ]
                            absolute_phase = reconstruct_absolute_phase_matrix(
                                *source_phase
                            )
                            try:
                                result = adaptive_transform_batch(
                                    absolute_phase,
                                    operator,
                                    block_size=moment_block_size,
                                    relative_rms_tolerance=relative_rms_tolerance,
                                    max_normalized_tolerance=max_normalized_tolerance,
                                    allow_unconverged_at_cap=allow_unconverged_at_cap,
                                )
                            except ConvergenceError as error:
                                wavelength_index = wavelength_start + error.cell_index
                                msg = (
                                    f"{model}, effective diameter {diameters[diameter_index]} "
                                    f"micron, wavelength {wavelengths_nm[wavelength_index]} nm: "
                                    f"{error}"
                                )
                                raise RuntimeError(msg) from error

                            used_moments = result.coefficients.shape[-1]
                            if used_moments > global_moments:
                                # Materialize zero padding for every cell. The normal NetCDF
                                # fill value is intentionally not used as a physical zero.
                                for name in COEFFICIENT_NAMES:
                                    variables[name][
                                        :, :, :, global_moments:used_moments
                                    ] = np.float32(0.0)
                                global_moments = used_moments

                            wave_slice = slice(wavelength_start, wavelength_stop)
                            for coefficient_index, name in enumerate(COEFFICIENT_NAMES):
                                variables[name][
                                    model_index,
                                    diameter_index,
                                    wave_slice,
                                    :used_moments,
                                ] = result.coefficients[:, coefficient_index].astype(
                                    np.float32
                                )
                            variables["required_moments"][
                                model_index, diameter_index, wave_slice
                            ] = result.required_moments
                            variables["phase_reconstruction_converged"][
                                model_index, diameter_index, wave_slice
                            ] = result.converged.astype(np.int8)
                            variables["phase_relative_rms_error"][
                                model_index, diameter_index, wave_slice
                            ] = result.rms_error.astype(np.float32)
                            variables["phase_max_normalized_error"][
                                model_index, diameter_index, wave_slice
                            ] = result.max_normalized_error.astype(np.float32)
                            variables["a1_normalization"][
                                model_index, diameter_index, wave_slice
                            ] = result.a1_normalization.astype(np.float32)
                            coefficient_asymmetry = result.coefficients[:, 0, 1] / 3.0
                            variables["coefficient_asymmetry_parameter"][
                                model_index, diameter_index, wave_slice
                            ] = coefficient_asymmetry.astype(np.float32)
                            source_asymmetry = np.asarray(
                                source["asymmetry_parameter"].isel(selection),
                                dtype=np.float64,
                            )
                            variables["asymmetry_parameter_mismatch"][
                                model_index, diameter_index, wave_slice
                            ] = (coefficient_asymmetry - source_asymmetry).astype(
                                np.float32
                            )

                            unconverged_cells += int(
                                np.count_nonzero(~result.converged)
                            )
                            rms_index = np.unravel_index(
                                np.argmax(result.rms_error), result.rms_error.shape
                            )
                            local_rms = float(result.rms_error[rms_index])
                            if local_rms > maximum_rms_error:
                                maximum_rms_error = local_rms
                                worst_rms_cell = (
                                    f"{model}, diameter={diameters[diameter_index]} um, "
                                    f"wavelength={wavelengths_nm[wavelength_start + rms_index[0]]} nm, "
                                    f"element={PHASE_MATRIX_ELEMENTS[rms_index[1]]}"
                                )
                            max_index = np.unravel_index(
                                np.argmax(result.max_normalized_error),
                                result.max_normalized_error.shape,
                            )
                            local_max = float(result.max_normalized_error[max_index])
                            if local_max > maximum_normalized_error:
                                maximum_normalized_error = local_max
                                worst_normalized_cell = (
                                    f"{model}, diameter={diameters[diameter_index]} um, "
                                    f"wavelength={wavelengths_nm[wavelength_start + max_index[0]]} nm, "
                                    f"element={PHASE_MATRIX_ELEMENTS[max_index[1]]}"
                                )
                        LOGGER.info(
                            "Finished %s diameter %.1f micron",
                            model,
                            diameters[diameter_index],
                        )
                        output.sync()

            variables["legendre"][:] = np.arange(global_moments, dtype=np.int32)
            output.maximum_required_moments = global_moments
            output.database_variant = "full"
            output.stored_moments = global_moments
            output.full_available_moments = global_moments
            output.standard_database_key = FULL_STANDARD_DATABASE_KEY
            output.phase_reconstruction_all_converged = int(unconverged_cells == 0)
            output.num_unconverged_cells = unconverged_cells
            output.phase_reconstruction_policy = (
                "cells that do not meet both tolerances by the safety cap retain "
                "the capped coefficients and are marked unconverged"
            )
            output.maximum_recorded_relative_rms_error = maximum_rms_error
            output.maximum_recorded_max_normalized_error = maximum_normalized_error
            output.worst_relative_rms_error_cell = worst_rms_cell
            output.worst_max_normalized_error_cell = worst_normalized_cell
        temporary_path.replace(output_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    LOGGER.info("Wrote %s with %d moments", output_path, global_moments)
    if unconverged_cells:
        LOGGER.warning(
            "Stored %d cells at the %d-moment cap without convergence; "
            "see phase_reconstruction_converged and global error attributes",
            unconverged_cells,
            global_moments,
        )
    return output_path


def _variable_creation_options(variable, dimension_sizes: dict[str, int]) -> dict:
    options = {}
    if "_FillValue" in variable.ncattrs():
        options["fill_value"] = variable.getncattr("_FillValue")

    chunking = variable.chunking()
    if isinstance(chunking, list):
        options["chunksizes"] = tuple(
            min(chunk, dimension_sizes[dimension])
            for chunk, dimension in zip(chunking, variable.dimensions, strict=True)
        )

    filters = variable.filters()
    if filters and filters.get("zlib", False):
        options.update(
            zlib=True,
            complevel=filters["complevel"],
            shuffle=filters["shuffle"],
            fletcher32=filters["fletcher32"],
        )
    return options


def _copy_subset_variable(source, target, num_moments: int) -> None:
    if source.ndim == 0:
        target.assignValue(source.getValue())
        return
    if "legendre" not in source.dimensions:
        target[:] = source[:]
        return

    legendre_axis = source.dimensions.index("legendre")
    selection = [slice(None)] * source.ndim
    selection[legendre_axis] = slice(0, num_moments)
    if legendre_axis != source.ndim - 1 or source.ndim < 3:
        target[:] = source[tuple(selection)]
        return

    # Coefficient variables are ordered (model, diameter, wavelength, moment).
    # Copy one model/diameter slab at a time to keep the subset operation's peak
    # memory independent of the size of the complete database.
    for prefix in np.ndindex(source.shape[:-2]):
        source_selection = (*prefix, slice(None), slice(0, num_moments))
        target_selection = (*prefix, slice(None), slice(None))
        target[target_selection] = source[source_selection]


def _verify_runtime_subset(
    full_database_path: Path, subset_path: Path, num_moments: int
) -> None:
    with Dataset(full_database_path) as source, Dataset(subset_path) as subset:
        source.set_auto_maskandscale(False)
        subset.set_auto_maskandscale(False)
        for name in COEFFICIENT_NAMES:
            source_variable = source.variables[name]
            subset_variable = subset.variables[name]
            for prefix in np.ndindex(source_variable.shape[:-2]):
                source_values = np.asarray(
                    source_variable[(*prefix, slice(None), slice(0, num_moments))]
                )
                subset_values = np.asarray(
                    subset_variable[(*prefix, slice(None), slice(None))]
                )
                if source_values.tobytes() != subset_values.tobytes():
                    msg = (
                        f"Runtime subset verification failed for {name} at "
                        f"index {prefix}"
                    )
                    raise RuntimeError(msg)


def create_runtime_subset(
    full_database_path: Path,
    output_path: Path,
    *,
    max_moments: int = DEFAULT_DATABASE_MOMENTS,
    full_database_key: str = FULL_STANDARD_DATABASE_KEY,
    overwrite: bool = False,
    verify: bool = True,
) -> Path:
    """Create the lightweight runtime database from a capped Baum database.

    The coefficient prefix is copied without numerical conversion. All scalar
    optical data and full-expansion diagnostics are retained; attributes identify
    that the coefficient expansion itself has been truncated.
    """
    full_database_path = Path(full_database_path)
    output_path = Path(output_path)
    if not full_database_path.exists():
        msg = f"Full Baum database does not exist: {full_database_path}"
        raise FileNotFoundError(msg)
    if full_database_path.resolve() == output_path.resolve():
        msg = "Full and runtime-subset database paths must be different"
        raise ValueError(msg)
    if isinstance(max_moments, (bool, np.bool_)) or not isinstance(
        max_moments, (int, np.integer)
    ):
        msg = "max_moments must be an integer"
        raise TypeError(msg)
    max_moments = int(max_moments)
    if max_moments <= 0:
        msg = "max_moments must be positive"
        raise ValueError(msg)
    if output_path.exists() and not overwrite:
        msg = f"Output already exists: {output_path}; pass overwrite=True to replace it"
        raise FileExistsError(msg)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent,
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)

    try:
        with Dataset(full_database_path) as source:
            if "legendre" not in source.dimensions:
                msg = "Full Baum database is missing the legendre dimension"
                raise ValueError(msg)
            available_moments = len(source.dimensions["legendre"])
            if max_moments > available_moments:
                msg = (
                    f"Requested {max_moments} moments, but the full Baum database "
                    f"only contains {available_moments}"
                )
                raise ValueError(msg)
            missing = [
                name for name in COEFFICIENT_NAMES if name not in source.variables
            ]
            if missing:
                msg = f"Full Baum database is missing coefficient variables: {missing}"
                raise ValueError(msg)

            source.set_auto_maskandscale(False)
            dimension_sizes = {
                name: max_moments if name == "legendre" else len(dimension)
                for name, dimension in source.dimensions.items()
            }
            with Dataset(temporary_path, "w", format="NETCDF4") as target:
                for name, dimension in source.dimensions.items():
                    if name == "legendre":
                        target.createDimension(name, max_moments)
                    else:
                        size = None if dimension.isunlimited() else len(dimension)
                        target.createDimension(name, size)

                target.setncatts(
                    {name: source.getncattr(name) for name in source.ncattrs()}
                )
                target.database_variant = "default"
                target.stored_moments = max_moments
                target.full_available_moments = available_moments
                target.standard_database_key = DEFAULT_STANDARD_DATABASE_KEY
                target.full_database_key = full_database_key
                target.coefficient_subset_note = (
                    f"Greek coefficients are the first {max_moments} moments of "
                    "the full database; required_moments and phase-error variables "
                    "describe the full capped expansion"
                )

                for name, source_variable in source.variables.items():
                    options = _variable_creation_options(
                        source_variable, dimension_sizes
                    )
                    target_variable = target.createVariable(
                        name,
                        source_variable.datatype,
                        source_variable.dimensions,
                        **options,
                    )
                    target_variable.setncatts(
                        {
                            attribute: source_variable.getncattr(attribute)
                            for attribute in source_variable.ncattrs()
                            if attribute != "_FillValue"
                        }
                    )
                    _copy_subset_variable(source_variable, target_variable, max_moments)

        if verify:
            _verify_runtime_subset(full_database_path, temporary_path, max_moments)
        temporary_path.replace(output_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    LOGGER.info(
        "Wrote %s with the first %d moments from %s",
        output_path,
        max_moments,
        full_database_path,
    )
    return output_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("/Volumes/SasktranFiles/BaumIceCrystals"),
    )
    parser.add_argument(
        "--output",
        "--full-output",
        dest="output",
        type=Path,
        required=True,
        help="Full capped database output path",
    )
    parser.add_argument(
        "--default-output",
        type=Path,
        help="Optional lightweight runtime database derived from --output",
    )
    parser.add_argument("--default-moments", type=int, default=DEFAULT_DATABASE_MOMENTS)
    parser.add_argument(
        "--subset-only",
        action="store_true",
        help="Skip generation and derive --default-output from the existing --output",
    )
    parser.add_argument("--relative-rms-tolerance", type=float, default=1e-4)
    parser.add_argument(
        "--max-normalized-tolerance",
        type=float,
        default=DEFAULT_MAX_NORMALIZED_TOLERANCE,
    )
    parser.add_argument("--moment-block-size", type=int, default=128)
    parser.add_argument("--safety-cap", type=int, default=16384)
    parser.add_argument("--quadrature-oversampling", type=float, default=2.0)
    parser.add_argument("--minimum-quadrature-order", type=int, default=4)
    parser.add_argument("--wavelength-batch-size", type=int, default=16)
    parser.add_argument(
        "--fail-on-nonconvergence",
        action="store_true",
        help="Fail instead of storing capped coefficients and convergence diagnostics",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    parser = _parser()
    args = parser.parse_args()
    if args.subset_only:
        if args.default_output is None:
            parser.error("--subset-only requires --default-output")
        full_database_path = args.output
    else:
        full_database_path = generate_database(
            args.input_dir,
            args.output,
            relative_rms_tolerance=args.relative_rms_tolerance,
            max_normalized_tolerance=args.max_normalized_tolerance,
            moment_block_size=args.moment_block_size,
            safety_cap=args.safety_cap,
            quadrature_oversampling=args.quadrature_oversampling,
            minimum_quadrature_order=args.minimum_quadrature_order,
            wavelength_batch_size=args.wavelength_batch_size,
            allow_unconverged_at_cap=not args.fail_on_nonconvergence,
            overwrite=args.overwrite,
        )
    if args.default_output is not None:
        create_runtime_subset(
            full_database_path,
            args.default_output,
            max_moments=args.default_moments,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    main()
