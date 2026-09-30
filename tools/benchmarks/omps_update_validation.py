"""Additional validation injected into a saved fixed-state OMPS runner.

The probe harness copies this helper and the original runner before execution.
Only NumPy copies and labeled radiance arrays are retained between evaluations;
native linearizations remain owned by the production retrieval cache.
"""

from __future__ import annotations

import inspect
import json
import time
from pathlib import Path

import numpy as np


class OmpsUpdateValidation:
    def __init__(self, objective, initial: np.ndarray, destination: Path) -> None:
        self.objective = objective
        self.initial = np.array(initial, copy=True)
        self.destination = destination
        self.cache = inspect.getclosurevars(objective).nonlocals["cache"]
        if type(self.cache).__name__ not in {
            "_LinearizedMeasurementCache",
            "_MatrixFreeStateCache",
        }:
            msg = f"Unexpected production objective cache: {type(self.cache).__name__}"
            raise TypeError(msg)
        self.direction = np.sin(0.013 * np.arange(initial.size) + 0.2)
        self.direction /= np.max(np.abs(self.direction))
        self.amplitude = 1e-3
        self.evaluations = {}
        self.adjoints = {}
        self.radiances = {}
        self.gradients = {}
        np.save(
            destination / "validation_direction.npy", self.direction, allow_pickle=False
        )
        np.save(
            destination / "validation_state_scale.npy",
            self.cache._state_scale,
            allow_pickle=False,
        )

    def evaluate(self, label: str, state: np.ndarray) -> tuple[float, np.ndarray]:
        forward_model = self.cache._forward_model
        original = forward_model.calculate_linearized_radiance
        radiances = {}

        def capture_radiances():
            result = original()
            for key, radiance in result.items():
                # The plain Dataset copy carries no native callbacks or engines.
                radiances[str(key)] = radiance.data[["radiance"]].copy(deep=True)
            return result

        forward_model.calculate_linearized_radiance = capture_radiances
        before = time.monotonic()
        try:
            value, gradient = self.objective(np.array(state, copy=True))
        finally:
            forward_model.calculate_linearized_radiance = original
        elapsed = time.monotonic() - before
        gradient = np.array(gradient, copy=True)
        if not radiances:
            msg = "Validation objective unexpectedly reused a cached evaluation"
            raise RuntimeError(msg)
        if not np.isfinite(value) or not np.all(np.isfinite(gradient)):
            msg = "Validation objective or gradient contains nonfinite values"
            raise ValueError(msg)
        self.radiances[label] = radiances
        self.gradients[label] = gradient
        np.save(
            self.destination / f"validation_{label}_state.npy",
            state,
            allow_pickle=False,
        )
        np.save(
            self.destination / f"validation_{label}_gradient.npy",
            gradient,
            allow_pickle=False,
        )
        radiance_files = []
        for index, (key, dataset) in enumerate(radiances.items()):
            filename = f"validation_{label}_radiance_{index}.nc"
            dataset.to_netcdf(self.destination / filename)
            radiance_files.append({"key": key, "file": filename})
        self.evaluations[label] = {
            "objective": float(value),
            "gradient_inf_norm": float(np.max(np.abs(gradient))),
            "evaluation_s": elapsed,
            "radiance_files": radiance_files,
            "state_file": f"validation_{label}_state.npy",
            "gradient_file": f"validation_{label}_gradient.npy",
        }
        return float(value), gradient

    def adjoint(self, label: str, state: np.ndarray) -> None:
        # The good-measurement operator already includes operational wavelength
        # combinations and normalization. Its input is in retrieval coordinates;
        # multiply by the optimizer's explicit state scale for this test.
        operator = self.cache.evaluate(state)["operator"]
        tangent = self.cache._state_scale * self.direction
        cotangent = np.cos(0.017 * np.arange(operator.shape[0]) + 0.1)
        before = time.monotonic()
        jvp = np.array(operator.matvec(tangent), copy=True)
        vjp = np.array(operator.rmatvec(cotangent), copy=True)
        elapsed = time.monotonic() - before
        if not np.all(np.isfinite(jvp)) or not np.all(np.isfinite(vjp)):
            msg = "Validation JVP or VJP contains nonfinite values"
            raise ValueError(msg)
        left = float(jvp @ cotangent)
        right = float(tangent @ vjp)
        scale = max(abs(left), abs(right), np.finfo(float).tiny)
        self.adjoints[label] = {
            "jvp_dot_cotangent": left,
            "tangent_dot_vjp": right,
            "absolute_difference": abs(left - right),
            "relative_difference": abs(left - right) / scale,
            "allclose_rtol1e_8_atol1e_10": bool(
                np.isclose(left, right, rtol=1e-8, atol=1e-10)
            ),
            "runtime_s": elapsed,
            "jvp_shape": list(jvp.shape),
            "vjp_shape": list(vjp.shape),
        }
        np.save(
            self.destination / f"validation_{label}_jvp.npy", jvp, allow_pickle=False
        )
        np.save(
            self.destination / f"validation_{label}_vjp.npy", vjp, allow_pickle=False
        )
        self.write_report()

    @staticmethod
    def differences(reference: np.ndarray, trial: np.ndarray) -> dict:
        delta = np.asarray(trial) - np.asarray(reference)
        norm = float(np.linalg.norm(reference))
        return {
            "count": int(delta.size),
            "bitwise_identical": bool(reference.tobytes() == trial.tobytes()),
            "max_absolute_difference": float(np.max(np.abs(delta))),
            "rms_difference": float(np.sqrt(np.mean(delta**2))),
            "relative_l2_difference": (
                float(np.linalg.norm(delta) / norm) if norm else None
            ),
            "allclose_rtol1e_10_atol1e_12": bool(
                np.allclose(trial, reference, rtol=1e-10, atol=1e-12)
            ),
        }

    def finish(self) -> None:
        # This runs after retrieval.nc and completed.json have been written.
        # The first perturbation forces the production cache to release its old
        # operator and update the atmosphere before constructing a new one.
        perturbed = self.initial + self.amplitude * self.direction
        self.evaluate("perturbed", perturbed)
        self.evaluate("restored", self.initial)
        self.adjoint("restored", self.initial)
        self.write_report()

    def write_report(self) -> None:
        report = {
            "cache_class": type(self.cache).__name__,
            "scope": (
                "Validation adds two objective/gradient evaluations plus initial and restored "
                "operational JVP/VJP products. Original nit=0/nfev=1 fixed-state outputs "
                "are saved before atmosphere-update evaluations and are not overwritten."
            ),
            "additional_objective_evaluations": max(len(self.evaluations) - 1, 0),
            "optimizer_coordinate_perturbation_amplitude": self.amplitude,
            "evaluations": self.evaluations,
            "operational_adjoint_checks": self.adjoints,
        }
        if "restored" in self.evaluations:
            report["restored_minus_initial_objective"] = (
                self.evaluations["restored"]["objective"]
                - self.evaluations["initial"]["objective"]
            )
            report["restored_gradient_comparison"] = self.differences(
                self.gradients["initial"], self.gradients["restored"]
            )
            report["restored_radiance_comparisons"] = {
                key: self.differences(
                    initial.radiance.values,
                    self.radiances["restored"][key].radiance.values,
                )
                for key, initial in self.radiances["initial"].items()
            }
            report["objective_directional_check"] = {
                "initial_gradient_dot_direction": float(
                    self.gradients["initial"] @ self.direction
                ),
                "one_sided_cost_difference_over_step": (
                    self.evaluations["perturbed"]["objective"]
                    - self.evaluations["initial"]["objective"]
                )
                / self.amplitude,
                "note": "One-sided difference includes curvature at a 1e-3 step; reported without an equality gate.",
            }
        (self.destination / "update_validation.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )
