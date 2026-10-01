from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

COMPARE = runpy.run_path(
    str(Path(__file__).resolve().parents[2] / "tools/benchmarks/omps_memory_compare.py")
)


def write_json(path, value):
    path.write_text(json.dumps(value))


def write_probe(folder, *, updates):
    folder.mkdir()
    state = np.arange(8, dtype=float)
    for name in ("state", "gradient"):
        np.save(folder / f"fixed_initial_{name}.npy", state)
    write_json(
        folder / "completed.json",
        {
            "experiment_fixed_ozone": 1,
            "ozone_minimizer_iterations": 0,
            "ozone_minimizer_function_evaluations": 1,
        },
    )
    write_json(folder / "parameters.json", {"kwargs": {"incoming": 110}})
    write_json(
        folder / "fixed_initial_evaluation.json",
        {"objective": 0.4, "evaluation_s": 1.0},
    )
    xr.Dataset(
        {
            "ozone_simulated_reflectance": ("ray", state),
            "ozone_vmr": ("altitude", state),
        },
        coords={"ray": state, "altitude": state},
    ).to_netcdf(folder / "retrieval.nc")
    xr.Dataset({"fixed_albedo": ("ray", state)}).to_netcdf(
        folder / "prescribed_scene.nc"
    )
    write_json(
        folder / "guard.summary.json",
        {
            "child_exit_code": 0,
            "stopped_for_memory": None,
            "monitor_error_or_interruption": None,
            "elapsed_s": 2.0,
            "peak_physical_footprint_gib": 1.0,
            "peak_rss_gib": 0.5,
            "sample_interval_s": 0.5,
        },
    )
    if not updates:
        return
    report = {"evaluations": {}, "internal_adjoint_passed": False}
    for label in ("initial", "perturbed", "restored"):
        files = {}
        for name in ("state", "gradient"):
            filename = f"validation_{label}_{name}.npy"
            np.save(folder / filename, state)
            files[f"{name}_file"] = filename
        radiance_file = f"validation_{label}_radiance_0.nc"
        xr.Dataset({"radiance": ("ray", state)}, coords={"ray": state}).to_netcdf(
            folder / radiance_file
        )
        report["evaluations"][label] = {
            **files,
            "radiance_files": [{"key": "0", "file": radiance_file}],
            "objective": 0.4,
            "evaluation_s": 1.0,
        }
    for label in ("initial", "restored"):
        for name in ("jvp", "vjp"):
            np.save(folder / f"validation_{label}_{name}.npy", state)
    write_json(folder / "update_validation.json", report)


def compare_probes(tmp_path, monkeypatch, *, require_updates=False):
    output = tmp_path / "comparison.json"
    arguments = [
        "omps_memory_compare.py",
        "--reference",
        str(tmp_path / "reference"),
        "--trial",
        str(tmp_path / "trial"),
        "--reference-summary",
        str(tmp_path / "reference/guard.summary.json"),
        "--trial-summary",
        str(tmp_path / "trial/guard.summary.json"),
        "--output",
        str(output),
    ]
    if require_updates:
        arguments.append("--require-updates")
    monkeypatch.setattr(sys, "argv", arguments)
    result = COMPARE["main"]()
    return result, json.loads(output.read_text())


@pytest.mark.parametrize("require_updates", [False, True])
def test_missing_update_artifacts_are_explicit(tmp_path, monkeypatch, require_updates):
    write_probe(tmp_path / "reference", updates=False)
    write_probe(tmp_path / "trial", updates=True)
    result, report = compare_probes(
        tmp_path, monkeypatch, require_updates=require_updates
    )
    assert result == int(require_updates)
    assert report["equivalent_at_saved_state"] is not require_updates
    assert report["compared_array_count"] == 4
    assert report["complete_update_validation_available"] is False
    assert report["equivalent_across_atmosphere_updates_and_products"] is False
    assert "comparison unavailable" in report["scope"]


def test_complete_coverage_preserves_internal_check_failures(tmp_path, monkeypatch):
    write_probe(tmp_path / "reference", updates=True)
    write_probe(tmp_path / "trial", updates=True)
    result, report = compare_probes(tmp_path, monkeypatch, require_updates=True)
    assert result == 0
    assert report["compared_array_count"] == 17
    assert report["all_compared_arrays_and_objectives_close"] is True
    assert report["equivalent_across_atmosphere_updates_and_products"] is True
    checks = report["atmosphere_update_comparisons"]
    assert checks["reference_internal_checks"]["internal_adjoint_passed"] is False
    assert checks["trial_internal_checks"]["internal_adjoint_passed"] is False


@pytest.mark.parametrize(
    "filename",
    [
        "fixed_initial_gradient.npy",
        "validation_perturbed_gradient.npy",
        "validation_restored_vjp.npy",
    ],
)
def test_every_array_participates_in_acceptance(tmp_path, monkeypatch, filename):
    write_probe(tmp_path / "reference", updates=True)
    write_probe(tmp_path / "trial", updates=True)
    changed = np.load(tmp_path / "trial" / filename)
    changed[-1] += 0.01
    np.save(tmp_path / "trial" / filename, changed)
    result, report = compare_probes(tmp_path, monkeypatch, require_updates=True)
    assert result == 1
    assert report["compared_array_count"] == 17
    assert report["all_compared_arrays_and_objectives_close"] is False
    assert report["equivalent_across_atmosphere_updates_and_products"] is False


def test_nonfinite_objective_is_rejected(tmp_path, monkeypatch):
    for name in ("reference", "trial"):
        write_probe(tmp_path / name, updates=False)
        write_json(
            tmp_path / name / "fixed_initial_evaluation.json",
            {"objective": float("inf"), "evaluation_s": 1.0},
        )
    with pytest.raises(ValueError, match="objective is nonfinite"):
        compare_probes(tmp_path, monkeypatch)


def test_nonfinite_guard_measurement_is_rejected(tmp_path):
    write_probe(tmp_path / "reference", updates=False)
    path = tmp_path / "reference/guard.summary.json"
    record = json.loads(path.read_text())
    record["peak_physical_footprint_gib"] = float("nan")
    write_json(path, record)
    with pytest.raises(ValueError, match="Invalid memory guard measurement"):
        COMPARE["memory_record"](path)


def test_guard_failure_rejects_complete_matching_arrays(tmp_path, monkeypatch):
    write_probe(tmp_path / "reference", updates=True)
    write_probe(tmp_path / "trial", updates=True)
    path = tmp_path / "trial/guard.summary.json"
    record = json.loads(path.read_text())
    record["stopped_for_memory"] = "physical footprint ceiling"
    write_json(path, record)
    result, report = compare_probes(tmp_path, monkeypatch, require_updates=True)
    assert result == 1
    assert report["all_compared_arrays_and_objectives_close"] is True
    assert report["memory"]["trial"]["successful"] is False
    assert report["equivalent_across_atmosphere_updates_and_products"] is False


def test_resolution_changes_reject_matching_arrays(tmp_path, monkeypatch):
    write_probe(tmp_path / "reference", updates=True)
    write_probe(tmp_path / "trial", updates=True)
    write_json(tmp_path / "trial/parameters.json", {"kwargs": {"incoming": 86}})
    result, report = compare_probes(tmp_path, monkeypatch, require_updates=True)
    assert result == 1
    assert report["all_compared_arrays_and_objectives_close"] is True
    assert report["same_numerical_settings"] is False
    assert report["equivalent_across_atmosphere_updates_and_products"] is False


@pytest.mark.parametrize("values", [np.array([]), np.array([float("nan")])])
def test_invalid_comparison_arrays_are_rejected(values):
    with pytest.raises(ValueError, match=r"empty|nonfinite"):
        COMPARE["array_comparison"](values, values, 1e-10, 1e-12)
