"""Compare complete saved OMPS reflectances, states, and ozone gradients.

Run with the pinned tomography Python. This script reads saved artifacts only;
it never instantiates an engine. Different source-column counts are reported
as different numerical resolutions, even when their arrays are close.

--initial-only compares the complete saved initial objective, state, gradient,
and captured radiances, including artifacts saved before additional validation
was interrupted. It does not certify successful completion of the whole probe.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def array_comparison(
    reference: np.ndarray, trial: np.ndarray, rtol: float, atol: float
) -> dict:
    reference = np.asarray(reference)
    trial = np.asarray(trial)
    if reference.shape != trial.shape:
        msg = f"Different array shapes: {reference.shape} versus {trial.shape}"
        raise ValueError(msg)
    if not np.all(np.isfinite(reference)) or not np.all(np.isfinite(trial)):
        msg = "A comparison array contains nonfinite values"
        raise ValueError(msg)
    delta = trial - reference
    reference_norm = float(np.linalg.norm(reference.ravel()))
    relative = np.divide(
        delta, np.abs(reference), out=np.zeros_like(delta), where=reference != 0
    )
    return {
        "shape": list(reference.shape),
        "count": int(reference.size),
        "reference_dtype": str(reference.dtype),
        "trial_dtype": str(trial.dtype),
        "bitwise_identical": bool(
            reference.dtype == trial.dtype
            and reference.tobytes(order="C") == trial.tobytes(order="C")
        ),
        "array_equal": bool(np.array_equal(reference, trial)),
        "allclose": bool(np.allclose(trial, reference, rtol=rtol, atol=atol)),
        "max_absolute_difference": float(np.max(np.abs(delta))),
        "rms_difference": float(np.sqrt(np.mean(delta**2))),
        "relative_l2_difference": (
            float(np.linalg.norm(delta.ravel()) / reference_norm)
            if reference_norm
            else None
        ),
        "max_relative_difference_where_reference_nonzero": float(
            np.max(np.abs(relative))
        ),
        "max_absolute_difference_where_reference_zero": float(
            np.max(np.abs(delta[reference == 0]), initial=0)
        ),
    }


def memory_record(path: Path) -> dict:
    record = read_json(path)
    return {
        "summary": str(path),
        "child_exit_code": record["child_exit_code"],
        "successful": bool(
            record["child_exit_code"] == 0
            and record["stopped_for_memory"] is None
            and record["monitor_error_or_interruption"] is None
        ),
        "elapsed_s": record["elapsed_s"],
        "peak_physical_footprint_decimal_gb": record["peak_physical_footprint_gib"]
        * 1024**3
        / 1e9,
        "peak_rss_decimal_gb": record["peak_rss_gib"] * 1024**3 / 1e9,
        "sample_interval_s": record["sample_interval_s"],
        "stopped_for_memory": record["stopped_for_memory"],
        "monitor_error_or_interruption": record["monitor_error_or_interruption"],
    }


def completed(path: Path) -> dict:
    record = read_json(path / "completed.json")
    if (
        int(record["experiment_fixed_ozone"]) != 1
        or int(record["ozone_minimizer_iterations"]) != 0
        or int(record["ozone_minimizer_function_evaluations"]) != 1
    ):
        msg = (
            f"Expected one completed fixed-state objective/gradient evaluation: {path}"
        )
        raise ValueError(msg)
    return record


def initial_radiance_comparison(
    reference: Path, trial: Path, rtol: float, atol: float
) -> dict:
    filename = "validation_initial_radiance_0.nc"
    for folder in (reference, trial):
        if len(list(folder.glob("validation_initial_radiance_*.nc"))) != 1:
            msg = f"Expected exactly one complete initial radiance dataset: {folder}"
            raise ValueError(msg)
    with (
        xr.open_dataset(reference / filename) as ref_ds,
        xr.open_dataset(trial / filename) as trial_ds,
    ):
        ref_radiance = ref_ds.radiance
        trial_radiance = trial_ds.radiance.transpose(*ref_radiance.dims)
        if ref_radiance.size != 111936 or trial_radiance.size != 111936:
            msg = "Initial-only comparison requires all 111,936 modeled radiances"
            raise ValueError(msg)
        xr.testing.assert_equal(
            ref_ds.coords.to_dataset(), trial_ds.coords.to_dataset()
        )
        coordinates = {}
        for name, ref_coordinate in ref_ds.coords.items():
            test_coordinate = trial_ds.coords[name].transpose(*ref_coordinate.dims)
            if ref_coordinate.dtype != test_coordinate.dtype:
                msg = f"Initial radiance coordinate dtype changed: {name}"
                raise ValueError(msg)
            coordinates[name] = {
                "dimensions": list(ref_coordinate.dims),
                "shape": list(ref_coordinate.shape),
                "dtype": str(ref_coordinate.dtype),
                "values_identical": True,
            }
        comparison = array_comparison(
            ref_radiance.values, trial_radiance.values, rtol, atol
        )
        comparison["dimensions"] = list(ref_radiance.dims)
        comparison["coordinates"] = coordinates
        comparison["artifact"] = filename
        comparison["variable"] = "radiance"
    return comparison


def initial_input_differences(reference: dict, trial: dict) -> dict:
    for parameters in (reference, trial):
        if parameters.get("fixed_ozone_evaluation") is not True:
            msg = "Initial-only comparison requires a fixed-ozone objective evaluation"
            raise ValueError(msg)
    names = (
        "case",
        "baseline",
        "l1g",
        "anc",
        "custom_scene",
        "scene_provenance",
        "grouping_origin",
        "omit_fixed_scene_derivatives",
    )
    return {
        name: [reference.get(name), trial.get(name)]
        for name in names
        if reference.get(name) != trial.get(name)
    }


def numerical_kwargs(parameters: dict) -> dict:
    result = dict(parameters["kwargs"])
    native = dict(result.pop("model_kwargs", None) or {})
    # A bounded cache changes retention and runtime, not physical sampling.
    native.pop("successive_orders_transport_cache_wavelengths", None)
    for field, uniform_name in (
        (
            "successive_orders_incoming_directions_by_altitude",
            "successive_orders_incoming",
        ),
        (
            "successive_orders_outgoing_directions_by_altitude",
            "successive_orders_outgoing",
        ),
    ):
        profile = native.get(field)
        if (
            profile is None
            or not profile
            or all(count == result.get(uniform_name) for count in profile)
        ):
            native.pop(field, None)
    if native:
        result["model_kwargs"] = native
    return result


def operational_signal_accuracy(reference: xr.Dataset, trial: xr.Dataset) -> dict:
    """Apply the production UV/visible transform to saved reflectances only."""
    # Ordinary comparisons require only NumPy/xarray; the optional study uses
    # the pinned processor's measurement transform without constructing an engine.
    from omps_tomography.ozone import ozone_measurement_vectors  # noqa: PLC0415
    from skretrieval.core.radianceformat import RadianceGridded  # noqa: PLC0415
    from skretrieval.retrieval.measvec import (  # noqa: PLC0415
        _grouped_altitude_triplet_plan,
    )

    image = reference.ozone_image.values
    height = reference.ozone_tangent_altitude.values
    observed = reference.ozone_measured_reflectance.transpose(
        "ozone_wavelength", "ozone_los"
    ).values
    noise = reference.ozone_reflectance_noise.transpose(
        "ozone_wavelength", "ozone_los"
    ).values
    radiance = RadianceGridded(
        xr.Dataset(
            {"radiance": (("wavelength", "los"), np.ones_like(observed))},
            coords={
                "wavelength": reference.ozone_wavelength.values,
                "los": np.arange(image.size),
                "image": ("los", image),
                "tangent_altitude": ("los", height),
            },
        )
    )
    mode = reference.attrs.get("ozone_measurement_mode", "v2_1_combined")
    vector = ozone_measurement_vectors(mode)[mode]
    plan = _grouped_altitude_triplet_plan(
        radiance,
        wavelength=vector._wavelength,
        weights=vector._weights,
        normalization_range=vector._normalization_range,
        altitude_weight_grid=vector._altitude_weight_grid,
        altitude_weight_values=vector._altitude_weight_values,
        altitude_range=vector._altitude_range,
        group_by="image",
        open_altitude_bounds=True,
    )
    _, first = np.unique(image, return_index=True)
    images = image[np.sort(first)]
    rows = np.concatenate(
        [
            np.flatnonzero((image == value) & (height > 5000) & (height < 59000))
            for value in images
        ]
    )
    if rows.size != plan.transform.shape[0]:
        msg = "Production operational measurement row order changed"
        raise ValueError(msg)
    with np.errstate(divide="ignore", invalid="ignore"):
        invalid = ~np.isfinite(np.log(observed.ravel()))
        relative_variance = (noise / observed).ravel() ** 2
        modeled = trial.ozone_simulated_reflectance.transpose(
            "ozone_wavelength", "ozone_los"
        ).values
        baseline = reference.ozone_simulated_reflectance.transpose(
            "ozone_wavelength", "ozone_los"
        ).values
        delta = np.asarray(plan.transform @ np.log(modeled / baseline).ravel()).ravel()
    valid = np.asarray(plan.validity_inputs @ invalid.astype(np.int8)).ravel() == 0
    variance = np.asarray(plan.legacy_variance_weights @ relative_variance).ravel()
    variance[variance <= 0] = 1.0
    if np.any(valid & ~np.isfinite(delta)):
        msg = "Nonfinite modeled operational signal on a valid measurement row"
        raise ValueError(msg)
    times = reference.ozone_time.values.astype("datetime64[ns]").astype("i8")
    latitude = np.interp(
        times,
        reference.ozone_reference_time.values.astype("datetime64[ns]").astype("i8"),
        reference.ozone_geodetic_latitude_deg.values,
    )
    scopes = {"all_operational_rows": valid}
    for lower, upper in ((0, 20), (5, 10)):
        scopes[f"lat{lower}_{upper}_height20.5_31.5km"] = (
            valid
            & (latitude[rows] >= lower)
            & (latitude[rows] < upper)
            & (height[rows] >= 20500)
            & (height[rows] <= 31500)
        )
    normalized = delta / np.sqrt(variance)
    report = {}
    for name, selection in scopes.items():
        values = normalized[selection]
        report[name] = {
            "count": int(values.size),
            "mean_assumed_error_units": float(np.mean(values)) if values.size else None,
            "rms_assumed_error_units": (
                float(np.sqrt(np.mean(values**2))) if values.size else None
            ),
            "max_absolute_assumed_error_units": (
                float(np.max(np.abs(values))) if values.size else None
            ),
        }
    return {
        "measurement_mode": mode,
        "scope": (
            "Saved reflectances only; production UV/visible combinations, image normalization "
            "and assumed-error propagation. Coarsening is an accuracy study and is not "
            "certified as numerical equivalence. Full gradient differences are reported separately."
        ),
        "signal_delta": report,
    }


def update_comparisons(reference: Path, trial: Path, rtol: float, atol: float) -> dict:
    paths = [folder / "update_validation.json" for folder in (reference, trial)]
    if not all(path.is_file() for path in paths):
        return {
            "available": False,
            "reference_has_update_validation": paths[0].is_file(),
            "trial_has_update_validation": paths[1].is_file(),
        }
    ref_report, trial_report = [read_json(path) for path in paths]
    report = {"available": True, "evaluations": {}, "operational_products": {}}
    passed = True
    bitwise = True
    for label in ("initial", "perturbed", "restored"):
        ref = ref_report["evaluations"][label]
        test = trial_report["evaluations"][label]
        result = {}
        for name in ("state", "gradient"):
            result[name] = array_comparison(
                np.load(reference / ref[f"{name}_file"], allow_pickle=False),
                np.load(trial / test[f"{name}_file"], allow_pickle=False),
                rtol,
                atol,
            )
            passed &= (
                result[name]["bitwise_identical"]
                if name == "state"
                else result[name]["allclose"]
            )
            bitwise &= result[name]["bitwise_identical"]
        ref_files = {item["key"]: item["file"] for item in ref["radiance_files"]}
        trial_files = {item["key"]: item["file"] for item in test["radiance_files"]}
        if set(ref_files) != set(trial_files):
            msg = f"Validation radiance keys changed in {label}"
            raise ValueError(msg)
        result["radiances"] = {}
        for key, filename in ref_files.items():
            with (
                xr.open_dataset(reference / filename) as ref_ds,
                xr.open_dataset(trial / trial_files[key]) as trial_ds,
            ):
                xr.testing.assert_equal(
                    ref_ds.coords.to_dataset(), trial_ds.coords.to_dataset()
                )
                comparison = array_comparison(
                    ref_ds.radiance.values,
                    trial_ds.radiance.transpose(*ref_ds.radiance.dims).values,
                    rtol,
                    atol,
                )
            result["radiances"][key] = comparison
            passed &= comparison["allclose"]
            bitwise &= comparison["bitwise_identical"]
        result["objective"] = {
            "reference": ref["objective"],
            "trial": test["objective"],
            "difference": test["objective"] - ref["objective"],
            "identical": test["objective"] == ref["objective"],
            "allclose": bool(
                np.isclose(test["objective"], ref["objective"], rtol=rtol, atol=atol)
            ),
        }
        result["runtime_s"] = [ref["evaluation_s"], test["evaluation_s"]]
        passed &= result["objective"]["allclose"]
        bitwise &= result["objective"]["identical"]
        report["evaluations"][label] = result
    for label in ("initial", "restored"):
        report["operational_products"][label] = {}
        for name in ("jvp", "vjp"):
            comparison = array_comparison(
                np.load(
                    reference / f"validation_{label}_{name}.npy", allow_pickle=False
                ),
                np.load(trial / f"validation_{label}_{name}.npy", allow_pickle=False),
                rtol,
                atol,
            )
            report["operational_products"][label][name] = comparison
            passed &= comparison["allclose"]
            bitwise &= comparison["bitwise_identical"]
    report["all_complete_arrays_and_objectives_close"] = bool(passed)
    report["all_complete_arrays_and_objectives_bitwise_identical"] = bool(bitwise)
    report["reference_internal_checks"] = ref_report
    report["trial_internal_checks"] = trial_report
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--trial", type=Path, required=True)
    parser.add_argument("--reference-summary", type=Path, required=True)
    parser.add_argument("--trial-summary", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument(
        "--accuracy-study",
        action="store_true",
        help="Report operational signal RMS in assumed-error units for a variable-grid study",
    )
    parser.add_argument(
        "--initial-only",
        action="store_true",
        help=(
            "Compare saved initial state/gradient/objective and complete captured "
            "radiances without requiring retrieval.nc or a completed whole probe"
        ),
    )
    args = parser.parse_args()
    if args.accuracy_study and args.initial_only:
        parser.error("--accuracy-study requires completed retrieval.nc datasets")
    if args.rtol < 0 or args.atol < 0 or not np.isfinite([args.rtol, args.atol]).all():
        parser.error("Comparison tolerances must be finite and nonnegative")
    reference = args.reference.resolve()
    trial = args.trial.resolve()
    if not args.initial_only:
        completed(reference)
        completed(trial)
    ref_evaluation = read_json(reference / "fixed_initial_evaluation.json")
    trial_evaluation = read_json(trial / "fixed_initial_evaluation.json")
    ref_parameters = read_json(reference / "parameters.json")
    trial_parameters = read_json(trial / "parameters.json")
    ref_numerical = numerical_kwargs(ref_parameters)
    trial_numerical = numerical_kwargs(trial_parameters)
    settings = {
        name: [ref_numerical.get(name), trial_numerical.get(name)]
        for name in sorted(set(ref_numerical) | set(trial_numerical))
        if ref_numerical.get(name) != trial_numerical.get(name)
    }
    comparisons = {}
    signal_accuracy = None
    for name in ("state", "gradient"):
        comparisons[name] = array_comparison(
            np.load(reference / f"fixed_initial_{name}.npy", allow_pickle=False),
            np.load(trial / f"fixed_initial_{name}.npy", allow_pickle=False),
            args.rtol,
            args.atol,
        )
    input_differences = {}
    if args.initial_only:
        input_differences = initial_input_differences(ref_parameters, trial_parameters)
        for name in ("state", "gradient"):
            if comparisons[name]["count"] != 15390:
                msg = f"Initial-only comparison requires all 15,390 {name} entries"
                raise ValueError(msg)
        comparisons["modeled_radiance"] = initial_radiance_comparison(
            reference, trial, args.rtol, args.atol
        )
        scene_identical = None
    else:
        with (
            xr.open_dataset(reference / "retrieval.nc") as ref_ds,
            xr.open_dataset(trial / "retrieval.nc") as trial_ds,
        ):
            name = "ozone_simulated_reflectance"
            ref_radiance = ref_ds[name]
            trial_radiance = trial_ds[name].transpose(*ref_radiance.dims)
            for coordinate in ref_radiance.coords:
                xr.testing.assert_equal(
                    ref_radiance[coordinate], trial_radiance[coordinate]
                )
            comparisons["modeled_reflectance"] = array_comparison(
                ref_radiance.values, trial_radiance.values, args.rtol, args.atol
            )
            comparisons["ozone_vmr"] = array_comparison(
                ref_ds.ozone_vmr.values, trial_ds.ozone_vmr.values, args.rtol, args.atol
            )
            if args.accuracy_study:
                signal_accuracy = operational_signal_accuracy(ref_ds, trial_ds)
        with (
            xr.open_dataset(reference / "prescribed_scene.nc") as ref_scene,
            xr.open_dataset(trial / "prescribed_scene.nc") as trial_scene,
        ):
            xr.testing.assert_identical(ref_scene, trial_scene)
        scene_identical = True
    objective_delta = trial_evaluation["objective"] - ref_evaluation["objective"]
    objective_close = bool(
        np.isclose(
            trial_evaluation["objective"],
            ref_evaluation["objective"],
            rtol=args.rtol,
            atol=args.atol,
        )
    )
    ref_memory = memory_record(args.reference_summary.resolve())
    trial_summary = args.trial_summary or trial.parent / "guard.summary.json"
    trial_memory = memory_record(trial_summary.resolve())
    updates = (
        {
            "available": False,
            "included_in_comparison": False,
            "reason": "Initial-only scope excludes additional update and adjoint checks",
        }
        if args.initial_only
        else update_comparisons(reference, trial, args.rtol, args.atol)
    )
    radiance_name = "modeled_radiance" if args.initial_only else "modeled_reflectance"
    equivalent = bool(
        not settings
        and not input_differences
        and comparisons["state"]["bitwise_identical"]
        and (args.initial_only or comparisons["ozone_vmr"]["bitwise_identical"])
        and comparisons["gradient"]["allclose"]
        and comparisons[radiance_name]["allclose"]
        and objective_close
        and (args.initial_only or ref_memory["successful"])
        and (args.initial_only or trial_memory["successful"])
        and (
            not updates["available"]
            or updates["all_complete_arrays_and_objectives_close"]
        )
    )
    bitwise_equivalent = bool(
        equivalent
        and comparisons["gradient"]["bitwise_identical"]
        and comparisons[radiance_name]["bitwise_identical"]
        and trial_evaluation["objective"] == ref_evaluation["objective"]
        and (
            not updates["available"]
            or updates["all_complete_arrays_and_objectives_bitwise_identical"]
        )
    )
    report = {
        "reference": str(reference),
        "trial": str(trial),
        "rtol": args.rtol,
        "atol": args.atol,
        "settings_differences_reference_then_trial": settings,
        "same_numerical_settings": not settings,
        "requested_native_config_overrides_reference_then_trial": [
            ref_parameters.get("candidate_native_config", {}).get(
                "requested_overrides", {}
            ),
            trial_parameters.get("candidate_native_config", {}).get(
                "requested_overrides", {}
            ),
        ],
        "operational_signal_accuracy_study": signal_accuracy,
        "prescribed_scene_identical": scene_identical,
        "complete_array_comparisons": comparisons,
        "atmosphere_update_comparisons": updates,
        "objective": {
            "reference": ref_evaluation["objective"],
            "trial": trial_evaluation["objective"],
            "difference": objective_delta,
            "identical": trial_evaluation["objective"] == ref_evaluation["objective"],
            "allclose": objective_close,
        },
        "evaluation_runtime_s": {
            "reference": ref_evaluation["evaluation_s"],
            "trial": trial_evaluation["evaluation_s"],
        },
        "memory": {"reference": ref_memory, "trial": trial_memory},
        "physical_footprint_saving_decimal_gb": ref_memory[
            "peak_physical_footprint_decimal_gb"
        ]
        - trial_memory["peak_physical_footprint_decimal_gb"],
        "trial_below_30_decimal_gb": trial_memory["peak_physical_footprint_decimal_gb"]
        < 30,
    }
    if args.initial_only:
        report.update(
            {
                "initial_input_metadata_differences_reference_then_trial": input_differences,
                "equivalent_initial_objective_gradient_and_radiances": equivalent,
                "bitwise_identical_initial_objective_gradient_and_radiances": bitwise_equivalent,
                "full_probe_validated": False,
                "scope": (
                    "Initial-only: complete saved initial objective, optimizer state, "
                    "ozone gradient and 111,936 captured modeled radiances with matching "
                    "coordinates. Whole-run guard outcomes and peaks include subsequent "
                    "validation, which may have been interrupted; numerical agreement "
                    "does not certify a successful full probe or full update/adjoint "
                    "validation. No retrieval.nc, completed.json or RT execution required."
                ),
            }
        )
    else:
        report["equivalent_at_saved_state"] = equivalent
        report["bitwise_equivalent_at_saved_state"] = bitwise_equivalent
        report["scope"] = (
            "Complete saved fixed-state reflectances and ozone gradient; no RT executed."
        )
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.resolve().write_text(text)
    sys.stdout.write(text)
    return 0 if equivalent or args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
