"""Dense-layer accuracy and timing harness for source integration methods."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import platform
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import sasktran2 as sk
import xarray as xr
from sasktran2.test_util.source_integration import (
    accuracy_metrics,
    benchmark_callable,
    minimum_layers_for_tolerance,
    observed_orders,
)

LOGGER = logging.getLogger("source-integration-harness")
TOP_OF_ATMOSPHERE_M = 60_000.0
WAVELENGTHS_NM = np.array([350.0, 500.0, 10_000.0])


@dataclass(frozen=True)
class Scenario:
    name: str
    solver: str
    geometry: str
    description: str
    reference_spacing_m: float
    confirmation_spacing_m: float


@dataclass
class ModelRun:
    engine: sk.Engine
    atmosphere: sk.Atmosphere
    result: xr.Dataset
    altitude_m: np.ndarray
    profiles: dict[str, np.ndarray]
    initial_run_seconds: float
    num_lines_of_sight: int
    viewing_labels: list[tuple[str, str]]


SCENARIOS = {
    scenario.name: scenario
    for scenario in [
        Scenario(
            "standard_emission",
            "standard_emission",
            "spherical_limb",
            "Line-of-sight integration of an atmospheric source function.",
            25.0,
            12.5,
        ),
        Scenario(
            "volume_emission_rate",
            "volume_emission_rate",
            "spherical_limb",
            "Line-of-sight integration of a volume emission rate.",
            25.0,
            12.5,
        ),
        Scenario(
            "single_scatter_exact",
            "single_scatter_exact",
            "spherical_limb",
            "Exact solar single-scatter integration through spherical layers.",
            25.0,
            12.5,
        ),
        Scenario(
            "discrete_ordinates_solar",
            "discrete_ordinates_solar",
            "plane_parallel",
            "Solar source interpolation from the discrete-ordinates solution.",
            100.0,
            50.0,
        ),
        Scenario(
            "discrete_ordinates_thermal",
            "discrete_ordinates_thermal",
            "plane_parallel",
            "Thermal source interpolation from the discrete-ordinates solution.",
            100.0,
            50.0,
        ),
        Scenario(
            "successive_orders",
            "successive_orders",
            "spherical_limb",
            "Successive-orders diffuse-source integration.",
            100.0,
            50.0,
        ),
        Scenario(
            "twostream_solar",
            "twostream_solar",
            "plane_parallel",
            "Existing C++ two-stream solar source integration.",
            100.0,
            50.0,
        ),
        Scenario(
            "twostream_thermal",
            "twostream_thermal",
            "plane_parallel",
            "Existing C++ two-stream thermal source integration.",
            100.0,
            50.0,
        ),
    ]
}


def _altitude_grid(spacing_m: float) -> np.ndarray:
    num_layers = round(TOP_OF_ATMOSPHERE_M / spacing_m)
    if spacing_m <= 0 or not np.isclose(num_layers * spacing_m, TOP_OF_ATMOSPHERE_M):
        msg = f"spacing {spacing_m:g} m must divide {TOP_OF_ATMOSPHERE_M:g} m"
        raise ValueError(msg)
    return np.linspace(0.0, TOP_OF_ATMOSPHERE_M, num_layers + 1)


def _configure(
    scenario: Scenario,
    num_threads: int,
    *,
    single_scatter_quadrature: bool,
    single_scatter_solar_transmission: str,
) -> sk.Config:
    config = sk.Config()
    config.num_threads = num_threads
    config.num_stokes = 1
    config.num_streams = 4
    config.num_singlescatter_moments = 8
    config.num_sza = 2
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource

    if scenario.solver == "standard_emission":
        config.emission_source = sk.EmissionSource.Standard
    elif scenario.solver == "volume_emission_rate":
        config.emission_source = sk.EmissionSource.VolumeEmissionRate
    elif scenario.solver == "single_scatter_exact":
        config.single_scatter_source = sk.SingleScatterSource.Exact
        config.single_scatter_source_quadrature = single_scatter_quadrature
        config.single_scatter_solar_transmission = {
            "exact": sk.SingleScatterSolarTransmission.Exact,
            "ray-table": sk.SingleScatterSolarTransmission.RayTable,
        }[single_scatter_solar_transmission]
    elif scenario.solver in {
        "discrete_ordinates_solar",
        "discrete_ordinates_thermal",
    }:
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        if scenario.solver == "discrete_ordinates_thermal":
            config.emission_source = sk.EmissionSource.DiscreteOrdinates
    elif scenario.solver == "successive_orders":
        config.single_scatter_source = sk.SingleScatterSource.Exact
        config.single_scatter_source_quadrature = True
        config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
        config.num_successive_orders_iterations = 3
        config.num_successive_orders_incoming = 26
        config.num_successive_orders_outgoing = 26
        config.init_successive_orders_with_discrete_ordinates = False
    elif scenario.solver == "twostream_solar":
        config.num_streams = 2
        config.multiple_scatter_source = sk.MultipleScatterSource.TwoStream
    elif scenario.solver == "twostream_thermal":
        config.num_streams = 2
        config.emission_source = sk.EmissionSource.TwoStream
    else:
        msg = f"unsupported solver {scenario.solver}"
        raise ValueError(msg)
    return config


def _viewing_geometry(scenario: Scenario) -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    if scenario.geometry == "spherical_limb":
        viewing.add_ray(sk.TangentAltitudeSolar(12_345.0, -0.4, 200_000.0, 0.2))
        viewing.add_ray(sk.TangentAltitudeSolar(27_123.0, 0.7, 200_000.0, 0.2))
        viewing.add_ray(sk.GroundViewingSolar(0.2, 0.35, 0.45, 200_000.0))
    else:
        viewing.add_ray(sk.GroundViewingSolar(0.2, -0.4, 0.25, 200_000.0))
        viewing.add_ray(sk.GroundViewingSolar(0.2, 0.7, 0.75, 200_000.0))
    return viewing


def _viewing_labels(scenario: Scenario) -> list[tuple[str, str]]:
    if scenario.geometry == "spherical_limb":
        return [
            ("limb", "tangent_12.345_km"),
            ("limb", "tangent_27.123_km"),
            ("ground", "ground_mu_0.45"),
        ]
    return [("ground", "ground_mu_0.25"), ("ground", "ground_mu_0.75")]


def _continuous_profiles(altitude_m: np.ndarray) -> dict[str, np.ndarray]:
    base_shape = np.exp(-altitude_m / 7_500.0)
    base_shape += 0.08 * np.exp(-0.5 * ((altitude_m - 18_000.0) / 2_000.0) ** 2)
    optical_depth = np.array([0.05, 0.5, 3.0])
    extinction = (
        base_shape[:, np.newaxis]
        * optical_depth[np.newaxis, :]
        / np.trapezoid(base_shape, altitude_m)
    )
    ssa = np.clip(
        0.88
        + 0.06 * np.exp(-altitude_m[:, np.newaxis] / 18_000.0)
        - np.array([0.02, 0.0, 0.04])[np.newaxis, :],
        0.0,
        1.0,
    )
    emission_shape = 2.0e-5 + 8.0e-5 * np.exp(
        -0.5 * ((altitude_m - 15_000.0) / 3_000.0) ** 2
    )
    emission = emission_shape[:, np.newaxis] * np.array([0.7, 1.0, 1.4])
    return {"extinction": extinction, "ssa": ssa, "emission": emission}


def _grid_linear_profiles(altitude_m: np.ndarray) -> dict[str, np.ndarray]:
    """Profiles represented identically on every candidate altitude grid."""
    optical_depth = np.array([0.05, 0.5, 3.0])
    extinction = (
        (1.5 - 0.5 * altitude_m[:, np.newaxis] / TOP_OF_ATMOSPHERE_M)
        * optical_depth[np.newaxis, :]
        / (1.25 * TOP_OF_ATMOSPHERE_M)
    )
    ssa = np.broadcast_to(
        np.array([0.90, 0.92, 0.88])[np.newaxis, :], extinction.shape
    ).copy()
    emission = (
        2.0e-5
        + 6.0e-5 * altitude_m[:, np.newaxis] / TOP_OF_ATMOSPHERE_M
    ) * np.array([0.7, 1.0, 1.4])[np.newaxis, :]
    return {"extinction": extinction, "ssa": ssa, "emission": emission}


def _build_run(
    scenario: Scenario,
    spacing_m: float,
    *,
    derivatives: bool,
    num_threads: int,
    profile_shape: str = "continuous",
    single_scatter_quadrature: bool = True,
    single_scatter_solar_transmission: str = "ray-table",
) -> ModelRun:
    start = perf_counter()
    altitude_m = _altitude_grid(spacing_m)
    config = _configure(
        scenario,
        num_threads,
        single_scatter_quadrature=single_scatter_quadrature,
        single_scatter_solar_transmission=single_scatter_solar_transmission,
    )
    geometry_type = (
        sk.GeometryType.Spherical
        if scenario.geometry == "spherical_limb"
        else sk.GeometryType.PlaneParallel
    )
    geometry = sk.Geometry1D(
        cos_sza=0.2,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitude_m,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=geometry_type,
    )
    viewing = _viewing_geometry(scenario)
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=derivatives,
        pressure_derivative=False,
        temperature_derivative=False,
        specific_humidity_derivative=False,
        legendre_derivative=False,
    )
    profiles = (
        _continuous_profiles(altitude_m)
        if profile_shape == "continuous"
        else _grid_linear_profiles(altitude_m)
    )
    atmosphere.storage.total_extinction[:] = profiles["extinction"]
    atmosphere.storage.ssa[:] = profiles["ssa"]
    atmosphere.storage.solar_irradiance[:] = 1.0
    if config.emission_source != sk.EmissionSource.NoSource:
        atmosphere.storage.emission_source[:] = profiles["emission"]
    atmosphere.leg_coeff.a1[0, :, :] = 1.0
    atmosphere.leg_coeff.a1[1, :, :] = 0.08
    atmosphere.leg_coeff.a1[2, :, :] = 0.45
    atmosphere.surface.albedo[:] = np.array([0.05, 0.15, 0.3])
    if config.emission_source != sk.EmissionSource.NoSource:
        atmosphere.surface.emission[:] = np.array([1.0e-5, 1.3e-5, 1.6e-5])

    if scenario.solver in {
        "standard_emission",
        "volume_emission_rate",
        "discrete_ordinates_thermal",
        "twostream_thermal",
    }:
        atmosphere.storage.solar_irradiance[:] = 0.0
    if scenario.solver in {"standard_emission", "volume_emission_rate"}:
        atmosphere.storage.ssa[:] = 0.0
        atmosphere.surface.albedo[:] = 0.0
        atmosphere.surface.emission[:] = 0.0
    if scenario.solver == "single_scatter_exact":
        atmosphere.surface.albedo[:] = 0.0

    engine = sk.Engine(config, geometry, viewing)
    result = engine.calculate_radiance(atmosphere)
    return ModelRun(
        engine=engine,
        atmosphere=atmosphere,
        result=result,
        altitude_m=altitude_m,
        profiles=profiles,
        initial_run_seconds=perf_counter() - start,
        num_lines_of_sight=len(result.los),
        viewing_labels=_viewing_labels(scenario),
    )


def _direction_modes(altitude_m: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "uniform": np.ones_like(altitude_m),
        "linear": 2.0 * altitude_m / TOP_OF_ATMOSPHERE_M - 1.0,
        "gaussian": np.exp(-0.5 * ((altitude_m - 15_000.0) / 3_000.0) ** 2),
    }


def _field_directions(run: ModelRun) -> dict[tuple[str, str], np.ndarray]:
    directions: dict[tuple[str, str], np.ndarray] = {}
    for mode_name, mode in _direction_modes(run.altitude_m).items():
        directions[("extinction", mode_name)] = (
            run.profiles["extinction"] * mode[:, np.newaxis]
        )
        directions[("ssa", mode_name)] = 0.01 * np.broadcast_to(
            mode[:, np.newaxis], run.profiles["ssa"].shape
        )
        directions[("emission", mode_name)] = (
            run.profiles["emission"] * mode[:, np.newaxis]
        )
    return directions


def _contract_derivatives(run: ModelRun) -> dict[tuple[str, str], np.ndarray]:
    contracted: dict[tuple[str, str], np.ndarray] = {}
    for (field, mode), direction in _field_directions(run).items():
        wf_name = f"wf_{field}"
        if wf_name not in run.result:
            continue
        weights = xr.DataArray(direction, dims=["altitude", "wavelength"])
        contracted[(field, mode)] = (
            (run.result[wf_name] * weights).sum("altitude").values
        )
    return contracted


def _finite_difference_rows(
    scenario: Scenario,
    run: ModelRun,
    *,
    spacing_m: float,
    step: float,
    floor_fraction: float,
) -> list[dict[str, Any]]:
    analytic = _contract_derivatives(run)
    arrays = {
        "extinction": run.atmosphere.storage.total_extinction,
        "ssa": run.atmosphere.storage.ssa,
    }
    if "wf_emission" in run.result:
        arrays["emission"] = run.atmosphere.storage.emission_source
    rows = []
    for key, direction in _field_directions(run).items():
        if key not in analytic:
            continue
        field, mode = key
        target = arrays[field]
        original = target.copy()
        plus_values = original + step * direction
        minus_values = original - step * direction
        if field == "ssa":
            plus_valid = bool(np.all((plus_values >= 0) & (plus_values <= 1)))
            minus_valid = bool(np.all((minus_values >= 0) & (minus_values <= 1)))
        else:
            plus_valid = bool(np.all(plus_values >= 0))
            minus_valid = bool(np.all(minus_values >= 0))
        if not plus_valid and not minus_valid:
            LOGGER.warning(
                "Skipping infeasible finite-difference direction %s/%s for %s",
                field,
                mode,
                scenario.name,
            )
            continue
        try:
            if plus_valid:
                target[:] = plus_values
                above = run.engine.calculate_radiance(
                    run.atmosphere
                ).radiance.values.copy()
            if minus_valid:
                target[:] = minus_values
                below = run.engine.calculate_radiance(
                    run.atmosphere
                ).radiance.values.copy()
        finally:
            target[:] = original
        if plus_valid and minus_valid:
            numerical = (above - below) / (2.0 * step)
            scheme = "central"
        elif plus_valid:
            numerical = (above - run.result.radiance.values) / step
            scheme = "forward"
        else:
            numerical = (run.result.radiance.values - below) / step
            scheme = "backward"
        metrics = accuracy_metrics(
            analytic[key], numerical, floor_fraction=floor_fraction
        )
        rows.append(
            {
                "scenario": scenario.name,
                "spacing_m": spacing_m,
                "field": field,
                "mode": mode,
                "scheme": scheme,
                **metrics.as_dict(),
            }
        )
    return rows


def _safe_ratio(value: float, baseline: float) -> float | None:
    if baseline == 0:
        return 1.0 if value == 0 else None
    return value / baseline


def _viewing_accuracy_rows(
    scenario: Scenario,
    run: ModelRun,
    reference: ModelRun,
    *,
    comparison: str,
    spacing_m: float,
    floor_fraction: float,
) -> list[dict[str, Any]]:
    rows = []
    for los_index, (viewing_kind, viewing_label) in enumerate(run.viewing_labels):
        metrics = accuracy_metrics(
            run.result.radiance.values[:, los_index, ...],
            reference.result.radiance.values[:, los_index, ...],
            floor_fraction=floor_fraction,
        )
        rows.append(
            {
                "scenario": scenario.name,
                "comparison": comparison,
                "spacing_m": spacing_m,
                "los_index": los_index,
                "viewing_kind": viewing_kind,
                "viewing_label": viewing_label,
                **metrics.as_dict(),
            }
        )
    return rows


def _package_version() -> str:
    try:
        return version("sasktran2")
    except PackageNotFoundError:
        return "unknown"


def _baseline_rows(path: Path | None) -> dict[tuple[str, float], dict[str, Any]]:
    if path is None:
        return {}
    with path.open() as stream:
        report = json.load(stream)
    return {
        (row["scenario"], float(row["spacing_m"])): row for row in report["accuracy"]
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0])
    fieldnames.extend(
        sorted({key for row in rows for key in row}.difference(fieldnames))
    )
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _parse_scenarios(value: str) -> list[Scenario]:
    names = list(SCENARIOS) if value == "all" else value.split(",")
    unknown = set(names).difference(SCENARIOS)
    if unknown:
        msg = f"unknown scenarios: {', '.join(sorted(unknown))}"
        raise ValueError(msg)
    return [SCENARIOS[name] for name in names]


def _parse_spacings(value: str) -> list[float]:
    spacings = sorted({float(item) for item in value.split(",")}, reverse=True)
    for spacing in spacings:
        _altitude_grid(spacing)
    return spacings


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenarios", default="all", help="comma-separated or all")
    parser.add_argument(
        "--spacings", default="2000,1000,500,250,125", help="coarse grids in m"
    )
    parser.add_argument("--output", type=Path, default=Path("build/source-integration"))
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument(
        "--profile-shape",
        choices=["continuous", "grid-linear"],
        default="continuous",
        help=(
            "continuous measures end-to-end altitude-grid convergence; "
            "grid-linear isolates source-integration convergence"
        ),
    )
    parser.add_argument(
        "--single-scatter-solar-transmission",
        choices=["exact", "ray-table"],
        default="ray-table",
        help="solar-path provider for candidate single-scatter runs",
    )
    parser.add_argument(
        "--single-scatter-candidate-integration",
        choices=["endpoint", "gauss8"],
        default="gauss8",
        help="within-cell rule for candidate single-scatter runs",
    )
    parser.add_argument(
        "--single-scatter-reference-integration",
        choices=["endpoint", "gauss8"],
        default="endpoint",
        help=(
            "independent dense-reference rule; endpoint avoids sharing the "
            "candidate's Gauss-8 implementation"
        ),
    )
    parser.add_argument(
        "--reference-spacing",
        type=float,
        help="override each scenario's primary dense spacing",
    )
    parser.add_argument(
        "--confirmation-spacing",
        type=float,
        help="override each scenario's reference-confirmation spacing",
    )
    parser.add_argument("--derivatives", action="store_true")
    parser.add_argument("--finite-difference", action="store_true")
    parser.add_argument("--finite-difference-step", type=float, default=1e-4)
    parser.add_argument("--timing-samples", type=int, default=0)
    parser.add_argument("--timing-warmups", type=int, default=3)
    parser.add_argument("--minimum-sample-seconds", type=float, default=0.05)
    parser.add_argument("--floor-fraction", type=float, default=1e-12)
    parser.add_argument("--target-tolerance", type=float, default=1e-3)
    parser.add_argument("--reference-tolerance", type=float, default=1e-4)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--max-error-ratio", type=float, default=1.10)
    parser.add_argument("--max-runtime-ratio", type=float, default=1.02)
    parser.add_argument("--fail-on-gate", action="store_true")
    parser.add_argument("--list-scenarios", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if args.list_scenarios:
        for scenario in SCENARIOS.values():
            LOGGER.info("%-28s %s", scenario.name, scenario.description)
        return 0

    scenarios = _parse_scenarios(args.scenarios)
    spacings = _parse_spacings(args.spacings)
    derivatives = args.derivatives or args.finite_difference
    baseline = _baseline_rows(args.baseline)
    accuracy_rows: list[dict[str, Any]] = []
    viewing_accuracy_rows: list[dict[str, Any]] = []
    directional_rows: list[dict[str, Any]] = []
    finite_difference_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    gate_failures: list[str] = []

    for scenario in scenarios:
        LOGGER.info("Running %s", scenario.name)
        reference_spacing_m = (
            scenario.reference_spacing_m
            if args.reference_spacing is None
            else args.reference_spacing
        )
        confirmation_spacing_m = (
            scenario.confirmation_spacing_m
            if args.confirmation_spacing is None
            else args.confirmation_spacing
        )
        _altitude_grid(reference_spacing_m)
        _altitude_grid(confirmation_spacing_m)
        if confirmation_spacing_m >= reference_spacing_m:
            msg = "confirmation spacing must be finer than reference spacing"
            raise ValueError(msg)
        confirmation = _build_run(
            scenario,
            confirmation_spacing_m,
            derivatives=derivatives,
            num_threads=args.num_threads,
            profile_shape=args.profile_shape,
            single_scatter_quadrature=(
                args.single_scatter_reference_integration == "gauss8"
            ),
            single_scatter_solar_transmission="exact",
        )
        reference = _build_run(
            scenario,
            reference_spacing_m,
            derivatives=derivatives,
            num_threads=args.num_threads,
            profile_shape=args.profile_shape,
            single_scatter_quadrature=(
                args.single_scatter_reference_integration == "gauss8"
            ),
            single_scatter_solar_transmission="exact",
        )
        reference_metrics = accuracy_metrics(
            reference.result.radiance.values,
            confirmation.result.radiance.values,
            floor_fraction=args.floor_fraction,
        )
        viewing_accuracy_rows.extend(
            _viewing_accuracy_rows(
                scenario,
                reference,
                confirmation,
                comparison="reference_confirmation",
                spacing_m=reference_spacing_m,
                floor_fraction=args.floor_fraction,
            )
        )
        if reference_metrics.max_normalized > args.reference_tolerance:
            gate_failures.append(
                f"{scenario.name}: dense-reference error "
                f"{reference_metrics.max_normalized:.3e} > "
                f"{args.reference_tolerance:.3e}"
            )

        confirmation_directional = (
            _contract_derivatives(confirmation) if derivatives else {}
        )
        reference_directional = _contract_derivatives(reference) if derivatives else {}
        for key, reference_value in reference_directional.items():
            if key not in confirmation_directional:
                continue
            directional_metrics = accuracy_metrics(
                reference_value,
                confirmation_directional[key],
                floor_fraction=args.floor_fraction,
            )
            directional_rows.append(
                {
                    "scenario": scenario.name,
                    "comparison": "reference_confirmation",
                    "spacing_m": reference_spacing_m,
                    "field": key[0],
                    "mode": key[1],
                    **directional_metrics.as_dict(),
                }
            )
            if directional_metrics.max_normalized > args.reference_tolerance:
                gate_failures.append(
                    f"{scenario.name} {key[0]}/{key[1]}: dense-reference "
                    f"derivative error {directional_metrics.max_normalized:.3e} > "
                    f"{args.reference_tolerance:.3e}"
                )
        scenario_rows = []
        for spacing_m in spacings:
            LOGGER.info("  spacing %g m", spacing_m)
            run = _build_run(
                scenario,
                spacing_m,
                derivatives=derivatives,
                num_threads=args.num_threads,
                profile_shape=args.profile_shape,
                single_scatter_quadrature=(
                    args.single_scatter_candidate_integration == "gauss8"
                ),
                single_scatter_solar_transmission=(
                    args.single_scatter_solar_transmission
                ),
            )
            metrics = accuracy_metrics(
                run.result.radiance.values,
                confirmation.result.radiance.values,
                floor_fraction=args.floor_fraction,
            )
            viewing_accuracy_rows.extend(
                _viewing_accuracy_rows(
                    scenario,
                    run,
                    confirmation,
                    comparison="candidate_confirmation",
                    spacing_m=spacing_m,
                    floor_fraction=args.floor_fraction,
                )
            )
            row: dict[str, Any] = {
                "scenario": scenario.name,
                "solver": scenario.solver,
                "geometry": scenario.geometry,
                "spacing_m": spacing_m,
                "layers": len(run.altitude_m) - 1,
                "reference_spacing_m": reference_spacing_m,
                "confirmation_spacing_m": confirmation_spacing_m,
                "initial_run_seconds": run.initial_run_seconds,
                **metrics.as_dict(),
            }
            if args.timing_samples > 0:
                timing = benchmark_callable(
                    lambda current_run=run: current_run.engine.calculate_radiance(
                        current_run.atmosphere
                    ),
                    warmups=args.timing_warmups,
                    samples=args.timing_samples,
                    minimum_sample_seconds=args.minimum_sample_seconds,
                )
                row.update(
                    {f"runtime_{key}": value for key, value in timing.as_dict().items()}
                )
                work_items = (
                    (len(run.altitude_m) - 1)
                    * len(WAVELENGTHS_NM)
                    * run.num_lines_of_sight
                )
                row["layer_wavelength_los_per_second"] = (
                    work_items / timing.median_seconds
                )

            baseline_row = baseline.get((scenario.name, spacing_m))
            if baseline_row is not None:
                row["error_ratio_to_baseline"] = _safe_ratio(
                    row["max_normalized"], baseline_row["max_normalized"]
                )
                if row["error_ratio_to_baseline"] is None or (
                    row["error_ratio_to_baseline"] > args.max_error_ratio
                ):
                    gate_failures.append(
                        f"{scenario.name} at {spacing_m:g} m: error ratio "
                        f"{row['error_ratio_to_baseline']}"
                    )
                if (
                    "runtime_median_seconds" in row
                    and "runtime_median_seconds" in baseline_row
                ):
                    row["runtime_ratio_to_baseline"] = _safe_ratio(
                        row["runtime_median_seconds"],
                        baseline_row["runtime_median_seconds"],
                    )
                    if row["runtime_ratio_to_baseline"] is None or (
                        row["runtime_ratio_to_baseline"] > args.max_runtime_ratio
                    ):
                        gate_failures.append(
                            f"{scenario.name} at {spacing_m:g} m: runtime ratio "
                            f"{row['runtime_ratio_to_baseline']}"
                        )

            if derivatives:
                for key, candidate in _contract_derivatives(run).items():
                    if key not in confirmation_directional:
                        continue
                    directional_metrics = accuracy_metrics(
                        candidate,
                        confirmation_directional[key],
                        floor_fraction=args.floor_fraction,
                    )
                    directional_rows.append(
                        {
                            "scenario": scenario.name,
                            "comparison": "candidate_confirmation",
                            "spacing_m": spacing_m,
                            "field": key[0],
                            "mode": key[1],
                            **directional_metrics.as_dict(),
                        }
                    )
            scenario_rows.append(row)

        error_values = [row["max_normalized"] for row in scenario_rows]
        ratios = [
            scenario_rows[index]["spacing_m"] / scenario_rows[index + 1]["spacing_m"]
            for index in range(len(scenario_rows) - 1)
        ]
        orders = (
            observed_orders(error_values, ratios)
            if len(error_values) > 1
            else np.array([])
        )
        for index, row in enumerate(scenario_rows):
            order = None if index == len(orders) else float(orders[index])
            row["observed_order_to_next_grid"] = (
                order if order is not None and np.isfinite(order) else None
            )
        accuracy_rows.extend(scenario_rows)
        summaries.append(
            {
                "scenario": scenario.name,
                "reference_max_normalized": reference_metrics.max_normalized,
                "reference_rms_normalized": reference_metrics.rms_normalized,
                "minimum_layers_for_target": minimum_layers_for_tolerance(
                    [row["layers"] for row in scenario_rows],
                    error_values,
                    args.target_tolerance,
                ),
            }
        )
        if args.finite_difference:
            finite_difference_rows.extend(
                _finite_difference_rows(
                    scenario,
                    reference,
                    spacing_m=reference_spacing_m,
                    step=args.finite_difference_step,
                    floor_fraction=args.floor_fraction,
                )
            )

    report = {
        "metadata": {
            "created_utc": datetime.now(UTC).isoformat(),
            "python": sys.version,
            "platform": platform.platform(),
            "sasktran2_version": _package_version(),
            "spacings_m": spacings,
            "num_threads": args.num_threads,
            "profile_shape": args.profile_shape,
            "derivatives": derivatives,
            "timing_samples": args.timing_samples,
            "single_scatter_solar_transmission": (
                args.single_scatter_solar_transmission
            ),
            "single_scatter_candidate_integration": (
                args.single_scatter_candidate_integration
            ),
            "single_scatter_reference_integration": (
                args.single_scatter_reference_integration
            ),
        },
        "summaries": summaries,
        "accuracy": accuracy_rows,
        "viewing_accuracy": viewing_accuracy_rows,
        "directional_derivatives": directional_rows,
        "finite_difference": finite_difference_rows,
        "gate_failures": gate_failures,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "report.json").open("w") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    _write_csv(args.output / "accuracy.csv", accuracy_rows)
    _write_csv(args.output / "viewing_accuracy.csv", viewing_accuracy_rows)
    _write_csv(args.output / "directional_derivatives.csv", directional_rows)
    _write_csv(args.output / "finite_difference.csv", finite_difference_rows)
    _write_csv(args.output / "summaries.csv", summaries)
    LOGGER.info("Wrote metrics to %s", args.output)
    for failure in gate_failures:
        LOGGER.warning("GATE: %s", failure)
    return 1 if args.fail_on_gate and gate_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
