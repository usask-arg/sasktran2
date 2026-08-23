# ruff: noqa: EM101, T201
"""Retrieve along-orbit Lambertian albedo from OMPS 350 nm limb radiances.

The objective is the mean squared relative Sun-normalized radiance residual
over the selected tangent-altitude range.  Each objective evaluation uses the
native orbital VJP for one gray albedo parameter per orbital grid position.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import psutil
import sasktran2 as sk
import xarray as xr
from omps_orbital_plane_2d import (
    DEFAULT_ANC,
    DEFAULT_L1G,
    load_omps_inputs,
    make_atmosphere,
    map_ancillary_state,
    select_viewing_geometry,
)
from scipy.optimize import OptimizeResult, minimize


class FixedManual(sk.constituent.Manual):
    """A native optical state that intentionally exposes no retrieval derivative."""

    def register_derivative(
        self, atmo: sk.Atmosphere, name: str
    ) -> None:  # noqa: ARG002
        return


@dataclass(frozen=True)
class Evaluation:
    loss: float
    data_loss: float
    smoothness_loss: float
    gradient: np.ndarray
    modeled: xr.DataArray
    elapsed_s: float
    rss_bytes: int


class AlbedoObjective:
    def __init__(
        self,
        engine: sk.OrbitalPlaneEngine,
        atmosphere: sk.Atmosphere,
        surface: sk.constituent.LambertianSurface2D,
        observed: np.ndarray,
        valid: np.ndarray,
        smoothness: float,
    ) -> None:
        self.engine = engine
        self.atmosphere = atmosphere
        self.surface = surface
        self.observed = observed
        self.valid = valid
        self.smoothness = smoothness
        self.history: list[Evaluation] = []
        self._last_x: np.ndarray | None = None
        self._last: Evaluation | None = None

    def evaluate(self, albedo: np.ndarray) -> Evaluation:
        albedo = np.asarray(albedo, dtype=np.float64)
        if self._last_x is not None and np.array_equal(albedo, self._last_x):
            return self._last

        started = time.perf_counter()
        self.surface.albedo = albedo
        linearization = self.engine.linearize(
            self.atmosphere, prepare_parameters=("surface_albedo",)
        )
        modeled = linearization.value.sel(stokes="I", drop=True)
        modeled_values = np.asarray(modeled.values[0], dtype=np.float64)

        residual = np.zeros_like(modeled_values)
        residual[self.valid] = (
            modeled_values[self.valid] - self.observed[self.valid]
        ) / self.observed[self.valid]
        data_loss = 0.5 * float(np.mean(residual[self.valid] ** 2))

        radiance_cotangent = np.zeros_like(modeled_values)
        radiance_cotangent[self.valid] = residual[self.valid] / (
            np.count_nonzero(self.valid) * self.observed[self.valid]
        )
        cotangent = xr.zeros_like(linearization.value)
        cotangent.loc[{"stokes": "I"}] = xr.DataArray(
            radiance_cotangent[np.newaxis, :],
            dims=modeled.dims,
            coords=modeled.coords,
        )
        gradient = np.asarray(
            linearization.vjp(
                cotangent, parameters=("surface_albedo",)
            ).surface_albedo.values,
            dtype=np.float64,
        )

        differences = np.diff(albedo)
        if self.smoothness > 0 and differences.size:
            smoothness_loss = 0.5 * self.smoothness * float(np.mean(differences**2))
            scale = self.smoothness / differences.size
            gradient[:-1] -= scale * differences
            gradient[1:] += scale * differences
        else:
            smoothness_loss = 0.0

        evaluation = Evaluation(
            loss=data_loss + smoothness_loss,
            data_loss=data_loss,
            smoothness_loss=smoothness_loss,
            gradient=gradient,
            modeled=modeled.copy(deep=True),
            elapsed_s=time.perf_counter() - started,
            rss_bytes=int(psutil.Process().memory_info().rss),
        )
        self.history.append(evaluation)
        self._last_x = albedo.copy()
        self._last = evaluation
        ratio = modeled_values[self.valid] / self.observed[self.valid]
        print(
            f"evaluation {len(self.history):3d}: loss={evaluation.loss:.8e}; "
            f"ratio median={np.median(ratio):.6f}; "
            f"albedo=[{albedo.min():.4f}, {albedo.max():.4f}]; "
            f"|gradient|inf={np.max(np.abs(gradient)):.3e}; "
            f"{evaluation.elapsed_s:.3f} s; "
            f"RSS={evaluation.rss_bytes / 1024**3:.3f} GiB",
            flush=True,
        )
        return evaluation

    def scipy_objective(self, albedo: np.ndarray) -> tuple[float, np.ndarray]:
        evaluation = self.evaluate(albedo)
        return evaluation.loss, evaluation.gradient


def derivative_atmosphere(
    geometry: sk.OrbitalPlaneGeometry,
    config: sk.Config,
    wavelengths_nm: np.ndarray,
    pressure_pa: np.ndarray,
    temperature_k: np.ndarray,
    ozone_vmr: np.ndarray,
    initial_albedo: np.ndarray,
) -> tuple[sk.Atmosphere, sk.constituent.LambertianSurface2D]:
    """Freeze atmospheric optics while retaining only the albedo derivative."""
    optical_atmosphere, _ = make_atmosphere(
        geometry,
        config,
        wavelengths_nm,
        pressure_pa,
        temperature_k,
        ozone_vmr,
        calculate_derivatives=False,
    )
    optical_atmosphere.internal_object()
    extinction = np.asarray(optical_atmosphere.storage.total_extinction).reshape(
        (*geometry.shape, len(wavelengths_nm))
    )
    ssa = np.asarray(optical_atmosphere.storage.ssa).reshape(
        (*geometry.shape, len(wavelengths_nm))
    )
    legendre = np.asarray(optical_atmosphere.storage.leg_coeff).reshape(
        (
            optical_atmosphere.storage.leg_coeff.shape[0],
            *geometry.shape,
            len(wavelengths_nm),
        )
    )

    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=wavelengths_nm,
        calculate_derivatives=True,
        pressure_derivative=False,
        temperature_derivative=False,
        specific_humidity_derivative=False,
        legendre_derivative=False,
    )
    atmosphere.pressure_pa = pressure_pa
    atmosphere.temperature_k = temperature_k
    atmosphere["fixed_optics"] = FixedManual(
        extinction.copy(), ssa.copy(), legendre.copy()
    )
    surface = sk.constituent.LambertianSurface2D(initial_albedo)
    atmosphere["surface"] = surface
    atmosphere.internal_object()
    return atmosphere, surface


def gradient_check(
    objective: AlbedoObjective,
    albedo: np.ndarray,
    relative_step: float = 1.0e-3,
) -> None:
    rng = np.random.default_rng(20260821)
    direction = rng.normal(size=albedo.size)
    direction /= np.linalg.norm(direction)
    room = np.minimum(albedo, 1.0 - albedo)
    step = min(relative_step, 0.25 * float(np.min(room)))
    center = objective.evaluate(albedo)
    plus = objective.evaluate(albedo + step * direction)
    minus = objective.evaluate(albedo - step * direction)
    finite_difference = (plus.loss - minus.loss) / (2.0 * step)
    analytic = float(center.gradient @ direction)
    relative_error = abs(finite_difference - analytic) / max(
        1.0e-12, abs(finite_difference), abs(analytic)
    )
    print(
        "gradient check: "
        f"finite difference={finite_difference:.9e}; "
        f"VJP={analytic:.9e}; relative error={relative_error:.3e}"
    )
    if relative_error > 2.0e-4:
        raise RuntimeError("Albedo VJP gradient check failed")
    objective.evaluate(albedo)


def make_output(
    result: OptimizeResult,
    objective: AlbedoObjective,
    initial_albedo: np.ndarray,
    initial: Evaluation,
    final: Evaluation,
    observed: np.ndarray,
    selected,
    geometry: sk.OrbitalPlaneGeometry,
    wavelengths_nm: np.ndarray,
    nearest_anc_scan: np.ndarray,
    l1g: Path,
    anc: Path,
) -> xr.Dataset:
    initial_modeled = np.asarray(initial.modeled.values[0])
    final_modeled = np.asarray(final.modeled.values[0])
    ratio_initial = initial_modeled / observed
    ratio_final = final_modeled / observed
    history = objective.history
    return xr.Dataset(
        data_vars={
            "retrieved_albedo": ("orbital_position", np.asarray(result.x)),
            "initial_albedo": ("orbital_position", initial_albedo),
            "ground_track_ecef_m": (
                ("orbital_position", "xyz"),
                geometry.ground_track_ecef_m,
            ),
            "surface_radius_m": ("orbital_position", geometry.surface_radii_m),
            "modeled_initial_sun_normalized_radiance": (
                ("wavelength", "los"),
                initial_modeled[np.newaxis, :],
            ),
            "modeled_retrieved_sun_normalized_radiance": (
                ("wavelength", "los"),
                final_modeled[np.newaxis, :],
            ),
            "observed_sun_normalized_radiance": (
                ("wavelength", "los"),
                observed[np.newaxis, :],
            ),
            "initial_model_to_observation_ratio": (
                ("wavelength", "los"),
                ratio_initial[np.newaxis, :],
            ),
            "retrieved_model_to_observation_ratio": (
                ("wavelength", "los"),
                ratio_final[np.newaxis, :],
            ),
            "objective": (
                "evaluation",
                np.asarray([evaluation.loss for evaluation in history]),
            ),
            "data_objective": (
                "evaluation",
                np.asarray([evaluation.data_loss for evaluation in history]),
            ),
            "smoothness_objective": (
                "evaluation",
                np.asarray([evaluation.smoothness_loss for evaluation in history]),
            ),
            "evaluation_time_s": (
                "evaluation",
                np.asarray([evaluation.elapsed_s for evaluation in history]),
            ),
            "evaluation_rss_bytes": (
                "evaluation",
                np.asarray([evaluation.rss_bytes for evaluation in history]),
            ),
        },
        coords={
            "orbital_position": geometry.orbital_positions,
            "along_track_angle_deg": (
                "orbital_position",
                np.rad2deg(geometry.cumulative_angles),
            ),
            "xyz": ["x", "y", "z"],
            "wavelength": wavelengths_nm,
            "los": np.arange(len(observed)),
            "time": ("los", selected.viewing.times),
            "scan_index": ("los", selected.scan_indices),
            "observed_tangent_altitude_m": (
                "los",
                selected.observed_tangent_altitude_m,
            ),
            "nearest_anc_scan": ("orbital_position", nearest_anc_scan),
            "evaluation": np.arange(len(history)),
        },
        attrs={
            "l1g_file": str(l1g),
            "ancillary_file": str(anc),
            "optimizer": "scipy.optimize.minimize L-BFGS-B",
            "optimizer_success": int(bool(result.success)),
            "optimizer_message": str(result.message),
            "optimizer_iterations": int(result.nit),
            "optimizer_function_evaluations": int(result.nfev),
            "initial_objective": initial.loss,
            "final_objective": final.loss,
            "albedo_smoothness_strength": objective.smoothness,
            "fitted_measurement": "Sun-normalized radiance L/E0",
        },
    )


def make_plot(dataset: xr.Dataset, output: Path) -> None:
    scan_indices = np.unique(dataset.scan_index.values)
    initial_scan_ratio = np.full(len(scan_indices), np.nan)
    final_scan_ratio = np.full(len(scan_indices), np.nan)
    for index, scan in enumerate(scan_indices):
        selected = dataset.scan_index.values == scan
        initial_scan_ratio[index] = np.nanmedian(
            dataset.initial_model_to_observation_ratio.values[0, selected]
        )
        final_scan_ratio[index] = np.nanmedian(
            dataset.retrieved_model_to_observation_ratio.values[0, selected]
        )

    figure, axes = plt.subplots(2, 1, figsize=(13, 8), constrained_layout=True)
    axes[0].plot(
        dataset.along_track_angle_deg,
        dataset.retrieved_albedo,
        color="tab:green",
        linewidth=1.8,
    )
    axes[0].set_ylabel("Retrieved Lambertian albedo")
    axes[0].set_xlabel("Atmosphere-grid along-track angle [deg]")
    axes[0].set_ylim(-0.02, 1.02)
    axes[0].grid(alpha=0.25)

    axes[1].plot(scan_indices, initial_scan_ratio, label="Initial", linewidth=1.5)
    axes[1].plot(scan_indices, final_scan_ratio, label="Retrieved", linewidth=1.8)
    axes[1].axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    axes[1].set_xlabel("OMPS image index along orbit")
    axes[1].set_ylabel("Median model / observation, 20-50 km")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    figure.suptitle("OMPS-LP 350 nm orbital Lambertian-albedo retrieval")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160)
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l1g", type=Path, default=DEFAULT_L1G)
    parser.add_argument("--anc", type=Path, default=DEFAULT_ANC)
    parser.add_argument("--slit", type=int, default=1)
    parser.add_argument("--scan-start", type=int, default=0)
    parser.add_argument("--num-scans", type=int)
    parser.add_argument("--scan-step", type=int, default=1)
    parser.add_argument("--wavelength-nm", type=float, default=350.0)
    parser.add_argument("--tangent-altitude-min-km", type=float, default=20.0)
    parser.add_argument("--tangent-altitude-max-km", type=float, default=50.0)
    parser.add_argument("--along-track-angle-deg", type=float, default=1.0)
    parser.add_argument("--path-padding-deg", type=float, default=10.0)
    parser.add_argument("--group-padding-deg", type=float, default=10.0)
    parser.add_argument("--time-group-duration-s", type=float, default=240.0)
    parser.add_argument("--num-sza", type=int, default=5)
    parser.add_argument("--source-altitude-points", type=int, default=25)
    parser.add_argument("--angular-points", type=int, default=110)
    parser.add_argument("--maximum-orders", type=int, default=50)
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument(
        "--max-refraction-tangent-altitude-km",
        type=float,
        default=25.0,
        help="Trace higher-tangent-altitude LOS without refraction; use inf for no cutoff",
    )
    parser.add_argument("--initial-albedo", type=float, default=0.1)
    parser.add_argument(
        "--smoothness",
        type=float,
        default=0.01,
        help="First-difference regularization strength; zero is unregularized",
    )
    parser.add_argument("--maximum-iterations", type=int, default=20)
    parser.add_argument("--function-tolerance", type=float, default=1.0e-9)
    parser.add_argument("--gradient-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--no-refraction", action="store_true")
    parser.add_argument("--gradient-check", action="store_true")
    parser.add_argument(
        "--derivative-execution",
        choices=("resident", "streaming"),
        default="resident",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plot", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 <= args.initial_albedo <= 1:
        raise ValueError("initial-albedo must lie in [0, 1]")
    if args.smoothness < 0:
        raise ValueError("smoothness must be non-negative")
    if args.maximum_iterations < 1:
        raise ValueError("maximum-iterations must be positive")
    if args.num_threads < 1:
        raise ValueError("num-threads must be positive")
    if (
        np.isnan(args.max_refraction_tangent_altitude_km)
        or args.max_refraction_tangent_altitude_km < 0
    ):
        raise ValueError(
            "max-refraction-tangent-altitude-km must be non-negative or inf"
        )

    data = load_omps_inputs(args.l1g, args.anc, args.slit)
    if args.num_scans is None:
        scan_indices = np.arange(args.scan_start, len(data.times), args.scan_step)
    else:
        scan_indices = args.scan_start + args.scan_step * np.arange(args.num_scans)
    if np.any(scan_indices < 0) or np.any(scan_indices >= len(data.times)):
        raise ValueError("Requested scan indices lie outside the OMPS orbit")

    wavelength_index = int(np.argmin(np.abs(data.wavelengths_nm - args.wavelength_nm)))
    wavelengths_nm = data.wavelengths_nm[[wavelength_index]]
    geoid = sk.WGS84()
    selected = select_viewing_geometry(
        data,
        scan_indices,
        args.tangent_altitude_min_km * 1.0e3,
        args.tangent_altitude_max_km * 1.0e3,
        np.asarray([wavelength_index]),
        geoid,
    )
    geometry = selected.viewing.construct_atmosphere_geometry(
        data.atmosphere_altitude_m,
        np.deg2rad(args.along_track_angle_deg),
        path_padding_angle=np.deg2rad(args.path_padding_deg),
    )
    pressure_pa, temperature_k, ozone_vmr, nearest_anc_scan = map_ancillary_state(
        data, geometry, geoid
    )

    config = sk.Config()
    config.num_threads = args.num_threads
    config.threading_model = sk.ThreadingModel.Wavelength
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
    config.occultation_source = sk.OccultationSource.NoSource
    config.los_refraction = not args.no_refraction
    config.los_refraction_max_tangent_altitude_m = (
        args.max_refraction_tangent_altitude_km * 1.0e3
    )
    config.num_sza = args.num_sza
    altitude_edges = np.linspace(
        data.atmosphere_altitude_m[0],
        data.atmosphere_altitude_m[-1],
        args.source_altitude_points + 1,
    )
    config.successive_orders_altitude_grid_m = 0.5 * (
        altitude_edges[:-1] + altitude_edges[1:]
    )
    config.num_successive_orders_incoming = args.angular_points
    config.num_successive_orders_outgoing = args.angular_points
    config.num_successive_orders_iterations = args.maximum_orders
    # Native successive-orders products differentiate the converged fixed
    # point. Use tighter primal convergence than an ordinary radiance run so
    # the optimization gradient and finite differences share that limit.
    config.successive_orders_relative_tolerance = 1.0e-10
    config.successive_orders_absolute_tolerance = 1.0e-13

    initial_albedo = np.full(geometry.shape[0], args.initial_albedo)
    atmosphere, surface = derivative_atmosphere(
        geometry,
        config,
        wavelengths_nm,
        pressure_pa,
        temperature_k,
        ozone_vmr,
        initial_albedo,
    )
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        selected.viewing,
        time_group_duration_s=args.time_group_duration_s,
        group_padding_angle=np.deg2rad(args.group_padding_deg),
        solar_handler=sk.solar.SolarGeometryHandlerAstropy(),
        derivative_execution=args.derivative_execution,
    )
    observed = np.asarray(selected.observed_sun_normalized_radiance[0])
    valid = np.isfinite(observed) & (observed > 0)
    objective = AlbedoObjective(
        engine, atmosphere, surface, observed, valid, args.smoothness
    )
    initial = objective.evaluate(initial_albedo)
    if args.gradient_check:
        gradient_check(objective, initial_albedo)

    started = time.perf_counter()
    result = minimize(
        objective.scipy_objective,
        initial_albedo,
        method="L-BFGS-B",
        jac=True,
        bounds=[(0.0, 1.0)] * initial_albedo.size,
        options={
            "maxiter": args.maximum_iterations,
            "ftol": args.function_tolerance,
            "gtol": args.gradient_tolerance,
            "maxls": 20,
        },
    )
    final = objective.evaluate(np.asarray(result.x))
    elapsed = time.perf_counter() - started

    initial_ratio = np.asarray(initial.modeled.values[0])[valid] / observed[valid]
    final_ratio = np.asarray(final.modeled.values[0])[valid] / observed[valid]
    print(f"Optimizer success: {result.success}; {result.message}")
    print(
        f"Iterations={result.nit}; evaluations={result.nfev}; "
        f"optimization elapsed={elapsed:.3f} s"
    )
    print(f"Composite group/wavelength threads: {args.num_threads}")
    print(
        "Maximum refracted tangent altitude: "
        f"{args.max_refraction_tangent_altitude_km:g} km"
    )
    print(
        "Model/observation ratio percentiles [5, 25, 50, 75, 95]:\n"
        f"  initial: {np.percentile(initial_ratio, [5, 25, 50, 75, 95]).tolist()}\n"
        f"  retrieved: {np.percentile(final_ratio, [5, 25, 50, 75, 95]).tolist()}"
    )
    print(
        "Retrieved albedo percentiles [0, 5, 25, 50, 75, 95, 100]: "
        f"{np.percentile(result.x, [0, 5, 25, 50, 75, 95, 100]).tolist()}"
    )

    dataset = make_output(
        result,
        objective,
        initial_albedo,
        initial,
        final,
        observed,
        selected,
        geometry,
        wavelengths_nm,
        nearest_anc_scan,
        args.l1g,
        args.anc,
    )
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        dataset.to_netcdf(args.output)
        print(f"Wrote {args.output}")
    if args.plot is not None:
        make_plot(dataset, args.plot)
        print(f"Wrote {args.plot}")


if __name__ == "__main__":
    main()
