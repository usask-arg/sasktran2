"""Successive-orders source-interpolation benchmark.

Compares the legacy interpolation (one globally oriented angular grid and
bilinear line-of-sight interpolation) with the default frame-aligned grids and
cubic line-of-sight interpolation, as a function of the number of horizontal
source columns. Reports, for each configuration:

- accuracy: maximum relative error of the multiple-scatter radiance against
  that scheme's own dense-column reference at the same angular resolution;
- wall time of engine construction plus a radiance calculation, with and
  without the full Jacobian;
- peak resident memory of a single-threaded run in a fresh process, with the
  diffuse transport nonzeros and source-weight bytes from
  ``SASKTRAN2_PROFILE_MEMORY``.

The script is self-contained: the scenes are synthetic Geometry2D limb scans and
a synthetic orbital-plane track.

Example::

    python tools/benchmarks/so_interpolation_benchmark.py --output build/so-interp
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import sasktran2 as sk

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.arange(0.0, 100_001.0, 1_000.0)
HORIZONTAL_RAD = np.deg2rad(np.arange(-20.0, 20.0 + 1e-9, 0.5))
WAVELENGTHS_NM = np.array([350.0, 525.0, 750.0])
TANGENTS_M = np.arange(10_000.0, 50_001.0, 5_000.0)
SOURCE_ALTITUDES_M = np.arange(1_000.0, 99_001.0, 4_000.0)
SCHEMES = ("legacy", "aligned")

# name: (solar zenith angle at the tangent point in degrees, solar azimuth)
CASES = {
    "fwd30": (30.0, 0.0),
    "fwd60": (60.0, 0.0),
    "back60": (60.0, np.pi),
    "side60": (60.0, np.pi / 2),
    "oblique70": (70.0, np.pi / 4),
    "fwd75": (75.0, 0.0),
    "back75": (75.0, np.pi),
    "side80": (80.0, np.pi / 2),
}


def make_config(
    scheme: str,
    num_columns: int,
    directions: int,
    threads: int,
    *,
    single_scatter: bool,
    multiple_scatter: bool = True,
) -> sk.Config:
    config = sk.Config()
    config.num_threads = threads
    config.num_stokes = 1
    config.num_streams = 16
    config.num_singlescatter_moments = 16
    config.single_scatter_source = (
        sk.SingleScatterSource.Exact
        if single_scatter
        else sk.SingleScatterSource.NoSource
    )
    config.multiple_scatter_source = (
        sk.MultipleScatterSource.SuccessiveOrders
        if multiple_scatter
        else sk.MultipleScatterSource.NoSource
    )
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    config.delta_m_scaling = False
    config.num_successive_orders_incoming = directions
    config.num_successive_orders_outgoing = directions
    config.num_successive_orders_iterations = 50
    config.successive_orders_relative_tolerance = 1e-7
    config.successive_orders_altitude_grid_m = SOURCE_ALTITUDES_M
    config.num_sza = num_columns
    config.successive_orders_legacy_interpolation = scheme == "legacy"
    return config


def make_geometry(case: str) -> sk.Geometry2D:
    sza, saa = CASES[case]
    return sk.Geometry2D(
        cos_sza=float(np.cos(np.deg2rad(sza))),
        solar_azimuth=float(saa),
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=HORIZONTAL_RAD,
    )


def make_viewing() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    for tangent in TANGENTS_M:
        viewing.add_ray(
            sk.TangentAltitude(
                tangent_altitude_m=float(tangent),
                observer_altitude_m=600_000.0,
                horizontal_angle_radians=0.0,
                viewing_azimuth_radians=0.0,
            )
        )
    return viewing


def fill_atmosphere(
    atmosphere: sk.Atmosphere, horizontal: np.ndarray, altitudes: np.ndarray
) -> None:
    """Rayleigh, ozone and stratospheric aerosol with albedo 0.3."""
    _, z = np.meshgrid(horizontal, altitudes, indexing="ij")
    z = z.ravel()
    wavelength = WAVELENGTHS_NM[np.newaxis, :]
    rayleigh = (2.55e25 * np.exp(-z / 7_500.0))[:, None] * (
        4.02e-32 * (wavelength / 1000.0) ** -4.04
    )
    ozone = (5e18 * np.exp(-((z - 22_000.0) ** 2) / (2 * 6_000.0**2)))[
        :, None
    ] * np.array([1e-27, 4.5e-25, 1e-25])[np.newaxis, :]
    aerosol = (1.5e-7 * np.exp(-((z - 18_000.0) ** 2) / (2 * 6_000.0**2)))[
        :, None
    ] * (wavelength / 750.0) ** -1.5
    extinction = rayleigh + ozone + aerosol
    scattering = rayleigh + 0.99 * aerosol
    atmosphere.storage.total_extinction[:] = extinction
    atmosphere.storage.ssa[:] = scattering / extinction
    num_legendre = atmosphere.leg_coeff.a1.shape[0]
    order = np.arange(num_legendre)
    henyey_greenstein = (2 * order + 1) * 0.7**order
    rayleigh_phase = np.zeros(num_legendre)
    rayleigh_phase[0] = 1.0
    rayleigh_phase[2] = 0.5
    atmosphere.leg_coeff.a1[:] = (
        rayleigh_phase[:, None, None] * (rayleigh / scattering)[None]
        + henyey_greenstein[:, None, None] * (0.99 * aerosol / scattering)[None]
    )
    atmosphere.surface.albedo[:] = 0.3


def run_2d(
    case: str,
    scheme: str,
    num_columns: int,
    directions: int,
    threads: int,
    *,
    single_scatter: bool = False,
    multiple_scatter: bool = True,
    derivatives: bool = False,
) -> tuple[np.ndarray, float]:
    config = make_config(
        scheme,
        num_columns,
        directions,
        threads,
        single_scatter=single_scatter,
        multiple_scatter=multiple_scatter,
    )
    geometry = make_geometry(case)
    start = time.perf_counter()
    engine = sk.Engine(config, geometry, make_viewing())
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=WAVELENGTHS_NM,
        calculate_derivatives=derivatives,
    )
    fill_atmosphere(atmosphere, HORIZONTAL_RAD, ALTITUDES_M)
    radiance = engine.calculate_radiance(atmosphere).radiance.values
    return (
        np.asarray(radiance).reshape(len(WAVELENGTHS_NM), -1),
        time.perf_counter() - start,
    )


def memory_probe(argv: list[str]) -> None:
    """Child process: run one configuration and print its peak memory."""
    case, scheme, num_columns, directions = argv[0], argv[1], int(argv[2]), int(argv[3])
    run_2d(case, scheme, num_columns, directions, 1, single_scatter=True)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    print(f"PEAK_RSS_BYTES {peak_bytes}", flush=True)


def measure_memory(case: str, scheme: str, num_columns: int, directions: int) -> dict:
    environment = dict(os.environ, SASKTRAN2_PROFILE_MEMORY="1")
    completed = subprocess.run(
        [
            sys.executable,
            __file__,
            "--memory-probe",
            case,
            scheme,
            str(num_columns),
            str(directions),
        ],
        capture_output=True,
        text=True,
        env=environment,
        check=True,
    )
    return parse_memory(completed)


def parse_memory(completed: subprocess.CompletedProcess) -> dict:
    """Peak RSS plus the diffuse and LOS interpolation records.

    Geometry construction compiles diffuse incoming rays first and observer
    lines of sight second, so their interpolation records appear in that order
    (once per local engine for orbital-plane runs; they are summed).
    """
    result = {
        "peak_rss_mb": None,
        "diffuse_transport_nonzeros": 0,
        "diffuse_source_weight_mb": 0.0,
        "los_transport_nonzeros": 0,
        "los_source_weight_mb": 0.0,
    }
    for line in completed.stdout.splitlines():
        if line.startswith("PEAK_RSS_BYTES"):
            result["peak_rss_mb"] = int(line.split()[1]) / 2**20
    records = [
        json.loads(line.split("SASKTRAN2_MEMORY ", 1)[1])
        for line in completed.stderr.splitlines()
        if '"kind":"interpolation_compaction"' in line
    ]
    for index, record in enumerate(records):
        prefix = "diffuse" if index % 2 == 0 else "los"
        result[f"{prefix}_transport_nonzeros"] += record["transport_nonzeros"]
        result[f"{prefix}_source_weight_mb"] += record["source_weight_bytes"] / 2**20
    return result


def benchmark_2d(args: argparse.Namespace) -> dict:
    results: dict = {}
    for case in args.cases:
        single, _ = run_2d(
            case,
            "legacy",
            1,
            args.directions,
            args.threads,
            single_scatter=True,
            multiple_scatter=False,
        )
        references = {
            scheme: run_2d(
                case, scheme, args.reference_columns, args.directions, args.threads
            )[0]
            for scheme in SCHEMES
        }
        case_result = {
            "reference_scheme_difference": float(
                np.max(np.abs(references["aligned"] / references["legacy"] - 1))
            ),
            "rows": [],
        }
        for scheme in SCHEMES:
            for num_columns in args.columns:
                multiple, _ = run_2d(
                    case, scheme, num_columns, args.directions, args.threads
                )
                error = np.abs(multiple - references[scheme])
                row = {
                    "scheme": scheme,
                    "columns": num_columns,
                    "max_relative_ms_error": float(np.max(error / references[scheme])),
                    "max_relative_total_error": float(
                        np.max(error / (references[scheme] + single))
                    ),
                }
                _, row["radiance_seconds"] = run_2d(
                    case,
                    scheme,
                    num_columns,
                    args.directions,
                    args.threads,
                    single_scatter=True,
                )
                if not args.skip_jacobian:
                    _, row["jacobian_seconds"] = run_2d(
                        case,
                        scheme,
                        num_columns,
                        args.directions,
                        args.threads,
                        single_scatter=True,
                        derivatives=True,
                    )
                if not args.skip_memory and case == args.cases[0]:
                    row.update(
                        measure_memory(case, scheme, num_columns, args.directions)
                    )
                case_result["rows"].append(row)
                print(case, json.dumps(row), flush=True)
        results[case] = case_result
    return results


def orbital_scene(
    num_columns: int, scheme: str, directions: int, threads: int, images: int = 16
):
    """Vertical limb scans along a synthetic orbit track in the x-z plane.

    Images are 60 s apart and grouped two per local engine; the solar zenith
    angle of each group ranges from 35 to 75 degrees with latitude.
    """
    image_angles = np.linspace(-0.35, 0.35, images)
    times, observers, tangents, slices = [], [], [], []
    start = np.datetime64("2026-01-01T00:00:00", "ns")
    for image, angle in enumerate(image_angles):
        coordinate = angle + 0.8
        up = np.array([np.sin(coordinate - 0.8), 0.0, np.cos(coordinate - 0.8)])
        forward = np.cross(np.array([0.0, 1.0, 0.0]), up)
        forward /= np.linalg.norm(forward)
        for tangent_altitude in TANGENTS_M:
            tangent_radius = EARTH_RADIUS_M + tangent_altitude
            distance = np.sqrt((EARTH_RADIUS_M + 700_000.0) ** 2 - tangent_radius**2)
            tangent = tangent_radius * up
            times.append(start + image * np.timedelta64(60, "s"))
            observers.append(tangent - distance * forward)
            tangents.append(tangent)
            slices.append(image)
    viewing = sk.OrbitalPlaneViewingGeometry.from_tangent_locations(
        np.asarray(times),
        np.asarray(observers),
        np.asarray(tangents),
        vertical_slice=np.asarray(slices),
    )
    geometry = viewing.construct_atmosphere_geometry(ALTITUDES_M, np.deg2rad(0.5))

    class SolarHandler:
        """Solar zenith varies from 35 to 75 degrees along the track."""

        def target_solar_angles(self, latitude, _longitude, _altitude, _time):
            return 35.0 + 40.0 * min(1.0, abs(float(latitude)) / 25.0), 40.0

    config = make_config(scheme, num_columns, directions, threads, single_scatter=False)
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=120,
        solar_handler=SolarHandler(),
    )
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=WAVELENGTHS_NM, calculate_derivatives=False
    )
    fill_atmosphere(atmosphere, np.arange(geometry.shape[0]), ALTITUDES_M)
    return engine, atmosphere


def orbital_memory_probe(argv: list[str]) -> None:
    """Child process: run one orbital-plane configuration single-threaded."""
    scheme, num_columns, directions, images = (
        argv[0],
        int(argv[1]),
        int(argv[2]),
        int(argv[3]),
    )
    engine, atmosphere = orbital_scene(num_columns, scheme, directions, 1, images)
    engine.calculate_radiance(atmosphere)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    print(f"PEAK_RSS_BYTES {peak_bytes}", flush=True)


def benchmark_orbital(args: argparse.Namespace) -> dict:
    results = {"images": args.orbital_images, "rows": []}
    references = {}
    for scheme in SCHEMES:
        engine, atmosphere = orbital_scene(
            args.orbital_reference_columns,
            scheme,
            args.directions,
            args.threads,
            args.orbital_images,
        )
        references[scheme] = engine.calculate_radiance(atmosphere).radiance.values
    for scheme in SCHEMES:
        for num_columns in args.orbital_columns:
            start = time.perf_counter()
            engine, atmosphere = orbital_scene(
                num_columns, scheme, args.directions, args.threads, args.orbital_images
            )
            radiance = engine.calculate_radiance(atmosphere).radiance.values
            row = {
                "scheme": scheme,
                "columns": num_columns,
                "seconds": time.perf_counter() - start,
                "max_relative_ms_error": float(
                    np.max(np.abs(radiance / references[scheme] - 1))
                ),
            }
            if not args.skip_memory:
                completed = subprocess.run(
                    [
                        sys.executable,
                        __file__,
                        "--orbital-memory-probe",
                        scheme,
                        str(num_columns),
                        str(args.directions),
                        str(args.orbital_images),
                    ],
                    capture_output=True,
                    text=True,
                    env=dict(os.environ, SASKTRAN2_PROFILE_MEMORY="1"),
                    check=True,
                )
                row.update(parse_memory(completed))
            results["rows"].append(row)
            print("orbital", json.dumps(row), flush=True)
    return results


def format_optional(value, spec: str) -> str:
    return "" if value is None else format(value, spec)


def write_markdown(results: dict, path: Path) -> None:
    lines = []
    for case, case_result in results.get("standard_2d", {}).items():
        lines.append(f"### {case}")
        lines.append("")
        lines.append(
            f"Dense-reference difference between schemes: {case_result['reference_scheme_difference']:.1e}"
        )
        lines.append("")
        lines.append(
            "| Scheme | Columns | Max MS error | Max total error | Radiance s | Jacobian s "
            "| Peak RSS MB | Diffuse weights MB | LOS weights MB |"
        )
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for row in case_result["rows"]:
            lines.append(
                f"| {row['scheme']} | {row['columns']} | {row['max_relative_ms_error']:.1e} | "
                f"{row['max_relative_total_error']:.1e} | {row['radiance_seconds']:.2f} | "
                f"{format_optional(row.get('jacobian_seconds'), '.2f')} | "
                f"{format_optional(row.get('peak_rss_mb'), '.0f')} | "
                f"{format_optional(row.get('diffuse_source_weight_mb'), '.1f')} | "
                f"{format_optional(row.get('los_source_weight_mb'), '.2f')} |"
            )
        lines.append("")
    if "orbital" in results:
        lines.append("### Orbital plane")
        lines.append("")
        lines.append(
            f"{results['orbital']['images']} images of {len(TANGENTS_M)} tangent altitudes."
        )
        lines.append("")
        lines.append(
            "| Scheme | Columns | Max MS error | Seconds | Peak RSS MB | Diffuse weights MB | LOS weights MB |"
        )
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for row in results["orbital"]["rows"]:
            lines.append(
                f"| {row['scheme']} | {row['columns']} | {row['max_relative_ms_error']:.1e} | {row['seconds']:.1f} | "
                f"{format_optional(row.get('peak_rss_mb'), '.0f')} | "
                f"{format_optional(row.get('diffuse_source_weight_mb'), '.1f')} | "
                f"{format_optional(row.get('los_source_weight_mb'), '.2f')} |"
            )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    if len(sys.argv) > 1 and sys.argv[1] == "--memory-probe":
        memory_probe(sys.argv[2:])
        return
    if len(sys.argv) > 1 and sys.argv[1] == "--orbital-memory-probe":
        orbital_memory_probe(sys.argv[2:])
        return
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--cases", nargs="*", default=list(CASES))
    parser.add_argument("--columns", nargs="*", type=int, default=[5, 7, 9, 11])
    parser.add_argument("--reference-columns", type=int, default=81)
    parser.add_argument("--directions", type=int, default=110)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-jacobian", action="store_true")
    parser.add_argument("--skip-memory", action="store_true")
    parser.add_argument("--skip-orbital", action="store_true")
    parser.add_argument("--orbital-columns", nargs="*", type=int, default=[5, 7, 11])
    parser.add_argument("--orbital-reference-columns", type=int, default=31)
    parser.add_argument("--orbital-images", type=int, default=16)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    results = {"directions": args.directions, "standard_2d": benchmark_2d(args)}
    if not args.skip_orbital:
        results["orbital"] = benchmark_orbital(args)
    (args.output / "results.json").write_text(json.dumps(results, indent=1))
    write_markdown(results, args.output / "results.md")


if __name__ == "__main__":
    main()
