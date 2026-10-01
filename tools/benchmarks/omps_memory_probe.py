"""Run the saved resident OMPS objective under its physical-footprint guard.

Build and extract a candidate wheel separately, then point --candidate-package
at its extraction directory. The pinned processor environment is never changed.
--output must name a fresh directory beneath build/omps-memory; calculation
artifacts are written to OUTPUT/LABEL. Only one calculation runs at a time,
using the existing analysis repository's lock and a 27 GiB (28.991 GB) ceiling.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
ANALYSIS = REPOSITORY.parent / "data-analysis"
REFERENCE_SUMMARY = (
    ANALYSIS
    / "outputs/omps_paper_5113/investigation/memory_logs"
    / "source_memory_ozone_ms5_in110_out110.summary.json"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def option(command: list[str], name: str) -> str:
    if command.count(name) != 1:
        msg = f"Expected exactly one {name} in saved command"
        raise ValueError(msg)
    return command[command.index(name) + 1]


def replace_option(command: list[str], name: str, value: str) -> None:
    option(command, name)
    command[command.index(name) + 1] = value


def absolute(path: str, root: Path) -> Path:
    return (root / path).resolve()


def checked_output(command: list[str], **kwargs) -> str:
    return subprocess.check_output(command, text=True, **kwargs).strip()


def runtime_identity(python: Path, environment: dict[str, str], root: Path) -> dict:
    # Importing the packages establishes which binary is used, without building
    # an engine or evaluating radiative transfer.
    code = (
        "import json, sys, numpy, sasktran2, sasktran2._core_rust as core; "
        "print(json.dumps(dict(python=sys.executable, "
        "python_version=sys.version, numpy_version=numpy.__version__, "
        "package=sasktran2.__file__, binary=core.__file__, "
        "version=getattr(sasktran2, '__version__', None))))"
    )
    identity = json.loads(
        checked_output([str(python), "-c", code], cwd=root, env=environment)
    )
    identity["binary_sha256"] = sha256(Path(identity["binary"]))
    identity["package_init_sha256"] = sha256(Path(identity["package"]))
    return identity


def package_manifest(package: Path) -> dict[str, str]:
    return {
        str(path.relative_to(package)): sha256(path)
        for path in sorted((package / "sasktran2").rglob("*"))
        if path.is_file() and "__pycache__" not in path.parts
    }


def validate_candidate_provenance(
    identity: dict,
    candidate: Path | None,
    build_provenance: dict | None,
    wheel: Path | None,
) -> dict:
    """Bind supplied build and wheel evidence to the package actually imported."""
    verified = {
        "build_binary_matches_runtime": None,
        "wheel_binary_matches_runtime": None,
        "wheel_package_matches_extracted_package": None,
    }
    if build_provenance is not None:
        if not build_provenance.get("build_completed_utc"):
            msg = "Supplied build provenance does not record a completed build"
            raise ValueError(msg)
        if build_provenance.get("binary_sha256") != identity["binary_sha256"]:
            msg = "Supplied build provenance does not match the imported runtime binary"
            raise ValueError(msg)
        verified["build_binary_matches_runtime"] = True
    if wheel is not None:
        if build_provenance is not None and build_provenance.get(
            "wheel_sha256"
        ) != sha256(wheel):
            msg = "Supplied wheel does not match the build provenance"
            raise ValueError(msg)
        with zipfile.ZipFile(wheel) as archive:
            entries = [
                entry
                for entry in archive.infolist()
                if entry.filename.startswith("sasktran2/") and not entry.is_dir()
            ]
            if len({entry.filename for entry in entries}) != len(entries):
                msg = "Supplied wheel contains duplicate package entries"
                raise ValueError(msg)
            wheel_files = {
                entry.filename: hashlib.sha256(archive.read(entry)).hexdigest()
                for entry in entries
            }
        if wheel_files.get("sasktran2/_core_rust.abi3.so") != identity["binary_sha256"]:
            msg = "Supplied wheel does not contain the imported runtime binary"
            raise ValueError(msg)
        verified["wheel_binary_matches_runtime"] = True
        if candidate is not None:
            if package_manifest(candidate) != wheel_files:
                msg = "Extracted candidate package differs from the supplied wheel"
                raise ValueError(msg)
            verified["wheel_package_matches_extracted_package"] = True
    return verified


def source_manifest() -> dict[str, str | None]:
    """Hash tracked source/configuration files without traversing build outputs."""
    tracked = subprocess.check_output(
        [
            "git",
            "ls-files",
            "-z",
            "--",
            "cpp",
            "rust",
            "src",
            "ci",
            ".cargo",
            "Cargo.toml",
            "Cargo.lock",
            "CMakeLists.txt",
            "pyproject.toml",
            "rust-toolchain.toml",
            "tools/openblas_support.py",
        ],
        cwd=REPOSITORY,
    )
    return {
        name: sha256(REPOSITORY / name) if (REPOSITORY / name).is_file() else None
        for name in sorted(tracked.decode().strip("\0").split("\0"))
        if name
    }


def inject_update_validation(source: str, helper: Path) -> str:
    replacements = {
        "from scipy.optimize import OptimizeResult\n": (
            "from scipy.optimize import OptimizeResult\n"
            "import runpy\n"
            f"OmpsUpdateValidation = runpy.run_path({str(helper)!r})['OmpsUpdateValidation']\n"
        ),
        "        value, gradient = fun(x0)\n        gradient = np.asarray(gradient)\n": (
            "        validation = OmpsUpdateValidation(fun, x0, destination)\n"
            "        captured['update_validation'] = validation\n"
            "        value, gradient = validation.evaluate('initial', x0)\n"
            "        gradient = np.asarray(gradient)\n"
        ),
        "        return OptimizeResult(\n": (
            "        validation.adjoint('initial', x0)\n"
            "        return OptimizeResult(\n"
        ),
        "    print(json.dumps(summary, indent=2), flush=True)\n": (
            "    print(json.dumps(summary, indent=2), flush=True)\n"
            "    if 'update_validation' in captured:\n"
            "        captured['update_validation'].finish()\n"
        ),
    }
    for marker, replacement in replacements.items():
        if source.count(marker) != 1:
            msg = f"Runner changed: expected exactly one injection marker {marker!r}"
            raise ValueError(msg)
        source = source.replace(marker, replacement)
    compile(source, "runner_snapshot.py", "exec")
    return source


def run_guard(command: list[str], root: Path, environment: dict, log: Path) -> int:
    with log.open("w") as stream:
        child = subprocess.Popen(
            command,
            cwd=root,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        try:
            for line in child.stdout:
                stream.write(line)
                stream.flush()
                sys.stdout.write(line)
                sys.stdout.flush()
            return child.wait()
        except BaseException:
            # The guard owns and terminates its calculation's process group.
            child.send_signal(signal.SIGTERM)
            child.wait(timeout=20)
            raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-summary", type=Path, default=REFERENCE_SUMMARY)
    parser.add_argument("--candidate-package", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--columns", type=int, default=5)
    parser.add_argument("--check-updates", action="store_true")
    parser.add_argument(
        "--profile-geometry",
        action="store_true",
        help="Log native memory accounting for geometry, caches and temporary workspace",
    )
    parser.add_argument("--wheel", type=Path)
    parser.add_argument(
        "--build-command", help="Command used to produce the candidate wheel"
    )
    parser.add_argument(
        "--build-provenance",
        type=Path,
        help="JSON recorded at build time with source/compiler configuration; copied into the run",
    )
    parser.add_argument("--analysis-root", type=Path, default=ANALYSIS)
    args = parser.parse_args()
    output = args.output.resolve()
    root = args.analysis_root.resolve()
    summary_path = args.reference_summary.resolve()
    if output.exists():
        parser.error("--output must not already exist")
    if not output.is_relative_to(REPOSITORY / "build/omps-memory"):
        parser.error("--output must be beneath this repository's build/omps-memory")
    if (
        not args.label
        or Path(args.label).name != args.label
        or args.label in {".", ".."}
    ):
        parser.error("--label must be a single directory name")
    if args.columns < 2:
        parser.error("--columns must be at least two")
    build_provenance_path = (
        args.build_provenance.resolve() if args.build_provenance else None
    )
    build_provenance_bytes = (
        build_provenance_path.read_bytes() if build_provenance_path else None
    )
    build_provenance = (
        json.loads(build_provenance_bytes)
        if build_provenance_bytes is not None
        else None
    )
    if build_provenance is not None and not isinstance(build_provenance, dict):
        parser.error("--build-provenance must contain a JSON object")
    saved = json.loads(summary_path.read_text())
    command = list(saved["command"])
    for flag in ("--resident", "--fixed-ozone", "--match-aerosol-time-origin"):
        if flag not in command:
            parser.error(f"Saved command must contain {flag}")
    for name, value in {
        "--threads": "1",
        "--forward-groups": "0",
        "--source-incoming": "110",
        "--source-outgoing": "110",
        "--group-seconds": "120",
    }.items():
        if option(command, name) != value:
            parser.error(f"Saved command must set {name} {value}")
    if "--low-memory" in command or "--resident-forward" in command:
        parser.error("Saved command must retain all forward and derivative groups")
    # A venv executable is often a symlink to the base interpreter. Resolving
    # that symlink would discard the pinned environment's site-packages.
    python = (root / command[0]).absolute()
    command[0] = str(python)
    runner_index = next(i for i, value in enumerate(command) if value.endswith(".py"))
    original_runner = absolute(command[runner_index], root)
    runner = output / "runner_snapshot.py"
    helper = REPOSITORY / "tools/benchmarks/omps_update_validation.py"
    helper_snapshot = output / "update_validation_snapshot.py"
    runner_source = original_runner.read_text()
    if args.check_updates:
        runner_source = inject_update_validation(runner_source, helper_snapshot)
    command[runner_index] = str(runner)
    input_files = {}
    for name in ("--baseline", "--l1g", "--anc", "--custom-scene"):
        path = absolute(option(command, name), root)
        replace_option(command, name, str(path))
        input_files[str(path)] = sha256(path)
    replace_option(command, "--output", str(output))
    replace_option(command, "--label", args.label)
    replace_option(command, "--source-columns", str(args.columns))
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    environment.pop("SASKTRAN2_PROFILE_MEMORY", None)
    if args.profile_geometry:
        environment["SASKTRAN2_PROFILE_MEMORY"] = "1"
    candidate = args.candidate_package.resolve() if args.candidate_package else None
    if candidate:
        if not (candidate / "sasktran2/__init__.py").is_file():
            parser.error("--candidate-package must contain extracted sasktran2 package")
        environment["PYTHONPATH"] = str(candidate)
    identity = runtime_identity(python, environment, root)
    if candidate and any(
        not Path(identity[name]).resolve().is_relative_to(candidate)
        for name in ("package", "binary")
    ):
        parser.error(
            "Pinned Python did not import the requested candidate package and binary"
        )
    try:
        verified_provenance = validate_candidate_provenance(
            identity,
            candidate,
            build_provenance,
            args.wheel.resolve() if args.wheel else None,
        )
    except ValueError as error:
        parser.error(str(error))
    guard = root / "scripts/run_omps_memory_guard.py"
    lock = root / "outputs/omps_paper_5113/investigation/retrieval.lock"
    guarded = [
        str(python),
        "-u",
        str(guard),
        "--max-rss-gib",
        "27",
        "--max-footprint-gib",
        "27",
        "--sample-interval-s",
        "0.5",
        "--min-available-percent",
        "25",
        "--lock",
        str(lock),
        "--log",
        str(output / "guard.jsonl"),
        "--",
        *command,
    ]
    diff = subprocess.check_output(["git", "diff", "--binary", "HEAD"], cwd=REPOSITORY)
    sources = source_manifest()
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_repository": str(REPOSITORY),
        "source_git_sha": checked_output(["git", "rev-parse", "HEAD"], cwd=REPOSITORY),
        "source_git_status": checked_output(
            ["git", "status", "--short"], cwd=REPOSITORY
        ),
        "source_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "source_diff_scope": "git diff --binary HEAD (tracked changes)",
        "source_manifest_at_probe": sources,
        "source_manifest_at_probe_sha256": hashlib.sha256(
            json.dumps(sources, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "source_manifest_scope": (
            "Tracked cpp/rust/src/ci files and build configuration, hashed at probe launch. "
            "Build-time source/compiler identity comes from the supplied build provenance."
        ),
        "build_command": args.build_command,
        "build_provenance": build_provenance,
        "build_provenance_path": (
            str(build_provenance_path) if build_provenance_path else None
        ),
        "build_provenance_sha256": (
            hashlib.sha256(build_provenance_bytes).hexdigest()
            if build_provenance_bytes is not None
            else None
        ),
        "wheel": str(args.wheel.resolve()) if args.wheel else None,
        "wheel_sha256": sha256(args.wheel.resolve()) if args.wheel else None,
        "candidate_package": str(candidate) if candidate else None,
        "candidate_package_files_sha256": (
            package_manifest(candidate) if candidate else None
        ),
        "runtime": identity,
        "verified_candidate_provenance": verified_provenance,
        "reference_summary": str(summary_path),
        "reference_summary_sha256": sha256(summary_path),
        "runner": str(runner),
        "runner_sha256": hashlib.sha256(runner_source.encode()).hexdigest(),
        "original_runner": str(original_runner),
        "original_runner_sha256": sha256(original_runner),
        "check_updates": args.check_updates,
        "update_validation_helper_sha256": (
            sha256(helper) if args.check_updates else None
        ),
        "guard": str(guard),
        "guard_sha256": sha256(guard),
        "inputs_sha256": input_files,
        "cwd": str(root),
        "command": guarded,
        "environment_overrides": {
            name: environment.get(name)
            for name in ("PYTHONPATH", "SASKTRAN2_PROFILE_MEMORY")
        },
        "physical_footprint_ceiling_decimal_gb": 27 * 1024**3 / 1e9,
        "calculation_directory": str(output / args.label),
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "source.patch").write_bytes(diff)
    runner.write_text(runner_source)
    if build_provenance_bytes is not None:
        (output / "build-provenance_snapshot.json").write_bytes(build_provenance_bytes)
    if args.check_updates:
        helper_snapshot.write_bytes(helper.read_bytes())
    provenance_path = output / "provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return_code = run_guard(guarded, root, environment, output / "stdout-stderr.log")
    provenance["guard_return_code"] = return_code
    provenance["finished_utc"] = datetime.now(timezone.utc).isoformat()
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
