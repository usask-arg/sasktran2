# Copilot instructions for this repo

This repository is a hybrid Python/C++/Rust radiative transfer framework (SASKTRAN2). These instructions capture the essentials so an AI coding agent can be productive fast.

## Big picture architecture
- Core C++ engine lives in `cpp/` (e.g., `cpp/include/sasktran2/engine.h`, `cpp/lib/…`) and is templated by NSTOKES for performance.
- Rust crates in `rust/` provide bindings and performance kernels:
  - `sasktran2-sys`: FFI bridge to the C++ core
  - `sasktran2-rs`: reusable Rust utilities (e.g., optical models)
  - `sasktran2-py-ext`: Python extension built with Maturin, exposes the engine to Python
- Python API in `src/sasktran2/` is the user-facing layer (e.g., `engine.py`, `atmosphere.py`, `config.py`, `viewinggeo/`, `constituent/`, `optical/`). Results are typically returned as labeled xarray DataArrays with full linearization/derivative mapping.
- Data flow: Python builds config + atmosphere ➜ calls Rust Python extension ➜ Rust FFI into C++ engine ➜ results/derivatives mapped back to Python.

## Build, test, and docs (Pixi tasks)
- Build everything (Python extension + C++/Rust): `pixi run build`
- Run Python tests (pytest): `pixi run test`
- Lint/format (Ruff via pre-commit): `pixi run pre-commit`
- Build docs (Sphinx + Doxygen): `pixi run docs`
Notes:
- Rebuild needed after C++ or Rust changes. Python-only edits do not require rebuild.
- Dev env requirements: Rust nightly (see `rust-toolchain.toml`), C++17, BLAS/LAPACK (OpenBLAS/MKL/Accelerate), CMake ≥ 3.20, Python ≥ 3.10.

## Project conventions and patterns
- C++ engine is template-heavy; avoid breaking the NSTOKES template and public headers (`cpp/include/sasktran2/**`).
- Builder-style configuration objects in Python (`src/sasktran2/config.py`).
- Optical property models split across Python (`src/sasktran2/optical/`) and Rust (`rust/sasktran2-rs/src/optical/`).
- New atmospheric constituents typically live under `src/sasktran2/constituent/` and integrate with the engine via the Python API.
- SKTRAN_DISCO (discrete ordinates) components are in `cpp/lib/sktran_disco/` and related headers.
- Results are xarray-first; keep coordinates/attrs consistent when adding outputs or derivatives.

## Testing layout and examples
- Python tests live in `tests/**` (organized by domain: `atmosphere/`, `engine/`, `mie/`, `viewing_geometry/`, etc.). Use pytest.
- C++ tests live in `cpp/lib/tests/` with a custom framework.
- Many development notebooks/scripts in `local_scripts/` demonstrate usage and performance benchmarks.

## Common change recipes
- Add a constituent: implement in `src/sasktran2/constituent/`, wire opticals in `src/sasktran2/optical/` (and optionally `rust/sasktran2-rs/src/optical/`), then validate via `pixi run test`.
- Modify optical properties: update models/data in the Python/Rust optical folders; ensure C++ interfaces still match via `sasktran2-sys`.
- Engine changes: edit `cpp/include/sasktran2/engine.h` (and related sources), then `pixi run build` to refresh bindings.

## Integration points and external deps
- BLAS/LAPACK must be discoverable by CMake; on macOS, Apple Accelerate works out of the box.
- Maturin builds the Rust-Python extension; Cargo workspaces are under `rust/`.
- xarray is expected at the Python boundary; keep label semantics intact for downstream tooling.

## Agent do's
- Keep patches minimal and public APIs stable; prefer incremental changes.
- After Rust/C++ edits: run `pixi run build`; before committing: run `pixi run pre-commit`.
- Reference: `CLAUDE.md` for a fuller overview; `README.md` for install/usage; `cpp/CMakeLists.txt`, `pyproject.toml`, and `rust/*/Cargo.toml` for build config.
