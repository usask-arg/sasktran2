# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Repository Overview

SASKTRAN2 is a radiative transfer framework developed at the University of Saskatchewan, designed for atmospheric modeling applications. It's a complete reimplementation of the original SASKTRAN with significant performance improvements, full linearizations, and an improved Python interface.

## Build System and Development Commands

This is a hybrid Python/C++/Rust project that uses:
- **Pixi** for environment and dependency management
- **Maturin** for building the Rust-Python extension
- **CMake** for the C++ components
- **Ruff** for Python linting

### Essential Commands

```bash
# Build the entire project (Python extension with C++/Rust components)
pixi run build

# Run tests
pixi run test

# Run linting
pixi run pre-commit

# Run Rust clippy checks
cargo clippy --all-targets --all-features

# Build documentation
pixi run docs
```

### Development Workflow

1. The project uses Pixi environments - ensure you're in the correct environment
2. Changes to Rust code require rebuilding with `pixi run build`
3. C++ changes also require rebuilding
4. Python-only changes don't require rebuilding
5. Always run `pixi run pre-commit` before committing to check linting and formatting
6. Run `cargo clippy --all-targets --all-features` before opening a PR when Rust code changes

## Architecture Overview

SASKTRAN2 follows a layered architecture with three main language components:

### Core C++ Engine (`cpp/`)
- **Engine** (`engine.h`): Main radiative transfer calculation engine, templated for different Stokes vector sizes
- **Atmosphere** (`atmosphere/`): Atmospheric property storage and management
- **Geometry** (`geometry/`): 1D spherical shell geometry handling
- **Raytracing** (`raytracing/`): Ray tracing through atmospheric layers
- **Source Integration** (`source_integrator.h`): Integration of various radiation sources
- **SKTRAN_DISCO** (`sktran_disco/`): Discrete ordinates method implementation

### Rust Bindings and Extensions (`rust/`)
- **sasktran2-rs**: Core Rust library with atmospheric modeling utilities
- **sasktran2-py-ext**: Python extension module providing Python bindings
- **sasktran2-sys**: FFI bindings to C++ components

### Python Interface (`src/sasktran2/`)
- **Engine** (`engine.py`): Main user-facing calculation interface
- **Atmosphere** (`atmosphere.py`): Atmospheric state specification
- **Config** (`config.py`): Configuration and settings
- **Viewing Geometry** (`viewinggeo/`): Observation geometry definitions
- **Constituents** (`constituent/`): Atmospheric species (gases, aerosols, surfaces)
- **Optical** (`optical/`): Optical property databases and models
- **Database** (`database/`): Data storage and retrieval systems

## Key Design Patterns

1. **Template-based C++ Core**: The C++ engine is templated on NSTOKES (number of Stokes parameters) for performance
2. **Builder Pattern**: Configuration objects are built up progressively
3. **Rust for Performance**: Critical computational kernels implemented in Rust
4. **xarray Integration**: Results returned as labeled xarray DataArrays
5. **Derivative Mapping**: Full linearization support for atmospheric parameters

## Testing Structure

- **C++ tests**: `cpp/lib/tests/` - Uses custom test framework
- **Python tests**: `tests/` - Uses pytest
- **Validation tests**: Include comparisons with reference models and literature

## Key Files for Common Tasks

- **Adding new atmospheric constituents**: `src/sasktran2/constituent/`
- **Modifying optical properties**: `src/sasktran2/optical/` and `rust/sasktran2-rs/src/optical/`
- **Core radiative transfer**: `cpp/lib/engine/` and `cpp/include/sasktran2/engine.h`
- **Python API changes**: `src/sasktran2/` and `rust/sasktran2-py-ext/src/`
- **Build configuration**: `pyproject.toml`, `cpp/CMakeLists.txt`, `rust/*/Cargo.toml`

## Development Environment

The project requires:
- Rust nightly toolchain (specified in `rust-toolchain.toml`)
- C++17 compiler
- BLAS/LAPACK implementation (OpenBLAS, MKL, or Apple Accelerate)
- CMake 3.20+
- Python 3.10+

All dependencies are managed through Pixi for consistent development environments.
