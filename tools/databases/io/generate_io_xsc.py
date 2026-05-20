from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np


def parse_metadata(path: Path) -> tuple[float, float]:
    """Return (temperature_k, pressure_pa) from header metadata."""
    temperature_k = None
    pressure_pa = None

    with path.open() as f:
        for line in f:
            if line.startswith("#IO"):
                parts = line[1:].split()
                # Expected: IO temp pressure resol ...
                if len(parts) >= 3:
                    temperature_k = float(parts[1])
                    pressure_pa = float(parts[2])
                break

    if temperature_k is None:
        msg = f"Could not determine temperature from {path}"
        raise ValueError(msg)

    if pressure_pa is None:
        pressure_pa = 101325.0

    return temperature_k, pressure_pa


def load_io_file(path: Path) -> tuple[np.ndarray, np.ndarray, float, float]:
    temperature_k, pressure_pa = parse_metadata(path)
    data = np.loadtxt(path, comments="#")

    wavenumber_cminv = data[:, 0]
    xs_cm2 = data[:, 1]

    wavelength_nm = 1.0e7 / wavenumber_cminv
    order = np.argsort(wavelength_nm)

    return wavelength_nm[order], xs_cm2[order], temperature_k, pressure_pa


def dedupe_sorted_grid(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    unique_x, inv = np.unique(x, return_inverse=True)
    sums = np.zeros_like(unique_x)
    counts = np.zeros_like(unique_x)

    np.add.at(sums, inv, y)
    np.add.at(counts, inv, 1.0)

    return unique_x, sums / counts


def average_to_0p1nm(
    wavelength_nm: np.ndarray, xs_cm2: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    wavelength_nm, xs_cm2 = dedupe_sorted_grid(wavelength_nm, xs_cm2)

    start = math.ceil(np.min(wavelength_nm) * 10.0) / 10.0
    end = math.floor(np.max(wavelength_nm) * 10.0) / 10.0

    centers = np.arange(start, end + 1e-12, 0.1)
    if len(centers) == 0:
        msg = "No 0.1 nm bins generated"
        raise ValueError(msg)

    edges = np.concatenate(([centers[0] - 0.05], centers + 0.05))

    hires = np.arange(edges[0], edges[-1] + 5e-4, 1e-3)
    xs_interp = np.interp(hires, wavelength_nm, xs_cm2, left=0.0, right=0.0)

    idx = np.digitize(hires, edges) - 1
    valid = (idx >= 0) & (idx < len(centers))

    binned_sum = np.zeros(len(centers))
    binned_count = np.zeros(len(centers))

    np.add.at(binned_sum, idx[valid], xs_interp[valid])
    np.add.at(binned_count, idx[valid], 1.0)

    binned_count[binned_count == 0.0] = 1.0
    return centers, binned_sum / binned_count


def write_hitran_like_block(
    f,
    molecule_name: str,
    wavelength_nm: np.ndarray,
    sigma_cm2: np.ndarray,
    temperature_k: float,
    pressure_pa: float,
) -> None:
    wavenumber_target = np.sort(1.0e7 / wavelength_nm)

    num_points = len(wavenumber_target)
    wv_start = float(wavenumber_target[0])
    wv_end = float(wavenumber_target[-1])

    # Rust parser assumes linear spacing in wavenumber from start/end/num_points.
    wv_uniform = np.linspace(wv_start, wv_end, num_points)

    wl_ascending = np.sort(wavelength_nm)
    xs_ascending = sigma_cm2[np.argsort(wavelength_nm)]
    sigma_uniform_cm2 = np.interp(1.0e7 / wv_uniform, wl_ascending, xs_ascending)

    header = (
        f"{molecule_name:>8s} {wv_start:10.4f} {wv_end:10.4f}"
        f" {num_points:8d} {temperature_k:8.1f} {pressure_pa:8.2f}"
        "  0.000E+00  0  GENERATED\n"
    )
    f.write(header)

    for i in range(0, num_points, 10):
        row = sigma_uniform_cm2[i : i + 10]
        f.write("".join(f"{v:10.3E}" for v in row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate IO fixed-width xsc file from source .asc files, "
            "averaged to 0.1 nm."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("tools/databases/io/data"),
        help="Folder containing IO source .asc files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cross-sections/xs/IO"),
        help="Output fixed-width xsc file path",
    )
    parser.add_argument(
        "--molecule-name",
        type=str,
        default="IO",
        help="Molecule name in output header",
    )
    args = parser.parse_args()

    src_files = sorted(args.data_dir.glob("IO_Unknown_*.asc"))
    if len(src_files) == 0:
        msg = f"No IO files found in {args.data_dir}"
        raise FileNotFoundError(msg)

    entries = []
    for path in src_files:
        wavelength_nm, xs_cm2, temperature_k, pressure_pa = load_io_file(path)
        wv_red, xs_red = average_to_0p1nm(wavelength_nm, xs_cm2)
        entries.append((temperature_k, pressure_pa, wv_red, xs_red, path.name))

    entries.sort(key=lambda x: x[0])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as f:
        for temperature_k, pressure_pa, wavelength_nm, xs_cm2, _ in entries:
            write_hitran_like_block(
                f,
                args.molecule_name,
                wavelength_nm,
                xs_cm2,
                temperature_k,
                pressure_pa,
            )

    print(f"Wrote {args.output} with {len(entries)} temperature blocks")


if __name__ == "__main__":
    main()
