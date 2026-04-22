from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np


def load_band_file(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a Schumann-Runge band file with columns WNO(cm^-1), SIGMA(cm^2)."""
    data = np.loadtxt(path, skiprows=2)
    wavenumber_cminv = data[:, 0]
    sigma_cm2 = data[:, 1]
    wavelength_nm = 1.0e7 / wavenumber_cminv
    return wavelength_nm, sigma_cm2


def load_continuum_file(path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    """Load a continuum file with columns wavelength(nm), sigma(cm^2)."""
    lines = path.read_text().splitlines()
    temp_k = infer_temperature_from_header(lines)

    values: list[tuple[float, float]] = []
    for line in lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            wav = float(parts[0])
            xs = float(parts[1])
        except ValueError:
            continue
        values.append((wav, xs))

    if not values:
        msg = f"No numeric continuum data found in {path}"
        raise ValueError(msg)

    arr = np.array(values)
    return arr[:, 0], arr[:, 1], temp_k


def infer_temperature_from_header(lines: list[str]) -> float:
    for line in lines[:5]:
        m = re.search(r"at\s+([0-9]+)K", line, flags=re.IGNORECASE)
        if m:
            return float(m.group(1))
    return 300.0


def dedupe_sorted_grid(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sort by x and merge duplicate x values by averaging y."""
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]

    unique_x, inv = np.unique(x_sorted, return_inverse=True)
    sums = np.zeros_like(unique_x)
    counts = np.zeros_like(unique_x)

    np.add.at(sums, inv, y_sorted)
    np.add.at(counts, inv, 1.0)

    return unique_x, sums / counts


def average_to_0p1nm(
    band_wavelength_nm: np.ndarray,
    band_xs_cm2: np.ndarray,
    cont_wavelength_nm: np.ndarray,
    cont_xs_cm2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Average total (band + continuum) cross section to 0.1 nm bins."""
    band_wavelength_nm, band_xs_cm2 = dedupe_sorted_grid(band_wavelength_nm, band_xs_cm2)
    cont_wavelength_nm, cont_xs_cm2 = dedupe_sorted_grid(cont_wavelength_nm, cont_xs_cm2)

    min_w = min(np.min(band_wavelength_nm), np.min(cont_wavelength_nm))
    max_w = max(np.max(band_wavelength_nm), np.max(cont_wavelength_nm))

    start = math.ceil(min_w * 10.0) / 10.0
    end = math.floor(max_w * 10.0) / 10.0
    centers = np.arange(start, end + 1e-12, 0.1)
    if len(centers) == 0:
        msg = "No 0.1 nm bins were generated from input data ranges"
        raise ValueError(msg)
    edges = np.concatenate(([centers[0] - 0.05], centers + 0.05))

    # Oversample once, then average by wavelength bins for stable results.
    hires = np.arange(edges[0], edges[-1] + 5e-4, 1e-3)

    band_interp = np.interp(hires, band_wavelength_nm, band_xs_cm2, left=0.0, right=0.0)
    cont_interp = np.interp(hires, cont_wavelength_nm, cont_xs_cm2, left=0.0, right=0.0)
    total_interp = band_interp + cont_interp

    bin_idx = np.digitize(hires, edges) - 1

    binned_sum = np.zeros(len(centers))
    binned_count = np.zeros(len(centers))
    valid = (bin_idx >= 0) & (bin_idx < len(centers))
    np.add.at(binned_sum, bin_idx[valid], total_interp[valid])
    np.add.at(binned_count, bin_idx[valid], 1.0)

    # Should not happen with the chosen oversampling, but keep this safe.
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
    """Write one fixed-width xsc block that the Rust fixed-width loader can read."""
    # The Rust loader assumes a linearly spaced wavenumber axis in the header.
    wavenumber_target = 1.0e7 / wavelength_nm
    wavenumber_target = np.sort(wavenumber_target)

    num_points = len(wavenumber_target)
    wv_start = float(wavenumber_target[0])
    wv_end = float(wavenumber_target[-1])

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
            "Generate a Schumann-Runge O2 fixed-width xsc file with two temperatures "
            "(90 K and 300 K) and 0.1 nm-averaged cross sections."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("tools/databases/sch/data"),
        help="Folder containing o2wb*x0.xsc and srcnt*.xsc input files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cross-sections/xs/O2SCHRUNG"),
        help="Output fixed-width xsc file path",
    )
    parser.add_argument(
        "--pressure-pa",
        type=float,
        default=101325.0,
        help="Pressure metadata written to each header",
    )
    parser.add_argument(
        "--molecule-name",
        type=str,
        default="O2-SR",
        help="Molecule name written into each header",
    )

    args = parser.parse_args()

    band_files = sorted(args.data_dir.glob("o2wb*x0.xsc"))
    if len(band_files) == 0:
        msg = f"No Schumann-Runge band files found in {args.data_dir}"
        raise FileNotFoundError(msg)

    band_w = []
    band_xs = []
    for band_file in band_files:
        w, xs = load_band_file(band_file)
        band_w.append(w)
        band_xs.append(xs)

    all_band_w = np.concatenate(band_w)
    all_band_xs = np.concatenate(band_xs)
    all_band_w, all_band_xs = dedupe_sorted_grid(all_band_w, all_band_xs)

    cont90_w, cont90_xs, temp90 = load_continuum_file(args.data_dir / "srcnt90.xsc")
    contrt_w, contrt_xs, temprt = load_continuum_file(args.data_dir / "srcntrt.xsc")

    w90, xs90 = average_to_0p1nm(all_band_w, all_band_xs, cont90_w, cont90_xs)
    wrt, xsrt = average_to_0p1nm(all_band_w, all_band_xs, contrt_w, contrt_xs)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as f:
        write_hitran_like_block(
            f,
            args.molecule_name,
            w90,
            xs90,
            temperature_k=temp90,
            pressure_pa=args.pressure_pa,
        )
        write_hitran_like_block(
            f,
            args.molecule_name,
            wrt,
            xsrt,
            temperature_k=300.0 if temprt == 295.0 else temprt,
            pressure_pa=args.pressure_pa,
        )

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
