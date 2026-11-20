"""Build the O2 ultraviolet cross sections used by ``sasktran2.optical.O2UV``.

Run from the repository root with the project Python environment::

    python tools/spectroscopy/build_o2_uv.py --source-dir /path/to/sources --output o2_uv.nc

Missing source files are downloaded into ``--source-dir``. The product covers
130 nm to the 242.4 nm dissociation threshold, as ``xs`` [m^2] on
``temperature_k`` (130-500 K in 5 K steps) and ``wavenumber_cminv``:

Schumann-Runge continuum, below 175.44 nm
    The CfA table already in the sasktran2 standard database
    (``cross_sections/o2/O2SCHRUNG``), measured at 90 K and 300 K, interpolated
    linearly in temperature and held constant outside that range.

Schumann-Runge bands, 175.44-204.08 nm (49000.5-57000 cm^-1)
    Minschwaner, K., G. P. Anderson, L. A. Hall and K. Yoshino (1992),
    Polynomial coefficients for calculating O2 Schumann-Runge cross sections
    at 0.5 cm^-1 resolution, J. Geophys. Res., 97(D9), 10103-10108. Coefficient
    sets for 130-190, 190-280 and 280-500 K give, in 1e-20 cm^2,
    ``a0 X^2 + a1 X + a2`` with ``X = ((T - 100) / 10)^2``. The readme the source
    page links to is no longer online; this coefficient order is the one that
    matches the CfA table to within the temperature difference and makes the
    three sets agree to 0.2% where their ranges meet.

Herzberg continuum, 205-240 nm
    Yoshino, K., A. S. C. Cheung, J. R. Esmond, W. H. Parkinson, D. E.
    Freeman, S. L. Guberman, A. Jenouvrier, B. Coquart and M. F. Merienne
    (1988), Improved absorption cross sections of oxygen in the wavelength
    region 205-240 nm of the Herzberg continuum, Planet. Space Sci., 36,
    1469-1475, from the MPI-Mainz UV/VIS Spectral Atlas (Keller-Rudek et al.,
    2013, doi:10.5194/essd-5-365-2013). Temperature independent. The continuum
    is added under the bands at its 205 nm value down to 194 nm, where the
    bands dominate, and decreases linearly to zero at 242.4 nm. The
    pressure-induced part, which matters only in the lower stratosphere and
    troposphere, is not included.

Lyman-alpha is not included; see ``sasktran2.optical.O2LymanAlpha``.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr
from sasktran2.database.web import StandardDatabase

MINSCHWANER_URL = "https://kestrel.nmt.edu/~krm/RESEARCH/SRBANDS/"
#: Coefficient file and its upper temperature limit [K].
MINSCHWANER_FILES = {
    "fitcoef_cold.txt": 190.0,
    "fitcoef_mid.txt": 280.0,
    "fitcoef_hot.txt": np.inf,
}
YOSHINO_URL = (
    "https://www.uv-vis-spectral-atlas-mainz.org/uvvis_data/cross_sections/Oxygen/"
    "O2_Yoshino(1988)_298K_205-240nm(rec).txt"
)
YOSHINO_FILE = "O2_Yoshino(1988)_298K_205-240nm(rec).txt"

TEMPERATURES_K = np.arange(130.0, 500.1, 5.0)
SCHUMANN_RUNGE_BANDS_CMINV = (49000.5, 57000.0)
HERZBERG_UNDER_BANDS_NM = 194.0
O2_DISSOCIATION_THRESHOLD_NM = 242.4
CM2_TO_M2 = 1.0e-4


def download(url: str, destination: Path) -> Path:
    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(destination.suffix + ".part")
        logging.info("Downloading %s", url)
        request = urllib.request.Request(
            url, headers={"User-Agent": "SASKTRAN2 spectroscopy data preparation"}
        )
        with (
            urllib.request.urlopen(request, timeout=300) as response,
            temporary.open("wb") as output,
        ):
            while block := response.read(1024 * 1024):
                output.write(block)
        temporary.replace(destination)
    return destination


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_minschwaner(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Wavenumbers [cm^-1] and (a0, a1, a2) columns of one coefficient file."""
    data = np.loadtxt(path, skiprows=3)
    return data[:, 0], data[:, 1:4]


def minschwaner_cross_section_m2(
    coefficients: list[tuple[float, np.ndarray]], temperature_k: float
) -> np.ndarray:
    """Schumann-Runge band cross section [m^2] at one temperature.

    ``coefficients`` holds ``(upper temperature limit, (n, 3) coefficients)``
    in increasing temperature order.
    """
    coeffs = next(c for upper, c in coefficients if temperature_k < upper)
    x = ((temperature_k - 100.0) / 10.0) ** 2
    a0, a1, a2 = coeffs.T
    return np.maximum(a0 * x**2 + a1 * x + a2, 0.0) * 1.0e-20 * CM2_TO_M2


def read_yoshino(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Wavelengths [nm] and Herzberg continuum cross sections [m^2]."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1] * CM2_TO_M2


def herzberg_continuum_m2(
    wavenumber_cminv: np.ndarray, table_nm: np.ndarray, table_m2: np.ndarray
) -> np.ndarray:
    """Herzberg continuum on a wavenumber grid, extended as described above."""
    wavelength = 1.0e7 / wavenumber_cminv
    nodes_nm = np.concatenate(
        [[HERZBERG_UNDER_BANDS_NM], table_nm, [O2_DISSOCIATION_THRESHOLD_NM]]
    )
    nodes_m2 = np.concatenate([[table_m2[0]], table_m2, [0.0]])
    return np.interp(wavelength, nodes_nm, nodes_m2, left=0.0, right=0.0)


def read_cfa_xsc(path: Path) -> dict[float, tuple[np.ndarray, np.ndarray]]:
    """HITRAN-format cross-section blocks keyed by temperature [K]: (wavenumber, m^2)."""

    def is_header(line: str) -> bool:
        tokens = line.split()
        if not tokens:
            return False
        try:
            float(tokens[0])
        except ValueError:
            return True
        return False

    lines = Path(path).read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if is_header(line)]
    blocks = {}
    for k, start in enumerate(starts):
        header = lines[start].split()
        numin, numax, npts, temperature = (
            float(header[1]),
            float(header[2]),
            int(header[3]),
            float(header[4]),
        )
        end = starts[k + 1] if k + 1 < len(starts) else len(lines)
        values = np.array(" ".join(lines[start + 1 : end]).split(), dtype=float)[:npts]
        blocks[temperature] = (np.linspace(numin, numax, npts), values * CM2_TO_M2)
    return blocks


def cfa_continuum_m2(
    blocks: dict[float, tuple[np.ndarray, np.ndarray]], temperature_k: float
) -> tuple[np.ndarray, np.ndarray]:
    """Continuum points above the bands, interpolated linearly in temperature."""
    (t_lo, (nu, lo)), (t_hi, (_, hi)) = sorted(blocks.items())
    keep = nu > SCHUMANN_RUNGE_BANDS_CMINV[1]
    weight = np.clip((temperature_k - t_lo) / (t_hi - t_lo), 0.0, 1.0)
    return nu[keep], (1.0 - weight) * lo[keep] + weight * hi[keep]


def build(source_dir: Path, cfa_path: Path) -> xr.Dataset:
    bands = []
    for name, upper in MINSCHWANER_FILES.items():
        band_nu, coeffs = read_minschwaner(source_dir / name)
        bands.append((upper, coeffs))
    yoshino_nm, yoshino_m2 = read_yoshino(source_dir / YOSHINO_FILE)
    cfa = read_cfa_xsc(cfa_path)

    herzberg_nu = 1.0e7 / np.concatenate([yoshino_nm, [O2_DISSOCIATION_THRESHOLD_NM]])
    continuum_nu, _ = cfa_continuum_m2(cfa, TEMPERATURES_K[0])
    wavenumber = np.concatenate([np.sort(herzberg_nu), band_nu, continuum_nu])

    herzberg = herzberg_continuum_m2(wavenumber, yoshino_nm, yoshino_m2)
    in_bands = (wavenumber >= SCHUMANN_RUNGE_BANDS_CMINV[0]) & (
        wavenumber <= SCHUMANN_RUNGE_BANDS_CMINV[1]
    )
    above_bands = wavenumber > SCHUMANN_RUNGE_BANDS_CMINV[1]

    xs = np.zeros((TEMPERATURES_K.size, wavenumber.size))
    for i, temperature in enumerate(TEMPERATURES_K):
        xs[i] = herzberg
        xs[i, in_bands] += minschwaner_cross_section_m2(bands, temperature)
        xs[i, above_bands] = cfa_continuum_m2(cfa, temperature)[1]

    sources = {
        name: sha256(source_dir / name) for name in [*MINSCHWANER_FILES, YOSHINO_FILE]
    }
    sources[cfa_path.name] = sha256(cfa_path)
    return xr.Dataset(
        {"xs": (("temperature_k", "wavenumber_cminv"), xs)},
        coords={"temperature_k": TEMPERATURES_K, "wavenumber_cminv": wavenumber},
        attrs={
            "title": "O2 ultraviolet absorption cross sections, 130-242.4 nm",
            "units": "m^2",
            "schumann_runge_continuum": "CfA O2SCHRUNG (sasktran2 standard database), 90 K and 300 K",
            "schumann_runge_bands": "Minschwaner et al. (1992), J. Geophys. Res., 97, 10103",
            "herzberg_continuum": "Yoshino et al. (1988), Planet. Space Sci., 36, 1469, via the MPI-Mainz UV/VIS Spectral Atlas",
            "notes": (
                "Herzberg continuum held at its 205 nm value under the bands down to "
                "194 nm, linear to zero at 242.4 nm; pressure-induced absorption and "
                "Lyman-alpha not included"
            ),
            "source_sha256": "; ".join(f"{k}={v}" for k, v in sources.items()),
            "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "builder": "tools/spectroscopy/build_o2_uv.py",
        },
    )


def write_netcdf(dataset: xr.Dataset, destination: Path) -> None:
    partial = destination.with_suffix(destination.suffix + ".part")
    try:
        dataset.to_netcdf(
            partial,
            engine="netcdf4",
            encoding={"xs": {"zlib": True, "complevel": 4}},
        )
        partial.replace(destination)
    finally:
        partial.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    for name in MINSCHWANER_FILES:
        download(MINSCHWANER_URL + name, args.source_dir / name)
    download(YOSHINO_URL, args.source_dir / YOSHINO_FILE)

    cfa_path = StandardDatabase().path("cross_sections/o2/O2SCHRUNG")
    write_netcdf(build(args.source_dir, cfa_path), args.output)
    logging.info("Wrote %s", args.output)


if __name__ == "__main__":
    main()
