"""Build the vacuum-ultraviolet tables that extend sasktran2.photolysis to 80 nm.

    python tools/spectroscopy/build_vuv.py --source-dir /Volumes/T9/data/vuv_solar_o2 --output-dir out/

Products, for the sasktran2 standard database:

``cross_sections/vuv/o2_vuv.nc``
    O2 photoabsorption, 80-130.1 nm, from the Leiden photodissociation
    database (Heays, Bosman and van Dishoeck 2017, A&A 602, A105; Holland et
    al. 1993 below 103 nm, Watanabe and Marmo 1956, Wu 2005, Ogawa and Ogawa
    1975 and Dose 1975 above), about 295 K. Zero within 0.05 nm of
    Lyman-alpha, which ``O2LymanAlpha`` and the Chabrillat and Kockarts
    parameterisation cover, and above 130.1 nm, where ``O2UV`` starts.
``cross_sections/vuv/n2_vuv.nc``
    N2 photoabsorption, 80-100 nm, Leiden database (Heays et al. 2014 and
    Shaw et al. 1992; 100 K, b = 1 km/s), averaged into 0.005 nm bins. The
    averaging under-represents saturation in the strongest bands.
``cross_sections/vuv/o_vuv.nc``
    Atomic O photoabsorption from the Leiden database, 80-100 nm.
``solar/solar_irradiance_whi2008_hsrs_composite.nc``
    TSIS-1 HSRS from 116 nm, and below it the WHI 2008 quiet-Sun reference
    spectrum (Woods et al. 2009, GRL 36, L01101; 10-16 April 2008,
    F10.7 = 68.9), 0.1 nm bins. They agree to 0.5% over 116-120 nm.
"""

from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

CM2_TO_M2 = 1.0e-4
LYMAN_ALPHA_NM = 121.567
HSRS_SPLICE_NM = 116.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_leiden(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Wavelength [nm] and photoabsorption cross section [cm^2]."""
    data = np.loadtxt(path, comments="#", usecols=(0, 1))
    return data[:, 0], data[:, 1]


def absorber(wavelength_nm, xs_cm2, title, source, attrs=None) -> xr.Dataset:
    order = np.argsort(1.0e7 / wavelength_nm)
    ds = xr.Dataset(
        {"xs": (("wavenumber_cminv",), xs_cm2[order] * CM2_TO_M2)},
        coords={"wavenumber_cminv": (1.0e7 / wavelength_nm)[order]},
        attrs={
            "title": title,
            "units": "m^2",
            "source_sha256": f"{source.name}={sha256(source)}",
            "history": f"{datetime.now(timezone.utc).isoformat()} tools/spectroscopy/build_vuv.py",
            **(attrs or {}),
        },
    )
    ds["xs"].attrs["units"] = "m^2"
    return ds


def binned(wavelength, values, lo, hi, step):
    edges = np.arange(lo, hi + step / 2, step)
    index = np.digitize(wavelength, edges) - 1
    keep = (index >= 0) & (index < edges.size - 1)
    total = np.bincount(index[keep], weights=values[keep], minlength=edges.size - 1)
    count = np.bincount(index[keep], minlength=edges.size - 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    ok = count > 0
    return centres[ok], total[ok] / count[ok]


def main(source_dir: Path, output_dir: Path) -> None:
    leiden = source_dir / "leiden_photodissociation"
    (output_dir / "cross_sections" / "vuv").mkdir(parents=True, exist_ok=True)
    (output_dir / "solar").mkdir(parents=True, exist_ok=True)
    citation = (
        "Heays, A. N., A. D. Bosman and E. F. van Dishoeck (2017), A&A 602, A105, "
        "and the original sources listed in the Leiden database file header"
    )

    w, xs = read_leiden(leiden / "O2.txt")
    keep = (w >= 80.0) & (w <= 130.1)
    w, xs = w[keep], xs[keep]
    xs = np.where(np.abs(w - LYMAN_ALPHA_NM) < 0.05, 0.0, xs)
    absorber(
        w,
        xs,
        "O2 photoabsorption 80-130.1 nm",
        leiden / "O2.txt",
        {"references": citation, "notes": "zero within 0.05 nm of Lyman-alpha"},
    ).to_netcdf(output_dir / "cross_sections" / "vuv" / "o2_vuv.nc")

    w, xs = read_leiden(leiden / "N2.txt")
    w, xs = binned(w, xs, 80.0, 100.0, 0.005)
    absorber(
        w,
        xs,
        "N2 photoabsorption 80-100 nm, 0.005 nm bin averages",
        leiden / "N2.txt",
        {"references": citation},
    ).to_netcdf(output_dir / "cross_sections" / "vuv" / "n2_vuv.nc")

    w, xs = read_leiden(leiden / "O.txt")
    keep = (w >= 80.0) & (w <= 100.0)
    absorber(
        w[keep],
        xs[keep],
        "Atomic O photoabsorption 80-100 nm",
        leiden / "O.txt",
        {"references": citation},
    ).to_netcdf(output_dir / "cross_sections" / "vuv" / "o_vuv.nc")

    whi_path = source_dir / "whi2008" / "ref_solar_irradiance_whi-2008_ver2.dat"
    rows = []
    for line in whi_path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 4 and not line.lstrip().startswith(";"):
            try:
                rows.append((float(parts[0]), float(parts[3])))
            except ValueError:
                continue
    whi = np.array(rows)
    hsrs_path = Path(
        "~/Library/Application Support/sasktran2/database/solar/"
        "solar_irradiance_hsrs_2022_11_30_extended.nc"
    ).expanduser()
    hsrs = xr.load_dataset(hsrs_path)
    low = whi[(whi[:, 0] < HSRS_SPLICE_NM) & (whi[:, 0] >= 0.5)]
    high = hsrs.where(hsrs["wavelength"] >= HSRS_SPLICE_NM, drop=True)
    wavelength = np.concatenate([low[:, 0], high["wavelength"].to_numpy()])
    irradiance = np.concatenate([low[:, 1], high["irradiance"].to_numpy()])
    xr.Dataset(
        {"irradiance": (("wavelength",), irradiance)},
        coords={"wavelength": wavelength},
        attrs={
            "title": "Solar spectral irradiance: WHI 2008 quiet Sun below 116 nm, TSIS-1 HSRS above",
            "references": (
                "Woods, T. N., et al. (2009), Geophys. Res. Lett. 36, L01101 (WHI 2008, "
                "10-16 Apr 2008); Coddington, O. M., et al. (2023), Earth Space Sci. 10, "
                "e2022EA002637 (TSIS-1 HSRS)"
            ),
            "source_sha256": f"{whi_path.name}={sha256(whi_path)}, {hsrs_path.name}={sha256(hsrs_path)}",
            "history": f"{datetime.now(timezone.utc).isoformat()} tools/spectroscopy/build_vuv.py",
        },
    ).assign(
        irradiance=lambda d: d["irradiance"].assign_attrs(units="W m-2 nm-1")
    ).to_netcdf(
        output_dir / "solar" / "solar_irradiance_whi2008_hsrs_composite.nc"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    main(args.source_dir, args.output_dir)
