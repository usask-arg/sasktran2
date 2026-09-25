"""Prepare local molecular resonance line data from published ExoMol lists.

Run from the repository root with the project Python environment, for example::

    python tools/spectroscopy/build_exomol_metals.py --species AlO MgO CaO TiO

This downloads the original sources into the database, then stores the subset
whose *vacuum* wavelengths fall between 270 and 820 nm. Radiative decay sums
are calculated from the COMPLETE transition list, before the wavelength cut.
The full spectral subset is archived under ``full/``; the operational files
retain lower-state energies <= 5000 cm-1 and declare a 500 K temperature limit.
No intensity pruning is performed. The omitted LTE population fraction at the
maximum temperature is recorded, but is not an opacity error bound.

The result describes a single isotopologue at unit isotopologue abundance.
It is not an emission spectrum, an instrument-convolved opacity, or a claim
that resonance fluorescence is elastic. The source data and derived tables
are CC BY-SA 4.0; this license applies to the data, not the surrounding code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.constants import c, e, epsilon_0, h, k, m_e

DATASETS = {
    "AlO": {
        "isotopologue": "27Al-16O",
        "dataset": "ATP",
        "references": [
            "https://doi.org/10.1093/mnras/stv507",
            "https://doi.org/10.1093/mnras/stab2525",
        ],
    },
    "MgO": {
        "isotopologue": "24Mg-16O",
        "dataset": "LiTY",
        "references": [
            "https://doi.org/10.1093/mnras/stz912",
            "https://doi.org/10.1093/rasti/rzae037",
        ],
    },
    "CaO": {
        "isotopologue": "40Ca-16O",
        "dataset": "VBATHY",
        "references": ["https://doi.org/10.1093/mnras/stv2858"],
    },
    "TiO": {
        "isotopologue": "48Ti-16O",
        "dataset": "Toto",
        "references": [
            "https://doi.org/10.1093/mnras/stz1818",
            "https://doi.org/10.1093/rasti/rzae037",
        ],
    },
}

DATA_LICENSE = """ExoMol molecular data and SASKTRAN2-derived line subsets

These DATA are licensed under Creative Commons Attribution-ShareAlike 4.0
International (CC BY-SA 4.0):
https://creativecommons.org/licenses/by-sa/4.0/
Source licensing statement: https://exomol.com/data/licence/

Attribution: ExoMol (https://exomol.com/), and the authors of the dataset-specific
papers cited in each NetCDF references attribute and the preparation script.
This license concerns the source and derived DATA, not SASKTRAN2 source code.

Modifications: source states and transitions were joined, wavelengths recomputed
from empirical state energies, oscillator strengths and radiative decay sums
calculated, wavelengths selected, and an explicit mesospheric lower-state energy
cut applied to the operational tables. No isotope-abundance weighting is applied.

The molecular directory contains operational files (T <= 500 K by default), full/
contains untruncated spectral subsets, and raw/exomol contains original downloads.
Per-dataset provenance JSON files record original URLs and SHA-256 hashes. Retain
this attribution and license, and the provenance metadata, when redistributing.
"""


def download(url: str, destination: Path) -> Path:
    """Cache a complete source file, avoiding partially downloaded cache hits."""
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
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def write_netcdf(dataset: xr.Dataset, destination: Path) -> None:
    """Keep a valid cached product intact if serialization is interrupted."""
    partial = destination.with_suffix(destination.suffix + ".part")
    try:
        dataset.to_netcdf(
            partial,
            engine="netcdf4",
            encoding={key: {"zlib": True, "complevel": 4} for key in dataset.data_vars},
        )
        partial.replace(destination)
    finally:
        partial.unlink(missing_ok=True)


def build(
    species: str,
    output_root: Path,
    source_root: Path,
    limits: tuple[float, float],
    lower_energy_max: float = 5000,
    temperature_max: float = 500,
) -> Path:
    info = DATASETS[species]
    iso, name = info["isotopologue"], info["dataset"]
    base = f"https://exomol.com/db/{species}/{iso}/{name}/"
    stem = f"{iso}__{name}"
    cache = source_root / species
    definitions = download(base + stem + ".def", cache / (stem + ".def"))
    states_path = download(base + stem + ".states.bz2", cache / (stem + ".states.bz2"))
    partition_path = download(base + stem + ".pf", cache / (stem + ".pf"))
    transitions_path = download(
        base + stem + ".trans.bz2", cache / (stem + ".trans.bz2")
    )

    definition_lines = definitions.read_text().splitlines()
    mass_line = next(line for line in definition_lines if "Isotopologue mass" in line)
    mass = float(mass_line.split()[0])
    version_line = next(line for line in definition_lines if "Version number" in line)
    version = version_line.split()[0]
    maximum_wavenumber = float(
        next(line for line in definition_lines if "Maximum wavenumber" in line).split()[
            0
        ]
    )
    source_temperature_max = float(
        next(
            line
            for line in definition_lines
            if "Maximum temperature of linelist" in line
        ).split()[0]
    )
    if not 0 < temperature_max <= source_temperature_max:
        msg = f"Requested temperature limit exceeds {species} source validity ({source_temperature_max:g} K)"
        raise ValueError(msg)
    # The first four fields are fixed by the ExoMol format, independent of
    # optional uncertainty/lifetime columns and molecule-specific quantum labels.
    states = pd.read_csv(
        states_path,
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["id", "energy", "g", "j"],
    )
    size = int(states["id"].max()) + 1
    energy = np.full(size, np.nan)
    degeneracy = np.full(size, np.nan)
    angular_momentum = np.full(size, np.nan)
    indices = states["id"].to_numpy(dtype=np.int64)
    energy[indices] = states["energy"]
    degeneracy[indices] = states["g"]
    angular_momentum[indices] = states["j"]
    decay_rate = np.zeros(size)
    retained = []
    count = 0
    # Ignore the optional transition-file frequency: updated empirical state
    # energies are authoritative and can differ from the original frequency.
    for frame in pd.read_csv(
        transitions_path,
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 2],
        names=["upper", "lower", "a"],
        chunksize=500_000,
    ):
        upper = frame["upper"].to_numpy(dtype=np.int64)
        lower = frame["lower"].to_numpy(dtype=np.int64)
        a = frame["a"].to_numpy(dtype=np.float64)
        if (
            np.any(upper <= 0)
            or np.any(lower <= 0)
            or np.any(upper >= size)
            or np.any(lower >= size)
        ):
            msg = "Transition references unknown state ID"
            raise ValueError(msg)
        if not np.all(np.isfinite(energy[upper])) or not np.all(
            np.isfinite(energy[lower])
        ):
            msg = "Transition references missing state ID"
            raise ValueError(msg)
        if not np.all(np.isfinite(a)) or np.any(a < 0):
            msg = "Invalid Einstein A coefficient"
            raise ValueError(msg)
        decay_rate += np.bincount(upper, weights=a, minlength=size)
        wavenumber = energy[upper] - energy[lower]
        good = (
            (wavenumber >= 1e7 / limits[1]) & (wavenumber <= 1e7 / limits[0]) & (a > 0)
        )
        if np.any(good):
            retained.append((upper[good], lower[good], a[good]))
        count += len(frame)
    if not retained:
        msg = f"No {species} lines in the requested wavelength interval"
        raise ValueError(msg)
    upper, lower, a = (np.concatenate([row[i] for row in retained]) for i in range(3))
    wavelength = 1e7 / (energy[upper] - energy[lower])
    order = np.argsort(wavelength)
    wavelength, upper, lower, a = (
        wavelength[order],
        upper[order],
        lower[order],
        a[order],
    )
    oscillator_strength = (
        m_e
        * epsilon_0
        * c
        / (2 * np.pi * e**2)
        * (wavelength * 1e-9) ** 2
        * a
        * degeneracy[upper]
        / degeneracy[lower]
    )
    pf = np.loadtxt(partition_path)
    ds = xr.Dataset(
        data_vars={
            "wavelength_nm": ("line", wavelength, {"units": "nm", "medium": "vacuum"}),
            "oscillator_strength": ("line", oscillator_strength, {"units": "1"}),
            "lower_energy_cminv": ("line", energy[lower], {"units": "cm-1"}),
            "lower_j": ("line", angular_momentum[lower]),
            "upper_j": ("line", angular_momentum[upper]),
            "lower_statistical_weight": ("line", degeneracy[lower]),
            "upper_statistical_weight": ("line", degeneracy[upper]),
            "einstein_a_s": ("line", a, {"units": "s-1"}),
            "upper_total_a_s": ("line", decay_rate[upper], {"units": "s-1"}),
            "lower_total_a_s": ("line", decay_rate[lower], {"units": "s-1"}),
            "lower_state_id": ("line", lower),
            "upper_state_id": ("line", upper),
            "energy_cminv": ("state", energy[indices], {"units": "cm-1"}),
            "statistical_weight": ("state", degeneracy[indices]),
            "state_j": ("state", angular_momentum[indices]),
            "partition_function": ("partition_temperature_k", pf[:, 1]),
        },
        coords={
            "partition_temperature_k": (
                "partition_temperature_k",
                pf[:, 0],
                {"units": "K"},
            ),
            "state": indices,
        },
        attrs={
            "schema_version": 1,
            "species": species,
            "isotopologue": iso,
            "dataset": name,
            "source_version": version,
            "mass_amu": mass,
            "source_url": f"https://exomol.com/data/molecules/{species}/{iso}/{name}/",
            "references": json.dumps(info["references"]),
            "license": "CC BY-SA 4.0",
            "license_url": "https://creativecommons.org/licenses/by-sa/4.0/",
            "source_license_url": "https://exomol.com/data/licence/",
            "prepared_utc": datetime.now(timezone.utc).isoformat(),
            "wavelength_min_nm": max(limits[0], 1e7 / maximum_wavenumber),
            "wavelength_max_nm": limits[1],
            "wavelength_medium": "vacuum",
            "source_maximum_wavenumber_cminv": maximum_wavenumber,
            "temperature_max_k": source_temperature_max,
            "source_transition_count": count,
            "isotopologue_abundance": 1.0,
            "population_model": "LTE; state weights and partition sum include nuclear-spin degeneracy",
            "spectroscopic_scope": "Published low-electronic-state model only; a wavelength within the file range is not proof that every physical molecular band is included. AlO ATP includes X, A, B states.",
            "radiative_sums": "All positive source Einstein A values, before wavelength filtering; no collisional or predissociation rates",
            "phase_model": "No molecular phase measurement supplied; isolated unpolarized E1 J-level model requires explicit approximation",
            "redistribution": "Only return to the identical lower state is elastic; other radiative branches change wavelength",
            "modifications": "Joined states and transitions; frequencies from updated energies; computed oscillator strengths and full radiative sums; selected vacuum wavelength interval",
        },
    )
    output_root.mkdir(parents=True, exist_ok=True)
    full_root = output_root / "full"
    full_root.mkdir(parents=True, exist_ok=True)
    output = full_root / f"{species}.nc"
    write_netcdf(ds, output)
    sources = [definitions, states_path, partition_path, transitions_path]
    provenance = {
        "species": species,
        "dataset": name,
        "isotopologue": iso,
        "source_transition_count": count,
        "retained_line_count": len(wavelength),
        "output_file": output.name,
        "output_sha256": sha256(output),
        "files": [
            {
                "filename": path.name,
                "url": base + path.name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in sources
        ],
        "license": "CC BY-SA 4.0",
        "license_url": "https://creativecommons.org/licenses/by-sa/4.0/",
    }
    (full_root / f"{species}.provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    logging.info(
        "%s: retained %s/%s lines -> %s (%.1f MB)",
        species,
        f"{len(wavelength):,}",
        f"{count:,}",
        output,
        output.stat().st_size / 1e6,
    )
    return reduce_for_mesosphere(output, output_root, lower_energy_max, temperature_max)


def reduce_for_mesosphere(
    source: Path, output_root: Path, lower_energy_max: float, temperature_max: float
) -> Path:
    """Write an explicitly temperature-limited working subset of a full file."""
    with xr.open_dataset(source, cache=False) as full:
        if not np.isfinite(lower_energy_max) or lower_energy_max <= 0:
            msg = "Lower-state energy limit must be finite and positive"
            raise ValueError(msg)
        if not 0 < temperature_max <= full.attrs["temperature_max_k"]:
            msg = "Requested temperature limit exceeds the full source dataset validity"
            raise ValueError(msg)
        selection = np.flatnonzero(full.lower_energy_cminv.values <= lower_energy_max)
        # Load contiguous variables before masking: random NetCDF selection of
        # hundreds of thousands of sorted wavelengths otherwise triggers many
        # tiny reads and is much slower than a sequential read.
        # Disable the input cache to avoid retaining the entire hot line list.
        reduced = xr.Dataset(
            {
                key: (
                    variable.dims,
                    (
                        variable.values[selection]
                        if variable.dims == ("line",)
                        else variable.values
                    ),
                    variable.attrs,
                )
                for key, variable in full.data_vars.items()
            },
            coords={
                key: variable.load().copy() for key, variable in full.coords.items()
            },
            attrs=full.attrs,
        )
        energies = full.energy_cminv.values
        populations = full.statistical_weight.values * np.exp(
            -h * c * 100 / k * energies / temperature_max
        )
        omitted_population = float(
            populations[energies > lower_energy_max].sum() / populations.sum()
        )
    reduced.attrs.update(
        temperature_max_k=temperature_max,
        lower_energy_max_cminv=lower_energy_max,
        omitted_lte_population_fraction_at_temperature_max=omitted_population,
        truncation_note="Operational mesospheric subset; excited-state lines omitted. The omitted population fraction is not an opacity error bound. Full spectral subset is archived under full/.",
        full_dataset_sha256=sha256(source),
    )
    output = output_root / source.name
    write_netcdf(reduced, output)
    provenance_path = source.with_suffix(".provenance.json")
    provenance = (
        json.loads(provenance_path.read_text()) if provenance_path.exists() else {}
    )
    provenance.update(
        full_dataset=str(source.relative_to(output_root)),
        lower_energy_max_cminv=lower_energy_max,
        temperature_max_k=temperature_max,
        omitted_lte_population_fraction_at_temperature_max=omitted_population,
        retained_line_count=len(selection),
        output_sha256=sha256(output),
    )
    output.with_suffix(".provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    logging.info(
        "%s: %s operational lines; omitted LTE population fraction %.3g at %.1f K",
        source.stem,
        f"{len(selection):,}",
        omitted_population,
        temperature_max,
    )
    return output


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", nargs="+", choices=DATASETS, default=["AlO"])
    parser.add_argument("--database-root", type=Path)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--wavelength-min", type=float, default=270.0)
    parser.add_argument("--wavelength-max", type=float, default=820.0)
    parser.add_argument("--lower-energy-max", type=float, default=5000.0)
    parser.add_argument("--temperature-max", type=float, default=500.0)
    args = parser.parse_args()
    if not 0 < args.wavelength_min < args.wavelength_max:
        parser.error("Require 0 < wavelength-min < wavelength-max")
    if args.lower_energy_max <= 0 or args.temperature_max <= 0:
        parser.error("Energy and temperature limits must be positive")
    if args.database_root is None:
        from sasktran2.appconfig import database_root  # noqa: PLC0415

        args.database_root = database_root()
    output = args.database_root / "spectroscopy/metals/molecular"
    output.mkdir(parents=True, exist_ok=True)
    (output / "LICENSE.txt").write_text(DATA_LICENSE)
    source = args.source_root or output / "raw/exomol"
    for species in args.species:
        build(
            species,
            output,
            source,
            (args.wavelength_min, args.wavelength_max),
            args.lower_energy_max,
            args.temperature_max,
        )


if __name__ == "__main__":
    main()
