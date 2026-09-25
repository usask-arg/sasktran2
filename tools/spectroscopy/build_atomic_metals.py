"""Acquire an auditable NIST ASD catalog for cold atmospheric metal spectroscopy.

Run with the project Python (numpy, xarray, netCDF4 installed). Requires curl.
Snapshots, request URLs and checksums are retained beside the NetCDF files.
The default search is deliberately a candidate screen, not a detectability claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import logging
import math
import subprocess
import time
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import xarray as xr

# Representative masses for unresolved terrestrial isotope mixtures. Interval
# elements use conventional values; ion electron-mass corrections are negligible
# compared with omitted isotope shifts. CIAAW: https://ciaaw.org/atomic-weights.htm
MASSES = {
    "Li": 6.94,
    "Be": 9.0121831,
    "Na": 22.98976928,
    "Mg": 24.305,
    "Al": 26.9815384,
    "Si": 28.085,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955907,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938043,
    "Fe": 55.845,
    "Co": 58.933194,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.630,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.905838,
    "Zr": 91.222,
    "Nb": 92.90637,
    "Mo": 95.95,
    "Ag": 107.8682,
    "Cd": 112.414,
    "In": 114.818,
    "Sn": 118.710,
    "Cs": 132.90545196,
    "Ba": 137.327,
    "Hf": 178.486,
    "Ta": 180.94788,
    "W": 183.84,
    "Pb": 207.2,
    "Bi": 208.98040,
}
BASE = "https://physics.nist.gov/cgi-bin/ASD/"
LOWER_LIMIT = 5000.0
UPPER_LIMIT = 50000.0
PRINCIPAL_WINDOWS = {
    "Li_I": (670, 672),
    "Na_I": (588, 591),
    "K_I": (766, 771),
    "Rb_I": (779, 796),
    "Cs_I": (850, 895),
    "Mg_II": (279, 281),
    "Be_II": (313, 314),
    "Mg_I": (285, 286),
}
RECONCILIATIONS_PATH = Path(__file__).with_name("atomic_strength_reconciliations.json")
STRENGTH_RECONCILIATIONS = {
    (entry["lower_level_id"], entry["upper_level_id"]): entry
    for entry in json.loads(RECONCILIATIONS_PATH.read_text())
}


def screening_category(species):
    """Keep atomic-data availability separate from ordinary-layer detectability."""
    if species in {"Na_I", "K_I"}:
        return "established_OSIRIS_atomic_layer"
    if species in {"Li_I", "Mg_I", "Mg_II", "Ca_I", "Ca_II", "Fe_I", "Fe_II", "Ni_I"}:
        return "natural_layer_observed_other_instruments"
    if species.split("_")[0] in {"Al", "Si", "Cr", "Mn", "Ti", "Co", "Cu", "Zn", "V"}:
        return "natural_meteoric_ablation_candidate_no_OSIRIS_detection_claim"
    if species.split("_")[0] in {"Rb", "Sr", "Ba"}:
        return "trace_meteoric_candidate_unestablished_detectability"
    return "atomic_database_screen_only_no_ambient_detectability_claim"


def clean(value):
    """Strip ASD spreadsheet wrappers, but retain spectroscopic annotations."""
    value = (value or "").strip()
    return value[2:-1] if value.startswith('="') and value.endswith('"') else value


def number(value):
    """Accept explicit numeric values (including bracketed estimated energies)."""
    try:
        return float(Fraction(clean(value).strip("[]()")))
    except (ValueError, ZeroDivisionError):
        return math.nan


def rows(path):
    raw = path.read_text(encoding="utf-8-sig")
    if "<html" in raw.lower():
        if "No lines" in raw or "No energy levels" in raw:
            return []
        msg = f"NIST returned an HTML error instead of CSV: {path}"
        raise ValueError(msg)
    return [
        {key: clean(value) for key, value in row.items() if key}
        for row in csv.DictReader(io.StringIO(raw))
    ]


def acquire(url, path, refresh=False):
    if refresh or not path.exists():
        partial = path.with_suffix(path.suffix + ".part")
        try:
            subprocess.run(
                [
                    "curl",
                    "--fail",
                    "--location",
                    "--silent",
                    "--show-error",
                    "--retry",
                    "2",
                    "--max-time",
                    "120",
                    "--user-agent",
                    "Mozilla/5.0",
                    url,
                    "--output",
                    str(partial),
                ],
                check=True,
            )
            rows(partial)  # Validate the response before replacing a cached file.
            partial.replace(path)
        finally:
            partial.unlink(missing_ok=True)
        time.sleep(0.25)
    return {
        "url": url,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "file": str(path.name),
        "retrieved_utc": datetime.fromtimestamp(
            path.stat().st_mtime, timezone.utc
        ).isoformat(),
    }


def wigner6j(a, b, c, d, e, f):
    """Racah factorial formula; arguments here are small integer/half integers."""

    def fact(x):
        if x < 0 or abs(x - round(x)) > 1e-8:
            return math.inf
        return math.factorial(round(x))

    def delta(x, y, z):
        if min(x + y - z, x - y + z, -x + y + z) < 0:
            return 0.0
        return math.sqrt(
            fact(x + y - z) * fact(x - y + z) * fact(-x + y + z) / fact(x + y + z + 1)
        )

    prefactor = delta(a, b, c) * delta(a, e, f) * delta(d, b, f) * delta(d, e, c)
    if not prefactor:
        return 0.0
    aa = [a + b + c, a + e + f, d + b + f, d + e + c]
    bb = [a + b + d + e, a + c + d + f, b + c + e + f]
    return prefactor * sum(
        (-1) ** z
        * fact(z + 1)
        / (math.prod(fact(z - v) for v in aa) * math.prod(fact(v - z) for v in bb))
        for z in range(math.ceil(max(aa)), math.floor(min(bb)) + 1)
    )


def electronic_w2(jl, ju):
    return 3 * (2 * ju + 1) * wigner6j(1, 1, 2, ju, ju, jl) ** 2


def level_key(row, side):
    # ASD's stable level IDs avoid rounding-related errors in matching branches.
    value = row.get(f"ID_{side}", "")
    if value:
        return value
    return (
        row.get(f"conf_{side}"),
        row.get(f"term_{side}"),
        row.get(f"J_{side}"),
        row.get("Ei(cm-1)" if side == "i" else "Ek(cm-1)"),
    )


def strength_reconciliation(row):
    """Return an audited resolution only while its exact source values match."""
    entry = STRENGTH_RECONCILIATIONS.get((row.get("ID_i"), row.get("ID_k")))
    if entry is None:
        return None
    for field, expected in (
        ("Aki(s^-1)", "nist_a_s"),
        ("fik", "nist_f"),
        ("J_i", "lower_j"),
        ("J_k", "upper_j"),
    ):
        if not math.isclose(number(row.get(field)), entry[expected], rel_tol=1e-10):
            return None
    return entry


def inconsistent_decay_strength(row):
    """Quarantine a gross A/f conflict from lifetime sums as well as opacity."""
    if row.get("Type", "") not in ("", "E1"):
        return False
    a, f = number(row.get("Aki(s^-1)")), number(row.get("fik"))
    jl, ju = number(row.get("J_i")), number(row.get("J_k"))
    wl = number(row.get("ritz_wl_vac(nm)"))
    if not math.isfinite(wl):
        wl = number(row.get("obs_wl_vac(nm)"))
    if not all(math.isfinite(x) for x in (a, f, jl, ju, wl)) or min(a, f, wl) <= 0:
        return False
    expected_f = 1.49919e-16 * a * (2 * ju + 1) / (2 * jl + 1) * (10 * wl) ** 2
    return abs(expected_f / f - 1) > 0.25 and strength_reconciliation(row) is None


def build_species(species, folder, refresh=False):
    raw_dir = folder / "raw" / "nist"
    raw_dir.mkdir(parents=True, exist_ok=True)
    line_params = {
        "spectra": species.replace("_", " "),
        "unit": 1,
        "format": 2,
        "line_out": 0,
        "en_unit": 0,
        "output": 0,
        "bibrefs": 1,
        "show_obs_wl": 1,
        "show_calc_wl": 1,
        "show_av": 3,
        "A_out": 0,
        "f_out": "on",
        "allowed_out": 1,
        "forbid_out": 1,
        "conf_out": "on",
        "term_out": "on",
        "enrg_out": "on",
        "J_out": "on",
        "max_upp_enrg": UPPER_LIMIT,
        "no_spaces": "on",
        "remove_js": "on",
        "ids_out": 1,
    }
    levels_params = {
        "spectrum": species.replace("_", " "),
        "units": 0,
        "format": 2,
        "output": 0,
        "conf_out": "on",
        "term_out": "on",
        "level_out": "on",
        "j_out": "on",
        "biblio": "on",
        "page_size": 15,
        "ids_out": 1,
    }
    provenance = {}
    for kind, endpoint, params in [
        ("lines", "lines1.pl", line_params),
        ("levels", "energy1.pl", levels_params),
    ]:
        path = raw_dir / f"{species}_{kind}.csv"
        provenance[kind] = acquire(
            BASE + endpoint + "?" + urlencode(params), path, refresh
        )
    lines = rows(raw_dir / f"{species}_lines.csv")
    levels = rows(raw_dir / f"{species}_levels.csv")
    if lines and "Ei(cm-1)" not in lines[0]:
        msg = f"Unexpected line columns: {list(lines[0])}"
        raise ValueError(msg)

    # A measured decay rate is not equivalent to a complete upper-state lifetime.
    # Sum distinct branches across ALL wavelengths, including forbidden branches.
    branches = defaultdict(dict)
    missing_a = defaultdict(int)
    inconsistent_decay_count = 0
    for row in lines:
        key = level_key(row, "k")
        if not key:
            continue
        a = number(row.get("Aki(s^-1)"))
        if inconsistent_decay_strength(row):
            inconsistent_decay_count += 1
            missing_a[key] += 1
            continue
        if not math.isfinite(a):
            missing_a[key] += 1
            continue
        pair = (level_key(row, "i"), row.get("Type", ""))
        if pair in branches[key] and not np.isclose(branches[key][pair], a, rtol=0.05):
            msg = f"Conflicting A values for {species} {key} -> {pair}"
            raise ValueError(msg)
        branches[key][pair] = a
    total_a = {key: sum(values.values()) for key, values in branches.items()}

    selected = []
    seen = set()
    rejected = defaultdict(int)
    inconsistent_strengths = []
    reconciled_strengths = []
    for row in lines:
        el, eu = number(row.get("Ei(cm-1)")), number(row.get("Ek(cm-1)"))
        wl = number(row.get("ritz_wl_vac(nm)"))
        if not math.isfinite(wl):
            wl = number(row.get("obs_wl_vac(nm)"))
        if not (274 <= wl <= 810 and 0 <= el <= LOWER_LIMIT):
            continue
        if row.get("Type", "") not in ("", "E1"):
            rejected["non_electric_dipole"] += 1
            continue
        a, f = number(row.get("Aki(s^-1)")), number(row.get("fik"))
        jl, ju = number(row.get("J_i")), number(row.get("J_k"))
        if not all(math.isfinite(x) for x in (el, eu, a, jl, ju)) or a <= 0:
            rejected["missing_strength_or_resolved_J"] += 1
            continue
        if abs(ju - jl) > 1 or (ju == 0 and jl == 0):
            rejected["invalid_E1_J"] += 1
            continue
        # ASD f is derived from A, not an independent measurement. Derive f from
        # canonical A and resolved J; flag discrepant displayed values and retain
        # gross conflicts only when individually reconciled against sources.
        expected_f = 1.49919e-16 * a * (2 * ju + 1) / (2 * jl + 1) * (10 * wl) ** 2
        strength_error = (
            abs(expected_f / f - 1) if math.isfinite(f) and f > 0 else math.nan
        )
        reconciliation = strength_reconciliation(row)
        if strength_error > 0.25 and reconciliation is None:
            rejected["inconsistent_A_and_f"] += 1
            inconsistent_strengths.append(
                {
                    "wavelength_nm": wl,
                    "nist_f": f,
                    "nist_a_s": a,
                    "f_implied_by_a": expected_f,
                    "tp_ref": row.get("tp_ref", ""),
                }
            )
            continue
        if reconciliation is not None:
            reconciled_strengths.append(reconciliation)
        kl, ku = level_key(row, "i"), level_key(row, "k")
        key = (kl, ku)
        if key in seen:
            rejected["duplicate_level_pair"] += 1
            continue
        seen.add(key)
        # These principal D levels (and Mg singlet P) have only the ground term
        # below them with opposite parity: E1-closed. The wavelength restriction
        # prevents assuming closure for a higher member with missing ASD rates.
        window = PRINCIPAL_WINDOWS.get(species, (math.inf, math.inf))
        closed = (
            window[0] <= wl <= window[1]
            and el == 0
            and (row.get("term_i"), row.get("term_k")) in (("2S", "2P*"), ("1S", "1P*"))
            and all(
                number(other.get("Ei(cm-1)")) == 0
                for other in lines
                if level_key(other, "k") == ku and other.get("Type", "") in ("", "E1")
            )
        )
        selected.append(
            {
                "wavelength_nm": wl,
                "oscillator_strength": expected_f,
                "nist_oscillator_strength": f,
                "a_f_relative_consistency_error": strength_error,
                "strength_reconciled": int(reconciliation is not None),
                "strength_reconciliation_references": json.dumps(
                    reconciliation["references"] if reconciliation is not None else []
                ),
                "lower_energy_cminv": el,
                "upper_energy_cminv": eu,
                "lower_j": jl,
                "upper_j": ju,
                "lower_statistical_weight": 2 * jl + 1,
                "einstein_a_s": a,
                "upper_total_a_s": total_a.get(ku, a),
                "lower_total_a_s": total_a.get(kl, 0),
                "upper_decay_data_complete": int(closed),
                "upper_known_missing_a_count": missing_a.get(ku, 0),
                "branch_reference_fraction": a / total_a.get(ku, a),
                "w2_electronic": electronic_w2(jl, ju),
                "lower_configuration": row.get("conf_i", ""),
                "upper_configuration": row.get("conf_k", ""),
                "lower_term": row.get("term_i", ""),
                "upper_term": row.get("term_k", ""),
                "lower_level_id": str(kl),
                "upper_level_id": str(ku),
                "nist_accuracy": row.get("Acc", ""),
                "nist_tp_ref": row.get("tp_ref", ""),
                "nist_line_ref": row.get("line_ref", ""),
            }
        )

    states_by_id = {}
    for row in levels:
        energy = number(row.get("Level (cm-1)"))
        j = number(row.get("J"))
        if 0 <= energy <= LOWER_LIMIT and math.isfinite(j):
            # Equal printed energies and J can belong to physically distinct
            # levels; retain both contributions to the partition function.
            key = row.get("Level ID") or (
                row.get("Configuration"),
                row.get("Term"),
                row.get("J"),
                energy,
            )
            states_by_id[key] = (energy, 2 * j + 1)
    states = sorted(states_by_id.values())
    record = {
        "species": species,
        "screening_category": screening_category(species),
        "raw_transition_count": len(lines),
        "lines": len(selected),
        "states": len(states),
        "rejected": dict(rejected),
        "inconsistent_strengths": inconsistent_strengths,
        "reconciled_strengths": reconciled_strengths,
        "quarantined_decay_rate_count": inconsistent_decay_count,
        "reconciled_decay_rate_count": sum(
            strength_reconciliation(row) is not None for row in lines
        ),
        "source": provenance,
    }
    if not selected:
        record["status"] = "no_usable_E1_lines_in_screen"
        # Remove only this builder's obsolete derived file after a new screen
        # finds no usable lines; otherwise an old result could silently survive.
        (folder / f"{species}.nc").unlink(missing_ok=True)
        return record
    if not states or states[0][0] != 0:
        msg = f"No classified ground level for {species}"
        raise ValueError(msg)
    selected.sort(key=lambda line: line["wavelength_nm"])
    ds = xr.Dataset(
        {key: ("line", [row[key] for row in selected]) for key in selected[0]}
    )
    ds["energy_cminv"] = ("state", [row[0] for row in states])
    ds["statistical_weight"] = ("state", [row[1] for row in states])
    ds.attrs.update(
        schema_version=1,
        species=species,
        mass_amu=MASSES[species.split("_")[0]],
        source_name="NIST Atomic Spectra Database",
        source_version="5.12 (November 2024); version checked 2026-09-24",
        source_url="https://physics.nist.gov/asd",
        source_citation="Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024), NIST Atomic Spectra Database (version 5.12), DOI 10.18434/T4W30F",
        source_copyright="Copyright US Department of Commerce on behalf of the United States; all rights reserved. https://www.nist.gov/pml/atomic-spectra-database",
        redistribution_status="Local research cache; website redistribution permission not established. Review NIST SRD terms: https://shop.nist.gov/ccrz__CCPage?pageKey=SRDTC",
        retrieved_utc=provenance["lines"]["retrieved_utc"],
        screening_category=screening_category(species),
        wavelength_medium="vacuum",
        wavelength_min_nm=274.0,
        wavelength_max_nm=810.0,
        lower_energy_max_cminv=LOWER_LIMIT,
        decay_query_upper_energy_max_cminv=UPPER_LIMIT,
        temperature_max_k=1000.0,
        partition_caveat="Only classified states <=5000 cm-1 included; cold-atmosphere approximation, not a complete high-temperature partition function.",
        decay_caveat="upper_total_a_s and lower_total_a_s sum available distinct ASD branches at all wavelengths, including forbidden types and excluding quarantined A/f conflicts. Missing rates make them LOWER BOUNDS, not measured lifetimes. upper_decay_data_complete marks only audited E1-closed principal doublets and Mg I singlet P; negligible forbidden channels are not certified.",
        phase_caveat="Isolated electronic E1 transitions, unpolarized lower state, no magnetic fields, hyperfine/isotope structure or J-state interference. W2=3(2Ju+1){1 1 2; Ju Ju Jl}^2.",
        strength_consistency_caveat="Oscillator strength is derived from NIST A, resolved statistical weights and vacuum wavelength. ASD displayed f is retained separately; conflicts >25 percent are quarantined unless individually corroborated in atomic_strength_reconciliations.json. Raw source data are unchanged.",
        strength_reconciliations_sha256=hashlib.sha256(
            RECONCILIATIONS_PATH.read_bytes()
        ).hexdigest(),
        mass_source="https://ciaaw.org/atomic-weights.htm; representative terrestrial mixture, conventional masses for interval elements",
        source_requests_json=json.dumps(provenance, sort_keys=True),
    )
    for name in ("wavelength_nm",):
        ds[name].attrs["units"] = "nm"
    for name in ("lower_energy_cminv", "upper_energy_cminv", "energy_cminv"):
        ds[name].attrs["units"] = "cm-1"
    for name in ("einstein_a_s", "upper_total_a_s", "lower_total_a_s"):
        ds[name].attrs["units"] = "s-1"
    dest = folder / f"{species}.nc"
    partial = dest.with_suffix(".nc.part")
    try:
        ds.to_netcdf(partial)
        partial.replace(dest)
    finally:
        partial.unlink(missing_ok=True)
    record.update(
        status="available",
        file=dest.name,
        sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
        max_oscillator_strength=max(row["oscillator_strength"] for row in selected),
        principal_wavelength_nm=max(
            selected, key=lambda row: row["oscillator_strength"]
        )["wavelength_nm"],
        strongest_lte_200k_wavelength_nm=max(
            selected,
            key=lambda row: (
                row["oscillator_strength"]
                * row["lower_statistical_weight"]
                * math.exp(-row["lower_energy_cminv"] * 1.438776877 / 200)
            ),
        )["wavelength_nm"],
    )
    return record


def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--species", nargs="+", help="ASD keys, e.g. Na_I Mg_II")
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / RECONCILIATIONS_PATH.name).write_bytes(
        RECONCILIATIONS_PATH.read_bytes()
    )
    species = args.species or [
        f"{symbol}_{stage}" for symbol in MASSES for stage in ("I", "II")
    ]
    catalog_path = args.output / "catalog.json"
    existing = json.loads(catalog_path.read_text()) if catalog_path.exists() else []
    catalog = {record["species"]: record for record in existing}
    failed = False
    for name in species:
        try:
            record = build_species(name, args.output, args.refresh)
        except (ValueError, subprocess.CalledProcessError) as error:
            record = {"species": name, "status": "error", "error": str(error)}
        catalog[name] = record
        failed |= record["status"] == "error"
        logging.info(
            "%s: %s (%s lines)", name, record["status"], record.get("lines", 0)
        )
        catalog_path.write_text(json.dumps(list(catalog.values()), indent=2) + "\n")
    if failed:
        msg = "Some NIST acquisitions failed; consult catalog.json"
        raise SystemExit(msg)


if __name__ == "__main__":
    main()
