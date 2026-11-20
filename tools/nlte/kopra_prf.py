"""Readers for KOPRA ``.prf`` profile files, as used by the CAIRT ERS archive.

The ERS archive (Zenodo 8256025) contains GRANADA non-LTE population ratios,
p/T, VMR and photolysis-rate profiles in KOPRA's text format. These readers
turn them into xarray Datasets for validating ``sasktran2-nlte``.

A ``.prf`` file is a sequence of blocks, each opened by a line starting with
``$``. The first numeric block holds the number of levels, the second the
altitude grid in km. Every later block with at least that many numbers is a
profile: the last ``n_levels`` numbers are the values, and any tokens before
them (species id, state index, parameter number and name) identify it.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import xarray as xr


@dataclass
class Profile:
    header: str
    """Text after the ``$`` that opened the block."""
    ids: list[str]
    """Tokens between the header and the values, e.g. ``["11", "1"]``."""
    values: np.ndarray


@dataclass
class PrfFile:
    altitude_km: np.ndarray
    profiles: list[Profile] = field(default_factory=list)


_FORTRAN_3_DIGIT_EXPONENT = re.compile(r"^([-+]?[0-9.]+)([-+][0-9]{3})$")


def _to_float(token: str) -> float:
    """Parse a number, including Fortran E-format overflow such as ``2.487055+100``.

    Fortran drops the ``E`` when a three-digit exponent does not fit the field;
    the ERS ratio files use this for very large O2(X, v>=29) ratios.
    """
    match = _FORTRAN_3_DIGIT_EXPONENT.match(token)
    if match:
        return float(f"{match[1]}e{match[2]}")
    return float(token)


def _is_number(token: str) -> bool:
    try:
        _to_float(token)
    except ValueError:
        return False
    return True


def read_prf(path: str | Path) -> PrfFile:
    blocks: list[tuple[str, list[str]]] = []
    # Files written by KOPRA's Mixer put a text label such as "pressure [hPa]"
    # on the line before the "$"; it describes the next block, not this one.
    pending_label: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("$"):
            header = " ".join([*pending_label, line[1:].strip()]).strip()
            blocks.append((header, []))
            pending_label = []
        elif not any(map(_is_number, line.split())):
            if line.strip():
                pending_label.append(line.strip())
        elif blocks:
            blocks[-1][1].extend(line.split())

    numeric_blocks = [
        (header, tokens) for header, tokens in blocks if any(map(_is_number, tokens))
    ]
    if len(numeric_blocks) < 2:
        msg = f"{path}: expected level count and altitude blocks"
        raise ValueError(msg)

    n_levels = int(next(t for t in numeric_blocks[0][1] if _is_number(t)))
    altitude_tokens = [t for t in numeric_blocks[1][1] if _is_number(t)]
    if len(altitude_tokens) != n_levels:
        msg = f"{path}: {len(altitude_tokens)} altitudes for {n_levels} levels"
        raise ValueError(msg)

    prf = PrfFile(altitude_km=np.array([_to_float(t) for t in altitude_tokens]))
    for header, tokens in numeric_blocks[2:]:
        if sum(map(_is_number, tokens)) < n_levels:
            continue
        values = np.array([_to_float(t) for t in tokens[-n_levels:]])
        prf.profiles.append(
            Profile(header=header, ids=tokens[:-n_levels], values=values)
        )
    return prf


def _profile_name(header: str) -> str:
    """Species name from a header such as ``(ppmv) H2O (HITRAN)  source: ERS v7``."""
    name = header.removeprefix("(ppmv)").split("source:")[0]
    return name.split(" (")[0].strip()


def _source(header: str) -> str:
    return header.split("source:")[1].strip() if "source:" in header else ""


def read_pt(path: str | Path) -> xr.Dataset:
    """Pressure [hPa] and temperature [K]."""
    prf = read_prf(path)
    ds = xr.Dataset(coords={"altitude_km": prf.altitude_km})
    for profile in prf.profiles:
        header = profile.header.lower()
        if header.startswith("pressure"):
            ds["pressure_hpa"] = ("altitude_km", profile.values)
        elif header.startswith("temperature"):
            ds["temperature_k"] = ("altitude_km", profile.values)
    return ds


def read_vmr(path: str | Path) -> xr.Dataset:
    """Volume mixing ratios [ppmv] with dims (species, altitude_km)."""
    prf = read_prf(path)
    return xr.Dataset(
        {
            "vmr_ppmv": (
                ("species", "altitude_km"),
                np.array([p.values for p in prf.profiles]),
            )
        },
        coords={
            "altitude_km": prf.altitude_km,
            "species": [_profile_name(p.header) for p in prf.profiles],
            "species_id": ("species", [int(p.ids[0]) for p in prf.profiles]),
            "source": ("species", [_source(p.header) for p in prf.profiles]),
        },
    )


def read_npar(path: str | Path) -> xr.Dataset:
    """GRANADA parameter profiles (photolysis rates, O(1D) density) by name."""
    prf = read_prf(path)
    ds = xr.Dataset(coords={"altitude_km": prf.altitude_km})
    for profile in prf.profiles:
        number, name = profile.ids[-2], profile.ids[-1]
        ds[name] = ("altitude_km", profile.values)
        ds[name].attrs["parameter_number"] = int(number)
    return ds


def read_ratio(path: str | Path) -> xr.Dataset:
    """GRANADA population ratios r = n/n_LTE with dims (state, altitude_km).

    ``species_id`` is HITRAN molecule*10 + isotope and ``label`` is the state
    description from the file (quantum numbers and energy). GRANADA normalises
    n_LTE over only the modelled states, so convert before comparing with
    populations relative to a full partition function.
    """
    prf = read_prf(path)
    return xr.Dataset(
        {
            "ratio": (
                ("state", "altitude_km"),
                np.array([p.values for p in prf.profiles]),
            )
        },
        coords={
            "altitude_km": prf.altitude_km,
            "species_id": ("state", [int(p.ids[0]) for p in prf.profiles]),
            "state_index": ("state", [int(p.ids[1]) for p in prf.profiles]),
            "label": ("state", [p.header for p in prf.profiles]),
        },
    )
