"""Excited-state (non-LTE) populations from kinetic mechanisms.

A :class:`Mechanism` lists excited states, background species and the
processes that couple them (reactions, photolysis, radiative transitions).
:func:`solve` gives the steady-state state populations, process rates and
production/loss budgets on a column of levels. The kinetics run in the
``sasktran2-nlte`` Rust crate; the mechanism file format is described in
``docs/sphinx/source/developer/nlte_mechanism_format.md``.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import xarray as xr

from sasktran2._core_rust import PyMechanism

__all__ = [
    "PHOTOCHEMICAL_SPECIES",
    "Mechanism",
    "add_photochemical_species",
    "budget",
    "solve",
]


class Mechanism:
    """A validated kinetic mechanism.

    Build one with :meth:`bundled`, :meth:`from_file` or :meth:`from_toml`.
    """

    def __init__(self, inner: PyMechanism):
        self._inner = inner

    @classmethod
    def bundled(cls, name: str) -> Mechanism:
        """A mechanism shipped with sasktran2; see :meth:`bundled_names`."""
        return cls(PyMechanism.bundled(name))

    @classmethod
    def from_toml(cls, text: str) -> Mechanism:
        return cls(PyMechanism.from_toml(text))

    @classmethod
    def from_file(cls, path: str | Path) -> Mechanism:
        return cls.from_toml(Path(path).read_text())

    @staticmethod
    def bundled_names() -> list[str]:
        return PyMechanism.bundled_names()

    @property
    def name(self) -> str:
        return self._inner.name

    @property
    def version(self) -> str:
        return self._inner.version

    @property
    def description(self) -> str:
        return self._inner.description

    @property
    def states(self) -> list[str]:
        """Ids of the excited states that are solved for."""
        return self._inner.states

    @property
    def background(self) -> list[str]:
        """Ids of the species whose number densities are inputs."""
        return self._inner.background

    @property
    def rate_inputs(self) -> list[str]:
        """Names of the per-molecule rates [s^-1] that are inputs, e.g. photolysis rates."""
        return self._inner.rate_inputs

    @property
    def processes(self) -> list[str]:
        return self._inner.process_ids

    @property
    def references(self) -> dict[str, str]:
        return self._inner.references

    def __repr__(self) -> str:
        return (
            f"Mechanism({self.name!r} v{self.version}: {len(self.states)} states, "
            f"{len(self.background)} background species, {len(self.processes)} processes, "
            f"{len(self.rate_inputs)} rate inputs)"
        )


def solve(
    mechanism: Mechanism,
    atmosphere: xr.Dataset,
    rates: xr.Dataset | Mapping[str, np.ndarray] | None = None,
) -> xr.Dataset:
    """Steady-state populations of the mechanism's states.

    Parameters
    ----------
    mechanism
        The kinetic mechanism.
    atmosphere
        ``temperature_k`` [K] on a single dimension, one variable per
        background species named by its id holding the number density
        [m^-3], and optionally ``pressure_pa`` [Pa], used to derive the total
        number density ``M`` when the mechanism needs it and it is not given.
        Other variables are ignored.
    rates
        Per-molecule rates [s^-1] named by the mechanism's rate inputs, on the
        same grid as ``atmosphere``.

    Returns
    -------
    xr.Dataset
        ``density`` [m^-3] and ``production``/``loss`` [m^-3 s^-1] by state;
        ``process_rate`` [m^-3 s^-1] by process, which for radiative
        processes is the photon volume emission rate; ``photon_ver`` for the
        radiative processes only, by transition; and ``relative_residual``,
        the steady-state imbalance relative to the largest production.
    """
    temperature = atmosphere["temperature_k"]
    if temperature.ndim != 1:
        msg = "atmosphere.temperature_k must be one-dimensional"
        raise ValueError(msg)
    (dim,) = temperature.dims

    def profile(values) -> np.ndarray:
        return np.ascontiguousarray(np.asarray(values, dtype=np.float64))

    densities = {
        name: profile(atmosphere[name])
        for name in mechanism.background
        if name in atmosphere
    }
    rates = {} if rates is None else rates
    rate_inputs = {
        name: profile(rates[name]) for name in mechanism.rate_inputs if name in rates
    }
    pressure = (
        profile(atmosphere["pressure_pa"]) if "pressure_pa" in atmosphere else None
    )

    result = mechanism._inner.solve_steady_state(
        profile(temperature), densities, rate_inputs, pressure
    )

    inner = mechanism._inner
    transitions = inner.transitions
    radiative = [index for index, _, _, _ in transitions]
    process_ids = inner.process_ids

    coords = {
        "state": mechanism.states,
        "process": process_ids,
        "process_kind": ("process", inner.process_kinds),
        "process_reference": ("process", inner.process_references),
        "transition": [process_ids[i] for i in radiative],
        "transition_upper": ("transition", [upper for _, upper, _, _ in transitions]),
        "transition_lower": ("transition", [lower for _, _, lower, _ in transitions]),
        "transition_wavelength_nm": (
            "transition",
            [np.nan if wl is None else wl for _, _, _, wl in transitions],
        ),
    }
    if dim in atmosphere.coords:
        coords[dim] = atmosphere[dim]

    ds = xr.Dataset(
        {
            "density": (("state", dim), result["state_density_m3"]),
            "production": (("state", dim), result["production_m3_s"]),
            "loss": (("state", dim), result["loss_m3_s"]),
            "process_rate": (("process", dim), result["process_rate_m3_s"]),
            "photon_ver": (
                ("transition", dim),
                result["process_rate_m3_s"][radiative, :],
            ),
            "relative_residual": ((dim,), result["relative_residual"]),
        },
        coords=coords,
        attrs={"mechanism": mechanism.name, "mechanism_version": mechanism.version},
    )
    ds["density"].attrs["units"] = "m^-3"
    for name in ("production", "loss", "process_rate"):
        ds[name].attrs["units"] = "m^-3 s^-1"
    ds["photon_ver"].attrs["units"] = "photons m^-3 s^-1"
    return ds


def budget(mechanism: Mechanism, solution: xr.Dataset, state: str) -> xr.DataArray:
    """Contributions of each process to a state's population change [m^-3 s^-1].

    Production is positive and loss negative; processes that do not affect
    the state are dropped. At steady state the contributions sum to about zero.
    """
    states = mechanism.states
    if state not in states:
        msg = f"'{state}' is not a state of mechanism '{mechanism.name}'"
        raise ValueError(msg)
    target = states.index(state)

    processes, state_index, coefficients = mechanism._inner.state_stoichiometry()
    rows = [
        (p, c)
        for p, s, c in zip(processes, state_index, coefficients, strict=True)
        if s == target and c != 0.0
    ]
    ids = mechanism.processes
    coefficient = xr.DataArray(
        [c for _, c in rows],
        dims="process",
        coords={"process": [ids[p] for p, _ in rows]},
    )
    # Select by id, so the solution may hold other mechanisms' processes too.
    contributions = (
        solution["process_rate"].sel(process=coefficient["process"]) * coefficient
    )
    contributions.attrs["units"] = "m^-3 s^-1"
    contributions.name = f"budget {state}"
    return contributions


# Uses Mechanism and solve above.
from .photochemistry import (  # noqa: E402
    PHOTOCHEMICAL_SPECIES,
    add_photochemical_species,
)
