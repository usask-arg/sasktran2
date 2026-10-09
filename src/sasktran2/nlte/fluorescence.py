"""Solar resonance fluorescence of OH A-X.

Each A-X line pumps its upper level at the rate absorbed from the actinic
flux, and each upper level re-emits on all of its lines with Einstein-A
branching. Lower levels are the thermal rotational levels of OH(X, v=0).

Lines come from HITRAN (molecule 13, the Yousefi et al. 2018 A-X list, which
is the same as MoLLIST-OH). Not included: absorption from vibrationally
excited OH(X, v>0), collisional quenching of OH(A) (a few percent of its
radiative rate at 60 km, more below), predissociation of v'>=2, and
hyperfine structure.

Yousefi, M., P. F. Bernath, J. Hodges and T. Masseron (2018), A new line
list for the A2Sigma+-X2Pi electronic transition of OH, J. Quant. Spectrosc.
Radiat. Transfer 217, 416-424.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np
import xarray as xr
from scipy import constants

from sasktran2.database.hitran_line import HITRANLineDatabase
from sasktran2.database.web import StandardDatabase

C2_CM_K = 1.4387769
#: Molar mass of 16OH [g mol^-1].
OH_MOLAR_MASS = 17.00274


@dataclass(frozen=True)
class OHAXLines:
    """16OH A-X lines from HITRAN."""

    wavenumber_cminv: np.ndarray
    einstein_a_s: np.ndarray
    lower_energy_cminv: np.ndarray
    g_upper: np.ndarray
    g_lower: np.ndarray
    v_upper: np.ndarray
    v_lower: np.ndarray
    #: Index of each line's upper level.
    upper_level: np.ndarray

    @property
    def wavelength_nm(self) -> np.ndarray:
        return 1.0e7 / self.wavenumber_cminv


@functools.cache
def oh_ax_lines(db: HITRANLineDatabase | None = None) -> OHAXLines:
    """Read the 16OH A-X lines from the HITRAN OH file."""
    path = (db or HITRANLineDatabase()).path("OH") / "OH.data"
    rows = []
    with path.open() as f:
        for line in f:
            if line[2] != "1" or "A" not in line[67:82]:
                continue
            rows.append(
                (
                    float(line[3:15]),
                    float(line[25:35]),
                    float(line[45:55]),
                    float(line[146:153]),
                    float(line[153:160]),
                    int(line[67:82].split()[-1]),
                    int(line[82:97].split()[-1]),
                    line[67:82].split()[0],
                )
            )
    nu, a, elow, gu, gl, vu, vl, component = (
        np.array(c) for c in zip(*rows, strict=True)
    )
    # In A2Sigma+ each spin component (F1, F2), v' and J (from g') is one
    # level of definite parity.
    _, upper_level = np.unique(
        np.char.add(np.char.add(component, vu.astype(str)), gu.astype(str)),
        return_inverse=True,
    )
    return OHAXLines(nu, a, elow, gu, gl, vu, vl, upper_level.ravel())


def oh_ax_fluorescence(
    lines: OHAXLines,
    temperature_k,
    oh_density_m3,
    actinic_flux_at_lines,
):
    """Photon volume emission rates of OH A-X fluorescence.

    Parameters
    ----------
    lines
        From :func:`oh_ax_lines`.
    temperature_k, oh_density_m3
        Profiles (altitude,) for the OH(X, v=0) rotational distribution and
        total OH number density [m^-3].
    actinic_flux_at_lines
        Actinic flux [photons m^-2 s^-1 nm^-1] at each line centre,
        (altitude, line).

    Returns
    -------
    photon_ver : np.ndarray
        Total photon VER [photons m^-3 s^-1], (altitude,).
    weights : np.ndarray
        Fraction of it in each line, (altitude, line).
    """
    temperature_k = np.atleast_1d(np.asarray(temperature_k, dtype=float))
    oh = np.atleast_1d(np.asarray(oh_density_m3, dtype=float))
    flux = np.atleast_2d(np.asarray(actinic_flux_at_lines, dtype=float))

    # Thermal populations of OH(X, v=0) levels, from the lower levels of the
    # v''=0 lines.
    ground = lines.v_lower == 0
    levels, index = np.unique(
        np.stack([np.round(lines.lower_energy_cminv, 3), lines.g_lower]),
        axis=1,
        return_inverse=True,
    )
    index = index.ravel()
    level_used = np.zeros(levels.shape[1], dtype=bool)
    level_used[index[ground]] = True
    energy, g = levels[0][level_used], levels[1][level_used]
    boltzmann = g[np.newaxis, :] * np.exp(
        -C2_CM_K * (energy - energy.min())[np.newaxis, :] / temperature_k[:, np.newaxis]
    )
    fraction_by_level = np.zeros((temperature_k.size, levels.shape[1]))
    fraction_by_level[:, level_used] = boltzmann / boltzmann.sum(axis=1, keepdims=True)
    lower_fraction = np.where(ground, fraction_by_level[:, index], 0.0)

    # Integrated absorption cross section [m^2 nm] times flux per nm.
    wavelength_m = 1.0e-2 / lines.wavenumber_cminv
    sigma_int = (
        wavelength_m**4
        / (8.0 * np.pi * constants.c)
        * lines.g_upper
        / lines.g_lower
        * lines.einstein_a_s
        * 1.0e9
    )
    absorbed = oh[:, np.newaxis] * lower_fraction * sigma_int * np.maximum(flux, 0.0)

    n_upper = lines.upper_level.max() + 1
    production = np.zeros((temperature_k.size, n_upper))
    for i in range(temperature_k.size):
        production[i] = np.bincount(
            lines.upper_level, weights=absorbed[i], minlength=n_upper
        )
    total_a = np.bincount(
        lines.upper_level, weights=lines.einstein_a_s, minlength=n_upper
    )
    branching = lines.einstein_a_s / total_a[lines.upper_level]
    line_ver = production[:, lines.upper_level] * branching
    photon_ver = line_ver.sum(axis=1)
    weights = np.divide(
        line_ver,
        photon_ver[:, np.newaxis],
        out=np.zeros_like(line_ver),
        where=photon_ver[:, np.newaxis] > 0,
    )
    return photon_ver, weights


@functools.cache
def _hsrs_photons() -> tuple[np.ndarray, np.ndarray]:
    ds = xr.load_dataset(
        StandardDatabase().path("solar/solar_irradiance_hsrs_2022_11_30_extended.nc")
    )
    wavelength = ds["wavelength"].to_numpy()
    photons = ds["irradiance"].to_numpy() / (
        constants.h * constants.c / (wavelength * 1.0e-9)
    )
    return wavelength, photons


def actinic_flux_at_lines(flux, wavelength_nm, earth_sun_distance_au=None):
    """Actinic flux [photons m^-2 s^-1 nm^-1] at line centres, (altitude, line).

    The TSIS-1 HSRS solar spectrum at each line, times the atmospheric
    transmission (actinic over top-of-atmosphere flux) of ``flux``, an
    :meth:`ActinicFlux.calculate` result, interpolated to the line.
    """
    grid = flux["wavelength"].to_numpy()
    transmission = (flux["actinic_flux"] / flux["solar_flux"]).to_numpy()
    at_lines = np.array([np.interp(wavelength_nm, grid, t) for t in transmission.T])
    distance = (
        flux.attrs.get("earth_sun_distance_au", 1.0)
        if earth_sun_distance_au is None
        else earth_sun_distance_au
    )
    hsrs_wavelength, hsrs = _hsrs_photons()
    solar = np.interp(wavelength_nm, hsrs_wavelength, hsrs) / distance**2
    return at_lines * solar[np.newaxis, :]
