"""Solar resonance fluorescence of OH A-X.

Each A-X line pumps its upper level at the rate absorbed from the actinic
flux, and each upper level re-emits on all of its lines with Einstein-A
branching. Lower levels are the thermal rotational levels of OH(X, v=0).

Collisions with N2 and O2 quench OH(A) and move v' to v'-1. Quenching removes
about 8% of OH(A, v'=0) at 70 km, 21% at 60 km and 75% at 40 km; below about
60 km, N2 moves most of v'=1 into v'=0 before it radiates. Transferred
molecules enter v'-1 with a thermal rotational distribution.

Lines come from HITRAN (molecule 13, the Yousefi et al. 2018 A-X list, which
is the same as MoLLIST-OH). Not included: absorption from vibrationally
excited OH(X, v>0), rotational energy transfer within OH(A) (which moves the
line distribution toward the local temperature below about 60 km, but not
the band total), predissociation of v'>=2, prompt OH(A) from H2O photolysis
at Lyman-alpha, and hyperfine structure.

Yousefi, M., P. F. Bernath, J. Hodges and T. Masseron (2018), A new line
list for the A2Sigma+-X2Pi electronic transition of OH, J. Quant. Spectrosc.
Radiat. Transfer 217, 416-424.

Bailey, A. E., D. E. Heard, P. H. Paul and M. J. Pilling (1997), Collisional
quenching of OH (A2Sigma+, v'=0) by N2, O2 and CO2 between 204 and 294 K,
J. Chem. Soc. Faraday Trans. 93, 2915-2920.

Williams, L. R. and D. R. Crosley (1996), Collisional vibrational energy
transfer of OH (A2Sigma+, v'=1), J. Chem. Phys. 104, 6507-6514.
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

#: OH(A, v'=0) electronic quenching cross sections [m^2] at 204 and 294 K
#: (Bailey et al. 1997), interpolated linearly in temperature and held
#: constant outside that range.
_QUENCH_V0_TEMPERATURE_K = (204.0, 294.0)
_QUENCH_V0_CROSS_SECTION_M2 = {"N2": (6.8e-20, 4.3e-20), "O2": (23.5e-20, 18.6e-20)}
#: Molar masses of the colliders [g mol^-1].
_COLLIDER_MOLAR_MASS = {"N2": 28.0134, "O2": 31.9988}
#: OH(A, v'>=1) at room temperature [m^3 s^-1] (Williams and Crosley 1996):
#: transfer v' -> v'-1 (N2 30.1, O2 2.8 A^2), and electronic quenching (O2;
#: N2 removal is almost all transfer). Transfer is faster at lower
#: temperatures, so these are lower limits. v'>=2 use the v'=1 values.
_TRANSFER_V1_M3_S = {"N2": 2.33e-16, "O2": 2.1e-17}
_QUENCH_V1_M3_S = {"N2": 0.0, "O2": 1.85e-16}


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


def _quench_v0_m3_s(collider: str, temperature_k: np.ndarray) -> np.ndarray:
    """OH(A, v'=0) electronic quenching rate coefficient [m^3 s^-1]."""
    t = np.clip(temperature_k, *_QUENCH_V0_TEMPERATURE_K)
    sigma = np.interp(
        t, _QUENCH_V0_TEMPERATURE_K, _QUENCH_V0_CROSS_SECTION_M2[collider]
    )
    mu = (
        OH_MOLAR_MASS
        * _COLLIDER_MOLAR_MASS[collider]
        / (OH_MOLAR_MASS + _COLLIDER_MOLAR_MASS[collider])
        * 1.0e-3
        / constants.N_A
    )
    return sigma * np.sqrt(8.0 * constants.k * temperature_k / (np.pi * mu))


def oh_ax_fluorescence(
    lines: OHAXLines,
    temperature_k,
    oh_density_m3,
    actinic_flux_at_lines,
    *,
    n2_density_m3=None,
    o2_density_m3=None,
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
    n2_density_m3, o2_density_m3
        Collider profiles [m^-3] for quenching and vibrational transfer in
        OH(A); without them, OH(A) only radiates.

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
    population = _upper_populations(
        lines, production, total_a, temperature_k, n2_density_m3, o2_density_m3
    )
    line_ver = population[:, lines.upper_level] * lines.einstein_a_s
    photon_ver = line_ver.sum(axis=1)
    # Rows without emission get uniform placeholder weights: the VER there is
    # zero, and line-list emission requires every row to sum to one.
    weights = np.divide(
        line_ver,
        photon_ver[:, np.newaxis],
        out=np.full_like(line_ver, 1.0 / line_ver.shape[1]),
        where=photon_ver[:, np.newaxis] > 0,
    )
    return photon_ver, weights


def _upper_populations(
    lines, production, total_a, temperature_k, n2_density_m3, o2_density_m3
):
    """Steady-state OH(A) level populations [m^-3], (altitude, level)."""
    if n2_density_m3 is None and o2_density_m3 is None:
        return production / total_a
    shape = temperature_k.shape
    colliders = {
        "N2": np.broadcast_to(
            np.zeros(shape) if n2_density_m3 is None else n2_density_m3, shape
        ),
        "O2": np.broadcast_to(
            np.zeros(shape) if o2_density_m3 is None else o2_density_m3, shape
        ),
    }
    # Per level: vibrational number, degeneracy and energy above OH(X).
    first = np.unique(lines.upper_level, return_index=True)[1]
    level_v = lines.v_upper[first]
    level_g = lines.g_upper[first]
    level_energy = (lines.lower_energy_cminv + lines.wavenumber_cminv)[first]

    quench_v0 = sum(_quench_v0_m3_s(c, temperature_k) * n for c, n in colliders.items())
    quench_v1 = sum(_QUENCH_V1_M3_S[c] * n for c, n in colliders.items())
    transfer = sum(_TRANSFER_V1_M3_S[c] * n for c, n in colliders.items())

    population = np.zeros_like(production)
    inflow = production.copy()
    for v in range(int(level_v.max()), -1, -1):
        levels = level_v == v
        quench = quench_v0 if v == 0 else quench_v1
        out = 0.0 if v == 0 else transfer
        population[:, levels] = inflow[:, levels] / (
            total_a[levels] + (quench + out)[:, np.newaxis]
        )
        below = level_v == v - 1
        if v == 0 or not below.any():
            continue
        moved = population[:, levels].sum(axis=1) * out
        energy = level_energy[below]
        thermal = level_g[below] * np.exp(
            -C2_CM_K
            * (energy - energy.min())[np.newaxis, :]
            / temperature_k[:, np.newaxis]
        )
        inflow[:, below] += (
            moved[:, np.newaxis] * thermal / thermal.sum(axis=1, keepdims=True)
        )
    return population


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
