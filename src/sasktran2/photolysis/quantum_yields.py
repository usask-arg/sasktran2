"""Photolysis quantum yields."""

from __future__ import annotations

import numpy as np

# Matsumi et al. (2002), Table 3.
_MATSUMI_A = (0.8036, 8.9061, 0.1192)
_MATSUMI_X_NM = (304.225, 314.957, 310.737)
_MATSUMI_OMEGA_NM = (5.576, 6.601, 2.187)
_MATSUMI_NU2_CM = 825.518
_MATSUMI_C = 0.0765
# Boltzmann constant in cm^-1 K^-1.
_K_CM_PER_K = 0.695


def o3_o1d_matsumi2002(wavelength_nm, temperature_k) -> np.ndarray:
    """Quantum yield of O(1D) from O3 photolysis.

    The Matsumi et al. (2002) recommendation, as adopted by the JPL
    evaluations: 0.90 up to 305 nm, the parameterised fall-off between 305
    and 328 nm, 0.08 from 328 to 340 nm and zero beyond. The arguments
    broadcast against each other.

    Matsumi, Y., et al. (2002), Quantum yields for production of O(1D) in the
    ultraviolet photolysis of ozone: Recommendation based on evaluation of
    laboratory data, J. Geophys. Res., 107(D3), 4024, doi:10.1029/2001JD000510.
    """
    w, t = np.broadcast_arrays(
        np.asarray(wavelength_nm, dtype=float), np.asarray(temperature_k, dtype=float)
    )
    q2 = np.exp(-_MATSUMI_NU2_CM / (_K_CM_PER_K * t))
    a1, a2, a3 = _MATSUMI_A
    x1, x2, x3 = _MATSUMI_X_NM
    o1, o2, o3 = _MATSUMI_OMEGA_NM
    fall_off = (
        _MATSUMI_C
        + a1 / (1.0 + q2) * np.exp(-(((x1 - w) / o1) ** 4))
        + a2 * (t / 300.0) ** 2 * q2 / (1.0 + q2) * np.exp(-(((x2 - w) / o2) ** 2))
        + a3 * (t / 300.0) ** 1.5 * np.exp(-(((x3 - w) / o3) ** 2))
    )
    return np.select(
        [w <= 305.0, w <= 328.0, w <= 340.0],
        [0.90, fall_off, 0.08],
        default=0.0,
    )


def o3_o3p_matsumi2002(wavelength_nm, temperature_k) -> np.ndarray:
    """Quantum yield of O(3P) from O3 photolysis: one minus :func:`o3_o1d_matsumi2002`."""
    return 1.0 - o3_o1d_matsumi2002(wavelength_nm, temperature_k)
