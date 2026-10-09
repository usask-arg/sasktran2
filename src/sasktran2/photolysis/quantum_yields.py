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


# Yankovsky and Vorobeva (2020), Table 1: threshold wavelengths [nm] for
# O2(a, v=0-5), the parameter x at each threshold, and the constants C_v.
_O2A_THRESHOLD_NM = np.array([310.0, 296.0, 284.0, 273.0, 263.0, 254.0])
_O2A_THRESHOLD_X = np.array([0.937, 0.689, 0.576, 0.483, 0.407, 0.339])
_O2A_C = np.array([1.068, 1.233, 0.564, 0.375, 0.377, 0.473])


def o3_o2a_fractions(wavelength_nm) -> np.ndarray:
    """Fractions of the spin-allowed O(1D) channel of O3 photolysis giving O2(a, v=0-5).

    Yankovsky and Vorobeva (2020), Eqs. 3-4: F_0 = C_0 x and
    F_v = C_v x (1 - x - x^2/2 - ... - x^v/2^(v-1)), with level v open only
    shortward of its threshold. Their Eq. 2 for x(lambda) does not reproduce
    their Table 1, so x is interpolated in log space between the tabulated
    (threshold, x) pairs and extrapolated with the end slopes. The fractions
    are renormalised to sum to 1. Longward of the v=0 threshold (310 nm) there
    is no spin-allowed channel and all fractions are zero; O(1D) there comes
    with O2(X) (see :func:`o3_o1d_o2x_fraction`).

    Returns an array of shape ``(6, *wavelength_nm.shape)``.

    Yankovsky, V. A. and E. V. Vorobeva (2020), Model of daytime oxygen
    emissions in the mesopause region and above: a review and new results,
    Atmosphere 11, 116, doi:10.3390/atmos11010116.
    """
    w = np.asarray(wavelength_nm, dtype=float)
    lam = _O2A_THRESHOLD_NM[::-1]
    lnx = np.log(_O2A_THRESHOLD_X[::-1])
    slope_lo = (lnx[1] - lnx[0]) / (lam[1] - lam[0])
    slope_hi = (lnx[-1] - lnx[-2]) / (lam[-1] - lam[-2])
    ln_x = np.where(
        w < lam[0],
        lnx[0] + slope_lo * (w - lam[0]),
        np.where(
            w > lam[-1], lnx[-1] + slope_hi * (w - lam[-1]), np.interp(w, lam, lnx)
        ),
    )
    x = np.minimum(np.exp(ln_x), 1.0)

    fractions = np.zeros((6, *w.shape))
    tail = np.zeros_like(x)
    for v in range(6):
        if v > 0:
            tail = tail + x**v / 2.0 ** (v - 1)
        shape = x if v == 0 else x * (1.0 - tail)
        fractions[v] = np.where(w <= _O2A_THRESHOLD_NM[v], _O2A_C[v] * shape, 0.0)
    fractions = np.maximum(fractions, 0.0)
    total = fractions.sum(axis=0)
    return np.divide(fractions, total, out=np.zeros_like(fractions), where=total > 0)


def o3_o1d_o2x_fraction(wavelength_nm) -> np.ndarray:
    """Fraction of O(1D) from O3 photolysis that comes with O2(X): 1 beyond 310 nm."""
    return 1.0 - o3_o2a_fractions(wavelength_nm).sum(axis=0)


class O3O1DYield:
    """Matsumi (2002) O(1D) yield times a product fraction, as a quantum-yield callable.

    ``level`` 0-5 selects O(1D) + O2(a, v=level); ``"X"`` selects O(1D) +
    O2(X), the spin-forbidden channel beyond 310 nm.
    """

    def __init__(self, level):
        self.level = level

    def __call__(self, wavelength_nm, temperature_k) -> np.ndarray:
        phi = o3_o1d_matsumi2002(wavelength_nm, temperature_k)
        w = np.broadcast_to(np.asarray(wavelength_nm, dtype=float), phi.shape)
        if self.level == "X":
            return phi * o3_o1d_o2x_fraction(w)
        return phi * o3_o2a_fractions(w)[int(self.level)]

    def __eq__(self, other):
        return isinstance(other, O3O1DYield) and other.level == self.level

    def __hash__(self):
        return hash(("O3O1DYield", self.level))
