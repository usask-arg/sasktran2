from __future__ import annotations

import numpy as np
from scipy import interpolate
from scipy.special import roots_legendre

from sasktran2._core_rust import (
    compute_greek_coefficients as _compute_greek_coefficients,
)


def compute_greek_coefficients(
    p11: np.ndarray,
    p12: np.ndarray,
    p22: np.ndarray,
    p33: np.ndarray,
    p34: np.ndarray,
    p44: np.ndarray,
    angle_grid: np.array,
    num_coeff: int,
):
    """
    Calculates the greek coefficients a1, a2, a3, a4, b1, b2 given the Legendre expansion of the phase function
    elements.

    Parameters
    ----------
    p11 : np.array
        Phase function P11
    p12 : np.array
        Phase function P12
    p22 : np.array
        Phase function P22
    p33 : np.array
        Phase function P33
    p34 : np.array
        Phase function P34
    p44 : np.array
        Phase function P44
    angle_grid : np.array
        Angular grid the phase functions are specified on.  Should fully span 0 to 180
    num_coeff: int
        Maximum number of coefficients to return in the expansion.

    Returns
    -------
    lm_a1 : np.array
        Greek coefficients for a1
    lm_a2 : np.array
        Greek coefficients for a2
    lm_a3 : np.array
        Greek coefficients for a3
    lm_a4 : np.array
        Greek coefficients for a4
    lm_b1 : np.array
        Greek coefficients for b1
    lm_b2 : np.array
        Greek coefficients for b2
    """
    cos_theta = np.cos(np.deg2rad(angle_grid))[::-1]

    c = 0.995
    nodes, weights = roots_legendre(num_coeff)

    nodes_left = (c - (-1)) / 2 * nodes + (c + (-1)) / 2
    weights_left = (c - (-1)) / 2 * weights

    nodes_right = (1 - c) / 2 * nodes + (1 + c) / 2
    weights_right = (1 - c) / 2 * weights
    all_weights = np.concatenate([weights_left, weights_right])
    cos_angle_grid = np.concatenate([nodes_left, nodes_right])

    # Keep SciPy's shape-preserving interpolation and split quadrature convention.
    # The all-order Wigner recurrence and six projections share the Rust kernel
    # used by the Mie integrator.
    phase = np.stack(
        [
            interpolate.PchipInterpolator(cos_theta, np.asarray(p)[:, ::-1], axis=1)(
                cos_angle_grid
            )
            for p in (p11, p12, p22, p33, p34, p44)
        ],
        axis=1,
    )
    coefficients = _compute_greek_coefficients(
        phase, cos_angle_grid, all_weights, num_coeff
    )
    return tuple(coefficients[:, i, :] for i in range(6))
