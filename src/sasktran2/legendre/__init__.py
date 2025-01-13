from __future__ import annotations

import numpy as np
from scipy import interpolate
from scipy.special import roots_legendre

from sasktran2 import util


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

    wigner00 = util.WignerD(0, 0)
    wigner22 = util.WignerD(2, 2)
    wigner2m2 = util.WignerD(2, -2)
    wigner02 = util.WignerD(0, 2)

    wavelength_dim = p11.shape[0]
    lm_a1 = np.zeros((wavelength_dim, num_coeff))
    lm_a2 = np.zeros((wavelength_dim, num_coeff))
    lm_a3 = np.zeros((wavelength_dim, num_coeff))
    lm_a4 = np.zeros((wavelength_dim, num_coeff))
    lm_b1 = np.zeros((wavelength_dim, num_coeff))
    lm_b2 = np.zeros((wavelength_dim, num_coeff))

    c = 0.995
    nodes, weights = roots_legendre(num_coeff)

    nodes_left = (c - (-1)) / 2 * nodes + (c + (-1)) / 2
    weights_left = (c - (-1)) / 2 * weights

    nodes_right = (1 - c) / 2 * nodes + (1 + c) / 2
    weights_right = (1 - c) / 2 * weights
    all_weights = np.concatenate([weights_left, weights_right])
    cos_angle_grid = np.concatenate([nodes_left, nodes_right])

    lpoly_00 = np.zeros((num_coeff, 2 * num_coeff))
    lpoly_22 = np.zeros((num_coeff, 2 * num_coeff))
    lpoly_2m2 = np.zeros((num_coeff, 2 * num_coeff))
    lpoly_02 = np.zeros((num_coeff, 2 * num_coeff))

    for i in range(num_coeff):
        lpoly_00[i] = wigner00.d(np.arccos(cos_angle_grid), i)
        lpoly_00[i] *= all_weights / (2.0 / (2.0 * i + 1))
        lpoly_22[i] = wigner22.d(np.arccos(cos_angle_grid), i)
        lpoly_22[i] *= all_weights / (2.0 / (2.0 * i + 1))
        lpoly_2m2[i] = wigner2m2.d(np.arccos(cos_angle_grid), i)
        lpoly_2m2[i] *= all_weights / (2.0 / (2.0 * i + 1))
        lpoly_02[i] = wigner02.d(np.arccos(cos_angle_grid), i)
        lpoly_02[i] *= all_weights / (2.0 / (2.0 * i + 1))

    p11_interp = interpolate.PchipInterpolator(
        cos_theta, np.transpose(np.flip(p11, axis=1))
    )(cos_angle_grid)
    p12_interp = interpolate.PchipInterpolator(
        cos_theta, np.transpose(np.flip(p12, axis=1))
    )(cos_angle_grid)
    p22_interp = interpolate.PchipInterpolator(
        cos_theta, np.transpose(np.flip(p22, axis=1))
    )(cos_angle_grid)
    p33_interp = interpolate.PchipInterpolator(
        cos_theta, np.transpose(np.flip(p33, axis=1))
    )(cos_angle_grid)
    p34_interp = interpolate.PchipInterpolator(
        cos_theta, np.transpose(np.flip(p34, axis=1))
    )(cos_angle_grid)
    p44_interp = interpolate.PchipInterpolator(
        cos_theta, np.transpose(np.flip(p44, axis=1))
    )(cos_angle_grid)

    lm_a1 = lpoly_00 @ p11_interp
    lm_a4 = lpoly_00 @ p44_interp

    lm_b1 = lpoly_02 @ p12_interp
    lm_b2 = lpoly_02 @ p34_interp * -1

    temp_1 = lpoly_22 @ (p22_interp + p33_interp)
    temp_2 = lpoly_2m2 @ (p22_interp - p33_interp)
    lm_a2 = (temp_1 + temp_2) / 2
    lm_a3 = (temp_1 - temp_2) / 2

    return (
        np.transpose(lm_a1),
        np.transpose(lm_a2),
        np.transpose(lm_a3),
        np.transpose(lm_a4),
        np.transpose(lm_b1),
        np.transpose(lm_b2),
    )
