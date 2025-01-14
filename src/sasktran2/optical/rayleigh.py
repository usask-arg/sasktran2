from __future__ import annotations

import numpy as np

from sasktran2.optical import pressure_temperature_to_numberdensity
from sasktran2.units import celsius_to_kelvin


def _o2_refrac_bates(wavelengths_um: np.array) -> np.array:
    """
    Calculates the O2 refractive index using the parameterization suggested by Bates.

    Parameters
    ----------
    wavelengths_nm : np.array
       Wavelengths in [nm]

    Returns
    -------
    np.array
       O2 Refractive index, (n-1) * 10**8
    """
    coeff_ranges = [
        {"coeff": [23796.7, 168988.4], "range": [0, 0.221]},
        {"coeff": [22120.4, 203187.6], "range": [0.221, 0.288]},
        {"coeff": [20564.8, 248089.9], "range": [0.288, 0.546]},
        {"coeff": [21351.1, 218567.0], "range": [0.546, np.inf]},
    ]

    result = np.zeros_like(wavelengths_um)

    for ele in coeff_ranges:
        good = (wavelengths_um > ele["range"][0]) & (wavelengths_um <= ele["range"][1])

        result[good] = ele["coeff"][0] + ele["coeff"][1] / (
            40.9 - wavelengths_um[good] ** (-2)
        )

    return result


def _n2_refrac_bates(wavelengths_um: np.array) -> np.array:
    """
    Calculates the N2 refractive index using the parameterization suggested by Bates.

    Parameters
    ----------
    wavelengths_nm : np.array
       Wavelengths in [nm]

    Returns
    -------
    np.array
       N2 Refractive index, (n-1) * 10**8
    """
    coeff_ranges = [
        {"coeff": [6998.749, 3233582.0], "range": [0, 0.254]},
        {"coeff": [5989.242, 3363266.3], "range": [0.254, 0.468]},
        {"coeff": [6855.200, 3243157.0], "range": [0.468, np.inf]},
    ]

    result = np.zeros_like(wavelengths_um)

    for ele in coeff_ranges:
        good = (wavelengths_um > ele["range"][0]) & (wavelengths_um <= ele["range"][1])
        delta_lambda = 0.468 - wavelengths_um[good]

        result[good] = (
            ele["coeff"][0]
            + ele["coeff"][1] / (144 - wavelengths_um[good] ** (-2))
            + 2.27684009 * np.sign(delta_lambda) * np.exp(-np.abs(delta_lambda) / 0.003)
        )

    return result


def _ar_refrac_bates(wavelengths_um: np.array) -> np.array:
    """
    Calculates the Ar refractive index using the parameterization suggested by Bates.

    Parameters
    ----------
    wavelengths_nm : np.array
       Wavelengths in [nm]

    Returns
    -------
    np.array
       Ar Refractive index, (n-1) * 10**8
    """
    nsq_m_1 = 5.547e-4 * (
        1 + 5.15e-3 * wavelengths_um ** (-2) + 4.19e-5 * wavelengths_um ** (-4)
    )

    return (np.sqrt(nsq_m_1 + 1) - 1) * 1e8


def _co2_refrac_bates(wavelengths_um: np.array) -> np.array:
    """
    Calculates the CO2 refractive index using the parameterization suggested by Bates.

    Parameters
    ----------
    wavelengths_nm : np.array
       Wavelengths in [nm]

    Returns
    -------
    np.array
       CO2 Refractive index, (n-1) * 10**8
    """
    return (
        22822.1
        + 117.8 * wavelengths_um ** (-2)
        + 2406030 / (130 - wavelengths_um ** (-2))
        + 15997 / (38.9 - wavelengths_um ** (-2))
    )


def _o2_king_bates(wavelengths_um: np.array) -> np.array:
    """
    Calculates the O2 King factor using the parameterization suggested by Bates.

    Parameters
    ----------
    wavelengths_nm : np.array
       Wavelengths in [nm]

    Returns
    -------
    np.array
       O2 King factor
    """
    return 1.096 + 1.385e-3 * wavelengths_um ** (-2) + 1.448e-4 * wavelengths_um ** (-4)


def _n2_king_bates(wavelengths_um: np.array) -> np.array:
    """
    Calculates the N2 King factor using the parameterization suggested by Bates.

    Parameters
    ----------
    wavelengths_nm : np.array
       Wavelengths in [nm]

    Returns
    -------
    np.array
       N2 King factor
    """
    return 1.034 + 3.17e-4 * wavelengths_um ** (-2)


def _ar_king_bates(wavelengths_um: np.array) -> np.array:
    """
    Calculates the Ar King factor using the parameterization suggested by Bates.

    Parameters
    ----------
    wavelengths_nm : np.array
       Wavelengths in [nm]

    Returns
    -------
    np.array
       Ar King factor
    """
    return np.ones_like(wavelengths_um)


def _co2_king_bates(wavelengths_um: np.array) -> np.array:
    """
    Calculates the CO2 King factor using the parameterization suggested by Bates.

    Parameters
    ----------
    wavelengths_nm : np.array
       Wavelengths in [nm]

    Returns
    -------
    np.array
       CO2 King factor
    """
    return np.ones_like(wavelengths_um) * 1.15


def rayleigh_cross_section_bates(
    wavelengths_um: np.array,
    n2_percentage: float = 78.084,
    o2_percentage: float = 20.946,
    ar_percentage: float = 0.934,
    co2_percentage: float = 0.036,
) -> np.array:
    """
    Calculates the Rayleigh scattering cross section in units of [m**2] for a given refractive index according to the Bates
    approximation.  Cross section is returned back at (0C, 1013.25 hPa)


    Parameters
    ----------
    wavelengths_um : np.array
        Wavelengths in [um]
    n2_percentage : float, Optional
        Percentage of N2, default 78.084%
    o2_percentage : float, Optional
        Percentage of O2, default 20.946%
    ar_percentage : float, Optional
        Percentage of Ar, default 0.934%
    co2_percentage : float, Optional
        Percentaage of CO2, default 0.036%

    Returns
    -------
    np.array
        The rayleigh scattering cross section in units [m**2]
    np.array
        The effective King factor [unitless]
    """
    lorenz_factors = (
        o2_percentage
        / 100
        * _o2_refrac_bates(wavelengths_um) ** 2
        * _o2_king_bates(wavelengths_um)
        + n2_percentage
        / 100
        * _n2_refrac_bates(wavelengths_um) ** 2
        * _n2_king_bates(wavelengths_um)
        + ar_percentage
        / 100
        * _ar_refrac_bates(wavelengths_um) ** 2
        * _ar_king_bates(wavelengths_um)
        + co2_percentage
        / 100
        * _co2_refrac_bates(wavelengths_um) ** 2
        * _co2_king_bates(wavelengths_um)
    )

    eff_king = (
        o2_percentage / 100 * _o2_king_bates(wavelengths_um)
        + n2_percentage / 100 * _n2_king_bates(wavelengths_um)
        + ar_percentage / 100 * _ar_king_bates(wavelengths_um)
        + co2_percentage / 100 * _co2_king_bates(wavelengths_um)
    )

    num_dens = pressure_temperature_to_numberdensity(
        1013.25 * 100, celsius_to_kelvin(0)
    )

    return (
        32 * np.pi**3 / (3 * num_dens**2 * wavelengths_um**4) * (lorenz_factors) * 1e8,
        eff_king,
    )
