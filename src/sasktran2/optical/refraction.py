from __future__ import annotations

import numpy as np

from sasktran2.units import kelvin_to_celsius


def ciddor_index_of_refraction(
    temperatures_k: np.array,
    pressure_pa: np.array,
    specific_humidity: np.array,
    x_co2: np.array,
    wavelength_nm: float,
) -> np.array:
    """
    Calculates the index of refraction of air using the formula of Ciddor (1996).

    Parameters
    ----------
    temperatures_k : np.array
        Temperature in Kelvin
    pressure_pa : np.array
        Pressures in Pascals
    specific_humidity : np.array
        Specific humidity [kg/kg]
    x_co2 : np.array
        Molar fraction of CO2 in [ppm]
    wavelength_nm : float
        Wavelength in [nm]

    Returns
    -------
    np.array
        Array of index of refraction values at the given conditions
    """
    temperatures_c = kelvin_to_celsius(temperatures_k)
    # Convert specific humidity to VMR
    x_h2o = specific_humidity / (
        specific_humidity + (1 - specific_humidity) / (18.01528 / 28.9647)
    )  #  Molar mass of water vapor / Molar mass of air

    # Constants
    w0 = 295.235
    w1 = 2.6422
    w2 = -0.03238
    w3 = 0.004028

    k0 = 238.0185
    k1 = 5792105
    k2 = 57.362
    k3 = 167917

    a0 = 1.58123e-6
    a1 = -2.9331e-8
    a2 = 1.1043e-10

    b0 = 5.707e-6
    b1 = -2.051e-8

    c0 = 1.9898e-4
    c1 = -2.376e-6

    d = 1.83e-11
    e = -0.765e-8

    Pr1 = 101325
    TR1 = 288.15

    Za = 0.9995922115
    Pvs = 0.00985938

    R = 8.314472
    Mv = 0.018015

    S = 1 / ((wavelength_nm / 1e3) ** 2)  # Convert wavelength to micrometers

    ras = 1e-8 * (k1 / (k0 - S) + k3 / (k2 - S))

    rvs = 1.022e-8 * (w0 + w1 * S + w2 * S**2 + w3 * S**3)

    Ma = 0.0289635 + 1.2011e-8 * (x_co2 - 400)

    raxs = ras * (1 + 5.34e-7 * (x_co2 - 450))

    T = temperatures_c + 273.15

    Zm = (
        1
        - (pressure_pa / T)
        * (
            a0
            + a1 * temperatures_c
            + a2 * temperatures_c**2
            + (b0 + b1 * temperatures_c) * x_h2o
            + (c0 + c1 * temperatures_c) * x_h2o**2
        )
        + (pressure_pa / T) ** 2 * (d + e * x_h2o**2)
    )

    rho_axs = Pr1 * Ma / (Za * R * TR1)

    rho_v = x_h2o * pressure_pa * Mv / (Zm * R * T)

    rho_a = (1 - x_h2o) * pressure_pa * Ma / (Zm * R * T)

    return 1 + (rho_a / rho_axs) * raxs + (rho_v / Pvs) * rvs
