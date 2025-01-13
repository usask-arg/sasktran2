from __future__ import annotations

import numpy as np
from sasktran2.optical.refraction import ciddor_index_of_refraction


def test_ciddor():
    # Check with 0 specific humidity
    n1 = ciddor_index_of_refraction(
        temperatures_k=273.15 + 20,
        pressure_pa=101325,
        specific_humidity=0,
        x_co2=450,
        wavelength_nm=633,
    )
    np.testing.assert_almost_equal(n1, 1.0002718)

    # Adjust pressure
    n2 = ciddor_index_of_refraction(
        temperatures_k=273.15 + 20,
        pressure_pa=60000,
        specific_humidity=0,
        x_co2=450,
        wavelength_nm=633,
    )
    np.testing.assert_almost_equal(n2, 1.000160924)

    # Adjust temperature
    n3 = ciddor_index_of_refraction(
        temperatures_k=273.15 + 50,
        pressure_pa=100000,
        specific_humidity=0,
        x_co2=450,
        wavelength_nm=633,
    )
    np.testing.assert_almost_equal(n3, 1.000243285)

    # Adjust wavelength
    n4 = ciddor_index_of_refraction(
        temperatures_k=273.15 + 20,
        pressure_pa=101325,
        specific_humidity=0,
        x_co2=450,
        wavelength_nm=1700,
    )
    np.testing.assert_almost_equal(n4, 1.000268479)

    # Adjust specific humidity
    # RH of 100% gives saturation vapor pressure of 2339 Pa at 20 C
    # Enhancement factor
    p = 100e3
    f = 1.00062 + 3.14e-8 * p + 5.60e-7 * (20**2)
    xv = f * 2339 / p

    # Convert to specific humidity
    q = xv / (xv + (1 - xv) * (18.01528 / 28.9647))

    n5 = ciddor_index_of_refraction(
        temperatures_k=273.15 + 20,
        pressure_pa=p,
        specific_humidity=q,
        x_co2=450,
        wavelength_nm=633,
    )
    np.testing.assert_allclose(n5, 1.000267394)
