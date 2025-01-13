from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_vandaele():
    """
    Tests the NO2Vandaele optical property against a few hardcoded test results obtained from SASKTRAN legacy
    """
    vand = sk.optical.NO2Vandaele()

    config = sk.Config()

    altitude_grid = np.arange(0, 65001, 1000)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    # Old SASKTRAN test results
    test_cases = []
    test_cases.append(
        {
            "temperature": 240,
            "values": np.array(
                [1.82452103e-19, 3.14023410e-19, 3.91827004e-19, 4.61405350e-20]
            ),
        }
    )
    test_cases.append(
        {
            "temperature": 228,
            "values": np.array(
                [1.81900948e-19, 3.11019054e-19, 3.85870165e-19, 4.60555094e-20]
            ),
        }
    )

    # once again old sasktran is in air wavelengths
    wavel = sk.optical.air_wavelength_to_vacuum_wavelength(
        np.array([310, 330, 350, 600])
    )

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    for test_case in test_cases:
        atmosphere.temperature_k[:] = test_case["temperature"]

        result = vand.atmosphere_quantities(atmo=atmosphere).extinction[10, :]

        np.testing.assert_allclose(result, test_case["values"] * 1e-4, rtol=1e-5)


def test_no2_hitran_compare():
    vand = sk.optical.NO2Vandaele()
    hitran = sk.optical.HITRANUV("no2")

    config = sk.Config()

    altitude_grid = np.arange(0, 65001, 1000)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    wavel = np.arange(350, 500, 0.001)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere.temperature_k[:] = 240

    xs_vand = vand.atmosphere_quantities(atmo=atmosphere).extinction[10, :]
    xs_hitran = hitran.atmosphere_quantities(atmo=atmosphere).extinction[10, :]

    np.testing.assert_allclose(xs_hitran, xs_vand, rtol=1e-3)
