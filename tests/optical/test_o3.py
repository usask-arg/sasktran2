from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def test_dbm():
    """
    Tests the DBM optical property against a few hardcoded test results obtained from SASKTRAN legacy
    """
    dbm = sk.optical.O3DBM()

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
                [8.71858864e-20, 3.03768968e-21, 1.03810358e-22, 5.20901605e-21]
            ),
        }
    )
    test_cases.append(
        {
            "temperature": 228,
            "values": np.array(
                [8.47813970e-20, 2.78021858e-21, 7.71502986e-23, 5.22090221e-21]
            ),
        }
    )

    # Old SASKTRAN wavelengths are messed up
    wavel = sk.optical.air_wavelength_to_vacuum_wavelength(
        np.array([310, 330, 350, 600])
    )

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    for test_case in test_cases:
        atmosphere.temperature_k[:] = test_case["temperature"]

        result = dbm.atmosphere_quantities(atmo=atmosphere).extinction[10, :]

        np.testing.assert_allclose(result, test_case["values"] * 1e-4)


def test_o3_hitran():
    """
    Checks O3DBM against O3 HITRAN.  They have significant differences so we check that this does fail.
    """
    dbm = sk.optical.O3DBM()
    hitran = sk.optical.HITRANUV("o3")

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

    wavel = np.arange(250, 350, 0.001)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    xs_dbm = dbm.atmosphere_quantities(atmo=atmosphere).extinction[10, :]
    xs_hitran = hitran.atmosphere_quantities(atmo=atmosphere).extinction[10, :]

    with pytest.raises(AssertionError):
        np.testing.assert_allclose(xs_dbm, xs_hitran)
