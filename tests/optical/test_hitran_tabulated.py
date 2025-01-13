from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


@pytest.mark.skipif(
    not sk.appconfig.are_extended_db_downloaded(),
    reason="Extended database not downloaded",
)
def test_hitran_moderate():
    opt_prop = sk.optical.HITRANTabulated("H2O")

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

    # Old SASKTRAN wavelengths are messed up
    wavel = np.arange(280, 800, 0.1)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    _ = opt_prop.atmosphere_quantities(atmo=atmosphere)
