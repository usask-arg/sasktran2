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


def test_hitran_db():
    pytest.importorskip("hapi")

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

    wvnum = np.arange(7294, 7340, 0.01)

    h2o_opt_prop = sk.database.HITRANDatabase(
        molecule="H2O",
        start_wavenumber=7294,
        end_wavenumber=7340,
        wavenumber_resolution=0.01,
        reduction_factor=1,
        backend="sasktran2",
    )
    atmosphere = sk.Atmosphere(geometry, config, wavenumber_cminv=wvnum)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    _ = h2o_opt_prop.atmosphere_quantities(atmo=atmosphere)
