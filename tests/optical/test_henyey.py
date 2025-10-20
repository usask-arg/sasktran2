from __future__ import annotations

import numpy as np
import sasktran2 as sk
from sasktran2.optical import HenyeyGreenstein


def test_henyey_construction():
    wavel = np.arange(300, 801, 10)
    xsec = np.ones_like(wavel) * 1e-4
    ssa = np.ones_like(wavel) * 0.9
    g = np.ones_like(wavel) * 0.75

    hg = HenyeyGreenstein.from_parameters(
        wavelength_nm=wavel,
        xs_total=xsec,
        ssa=ssa,
        g=g,
    )

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

    wavel = np.arange(350, 500, 0.1)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    aq = hg.atmosphere_quantities(
        atmo=atmosphere,
    )

    np.testing.assert_allclose(aq.extinction, 1e-4, rtol=1e-6)
    np.testing.assert_allclose(aq.ssa, 0.9, rtol=1e-6)

    np.testing.assert_allclose(aq.leg_coeff[0], 1.0, rtol=1e-6)
    np.testing.assert_allclose(aq.leg_coeff[5], (0.75) ** 5 * (2 * 5 + 1), rtol=1e-6)
