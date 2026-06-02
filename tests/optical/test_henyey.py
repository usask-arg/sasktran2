from __future__ import annotations

import numpy as np
import sasktran2 as sk
import xarray as xr
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


def test_henyey_dim4_construction():
    wavelengths_nm = np.array([350.0, 500.0])
    param0 = np.array([0.1, 0.2])
    param1 = np.array([1.0, 2.0])
    param2 = np.array([10.0, 20.0])

    shape = (len(param0), len(param1), len(param2), len(wavelengths_nm))
    xs_total = np.full(shape, 1e-4)
    ssa = np.full(shape, 0.9)
    g = np.full(shape, 0.75)

    ds = xr.Dataset(
        {
            "xs_total": (("p0", "p1", "p2", "wavelength_nm"), xs_total),
            "ssa": (("p0", "p1", "p2", "wavelength_nm"), ssa),
            "asymmetry_parameter": (("p0", "p1", "p2", "wavelength_nm"), g),
        },
        coords={
            "p0": param0,
            "p1": param1,
            "p2": param2,
            "wavelength_nm": wavelengths_nm,
        },
    )

    hg = HenyeyGreenstein(db=ds, max_num_moments=16)

    altitudes = np.array([0.0, 1000.0, 2000.0])
    quants = hg.cross_sections(
        wavelengths_nm=wavelengths_nm,
        altitudes_m=altitudes,
        p0=np.full_like(altitudes, 0.15),
        p1=np.full_like(altitudes, 1.5),
        p2=np.full_like(altitudes, 15.0),
    )

    np.testing.assert_allclose(quants.extinction, 1e-4, rtol=1e-6)
    np.testing.assert_allclose(quants.ssa, 0.9, rtol=1e-6)
