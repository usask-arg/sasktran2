from __future__ import annotations

import numpy as np

import sasktran2 as sk


def default_pure_scattering_atmosphere(
    config: sk.Config, geometry: sk.Geometry1D, ssa=1, albedo=0
):
    new_alt_grid = geometry.altitudes()

    alt_grid = np.arange(0, 100001, 1000)

    extinction = np.array(
        [
            7.07906113e-05,
            6.46250950e-05,
            5.86431083e-05,
            5.29850715e-05,
            4.77339013e-05,
            4.29288557e-05,
            3.85773022e-05,
            3.46642865e-05,
            3.11600517e-05,
            2.80258050e-05,
            2.52180748e-05,
            2.26734428e-05,
            2.02816648e-05,
            1.79778464e-05,
            1.57467704e-05,
            1.36034281e-05,
            1.15882231e-05,
            9.77118267e-06,
            8.18898344e-06,
            6.84554061e-06,
            5.72584994e-06,
            4.80319926e-06,
            4.04164882e-06,
            3.41027519e-06,
            2.88467502e-06,
            2.44547520e-06,
            2.07720643e-06,
            1.76744635e-06,
            1.50616412e-06,
            1.28521627e-06,
            1.09795686e-06,
            9.38934589e-07,
            8.03656135e-07,
            6.88499846e-07,
            5.90874787e-07,
            5.08080278e-07,
            4.37762355e-07,
            3.77949584e-07,
            3.26990467e-07,
            2.83500936e-07,
            2.46320382e-07,
            2.14474894e-07,
            1.87146548e-07,
            1.63647786e-07,
            1.43400047e-07,
            1.25915952e-07,
            1.10782493e-07,
            9.76450474e-08,
            8.62059194e-08,
            7.62163724e-08,
            6.74679955e-08,
            5.97856954e-08,
            5.30219866e-08,
            4.70523226e-08,
            4.17712677e-08,
            3.70893448e-08,
            3.29295331e-08,
            2.92235209e-08,
            2.59140539e-08,
            2.29536982e-08,
            2.03028465e-08,
            1.79281189e-08,
            1.58010825e-08,
            1.38964913e-08,
            1.21932261e-08,
            1.06731597e-08,
            9.31932242e-09,
            8.11636527e-09,
            7.05027714e-09,
            6.10817608e-09,
            5.27815905e-09,
            4.54919569e-09,
            3.91105464e-09,
            3.35204107e-09,
            2.86870538e-09,
            2.45077836e-09,
            2.09082449e-09,
            1.78181885e-09,
            1.51726918e-09,
            1.29127795e-09,
            1.09856137e-09,
            9.35341764e-10,
            7.95620813e-10,
            6.76822327e-10,
            5.75848867e-10,
            4.90034206e-10,
            4.17093171e-10,
            3.55073944e-10,
            3.02313970e-10,
            2.57400129e-10,
            2.19133412e-10,
            1.86465279e-10,
            1.58426278e-10,
            1.34262300e-10,
            1.13404218e-10,
            9.54097598e-11,
            7.99225004e-11,
            6.66436434e-11,
            5.53133171e-11,
            4.56988383e-11,
            3.75879135e-11,
        ]
    )

    extinction = np.interp(new_alt_grid, alt_grid, extinction).reshape(
        len(new_alt_grid), 1
    )

    atmo = sk.Atmosphere(geometry, config, numwavel=1)

    atmo.storage.total_extinction[:] = extinction
    atmo.storage.ssa[:] = np.ones_like(extinction) * ssa

    atmo.leg_coeff.a1[0, :, 0] = 1
    atmo.leg_coeff.a1[2, :, 0] = 0.5

    if atmo.nstokes == 3:
        atmo.leg_coeff.a2[2] = 3
        atmo.leg_coeff.b1[2] = np.sqrt(6.0) / 2

    sk.climatology.us76.add_us76_standard_atmosphere(atmo)

    atmo.surface.albedo[:] = albedo

    return atmo


def test_aerosol_constituent(altitude_grid: np.array, extinction_space=True):
    alts = np.arange(0, 40000, 1000.0)

    ext = np.array(
        [
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            2.88027019e-07,
            3.69604997e-07,
            3.10548639e-07,
            2.75319733e-07,
            2.43405259e-07,
            2.64027971e-07,
            2.57261097e-07,
            2.50073674e-07,
            2.45212374e-07,
            2.09064034e-07,
            1.71399035e-07,
            1.57256087e-07,
            1.51544489e-07,
            1.49381653e-07,
            1.48309802e-07,
            1.28552798e-07,
            1.03291371e-07,
            8.03122894e-08,
            6.15167010e-08,
            3.81730206e-08,
            2.27081327e-08,
            7.19716081e-09,
            7.79190668e-09,
            6.16738043e-09,
            4.63396327e-09,
        ]
    )

    ext = np.interp(altitude_grid, alts, ext)

    ext_wavel = 525

    mie = sk.optical.database.OpticalDatabaseGenericScatterer(
        sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
    )
    radius = np.ones_like(altitude_grid) * 105
    const = sk.constituent.ExtinctionScatterer(
        mie, altitude_grid, ext, ext_wavel, lognormal_median_radius=radius
    )

    if extinction_space:
        return const

    return sk.constituent.NumberDensityScatterer(
        mie, altitude_grid, const.number_density, lognormal_median_radius=radius
    )
