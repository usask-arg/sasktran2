from __future__ import annotations

from copy import copy

import numpy as np
import sasktran2 as sk
import xarray as xr
from sasktran2.optical.database import (
    OpticalDatabaseGenericScatterer,
    OpticalDatabaseGenericScattererRust,
)


def _storage_delta(atmosphere: sk.Atmosphere, constituent):
    ssa = copy(atmosphere.storage.ssa)
    ext = copy(atmosphere.storage.total_extinction)
    legendre = copy(atmosphere.storage.leg_coeff)

    constituent.add_to_atmosphere(atmosphere)

    return (
        copy(atmosphere.storage.ssa) - ssa,
        copy(atmosphere.storage.total_extinction) - ext,
        copy(atmosphere.storage.leg_coeff) - legendre,
    )


def _test_scenarios():
    config = sk.Config()
    # config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 4

    altitude_grid = np.arange(0, 65001, 1000.0)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    for tan_alt in np.arange(10000, 60000, 2000):
        viewing_geo.add_ray(sk.TangentAltitudeSolar(tan_alt, 0, 600000, 0.6))

    wavel = np.array([310, 330, 350, 600])

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    vmr_altitude_grid = np.array([10000.0, 30000.0, 60000.0])
    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(), vmr_altitude_grid, np.ones_like(vmr_altitude_grid) * 1e-6
    )
    atmosphere.surface.albedo[:] = 0.3

    atmos = []
    atmos.append(atmosphere)

    scen = []

    for atmo in atmos:
        scen.append(
            {
                "config": config,
                "geometry": geometry,
                "viewing_geo": viewing_geo,
                "atmosphere": atmo,
            }
        )

    return scen


def test_scattering_db_construction():
    mie = sk.optical.database.OpticalDatabaseGenericScatterer(
        sk.appconfig.database_root().joinpath("cross_sections/mie/sulfate_test.nc")
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

    _ = mie.atmosphere_quantities(
        atmo=atmosphere,
        lognormal_median_radius=np.ones_like(atmosphere.model_geometry.altitudes())
        * 105,
    ).extinction[10, :]


def test_scattering_db_wf():
    """
    Tests that weighting functions with respect to the scattering DB auxillary variables (e.g. particle size)
    is working correctly

    The precision of these weighting functions can be quite poor
    """
    scens = _test_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

        atmosphere["strat_aerosol"] = sk.test_util.scenarios.test_aerosol_constituent(
            atmosphere.model_geometry.altitudes(), extinction_space=False
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        atmosphere["strat_aerosol"].lognormal_median_radius[:] = 107

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["strat_aerosol"].lognormal_median_radius,
            0.0001,
            engine,
            atmosphere,
            "wf_strat_aerosol_lognormal_median_radius",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_strat_aerosol_lognormal_median_radius"],
            radiance["wf_strat_aerosol_lognormal_median_radius_numeric"],
            wf_dim="strat_aerosol_altitude",
            decimal=3,
        )


def test_scattering_db_wf_extinction():
    """
    Tests that weighting functions with respect to the scattering DB auxillary variables (e.g. particle size)
    is working correctly

    The precision of these weighting functions can be quite poor
    """
    scens = _test_scenarios()

    for scen in scens:
        atmosphere = scen["atmosphere"]

        atmosphere["strat_aerosol"] = sk.test_util.scenarios.test_aerosol_constituent(
            atmosphere.model_geometry.altitudes(), extinction_space=True
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        atmosphere["strat_aerosol"].lognormal_median_radius[:] = 107

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["strat_aerosol"].lognormal_median_radius,
            0.000001,
            engine,
            atmosphere,
            "wf_strat_aerosol_lognormal_median_radius",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_strat_aerosol_lognormal_median_radius"],
            radiance["wf_strat_aerosol_lognormal_median_radius_numeric"],
            wf_dim="strat_aerosol_altitude",
            decimal=2,
        )


def test_scattering_db_rust_impl():
    """
    Tests that the scattering db rust implementation agrees with the Python reference implementation
    """
    scens = _test_scenarios()

    alt_grid = np.arange(0, 65001, 1000.0)
    ext = np.zeros(len(alt_grid))
    ext[:] = 1e-7

    wavel = np.array([310, 330, 350, 600, 745.0])

    r_constituents = []
    py_constituents = []

    db2 = sk.database.MieDatabase(
        sk.mie.distribution.LogNormalDistribution(),
        sk.mie.refractive.H2SO4(),
        wavel,
        median_radius=np.array([100, 200]),
        mode_width=np.array([1.5, 1.6, 1.7]),
    )

    r_constituents.append(
        sk.constituent.ExtinctionScatterer(
            db2,
            alt_grid,
            ext,
            745.0,
            median_radius=np.ones_like(ext) * 120.0,
            mode_width=np.ones_like(ext) * 1.55,
        )
    )

    py_constituents.append(
        sk.constituent.ExtinctionScatterer(
            OpticalDatabaseGenericScatterer(db2.path()),
            alt_grid,
            ext,
            745.0,
            median_radius=np.ones_like(ext) * 120.0,
            mode_width=np.ones_like(ext) * 1.55,
        )
    )

    db1 = sk.database.MieDatabase(
        sk.mie.distribution.LogNormalDistribution().freeze(mode_width=1.6),
        sk.mie.refractive.H2SO4(),
        wavel,
        median_radius=np.array([100, 200]),
    )

    r_constituents.append(
        sk.constituent.ExtinctionScatterer(
            db1,
            alt_grid,
            ext,
            745.0,
            median_radius=np.ones_like(ext) * 120.0,
        )
    )

    py_constituents.append(
        sk.constituent.ExtinctionScatterer(
            OpticalDatabaseGenericScatterer(db1.path()),
            alt_grid,
            ext,
            745.0,
            median_radius=np.ones_like(ext) * 120.0,
        )
    )

    db0 = sk.database.MieDatabase(
        sk.mie.distribution.LogNormalDistribution().freeze(
            mode_width=1.6, median_radius=120
        ),
        sk.mie.refractive.H2SO4(),
        wavel,
    )

    r_constituents.append(sk.constituent.ExtinctionScatterer(db0, alt_grid, ext, 745.0))

    py_constituents.append(
        sk.constituent.ExtinctionScatterer(
            OpticalDatabaseGenericScatterer(db0.path()), alt_grid, ext, 745.0
        )
    )

    for scen in scens:
        atmosphere = scen["atmosphere"]

        for r_const, py_const in zip(r_constituents, py_constituents, strict=False):

            r_ssa, r_k, r_leg = _storage_delta(atmosphere, r_const)
            py_ssa, py_k, py_leg = _storage_delta(atmosphere, py_const)

            np.testing.assert_allclose(
                r_ssa,
                py_ssa,
                rtol=1e-5,
                atol=1e-5,
                err_msg="SSA does not match between Rust and Python",
            )

            np.testing.assert_allclose(
                r_k,
                py_k,
                rtol=1e-5,
                atol=1e-5,
                err_msg="Extinction does not match between Rust and Python",
            )

            np.testing.assert_allclose(
                r_leg,
                py_leg,
                rtol=1e-5,
                atol=1e-5,
                err_msg="Legendre coefficients do not match between Rust and Python",
            )


def test_scattering_db_dim4_construction(tmp_path):
    wavelengths_nm = np.array([350.0, 500.0])
    legendre_order = np.array([0, 1])
    param0 = np.array([0.1, 0.2])
    param1 = np.array([1.0, 2.0])
    param2 = np.array([10.0, 20.0])
    coeffs = {
        "p0": 1e-6,
        "p1": 2e-7,
        "p2": 3e-8,
        "wavelength_nm": 1e-10,
    }

    shape = (
        len(param0),
        len(param1),
        len(param2),
        len(wavelengths_nm),
        len(legendre_order),
    )
    p0_grid, p1_grid, p2_grid, wavelength_grid = np.meshgrid(
        param0, param1, param2, wavelengths_nm, indexing="ij"
    )
    xs_total = (
        1e-4
        + coeffs["p0"] * p0_grid
        + coeffs["p1"] * p1_grid
        + coeffs["p2"] * p2_grid
        + coeffs["wavelength_nm"] * wavelength_grid
    )
    xs_scattering = 0.8 * xs_total
    lm_a1 = np.ones(shape)
    lm_a2 = np.zeros(shape)
    lm_a3 = np.zeros(shape)
    lm_a4 = np.zeros(shape)
    lm_b1 = np.zeros(shape)
    lm_b2 = np.zeros(shape)

    ds = xr.Dataset(
        {
            "xs_total": (("p0", "p1", "p2", "wavelength_nm"), xs_total),
            "xs_scattering": (("p0", "p1", "p2", "wavelength_nm"), xs_scattering),
            "lm_a1": (("p0", "p1", "p2", "wavelength_nm", "legendre"), lm_a1),
            "lm_a2": (("p0", "p1", "p2", "wavelength_nm", "legendre"), lm_a2),
            "lm_a3": (("p0", "p1", "p2", "wavelength_nm", "legendre"), lm_a3),
            "lm_a4": (("p0", "p1", "p2", "wavelength_nm", "legendre"), lm_a4),
            "lm_b1": (("p0", "p1", "p2", "wavelength_nm", "legendre"), lm_b1),
            "lm_b2": (("p0", "p1", "p2", "wavelength_nm", "legendre"), lm_b2),
        },
        coords={
            "p0": param0,
            "p1": param1,
            "p2": param2,
            "wavelength_nm": wavelengths_nm,
            "legendre": legendre_order,
        },
    )

    db_file = tmp_path / "dim4_scattering.nc"
    ds.to_netcdf(db_file)

    scatterer = OpticalDatabaseGenericScattererRust(db_file)

    altitudes = np.array([0.0, 1000.0, 2000.0])
    p0_profile = np.array([0.12, 0.15, 0.18])
    p1_profile = np.array([1.2, 1.5, 1.8])
    p2_profile = np.array([12.0, 15.0, 18.0])
    quants = scatterer.cross_sections(
        wavelengths_nm=wavelengths_nm,
        altitudes_m=altitudes,
        p0=p0_profile,
        p1=p1_profile,
        p2=p2_profile,
    )

    expected_extinction = (
        1e-4
        + coeffs["p0"] * p0_profile[:, np.newaxis]
        + coeffs["p1"] * p1_profile[:, np.newaxis]
        + coeffs["p2"] * p2_profile[:, np.newaxis]
        + coeffs["wavelength_nm"] * wavelengths_nm[np.newaxis, :]
    )
    np.testing.assert_allclose(quants.extinction, expected_extinction, rtol=1e-12)
    np.testing.assert_allclose(quants.ssa, 0.8 * expected_extinction, rtol=1e-12)

    derivs = scatterer.cross_section_derivatives(
        wavelengths_nm=wavelengths_nm,
        altitudes_m=altitudes,
        p0=p0_profile,
        p1=p1_profile,
        p2=p2_profile,
    )

    assert set(derivs.keys()) == {"p0", "p1", "p2"}

    for key in ["p0", "p1", "p2"]:
        np.testing.assert_allclose(
            derivs[key].reshape(len(altitudes), len(wavelengths_nm)),
            coeffs[key],
            rtol=1e-12,
            atol=1e-18,
        )


def test_scattering_db_dim3_second_parameter_derivative(tmp_path):
    wavelengths_nm = np.array([350.0, 500.0])
    legendre_order = np.array([0, 1])
    param0 = np.array([0.1, 0.2])
    param1 = np.array([1.0, 2.0])
    coeffs = {"p0": 1e-6, "p1": 2e-7, "wavelength_nm": 1e-10}

    p0_grid, p1_grid, wavelength_grid = np.meshgrid(
        param0, param1, wavelengths_nm, indexing="ij"
    )
    xs_total = (
        1e-4
        + coeffs["p0"] * p0_grid
        + coeffs["p1"] * p1_grid
        + coeffs["wavelength_nm"] * wavelength_grid
    )
    xs_scattering = 0.8 * xs_total
    shape = (*xs_total.shape, len(legendre_order))
    lm_a1 = np.ones(shape)
    lm_a2 = np.zeros(shape)
    lm_a3 = np.zeros(shape)
    lm_a4 = np.zeros(shape)
    lm_b1 = np.zeros(shape)
    lm_b2 = np.zeros(shape)

    ds = xr.Dataset(
        {
            "xs_total": (("p0", "p1", "wavelength_nm"), xs_total),
            "xs_scattering": (("p0", "p1", "wavelength_nm"), xs_scattering),
            "lm_a1": (("p0", "p1", "wavelength_nm", "legendre"), lm_a1),
            "lm_a2": (("p0", "p1", "wavelength_nm", "legendre"), lm_a2),
            "lm_a3": (("p0", "p1", "wavelength_nm", "legendre"), lm_a3),
            "lm_a4": (("p0", "p1", "wavelength_nm", "legendre"), lm_a4),
            "lm_b1": (("p0", "p1", "wavelength_nm", "legendre"), lm_b1),
            "lm_b2": (("p0", "p1", "wavelength_nm", "legendre"), lm_b2),
        },
        coords={
            "p0": param0,
            "p1": param1,
            "wavelength_nm": wavelengths_nm,
            "legendre": legendre_order,
        },
    )

    db_file = tmp_path / "dim3_scattering.nc"
    ds.to_netcdf(db_file)
    scatterer = OpticalDatabaseGenericScattererRust(db_file)

    altitudes = np.array([0.0, 1000.0, 2000.0])
    derivs = scatterer.cross_section_derivatives(
        wavelengths_nm=wavelengths_nm,
        altitudes_m=altitudes,
        p0=np.array([0.12, 0.15, 0.18]),
        p1=np.array([1.2, 1.5, 1.8]),
    )

    assert set(derivs.keys()) == {"p0", "p1"}

    for key in ["p0", "p1"]:
        np.testing.assert_allclose(
            derivs[key].reshape(len(altitudes), len(wavelengths_nm)),
            coeffs[key],
            rtol=1e-12,
            atol=1e-18,
        )
