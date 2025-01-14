from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
from sasktran2 import appconfig


def test_mie_construction():
    mie = sk.mie.LinearizedMie()
    size_param = np.array([3.0, 4.0, 5.0])
    cos_angles = np.linspace(-1, 1, 100)
    refractive_index = 1.5 + 0.0j

    output = mie.calculate(size_param, refractive_index, cos_angles, True)

    np.testing.assert_allclose(output.size_parameter, size_param)
    np.testing.assert_allclose(output.cos_angles, cos_angles)


def test_qext_q_Sca_s1_s2():
    mie = sk.mie.LinearizedMie()
    refractive_index = 0.75
    size_param = np.array(
        [0.101, 0.2]
    )  # not checking values for 0.2, just there to show that we can calculate with multiple size param and angles
    angles = np.arange(0, 181, 30)
    cos_angles = np.cos(angles * np.pi / 180)
    output = mie.calculate(size_param, refractive_index, cos_angles, True)

    np.testing.assert_allclose(output.values.Qext[0], 0.000008, atol=1e-06, rtol=0)
    np.testing.assert_allclose(output.values.Qsca[0], 0.000008, atol=1e-06, rtol=0)

    np.testing.assert_allclose(
        output.values.S1[0][0].real, 2.048754e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][0].imag, -1.756419e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][0].real, 2.048754e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][0].imag, -1.756419e-04, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][1].real, 2.048754e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][1].imag, -1.755965e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][1].real, 1.774273e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][1].imag, -1.520629e-04, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][2].real, 2.048753e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][2].imag, -1.754726e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][2].real, 1.024377e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][2].imag, -8.771198e-05, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][3].real, 2.048751e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][3].imag, -1.753033e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][3].real, 1.845057e-15, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][3].imag, -3.247270e-08, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][4].real, 2.048750e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][4].imag, -1.751341e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][4].real, -1.024375e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][4].imag, 8.759147e-05, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][5].real, 2.048749e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][5].imag, -1.750103e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][5].real, -1.774269e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][5].imag, 1.515715e-04, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][4].real, 2.048749e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][5].imag, -1.749650e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][6].real, -2.048749e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][6].imag, 1.749650e-04, atol=1e-06, rtol=0
    )


def test_mie_integration():
    mie = sk.mie.LinearizedMie()
    refractive_index = 1 - 0.1j
    size_param = np.array(
        [1, 2, 3]
    )  # not checking values for these, just there to show that we can calculate with multiple size param and angles
    angles = np.linspace(0, 180, 3)
    cos_angles = np.cos(angles * np.pi / 180)
    output = mie.calculate(size_param, refractive_index, cos_angles, True)
    np.testing.assert_allclose(
        output.values.Qext[0], 0.25925991210268584, atol=1e-06, rtol=0
    )
    from scipy.stats import rv_continuous

    class gaussian_gen(rv_continuous):
        "Gaussian distribution"

        def _pdf(self, x):
            mu = 500  # 500
            sig = 100  # 100
            return np.exp(-(((x - mu) / sig) ** 2) / 2.0) / np.sqrt(2.0 * np.pi) / sig

    prob_dist = gaussian_gen(name="gaussian")

    def refrac_index_fn(w):
        return 1.0 - 0.1j + 0 * w

    wavelengths = np.array([400, 500])
    my_mie = sk.mie.LinearizedMie()

    out = sk.mie.integrate_mie(
        my_mie,
        prob_dist,
        refrac_index_fn,
        wavelengths,
        num_angles=3,
        num_quad=20,
        compute_coeffs=True,
        num_coeffs=2,
    )
    np.testing.assert_allclose(out["p11"][0][0].data, 8.11553105e01, atol=1e-06, rtol=0)
    np.testing.assert_allclose(
        out["p11"][0][1].data, 1.16641980e-02, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        out["p11"][0][2].data, 5.64571171e-03, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(out["p12"][0][0].data, 0, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p12"][0][1].data, 0.01117573, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p12"][0][2].data, 0, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p33"][0][0].data, 8.11553105e01, atol=1e-06, rtol=0)
    np.testing.assert_allclose(
        out["p33"][0][1].data, -6.33399061e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        out["p33"][0][2].data, -5.64571171e-03, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(out["p34"][0][0].data, 0, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p34"][0][1].data, -0.00099857, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p34"][0][2].data, 0, atol=1e-06, rtol=0)

    np.testing.assert_allclose(
        out["xs_total"][0].data, 1077477.99015462, atol=1e-03, rtol=0
    )
    np.testing.assert_allclose(
        out["xs_scattering"][0].data, 348484.73075379, atol=1e-03, rtol=0
    )
    np.testing.assert_allclose(
        out["xs_absorption"][0].data, 728993.25940083, atol=1e-03, rtol=0
    )


def test_mie_database():
    _ = pytest.importorskip("sasktran")

    db_root_legacy = appconfig.database_root().joinpath("legacy")
    db_root_current = appconfig.database_root().joinpath("current")

    refrac = sk.mie.refractive.H2SO4()
    distribution = sk.mie.distribution.LogNormalDistribution()
    mie_db = sk.database.MieDatabase(
        distribution,
        refrac,
        wavelengths_nm=np.arange(270, 1000, 50.0),
        median_radius=[100, 150, 200],
        mode_width=[1.5, 1.7],
        db_root=db_root_legacy,
        max_legendre_moments=100,
    )
    mie_db.load_ds()

    mie_db_new = sk.database.MieDatabase(
        distribution,
        refrac,
        wavelengths_nm=np.arange(270, 1000, 50.0),
        median_radius=[100, 150, 200],
        mode_width=[1.5, 1.7],
        backend="sasktran2",
        db_root=db_root_current,
        max_legendre_moments=100,
    )

    # testing select legendre moments for both to make sure we get reasonably close results
    tolerance = 1e-7
    np.testing.assert_allclose(
        mie_db._database["lm_a1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        mie_db_new._database["lm_a1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        mie_db_new._database["lm_a1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        mie_db_new._database["lm_a1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        mie_db_new._database["lm_a1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        atol=tolerance,
        rtol=0,
    )

    np.testing.assert_allclose(
        mie_db._database["lm_a2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        mie_db_new._database["lm_a2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        mie_db_new._database["lm_a2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        mie_db_new._database["lm_a2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        mie_db_new._database["lm_a2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        atol=tolerance,
        rtol=0,
    )

    np.testing.assert_allclose(
        mie_db._database["lm_a3"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        mie_db_new._database["lm_a3"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a3"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        mie_db_new._database["lm_a3"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a3"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        mie_db_new._database["lm_a3"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a3"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        mie_db_new._database["lm_a3"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        atol=tolerance,
        rtol=0,
    )

    np.testing.assert_allclose(
        mie_db._database["lm_a4"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        mie_db_new._database["lm_a4"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a4"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        mie_db_new._database["lm_a4"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a4"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        mie_db_new._database["lm_a4"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_a4"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        mie_db_new._database["lm_a4"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        atol=tolerance,
        rtol=0,
    )

    np.testing.assert_allclose(
        mie_db._database["lm_b1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        mie_db_new._database["lm_b1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_b1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        mie_db_new._database["lm_b1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_b1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        mie_db_new._database["lm_b1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_b1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        mie_db_new._database["lm_b1"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        atol=tolerance,
        rtol=0,
    )

    np.testing.assert_allclose(
        mie_db._database["lm_b2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        mie_db_new._database["lm_b2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=0)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_b2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        mie_db_new._database["lm_b2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=20)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_b2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        mie_db_new._database["lm_b2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=40)
        .data,
        atol=tolerance,
        rtol=0,
    )
    np.testing.assert_allclose(
        mie_db._database["lm_b2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        mie_db_new._database["lm_b2"]
        .isel(median_radius=0, mode_width=0, wavelength_nm=0, legendre=60)
        .data,
        atol=tolerance,
        rtol=0,
    )


def test_mie_db_freezing():
    refrac = sk.mie.refractive.H2SO4()
    # try freezing
    distribution = sk.mie.distribution.LogNormalDistribution().freeze(
        median_radius=160, mode_width=1.6
    )
    mie_db_f = sk.database.MieDatabase(
        distribution,
        refrac,
        wavelengths_nm=np.arange(270, 1000, 50.0),
    )

    mie_db_f.load_ds()
