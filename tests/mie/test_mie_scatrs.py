from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import sasktran2.mie.distribution as distribution
import xarray as xr
from sasktran2.mie.distribution import integrate_mie_cpp, integrate_mie_scatrs
from scipy.stats import lognorm

COEFFS = ["lm_a1", "lm_a2", "lm_a3", "lm_a4", "lm_b1", "lm_b2"]


def _refractive_index(wavelength):
    return (1.45 - 0.01j) + 0 * wavelength


@pytest.mark.parametrize(
    ("dist", "dist_args"),
    [
        (distribution.LogNormalDistribution(), {"median_radius": 150, "mode_width": 1.6}),
        (distribution.GammaDistribution(), {"alpha": 4, "beta": 0.02}),
        (distribution.UniformDistribution(), {"min_radius": 50, "max_radius": 300}),
        (distribution.UniformDistribution(), {"min_radius": 0, "max_radius": 300}),
        (
            distribution.TriangularDistribution(),
            {"min_radius": 50, "max_radius": 400, "center_radius": 150},
        ),
        (
            distribution.TriangularDistribution(),
            {"min_radius": 0, "max_radius": 400, "center_radius": 400},
        ),
        (
            distribution.LogNormalDistribution().freeze(mode_width=1.4),
            {"median_radius": 120},
        ),
    ],
)
def test_scatrs_distributions_match_scipy(dist, dist_args):
    """
    Each scatrs distribution is the scipy one: the scatrs integrals agree with the fixed
    quadrature of integrate_mie_cpp over the scipy pdf, within the latter's accuracy
    """
    assert dist.supports_scatrs
    wavelengths = np.array([350.0, 600.0, 1000.0])
    reference = integrate_mie_cpp(
        [dist.distribution(**dist_args)],
        _refractive_index,
        wavelengths,
        num_coeffs=16,
        num_quad=64,
    )
    result = integrate_mie_scatrs(
        [dist.scatrs_distribution(**dist_args)],
        _refractive_index,
        wavelengths,
        num_coeffs=16,
        rtol=1e-7,
    )
    for key in ["xs_total", "xs_scattering", "xs_absorption"]:
        np.testing.assert_allclose(result[key], reference[key], rtol=2e-4)
    for key in COEFFS:
        np.testing.assert_allclose(result[key], reference[key], atol=2e-2)
    np.testing.assert_array_equal(result["lm_a1"].isel(legendre=0), 1.0)


def test_scatrs_threads_do_not_change_results():
    dists = [
        distribution.LogNormalDistribution().scatrs_distribution(
            median_radius=r, mode_width=1.5
        )
        for r in [80, 160]
    ]
    inputs = (dists, _refractive_index, np.array([400.0, 800.0]))
    serial = integrate_mie_scatrs(*inputs, num_coeffs=16, num_threads=1)
    threaded = integrate_mie_scatrs(*inputs, num_coeffs=16, num_threads=0)
    xr.testing.assert_identical(serial, threaded)


def test_scatrs_warns_when_not_converged():
    dists = [
        distribution.LogNormalDistribution().scatrs_distribution(
            median_radius=500, mode_width=1.8
        )
    ]
    with pytest.warns(UserWarning, match="did not meet rtol"):
        integrate_mie_scatrs(
            dists, _refractive_index, np.array([500.0]), rtol=1e-10, max_panels=1
        )


def test_scatrs_rejects_gain_medium():
    dists = [("lognormal", [100.0, 1.5])]
    with pytest.raises(ValueError, match="gain medium"):
        integrate_mie_scatrs(dists, lambda _: 1.5 + 0.01j, np.array([500.0]))


class _ScipyOnly(distribution.ParticleSizeDistribution):
    def __init__(self) -> None:
        super().__init__("scipy_only_lognormal")

    def distribution(self, **kwargs):
        return lognorm(np.log(1.5), scale=kwargs["median_radius"])

    def args(self):
        return ["median_radius"]


def test_mie_database_backends(tmp_path):
    refractive = sk.mie.refractive.RefractiveIndex(_refractive_index, "test_index")
    wavelengths = np.array([500.0, 800.0])

    db = sk.database.MieDatabase(
        distribution.LogNormalDistribution(),
        refractive,
        wavelengths,
        db_root=tmp_path,
        median_radius=[80.0, 120.0],
        mode_width=[1.4, 1.6],
    )
    assert db._backend == "scatrs"
    ds = db.load_ds()
    assert ds["xs_total"].dims == ("wavelength_nm", "median_radius", "mode_width")
    assert ds["lm_a1"].sizes["legendre"] == 64

    # The table holds each distribution as integrated on its own
    single = integrate_mie_scatrs(
        [
            distribution.LogNormalDistribution().scatrs_distribution(
                median_radius=120.0, mode_width=1.4
            )
        ],
        _refractive_index,
        wavelengths,
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        ds["xs_total"].sel(median_radius=120.0, mode_width=1.4),
        single["xs_total"].isel(distribution=0),
        rtol=3e-4,
    )

    # Distributions without a scatrs equivalent keep the sasktran2_cpp backend
    db = sk.database.MieDatabase(
        _ScipyOnly(), refractive, wavelengths, db_root=tmp_path, median_radius=[100.0]
    )
    assert db._backend == "sasktran2_cpp"
    with pytest.raises(ValueError, match="scatrs backend does not support"):
        sk.database.MieDatabase(
            _ScipyOnly(),
            refractive,
            wavelengths,
            db_root=tmp_path,
            backend="scatrs",
            median_radius=[100.0],
        )


def test_scatrs_phase_matrix_matches_sasktran2_cpp():
    """
    The optional phase matrix has the sasktran2 conventions: the normalization and signs of
    integrate_mie_cpp, at its cosines
    """
    dist = distribution.LogNormalDistribution()
    dist_args = {"median_radius": 150, "mode_width": 1.5}
    wavelengths = np.array([500.0])
    reference = integrate_mie_cpp(
        [dist.distribution(**dist_args)], _refractive_index, wavelengths, num_coeffs=64
    )
    cos_angles = reference["cos_angle"].to_numpy()
    result = integrate_mie_scatrs(
        [dist.scatrs_distribution(**dist_args)],
        _refractive_index,
        wavelengths,
        rtol=1e-7,
        phase_cos_angles=cos_angles,
    )
    np.testing.assert_array_equal(result["cos_angle"], cos_angles)
    scale = np.abs(reference["p11"]).max().item()
    for key in ["p11", "p12", "p33", "p34"]:
        np.testing.assert_allclose(result[key], reference[key], atol=1e-4 * scale)

    without = integrate_mie_scatrs(
        [dist.scatrs_distribution(**dist_args)], _refractive_index, wavelengths
    )
    assert "p11" not in without
    assert "cos_angle" not in without.coords


def test_mie_database_phase_matrix(tmp_path):
    refractive = sk.mie.refractive.RefractiveIndex(_refractive_index, "test_index")
    args = {
        "db_root": tmp_path,
        "median_radius": [80.0, 120.0],
    }
    dist = distribution.LogNormalDistribution().freeze(mode_width=1.5)
    wavelengths = np.array([500.0, 800.0])

    default = sk.database.MieDatabase(dist, refractive, wavelengths, **args)
    assert "p11" not in default.load_ds()

    db = sk.database.MieDatabase(
        dist, refractive, wavelengths, phase_matrix=True, **args
    )
    assert db.path() != default.path()
    ds = db.load_ds()
    assert ds["p11"].sizes["cos_angle"] == 181
    np.testing.assert_allclose(ds["cos_angle"], np.cos(np.deg2rad(np.arange(181))))
    # Normalized: half the integral of p11 over the cosine is 1
    p11 = ds["p11"].isel(wavelength_nm=1, median_radius=0).to_numpy()
    mu = ds["cos_angle"].to_numpy()
    np.testing.assert_allclose(-0.5 * np.trapezoid(p11, mu), 1.0, rtol=1e-3)

    with pytest.raises(ValueError, match="phase_matrix needs the scatrs backend"):
        sk.database.MieDatabase(
            dist,
            refractive,
            wavelengths,
            backend="sasktran2_cpp",
            phase_matrix=True,
            **args,
        )
