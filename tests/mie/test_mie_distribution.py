from __future__ import annotations

import pytest
import sasktran2.mie.distribution as distribution

test_dists = [
    (distribution.LogNormalDistribution, {"median_radius": 200, "mode_width": 1.6}),
    (distribution.GammaDistribution, {"alpha": 2, "beta": 1}),
    (distribution.UniformDistribution, {"min_radius": 0, "max_radius": 1}),
    (
        distribution.TriangularDistribution,
        {"min_radius": 0, "max_radius": 1, "center_radius": 0.5},
    ),
]


@pytest.mark.parametrize(("dist", "dist_args"), test_dists)
def test_dist_construction(dist, dist_args):  # noqa: ARG001
    """
    Verify that the different distribution classes can be constructed
    """
    dist()


@pytest.mark.parametrize(("dist", "dist_args"), test_dists)
def test_dist_args(dist, dist_args):
    """
    Verifies that the distribution arguments are correct, and that the scipy distribution can be
    constructed
    """
    prob_dist = dist()

    args = prob_dist.args()

    for arg in dist_args:
        assert arg in args

    prob_dist.distribution(**dist_args)
