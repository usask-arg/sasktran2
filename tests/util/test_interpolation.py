from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_linear_interpolating_matrix():
    """
    Test that the linear_interpolating_matrix function works as expected by comparing against the
    numpy interp function
    """
    from_grid = np.array([0.0, 1000.0, 3000.0, 5000.0, 10000.0])

    to_grid = np.array([-20.0, 500.0, 3000.0, 4000.0, 8000.0, 12000.0])

    rng = np.random.default_rng()
    vals = rng.random(len(from_grid))

    M = sk.util.interpolation.linear_interpolating_matrix(from_grid, to_grid, "extend")
    np.testing.assert_array_almost_equal(
        M @ vals, np.interp(to_grid, from_grid, vals, left=vals[0], right=vals[-1])
    )

    M = sk.util.interpolation.linear_interpolating_matrix(from_grid, to_grid, "zero")
    np.testing.assert_array_almost_equal(
        M @ vals, np.interp(to_grid, from_grid, vals, left=0, right=0)
    )
