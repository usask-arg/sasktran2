from __future__ import annotations

import numpy as np


def linear_interpolating_matrix(
    from_grid: np.array, to_grid: np.array, out_of_bounds_mode: str
) -> np.ndarray:
    """
    Constructs an interpolating matrix that can be used to interpolate from one grid to another.

    Parameters
    ----------
    from_grid : np.array
        Grid with the original values. Must be sorted in ascending order
    to_grid : np.array
        Grid we are going to interpolate to. Must be sorted in ascending order
    out_of_bounds_mode : str
        One of ["zero", "extend"], "zero" will set the result to 0 outside the bounds of the original grid,
        "extend" will extend the lowest and highest values

    Returns
    -------
    np.ndarray
        A matrix such that M @ (values on from_grid) = (values on to_grid)
    """
    M = np.zeros((len(to_grid), len(from_grid)))

    for idx, ele in enumerate(to_grid):
        if ele < from_grid[0] or ele > from_grid[-1]:
            if out_of_bounds_mode == "zero":
                continue
            if out_of_bounds_mode == "extend":
                if ele < from_grid[0]:
                    M[idx, 0] = 1
                else:
                    M[idx, -1] = 1
                continue
            msg = f"Invalid out_of_bounds_mode {out_of_bounds_mode}"
            raise ValueError(msg)

        # else we are inside the grid
        idx_above = np.nonzero(from_grid > ele)[0]
        idx_above = len(from_grid) if len(idx_above) == 0 else idx_above[0]

        interp_ele = ele
        if idx_above == 0:
            idx_above = 1
            interp_ele = from_grid[0]

        if idx_above == len(from_grid):
            M[idx, idx_above - 1] = 1
        else:
            w = (from_grid[idx_above] - interp_ele) / (
                from_grid[idx_above] - from_grid[idx_above - 1]
            )
            M[idx, idx_above] = 1 - w
            M[idx, idx_above - 1] = w

    return M
