from __future__ import annotations

import numpy as np

from sasktran2._core_rust import PyGrid

from .basis import Basis, Delta, Rectangle, Triangle


def _left_right_splits(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Splits an array into left and right boundaries for grid points.

    Parameters:
        x (np.ndarray): Array of grid points.

    Returns:
        tuple[np.ndarray, np.ndarray]: Left and right boundaries.
    """
    left = np.zeros_like(x)
    right = np.zeros_like(x)

    left[0] = x[0]
    right[-1] = x[-1]

    left[1:] = (x[:-1] + x[1:]) / 2
    right[:-1] = (x[:-1] + x[1:]) / 2

    return left, right


def _left_right_triangle_splits(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Splits an array into left and right boundaries for grid points.

    Parameters:
        x (np.ndarray): Array of grid points.

    Returns:
        tuple[np.ndarray, np.ndarray]: Left and right boundaries.
    """
    left = np.zeros_like(x)
    right = np.zeros_like(x)

    left[0] = x[0]
    right[-1] = x[-1]

    left[1:] = x[:-1]
    right[:-1] = x[1:]

    return left, right


class Grid:
    def __init__(self, basis_list: list[Basis]):
        self._grid = PyGrid([b._internal_object() for b in basis_list])

    @classmethod
    def from_rectangles(cls, grid_points: np.ndarray):
        gp = np.atleast_1d(grid_points)

        left, right = _left_right_splits(gp)

        basis_list = [Rectangle(le, r) for le, r in zip(left, right, strict=False)]

        return cls(basis_list)

    @classmethod
    def from_deltas(cls, grid_points: np.ndarray):
        gp = np.atleast_1d(grid_points)
        basis_list = [Delta(x) for x in gp]
        return cls(basis_list)

    @classmethod
    def from_triangles(cls, grid_points: np.ndarray):
        gp = np.atleast_1d(grid_points)
        left, right = _left_right_triangle_splits(gp)
        basis_list = [
            Triangle(le, r, c)
            for le, r, c in zip(left, right, grid_points, strict=False)
        ]

        return cls(basis_list)

    def _internal_object(self) -> PyGrid:
        return self._grid

    def mapping_to(self, grid: Grid, normalize=True) -> np.ndarray:
        mapping = self._grid.mapping_to(grid._internal_object())

        if normalize:
            mapping /= np.sum(mapping, axis=1, keepdims=True)

        return mapping
