from __future__ import annotations

import numpy as np
from sasktran2.basis import Delta, Grid, Rectangle


def test_grid_construction():
    rect1 = Rectangle(3.0, 4.0)
    rect2 = Rectangle(5.0, 6.0)

    _ = Grid([rect1, rect2])


def test_mapping():
    rect1 = Rectangle(3.0, 4.0)
    rect2 = Rectangle(5.0, 6.0)

    grid1 = Grid([rect1, rect2])

    delta = Delta(5.5)

    grid2 = Grid([delta])

    result = grid1.mapping_to(grid2)

    pass

def test_grid_from_rectangles():
    grid_points = np.arange(0, 10, 0.01)

    grid1 = Grid.from_rectangles(grid_points)

    np.testing.assert_array_equal(
        grid1.mapping_to(grid1),
        np.eye(len(grid_points), dtype=float)
    )