from __future__ import annotations

import datetime as dt

import numpy as np
import pytest
import sasktran2 as sk


def test_exact_grid_lookup():
    density = sk.climatology.atomic_oxygen.number_density([130000.0], 1, -80)

    np.testing.assert_allclose(density, [2.5e16])


def test_latitude_interpolation():
    density = sk.climatology.atomic_oxygen.number_density([130000.0], 1, -75)

    np.testing.assert_allclose(density, [2.65e16])


def test_altitude_interpolation():
    density = sk.climatology.atomic_oxygen.number_density([127500.0], 1, -80)

    np.testing.assert_allclose(density, [3.05e16])


def test_time_interpolation():
    density = sk.climatology.atomic_oxygen.number_density([130000.0], 1.5, -80)

    np.testing.assert_allclose(density, [2.6e16])


def test_date_time_coordinate():
    density = sk.climatology.atomic_oxygen.number_density(
        [130000.0], dt.date(2026, 1, 16), -80
    )

    np.testing.assert_allclose(density, [2.596774193548387e16])


def test_shape_preserved():
    altitudes = np.array([[130000.0, 125000.0], [120000.0, 115000.0]])

    density = sk.climatology.atomic_oxygen.number_density(altitudes, 1, -80)

    assert density.shape == altitudes.shape


def test_default_altitude_boundary_handling():
    density = sk.climatology.atomic_oxygen.number_density([39000.0, 131000.0], 1, -80)

    np.testing.assert_allclose(density, [0.0, 2.5e16])


def test_custom_altitude_boundary_handling():
    density = sk.climatology.atomic_oxygen.number_density(
        [39000.0, 131000.0],
        1,
        -80,
        lower_fill_value=None,
        upper_fill_value=1.0e15,
    )

    np.testing.assert_allclose(density, [2.9e14, 1.0e15])


def test_latitude_range_errors():
    with pytest.raises(ValueError, match="latitude_degrees"):
        sk.climatology.atomic_oxygen.number_density([130000.0], 1, -85)
