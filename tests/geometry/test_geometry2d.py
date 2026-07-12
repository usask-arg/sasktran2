from __future__ import annotations

import numpy as np
import pytest

import sasktran2 as sk


def geometry2d(
    interpolation_method: sk.InterpolationMethod = sk.InterpolationMethod.LinearInterpolation,
) -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.2,
        earth_radius_m=6_371_000.0,
        altitude_grid_m=np.array([0.0, 10_000.0, 20_000.0]),
        horizontal_angle_grid_radians=np.array([-0.3, 0.0, 0.4, 0.7]),
        interpolation_method=interpolation_method,
    )


def test_geometry2d_constructor_grids_shape_and_indexing():
    geometry = geometry2d()

    np.testing.assert_array_equal(
        geometry.altitudes(), np.array([0.0, 10_000.0, 20_000.0])
    )
    np.testing.assert_array_equal(
        geometry.horizontal_angles(), np.array([-0.3, 0.0, 0.4, 0.7])
    )
    assert geometry.shape == (4, 3)

    for horizontal_index in range(4):
        for altitude_index in range(3):
            assert geometry.location_index(altitude_index, horizontal_index) == (
                horizontal_index * 3 + altitude_index
            )

    with pytest.raises(RuntimeError):
        geometry.location_index(3, 0)
    with pytest.raises(RuntimeError):
        geometry.location_index(0, 4)


@pytest.mark.parametrize(
    "interpolation_method",
    [
        sk.InterpolationMethod.LinearInterpolation,
        sk.InterpolationMethod.ShellInterpolation,
        sk.InterpolationMethod.LowerInterpolation,
    ],
)
def test_geometry2d_supports_altitude_interpolation_modes(interpolation_method):
    geometry = geometry2d(interpolation_method)
    assert geometry.shape == (4, 3)


def test_geometry2d_owns_contiguous_copies_of_input_grids():
    altitude_storage = np.array([-1.0, 0.0, -1.0, 10_000.0, -1.0, 20_000.0])
    horizontal_storage = np.array([-0.3, 99.0, 0.0, 99.0, 0.4, 99.0])
    altitudes = altitude_storage[1::2]
    horizontal_angles = horizontal_storage[::2]
    assert not altitudes.flags.c_contiguous
    assert not horizontal_angles.flags.c_contiguous

    geometry = sk.Geometry2D(
        0.6,
        0.2,
        6_371_000.0,
        altitudes,
        horizontal_angles,
    )
    altitudes[:] = 100.0
    horizontal_angles[:] = 0.1

    np.testing.assert_array_equal(
        geometry.altitudes(), np.array([0.0, 10_000.0, 20_000.0])
    )
    np.testing.assert_array_equal(
        geometry.horizontal_angles(), np.array([-0.3, 0.0, 0.4])
    )

    returned_altitudes = geometry.altitudes()
    returned_altitudes[:] = -1.0
    np.testing.assert_array_equal(
        geometry.altitudes(), np.array([0.0, 10_000.0, 20_000.0])
    )


@pytest.mark.parametrize(
    ("altitudes", "horizontal_angles"),
    [
        ([0.0], [-0.2, 0.2]),
        ([10_000.0, 0.0], [-0.2, 0.2]),
        ([0.0, 10_000.0], [0.0]),
        ([0.0, 10_000.0], [0.2, -0.2]),
        ([0.0, 10_000.0], [-np.pi / 2, np.pi / 2]),
        ([0.0, np.nan], [-0.2, 0.2]),
        ([0.0, 10_000.0], [-0.2, np.nan]),
    ],
)
def test_geometry2d_rejects_invalid_grids(altitudes, horizontal_angles):
    with pytest.raises(RuntimeError, match="Failed to create Geometry2D"):
        sk.Geometry2D(
            0.6,
            0.2,
            6_371_000.0,
            np.asarray(altitudes),
            np.asarray(horizontal_angles),
        )
