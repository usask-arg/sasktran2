from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def test_config_creation():
    """Test basic config creation"""
    config = sk.Config()
    assert config is not None


def test_log_level_default():
    """Test that default log level is Warn"""
    config = sk.Config()
    assert config.log_level == sk.LogLevel.Warn


def test_log_level_set_get():
    """Test setting and getting different log levels"""
    config = sk.Config()

    # Test all log levels
    levels = [
        sk.LogLevel.Trace,
        sk.LogLevel.Debug,
        sk.LogLevel.Info,
        sk.LogLevel.Warn,
        sk.LogLevel.Error,
        sk.LogLevel.Critical,
        sk.LogLevel.Off,
    ]

    for level in levels:
        config.log_level = level
        assert config.log_level == level, f"Failed to set/get log level {level}"


def test_log_level_enum_values():
    """Test that log level enums have expected values"""
    # Test that enums exist and are different
    levels = [
        sk.LogLevel.Trace,
        sk.LogLevel.Debug,
        sk.LogLevel.Info,
        sk.LogLevel.Warn,
        sk.LogLevel.Error,
        sk.LogLevel.Critical,
        sk.LogLevel.Off,
    ]

    # Check all levels are distinct by comparing their string representations
    level_names = [str(level) for level in levels]
    assert len(set(level_names)) == len(
        level_names
    ), "Log levels should all be distinct"

    # Also check that they're not None
    for level in levels:
        assert level is not None


def test_config_multiple_instances():
    """Test that multiple config instances can have different log levels"""
    config1 = sk.Config()
    config2 = sk.Config()

    # Set different log levels
    config1.log_level = sk.LogLevel.Debug
    config2.log_level = sk.LogLevel.Error

    # Verify they're independent
    assert config1.log_level == sk.LogLevel.Debug
    assert config2.log_level == sk.LogLevel.Error


def test_multiple_scatter_source_preserves_public_integer_values():
    assert sk.MultipleScatterSource.DiscreteOrdinates == 0
    assert sk.MultipleScatterSource.SuccessiveOrdersLegacy == 1
    assert sk.MultipleScatterSource.TwoStream == 2
    assert sk.MultipleScatterSource.NoSource == 3
    assert sk.MultipleScatterSource.SuccessiveOrders == 4
    assert not hasattr(sk.MultipleScatterSource, "SuccessiveOrdersCpp")

    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersLegacy
    assert (
        config.multiple_scatter_source
        == sk.MultipleScatterSource.SuccessiveOrdersLegacy
    )

    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
    assert config.multiple_scatter_source == sk.MultipleScatterSource.SuccessiveOrders


def test_num_successive_order_points_round_trip():
    """The Python property maps to the pluralized native config field."""
    config = sk.Config()

    config.num_successive_order_points = 17

    assert config.num_successive_order_points == 17


def test_wavelength_batch_size_round_trip_and_validation():
    config = sk.Config()

    assert config.wavelength_batch_size == 1
    config.wavelength_batch_size = 16
    assert config.wavelength_batch_size == 16

    with pytest.raises(RuntimeError, match="at least 1"):
        config.wavelength_batch_size = 0


def test_los_refraction_tangent_altitude_limit_round_trip_and_validation():
    config = sk.Config()

    assert np.isinf(config.los_refraction_max_tangent_altitude_m)
    config.los_refraction_max_tangent_altitude_m = 25_000.0
    assert config.los_refraction_max_tangent_altitude_m == 25_000.0
    config.los_refraction_max_tangent_altitude_m = np.inf
    assert np.isinf(config.los_refraction_max_tangent_altitude_m)

    for invalid in (-1.0, np.nan):
        with pytest.raises(RuntimeError, match="non-negative or infinity"):
            config.los_refraction_max_tangent_altitude_m = invalid


def test_cpp_successive_orders_controls_round_trip_and_validation():
    config = sk.Config()

    assert config.num_successive_orders_iterations == 50
    assert config.num_successive_orders_incoming == 110
    assert config.num_successive_orders_outgoing == 110
    assert config.successive_orders_relative_tolerance == 1.0e-6
    assert config.successive_orders_absolute_tolerance == 1.0e-12
    assert config.successive_orders_anderson_depth == 3
    assert config.successive_orders_damping == 1.0

    config.num_successive_orders_iterations = 31
    config.num_successive_orders_incoming = 74
    config.num_successive_orders_outgoing = 86
    config.successive_orders_relative_tolerance = 2.0e-7
    config.successive_orders_absolute_tolerance = 3.0e-13
    config.successive_orders_anderson_depth = 4
    config.successive_orders_damping = 0.85

    assert config.num_successive_orders_iterations == 31
    assert config.num_successive_orders_incoming == 74
    assert config.num_successive_orders_outgoing == 86
    assert config.successive_orders_relative_tolerance == 2.0e-7
    assert config.successive_orders_absolute_tolerance == 3.0e-13
    assert config.successive_orders_anderson_depth == 4
    assert config.successive_orders_damping == 0.85

    for invalid in (-1.0, np.nan, np.inf):
        with pytest.raises(RuntimeError, match="finite and non-negative"):
            config.successive_orders_relative_tolerance = invalid
        with pytest.raises(RuntimeError, match="finite and non-negative"):
            config.successive_orders_absolute_tolerance = invalid

    for invalid in (0.0, -0.1, 1.1, np.nan, np.inf):
        with pytest.raises(RuntimeError, match="interval"):
            config.successive_orders_damping = invalid


def test_cpp_successive_orders_altitude_grid_round_trip_and_validation():
    config = sk.Config()
    assert config.successive_orders_altitude_grid_m is None

    expected = np.array([500.0, 2_500.0, 10_000.0])
    config.successive_orders_altitude_grid_m = expected
    np.testing.assert_array_equal(config.successive_orders_altitude_grid_m, expected)

    returned = config.successive_orders_altitude_grid_m
    returned[0] = -1.0
    np.testing.assert_array_equal(config.successive_orders_altitude_grid_m, expected)

    for invalid in (
        np.array([[1_000.0, 2_000.0]]),
        np.array([1_000.0, 1_000.0]),
        np.array([2_000.0, 1_000.0]),
        np.array([1_000.0, np.nan]),
    ):
        with pytest.raises(ValueError, match="successive_orders_altitude_grid_m"):
            config.successive_orders_altitude_grid_m = invalid

    config.successive_orders_altitude_grid_m = np.array([])
    assert config.successive_orders_altitude_grid_m is None
    config.successive_orders_altitude_grid_m = expected
    config.successive_orders_altitude_grid_m = None
    assert config.successive_orders_altitude_grid_m is None


def test_cpp_successive_orders_horizontal_angle_grid_round_trip_and_validation():
    config = sk.Config()
    assert config.successive_orders_horizontal_angle_grid_radians is None

    expected = np.array([-0.3, -0.05, 0.1, 0.4])
    config.successive_orders_horizontal_angle_grid_radians = expected
    np.testing.assert_array_equal(
        config.successive_orders_horizontal_angle_grid_radians, expected
    )

    returned = config.successive_orders_horizontal_angle_grid_radians
    returned[0] = -1.0
    np.testing.assert_array_equal(
        config.successive_orders_horizontal_angle_grid_radians, expected
    )

    for invalid in (
        np.array([[-0.2, 0.2]]),
        np.array([-0.1, -0.1]),
        np.array([0.2, -0.2]),
        np.array([0.0, np.nan]),
    ):
        with pytest.raises(
            ValueError, match="successive_orders_horizontal_angle_grid_radians"
        ):
            config.successive_orders_horizontal_angle_grid_radians = invalid

    config.successive_orders_horizontal_angle_grid_radians = np.array([])
    assert config.successive_orders_horizontal_angle_grid_radians is None
    config.successive_orders_horizontal_angle_grid_radians = expected
    config.successive_orders_horizontal_angle_grid_radians = None
    assert config.successive_orders_horizontal_angle_grid_radians is None
