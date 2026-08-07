from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def test_config_creation():
    """Test basic config creation"""
    config = sk.Config()
    assert config is not None


def test_two_stream_backend_roundtrip():
    config = sk.Config()
    assert config.two_stream_backend == sk.TwoStreamBackend.Rust
    config.two_stream_backend = sk.TwoStreamBackend.Cpp
    assert config.two_stream_backend == sk.TwoStreamBackend.Cpp


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
    assert sk.MultipleScatterSource.SuccessiveOrders == 1
    assert sk.MultipleScatterSource.TwoStream == 2
    assert sk.MultipleScatterSource.NoSource == 3
    assert sk.MultipleScatterSource.SuccessiveOrdersRust == 4


def test_num_successive_order_points_round_trip():
    config = sk.Config()
    assert config.num_successive_order_points == -1

    config.num_successive_order_points = 17
    assert config.num_successive_order_points == 17

    config.num_successive_order_points = -1
    assert config.num_successive_order_points == -1

    for invalid in (-2, 0, 1):
        with pytest.raises(ValueError, match="must be -1 or at least 2"):
            config.num_successive_order_points = invalid


def test_wavelength_batch_size_round_trip_and_validation():
    config = sk.Config()

    assert config.wavelength_batch_size == 1
    config.wavelength_batch_size = 16
    assert config.wavelength_batch_size == 16

    with pytest.raises(RuntimeError, match="at least 1"):
        config.wavelength_batch_size = 0


def test_rust_successive_orders_config_round_trip():
    config = sk.Config()

    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersRust
    config.successive_orders_max_iterations = 31
    config.successive_orders_relative_tolerance = 2.0e-7
    config.successive_orders_absolute_tolerance = 3.0e-13
    config.successive_orders_anderson_depth = 4
    config.successive_orders_damping = 0.85

    assert (
        config.multiple_scatter_source == sk.MultipleScatterSource.SuccessiveOrdersRust
    )
    assert config.successive_orders_max_iterations == 31
    assert config.successive_orders_relative_tolerance == 2.0e-7
    assert config.successive_orders_absolute_tolerance == 3.0e-13
    assert config.successive_orders_anderson_depth == 4
    assert config.successive_orders_damping == 0.85


def test_successive_orders_initialization_round_trip_and_legacy_alias():
    config = sk.Config()

    assert (
        config.successive_orders_initialization
        == sk.SuccessiveOrdersInitialization.NoInitialization
    )
    config.successive_orders_initialization = (
        sk.SuccessiveOrdersInitialization.TwoStream
    )
    assert (
        config.successive_orders_initialization
        == sk.SuccessiveOrdersInitialization.TwoStream
    )
    assert not config.init_successive_orders_with_discrete_ordinates

    config.init_successive_orders_with_discrete_ordinates = True
    assert (
        config.successive_orders_initialization
        == sk.SuccessiveOrdersInitialization.DiscreteOrdinates
    )
    config.init_successive_orders_with_discrete_ordinates = False
    assert (
        config.successive_orders_initialization
        == sk.SuccessiveOrdersInitialization.NoInitialization
    )


def test_successive_orders_altitude_grid_m_round_trip_and_validation():
    config = sk.Config()
    assert config.successive_orders_altitude_grid_m is None
    expected = np.array([500.0, 2_500.0, 10_000.0])
    config.successive_orders_altitude_grid_m = expected
    np.testing.assert_array_equal(config.successive_orders_altitude_grid_m, expected)

    returned = config.successive_orders_altitude_grid_m
    returned[0] = -1.0
    np.testing.assert_array_equal(config.successive_orders_altitude_grid_m, expected)

    config.successive_orders_altitude_grid_m = np.array([4_000.0])
    np.testing.assert_array_equal(
        config.successive_orders_altitude_grid_m, np.array([4_000.0])
    )

    for invalid in (
        np.array([[1_000.0, 2_000.0]]),
        np.array([1_000.0, 1_000.0]),
        np.array([2_000.0, 1_000.0]),
        np.array([1_000.0, np.nan]),
    ):
        with pytest.raises(ValueError, match="successive_orders_altitude_grid_m"):
            config.successive_orders_altitude_grid_m = invalid

    config.successive_orders_altitude_grid_m = None
    assert config.successive_orders_altitude_grid_m is None


def test_successive_orders_source_profiles_aliases_num_sza():
    config = sk.Config()
    config.num_successive_orders_source_profiles = 7
    assert config.num_successive_orders_source_profiles == 7
    assert config.num_sza == 7

    config.num_sza = 3
    assert config.num_successive_orders_source_profiles == 3
