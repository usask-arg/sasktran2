from __future__ import annotations

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
