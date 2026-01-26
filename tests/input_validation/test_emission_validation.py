from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _default_settings(geometry_type=sk.GeometryType.Spherical):
    """Create default settings for testing."""
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.num_streams = 2

    model_geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.arange(0, 65001, 1000),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=geometry_type,
    )

    viewing_geo = sk.ViewingGeometry()

    viewing_geo.add_ray(sk.GroundViewingSolar(0.6, 0.0, 1.0, 200000))

    return config, model_geometry, viewing_geo


def test_emission_do_requires_single_scatter_do():
    """Test that emission_source=DiscreteOrdinates requires single_scatter_source=DiscreteOrdinates.

    Note: Due to how exceptions are handled in the Rust FFI layer, the exception
    may not propagate to Python, but the validation logic logs a critical message.
    The capfd fixture captures C++ stdout/stderr at the file descriptor level.
    """
    config, model_geometry, viewing_geo = _default_settings()

    config.emission_source = sk.EmissionSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.Exact  # Wrong setting

    try:
        _ = sk.Engine(config, model_geometry, viewing_geo)
        msg = "Expected exception was not raised."
        raise AssertionError(msg)
    except RuntimeError:
        # Should fail
        pass


def test_emission_do_requires_multiple_scatter_do():
    """Test that emission_source=DiscreteOrdinates requires multiple_scatter_source=DiscreteOrdinates.

    Note: Due to how exceptions are handled in the Rust FFI layer, the exception
    may not propagate to Python, but the validation logic logs a critical message.
    The capfd fixture captures C++ stdout/stderr at the file descriptor level.
    """
    config, model_geometry, viewing_geo = _default_settings()

    config.emission_source = sk.EmissionSource.DiscreteOrdinates
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource  # Wrong setting

    try:
        _ = sk.Engine(config, model_geometry, viewing_geo)
        msg = "Expected exception was not raised."
        raise AssertionError(msg)
    except RuntimeError:
        # Should fail
        pass


def test_emission_do_valid_with_plane_parallel():
    """Test that emission_source=DiscreteOrdinates works with plane parallel geometry."""
    config, model_geometry, viewing_geo = _default_settings(
        geometry_type=sk.GeometryType.PlaneParallel
    )

    config.emission_source = sk.EmissionSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates

    # Should not raise an exception
    _ = sk.Engine(config, model_geometry, viewing_geo)


def test_emission_do_valid_with_spherical():
    """Test that emission_source=DiscreteOrdinates works with spherical geometry."""
    config, model_geometry, viewing_geo = _default_settings(
        geometry_type=sk.GeometryType.Spherical
    )

    config.emission_source = sk.EmissionSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates

    # Should not raise an exception
    _ = sk.Engine(config, model_geometry, viewing_geo)


def test_emission_do_valid_with_pseudospherical():
    """Test that emission_source=DiscreteOrdinates works with pseudospherical geometry."""
    config, model_geometry, viewing_geo = _default_settings(
        geometry_type=sk.GeometryType.PseudoSpherical
    )

    config.emission_source = sk.EmissionSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates

    # Should not raise an exception
    _ = sk.Engine(config, model_geometry, viewing_geo)
