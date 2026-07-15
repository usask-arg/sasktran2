from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_atmosphere_construction():
    config = sk.Config()

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    _ = sk.Atmosphere(geometry, config, numwavel=17)


def test_atmosphere_vector_construction():
    config = sk.Config()
    config.num_stokes = 3

    altitude_grid = np.arange(0, 100001, 1000)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(geometry, config, numwavel=17)

    _ = atmosphere.storage


def test_1d_native_access_preserves_existing_zero_copy_state_objects():
    config = sk.Config()
    geometry = sk.Geometry1D(
        0.6,
        0,
        6_327_000,
        np.arange(0.0, 4_001.0, 1_000.0),
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(geometry, config, numwavel=2)
    pressure = np.linspace(100_000.0, 60_000.0, 5)
    temperature = np.linspace(280.0, 240.0, 5)
    atmosphere.pressure_pa = pressure
    atmosphere.temperature_k = temperature

    assert atmosphere._native_state("pressure_pa") is pressure
    assert atmosphere._native_state("temperature_k") is temperature
    assert atmosphere._native_state_equation() is atmosphere.state_equation
