from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_us76_log_pressure_extrapolates_above_top_limit():
    altitudes_m = np.array([80_000.0, 90_000.0, 100_000.0])
    geometry = sk.Geometry1D(
        0.5,
        0.0,
        6371000.0,
        altitudes_m,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PlaneParallel,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        sk.Config(),
        np.array([500.0]),
        calculate_derivatives=False,
    )

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    log_pressure_table = np.log(sk.climatology.us76._PRESSURE * 1e4)
    top_slope = (log_pressure_table[-1] - log_pressure_table[-2]) / (
        sk.climatology.us76._ALTS[-1] - sk.climatology.us76._ALTS[-2]
    )
    expected_log_pressure = log_pressure_table[-1] + top_slope * (
        altitudes_m - sk.climatology.us76._ALTS[-1]
    )

    np.testing.assert_allclose(
        np.log(atmosphere.pressure_pa),
        expected_log_pressure,
    )
    assert np.all(np.diff(atmosphere.pressure_pa) < 0.0)
