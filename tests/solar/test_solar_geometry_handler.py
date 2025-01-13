from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import sasktran2 as sk


def test_forced_handler():
    forced_solar_zenith = 45.0
    forced_solar_azimuth = 180.0

    handler = sk.solar.SolarGeometryHandlerForced(
        forced_solar_zenith, forced_solar_azimuth
    )

    latitude = 0
    longitude = 0

    solar_zenith, solar_azimuth = handler.target_solar_angles(
        latitude, longitude, 0.0, pd.Timestamp.now()
    )

    np.testing.assert_allclose(solar_zenith, forced_solar_zenith)
    np.testing.assert_allclose(solar_azimuth, forced_solar_azimuth)


def test_astropy_handler():
    pytest.importorskip("astropy")

    handler = sk.solar.SolarGeometryHandlerAstropy()

    latitude = 39.74
    longitude = -104.99

    time = pd.Timestamp("2024-11-12 20:00:00")

    solar_zenith, solar_azimuth = handler.target_solar_angles(
        latitude, longitude, 0.0, time
    )

    # Values calculated from NOAA calculator, I assume astropy is more accurate
    np.testing.assert_allclose(solar_zenith, 90.0 - 29.65, atol=1e-1)
    np.testing.assert_allclose(solar_azimuth, 200.82, atol=1e-1)
