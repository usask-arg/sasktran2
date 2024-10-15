import numpy as np
import pandas as pd
import sasktran2 as sk


def test_forced_handler():
    forced_solar_zenith = 45.0
    forced_solar_azimuth = 180.0

    handler = sk.solar.SolarHandlerForced(forced_solar_zenith, forced_solar_azimuth)

    latitude = 0
    longitude = 0

    solar_zenith, solar_azimuth = handler.target_solar_angles(
        latitude, longitude, 0.0, pd.Timestamp.now()
    )

    np.testing.assert_allclose(solar_zenith, forced_solar_zenith)
    np.testing.assert_allclose(solar_azimuth, forced_solar_azimuth)
