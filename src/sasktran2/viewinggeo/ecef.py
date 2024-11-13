import numpy as np
import pandas as pd

from sasktran2 import (
    Geodetic,
    GroundViewingSolar,
    TangentAltitudeSolar,
    ViewingGeometryBase,
)
from sasktran2.geodetic import WGS84
from sasktran2.solar import SolarGeometryHandlerBase, SolarGeometryHandlerForced


def ecef_to_sasktran2_ray(
    observer: np.array,
    look_vector: np.array,
    time: pd.Timestamp,
    geoid: Geodetic | None = None,
    solar_handler: SolarGeometryHandlerBase | None = None,
    ground_elevation: float = 0.0,
) -> ViewingGeometryBase:
    """
    Converts an observer, look vector in ECEF coordinates and time to a sasktran2 viewing ray object

    Parameters
    ----------
    observer : np.array
        Observer position in ECEF (ITRF) coordinates
    look_vector : np.array
        Local lok vector in ECEF (ITRF) coordinates
    time: pd.Timestamp
        Time of the observation.  Used to lookup the sun position if applicable.
    geoid : sk.Geodetic
        Geodetic object to use for the internal calculations. Defaults to WGS84
    solar_handler : sk.solar.SolarGeometryHandlerBase
        Object used to lookup the solar angles for each ray.  Defaults to None in which case a SolarHandlerForced object is used with 0 solar zenith and 0 solar azimuth

    Returns
    -------
    sk.ViewingGeometryBase
    """
    if solar_handler is None:
        solar_handler = SolarGeometryHandlerForced(0, 0)

    if geoid is None:
        geoid = WGS84()

    # Start by determining the observer's geodetic coordinates
    geoid.from_xyz(observer)

    # obs_lat = geoid.latitude
    # obs_lon = geoid.longitude
    obs_alt = geoid.altitude

    # Now start by checking if this is a limb viewing ray
    geoid.from_tangent_point(observer, look_vector)

    if geoid.altitude > ground_elevation:
        # Limb viewing
        solar_zenith, solar_azimuth = solar_handler.target_solar_angles(
            geoid.latitude, geoid.longitude, geoid.altitude, time
        )

        # Get the viewing azimuth angle
        viewing_azimuth = -np.arctan2(
            np.dot(look_vector, geoid.local_west),
            -np.dot(look_vector, geoid.local_south),
        )

        return TangentAltitudeSolar(
            geoid.altitude,
            np.deg2rad(solar_azimuth - viewing_azimuth),
            obs_alt,
            np.cos(np.deg2rad(solar_zenith)),
        )
    # ground viewing
    geoid.from_xyz(
        geoid.altitude_intercepts(ground_elevation, observer, look_vector)[0]
    )
    solar_zenith, solar_azimuth = solar_handler.target_solar_angles(
        geoid.latitude, geoid.longitude, geoid.altitude, time
    )

    cos_viewing_zenith = np.dot(look_vector, geoid.local_up)

    if np.abs(cos_viewing_zenith) > (1 - 1e-8):
        # Basically nadir viewing and so we can't get an observer azimuth
        viewing_azimuth = 0.0
    else:
        # Get the viewing azimuth angle
        viewing_azimuth = -np.arctan2(
            np.dot(look_vector, geoid.local_west),
            -np.dot(look_vector, geoid.local_south),
        )

    return GroundViewingSolar(
        np.cos(np.deg2rad(solar_zenith)),
        np.deg2rad(solar_azimuth - viewing_azimuth),
        cos_viewing_zenith,
        obs_alt,
    )
