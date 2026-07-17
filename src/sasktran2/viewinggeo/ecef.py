from __future__ import annotations

import numpy as np
import pandas as pd

import sasktran2 as sk


def ecef_to_sasktran2_ray(
    observer: np.ndarray,
    look_vector: np.ndarray,
    time: pd.Timestamp,
    geoid: sk.Geodetic | None = None,
    solar_handler: sk.solar.SolarGeometryHandlerBase | None = None,
    ground_elevation: float = 0.0,
) -> sk.GroundViewingSolar | sk.TangentAltitudeSolar:
    """
    Converts an observer, look vector in ECEF coordinates and time to a sasktran2 viewing ray object

    Parameters
    ----------
    observer : np.ndarray
        Observer position in ECEF (ITRF) coordinates
    look_vector : np.ndarray
        Look vector from the observer in ECEF (ITRF) coordinates. The vector
        does not need to be normalized.
    time: pd.Timestamp
        Time of the observation.  Used to lookup the sun position if applicable.
    geoid : sk.Geodetic
        Geodetic object to use for the internal calculations. Defaults to WGS84
    solar_handler : sk.solar.SolarGeometryHandlerBase
        Object used to lookup the solar angles for each ray.  Defaults to None in which case a SolarHandlerForced object is used with 0 solar zenith and 0 solar azimuth

    Returns
    -------
    sk.GroundViewingSolar or sk.TangentAltitudeSolar
        A ground-viewing ray when the look vector intersects the requested
        ground elevation, otherwise a limb-viewing ray.
    """
    if solar_handler is None:
        solar_handler = sk.solar.SolarGeometryHandlerForced(0, 0)

    if geoid is None:
        geoid = sk.WGS84()

    observer = np.asarray(observer, dtype=float)
    look_vector = np.asarray(look_vector, dtype=float)
    look_vector_norm = np.linalg.norm(look_vector)
    if not np.isfinite(look_vector_norm) or look_vector_norm == 0:
        msg = "look_vector must be finite and non-zero"
        raise ValueError(msg)
    look_vector = look_vector / look_vector_norm

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
        viewing_azimuth = -np.rad2deg(
            np.arctan2(
                np.dot(look_vector, geoid.local_west),
                -np.dot(look_vector, geoid.local_south),
            )
        )

        return sk.TangentAltitudeSolar(
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

    # GroundViewingSolar defines the viewing cosine from the ground toward the
    # observer, opposite to the supplied observer-to-ground look vector.
    cos_viewing_zenith = -np.dot(look_vector, geoid.local_up)

    if np.abs(cos_viewing_zenith) > (1 - 1e-8):
        # Basically nadir viewing and so we can't get an observer azimuth
        viewing_azimuth = 0.0
    else:
        # Get the viewing azimuth angle
        viewing_azimuth = -np.rad2deg(
            np.arctan2(
                np.dot(look_vector, geoid.local_west),
                -np.dot(look_vector, geoid.local_south),
            )
        )

    return sk.GroundViewingSolar(
        np.cos(np.deg2rad(solar_zenith)),
        np.deg2rad(solar_azimuth - viewing_azimuth),
        cos_viewing_zenith,
        obs_alt,
    )
