from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
from scipy import interpolate

from sasktran2._core import (
    Geometry1D,
    GeometryType,
    InterpolationMethod,
    TangentAltitudeSolar,
    ViewingGeometry,
)
from sasktran2.geodetic import WGS84
from sasktran2.solar import SolarGeometryHandlerBase

from .base import ViewingGeometryContainer


class LimbVertical(ViewingGeometry, ViewingGeometryContainer):
    def __init__(
        self,
        solar_handler: SolarGeometryHandlerBase,
        tangent_altitudes: np.array,
        tangent_latitude: float | np.ndarray,
        tangent_longitude: float | np.ndarray,
        time: pd.Timestamp | np.ndarray,
        observer_altitude: float | np.ndarray,
        observer_latitude: float | np.ndarray,
        observer_longitude: float | np.ndarray,
        reference_altitude: float = 25000,
    ):
        """
        A Limb Vertical image is a geometry container where the observer is looking through the limb of the atmosphere, at a set
        of tangent altitudes.

        We recommend constructing this class through one of the various `from_` class methods.  No validation checking
        is done through this constructor.

        Parameters
        ----------
        solar_handler : sk.solar.SolarGeometryHandlerBase
            Solar geometry handler to use for calculating solar angles.
        tangent_altitudes : np.array
            Altitudes of the tangent points.
        tangent_latitude : float | np.ndarray
            Latitude of the tangent points. Can be a scalar or same size as `tangent_altitudes`.
        tangent_longitude : float | np.ndarray
            Longitude of the tangent points. Can be a scalar or same size as `tangent_altitudes`.
        time : pd.Timestamp | np.ndarray
            Time of the observation. Can be a scalar or same size as `tangent_altitudes`.
        observer_altitude : float | np.ndarray
            Altitude of the observer. Can be a scalar or same size as `tangent_altitudes`.
        observer_latitude : float | np.ndarray
            Latitude of the observer. Can be a scalar or same size as `tangent_altitudes`.
        observer_longitude : float | np.ndarray
            Longitude of the observer. Can be a scalar or same size as `tangent_altitudes`.
        reference_altitude : float, optional
            The tangent altitude where reference parameters such as earth radius
            are defined at, by default 25000.
        """
        self._tangent_altitudes = np.array(tangent_altitudes)

        # Determine the expected length
        n = len(self._tangent_altitudes)

        # Helper function to convert scalars to arrays
        def to_array(param):
            if np.isscalar(param) | isinstance(param, pd.Timestamp):
                return np.full(n, param)
            param = np.array(param)
            if len(param) != n:
                msg = "Parameter length mismatch."
                raise ValueError(msg)
            return param

        self._tangent_latitude = to_array(tangent_latitude)
        self._tangent_longitude = to_array(tangent_longitude)
        self._time = to_array(time)
        self._observer_altitude = to_array(observer_altitude)
        self._observer_latitude = to_array(observer_latitude)
        self._observer_longitude = to_array(observer_longitude)

        self._reference_altitude = reference_altitude

        ViewingGeometry.__init__(self)

        tangent_geo = WGS84()
        observer_geo = WGS84()

        viewing_zenith = np.zeros(n)
        viewing_azimuth = np.zeros(n)

        self._cos_sza = np.zeros(n)
        self._earth_radius = np.zeros(n)
        self._solar_azimuth = np.zeros(n)
        self._observer_azimuth = np.zeros(n)

        for i, (alt, lat, lon, t, obs_alt, obs_lat, obs_lon) in enumerate(
            zip(
                tangent_altitudes,
                self._tangent_latitude,
                self._tangent_longitude,
                self._time,
                self._observer_altitude,
                self._observer_latitude,
                self._observer_longitude,
                strict=False,
            )
        ):
            tangent_geo.from_lat_lon_alt(lat, lon, alt)
            observer_geo.from_lat_lon_alt(obs_lat, obs_lon, obs_alt)

            solar_zenith, solar_azimuth = solar_handler.target_solar_angles(
                lat, lon, alt, t
            )

            # Need to calculate the relative azimuth
            lv = tangent_geo.location - observer_geo.location
            lv /= np.linalg.norm(lv)

            # Need angle in the (north, west) plane measured clockwise
            observer_azimuth = -np.arctan2(
                np.dot(lv, tangent_geo.local_west), -np.dot(lv, tangent_geo.local_south)
            )

            # Now observer_azimuth is pointing away from the observer, and solar azimuth is pointing towards the sun
            # If they are equal this is forward scatter as expected,

            self.add_ray(
                TangentAltitudeSolar(
                    tangent_altitude_m=alt,
                    relative_azimuth=np.deg2rad(solar_azimuth - observer_azimuth),
                    observer_altitude_m=obs_alt,
                    cos_sza=np.cos(np.deg2rad(solar_zenith)),
                )
            )

            self._observer_azimuth[i] = observer_azimuth
            self._solar_azimuth[i] = solar_azimuth

            # Assign the viewing zenith and azimuthal angles
            viewing_zenith[i] = np.rad2deg(np.arccos(np.dot(lv, observer_geo.local_up)))

            viewing_azimuth[i] = -np.arctan2(
                np.dot(lv, observer_geo.local_west),
                -np.dot(lv, observer_geo.local_south),
            )

            # Store parameters
            self._cos_sza[i] = np.cos(np.deg2rad(solar_zenith))
            self._earth_radius[i] = np.linalg.norm(
                tangent_geo.location - alt * tangent_geo.local_up
            )

        geometry_ds = xr.Dataset(
            {
                "tangent_altitude": (["los"], tangent_altitudes),
                "tangent_longitude": (
                    ["los"],
                    np.ones_like(tangent_altitudes) * tangent_longitude,
                ),
                "tangent_latitude": (
                    ["los"],
                    np.ones_like(tangent_altitudes) * tangent_latitude,
                ),
                "time": (["los"], self._time.astype(np.datetime64)),
                "observer_altitude": (["los"], self._observer_altitude),
                "observer_latitude": (["los"], self._observer_latitude),
                "observer_longitude": (["los"], self._observer_longitude),
                "tangent_cos_sza": (["los"], self._cos_sza),
                "tangent_solar_azimuth": (["los"], self._solar_azimuth),
                "tangent_observer_azimuth": (["los"], self._observer_azimuth),
                "viewing_zenith": (["los"], viewing_zenith),
                "viewing_azimuth": (["los"], viewing_azimuth),
            },
        )

        ViewingGeometryContainer.__init__(self, geometry_ds)

    def recommended_cos_sza(self) -> float:
        """
        Returns the cosine of the solar zenith angle at the reference altitude.

        Returns
        -------
        float
        """
        return interpolate.interp1d(self._tangent_altitudes, self._cos_sza)(
            self._reference_altitude
        )

    def recommended_earth_radius(self) -> float:
        """
        Returns the earth radius at the reference altitude.

        Returns
        -------
        float
        """
        return interpolate.interp1d(self._tangent_altitudes, self._earth_radius)(
            self._reference_altitude
        )

    def model_geometry(self, altitude_grid_m: np.ndarray) -> Geometry1D:
        return Geometry1D(
            self.recommended_cos_sza(),
            0.0,
            self.recommended_earth_radius(),
            altitude_grid_m,
            InterpolationMethod.LinearInterpolation,
            GeometryType.Spherical,
        )

    @classmethod
    def from_tangent_parameters(
        cls,
        solar_handler: SolarGeometryHandlerBase,
        tangent_altitudes: np.ndarray,
        tangent_latitude: float,
        tangent_longitude: float,
        time: pd.Timestamp,
        observer_altitude: float,
        viewing_azimuth: float,
        reference_altitude: float = 25000,
        forced_constant_tangent: bool = False,
    ):
        """
        Initialize a LimbVertical object from a set of tangent parameters assuming a vertical image is taken. From
        a single observer location.

        Parameters
        ----------
        solar_handler : sk.solar.SolarGeometryHandlerBase
            Solar geometry handler to use for calculating solar angles.
        tangent_altitudes : np.ndarray
            The altitudes of the tangent points.
        tangent_latitude : float
            Latitude of the tangent point at reference_altitude
        tangent_longitude : float
            Longitude of the tangent piont at reference_altitude
        time : pd.Timestamp
            Time of the observation
        observer_altitude : float
            Altitude of the observer
        viewing_azimuth : float
            Azimuth of the viewing direction in degrees clockwise from north at the tangent point.
        reference_altitude : float, optional
            Altitude the tangent point parameters are defined at, by default 25000
        forced_constant_tangent : bool, optional
            If true, then variation in latitude/longitude is not taken into account as a function of
            tangent altitude, by default False
        """
        tangent_geo = WGS84()

        tangent_geo.from_lat_lon_alt(
            tangent_latitude, tangent_longitude, reference_altitude
        )

        lv = -tangent_geo.local_south * np.cos(
            np.deg2rad(viewing_azimuth)
        ) - tangent_geo.local_west * np.sin(np.deg2rad(viewing_azimuth))

        # Now we have to find the observer location
        observer_geo = WGS84()
        observer_geo.from_xyz(
            observer_geo.altitude_intercepts(
                observer_altitude, tangent_geo.location, lv
            )[0]
        )

        # Then we can calculate all of the actual tangent latitudes and longitudes
        actual_tangent_latitude = np.zeros(len(tangent_altitudes))
        actual_tangent_longitude = np.zeros(len(tangent_altitudes))

        if forced_constant_tangent:
            actual_tangent_latitude = tangent_latitude
            actual_tangent_longitude = tangent_longitude
        else:
            for i, alt in enumerate(tangent_altitudes):
                tangent_geo.from_tangent_altitude(alt, observer_geo.location, lv)
                actual_tangent_latitude[i] = tangent_geo.latitude
                actual_tangent_longitude[i] = tangent_geo.longitude

        return cls(
            solar_handler,
            tangent_altitudes,
            actual_tangent_latitude,
            actual_tangent_longitude,
            time,
            observer_altitude,
            observer_geo.latitude,
            observer_geo.longitude,
            reference_altitude,
        )
