from __future__ import annotations

import numpy as np

from sasktran2._core_rust import PyGeodetic


class Geodetic:
    _internal: PyGeodetic

    def __init__(self, radius: float, flattening: float):
        """
        A geodetic object that can be used to represent a location on the
        surface of the Earth. The geodetic object is based on the WGS84
        ellipsoid.

        Parameters
        ----------
        radius: float
            Radius of the ellipsoid in [m]
        flattening: float
            Flattening of the ellipsoid (1 - b/a)
        """
        self._internal = PyGeodetic(radius, flattening)

    @property
    def latitude(self) -> float:
        return self._internal.latitude

    @property
    def longitude(self) -> float:
        return self._internal.longitude

    @property
    def altitude(self) -> float:
        return self._internal.altitude

    @property
    def valid(self) -> bool:
        return self._internal.is_valid()

    def from_lat_lon_alt(self, latitude: float, longitude: float, altitude: float):
        self._internal.from_lat_lon_alt(
            float(latitude), float(longitude), float(altitude)
        )

    def from_xyz(self, xyz: np.ndarray):
        self._internal.from_xyz(np.atleast_1d(xyz).astype(np.float64))

    def from_tangent_altitude(
        self, altitude: float, observer: np.ndarray, boresight: np.ndarray
    ) -> np.ndarray:
        return self._internal.from_tangent_altitude(
            altitude,
            np.atleast_1d(observer).astype(np.float64),
            np.atleast_1d(boresight).astype(np.float64),
        )

    def from_tangent_point(self, observer, look_vector):
        return self._internal.from_tangent_point(
            np.atleast_1d(observer).astype(np.float64),
            np.atleast_1d(look_vector).astype(np.float64),
        )

    @property
    def location(self) -> np.ndarray:
        return self._internal.location

    @property
    def local_south(self) -> np.ndarray:
        return self._internal.local_south

    @property
    def local_up(self) -> np.ndarray:
        return self._internal.local_up

    @property
    def local_west(self) -> np.ndarray:
        return self._internal.local_west

    def altitude_intercepts(
        self, altitude: float, observer: np.ndarray, look_vector: np.ndarray
    ) -> (np.ndarray, np.ndarray):
        return self._internal.altitude_intercepts(
            altitude,
            np.atleast_1d(observer).astype(np.float64),
            np.atleast_1d(look_vector).astype(np.float64),
        )


class WGS84(Geodetic):
    def __init__(self):
        """
        A geodetic object based upon the standard WGS84 ellipsoid. See
        :py:class:`sasktran2.Geodetic` for more information and usage.
        """
        super().__init__(6378137.0, 1.0 / 298.257223563)

    def __repr__(self):
        if self.valid:
            return f"WGS84 Location:\nLatitude: {self.latitude}, Longitude: {self.longitude}, Altitude: {self.altitude}"
        else:  # noqa: RET505
            return "WGS84 Unitialized"


class SphericalGeoid(Geodetic):
    def __init__(self, radius: float):
        """
        A geoid that is represented as a perfect sphere. See
        :py:class:`sasktran2.Geodetic` for more information and usage.

        Parameters
        ----------
        radius: float
            Radius of the sphere in [m]
        """
        super().__init__(radius, 0.0)

    def __repr__(self):
        if self.valid:
            return f"Spherical Geoid Location:\nLatitude: {self.latitude}, Longitude: {self.longitude}, Altitude: {self.altitude}"
        else:  # noqa: RET505
            return "Spherical Geoid Unitialized"
