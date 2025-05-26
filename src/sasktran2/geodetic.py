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
        """
        Returns
        -------
        float
            Geodetic latitude in degrees
        """
        return self._internal.latitude

    @property
    def longitude(self) -> float:
        """
        Returns
        -------
        float
            Geodetic longitude in degrees
        """
        return self._internal.longitude

    @property
    def altitude(self) -> float:
        """
        Returns
        -------
        float
            Altitude in [m] above the surface of the ellipsoid.
        """
        return self._internal.altitude

    @property
    def valid(self) -> bool:
        """
        Returns
        -------
        bool
            True if the geodetic object has been initialized, False otherwise.
        """
        return self._internal.is_valid()

    def from_lat_lon_alt(self, latitude: float, longitude: float, altitude: float):
        """
        Initializes the Geodetic based on a specifiec latitude, longitude, and altitude.

        Parameters
        ----------
        latitude : float
            Latitude in degrees (-90 to 90)
        longitude : float
            Longitude in degrees (0 to 360 or -180 to 180)
        altitude : float
            Altitude above the geoid in metres

        Examples
        --------
        >>> import sasktran2 as sk
        >>> geodetic = sk.WGS84()
        >>> geodetic.from_lat_lon_alt(latitude=-15, longitude=-20, altitude=7342)
        >>> print(geodetic)
        WGS84 Location:
        Latitude: -15.0, Longitude: 340.0, Altitude: 7342.0
        """
        self._internal.from_lat_lon_alt(
            float(latitude), float(longitude), float(altitude)
        )

    def from_xyz(self, location: np.ndarray):
        """
        Initializes the Geodetic from a geocentric location

        Parameters
        ----------
        location : np.ndarray
            Three element vector containing a location in geocentric coordinates

        Examples
        --------
        >>> import sasktran2 as sk
        >>> geodetic = sk.WGS84()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic)
        WGS84 Location:
        Latitude: -14.999999973747736, Longitude: 340.00000000000006, Altitude: 7344.999610390202
        """
        self._internal.from_xyz(np.atleast_1d(location).astype(np.float64))

    def from_tangent_altitude(
        self, altitude: float, observer: np.ndarray, boresight: np.ndarray
    ) -> np.ndarray:
        """
        Initialized the Geodetic from a specified tangent altitude, obsever location, and bore sight plane.

        Parameters
        ----------
        altitude : float
            Tangent altitude in meters
        observer : np.ndarray
            Three element array containing the obsever position in geocentric coordinates
        boresight : np.ndarray
            Three element array containing a normalized look vector that is within the bore sight plane.

        Returns
        -------
        np.ndarray
            Three element array containing the normalized look vector to the tangent point.

        Examples
        --------
        >>> import sasktran2 as sk
        >>> geodetic = sk.WGS84()
        >>> look = geodetic.from_tangent_altitude(15322, [ 3.676013154788849600e+005, 1.009976313640051500e+006,\
                                                            -6.871601202127538600e+006], [0, 0, 1])
        >>> print(look)
        [0.28880556 0.79348676 0.53569591]
        >>> print(geodetic)
        WGS84 Location:
        Latitude: -57.60888188776806, Longitude: 70.00000000000001, Altitude: 15321.971935882739
        """
        return self._internal.from_tangent_altitude(
            altitude,
            np.atleast_1d(observer).astype(np.float64),
            np.atleast_1d(boresight).astype(np.float64),
        )

    def from_tangent_point(self, observer, look_vector):
        """
        Initializes  the Geodetic by calculating the tangent point from an observer position and look vector

        Parameters
        ----------
        observer : np.ndarray
            Three element array containing the observer position in geocentric coordinates
        look_vector : np.ndarray
            Three element array containing a normalized look vector

        Examples
        --------
        >>> import sasktran2 as sk
        >>> geodetic = sk.WGS84()
        >>> geodetic.from_tangent_point([ 3.676013154788849600e+005, 1.009976313640051500e+006,\
                                            -6.871601202127538600e+006], [ 2.884568631765662100e-001,\
                                            7.925287180643269000e-001,  5.372996083468238900e-001])
        >>> print(geodetic)
        WGS84 Location:
        Latitude: -57.500000192733594, Longitude: 70.0, Altitude: 10002.99586173162
        """
        return self._internal.from_tangent_point(
            np.atleast_1d(observer).astype(np.float64),
            np.atleast_1d(look_vector).astype(np.float64),
        )

    @property
    def location(self) -> np.ndarray:
        """
        Returns
        -------
        np.ndarray
            Geocentric location in cartesian coordinates
        """
        return self._internal.location

    @property
    def local_south(self) -> np.ndarray:
        """
        Returns
        -------
        np.ndarray
            A unit vector pointing in the local south direction
        """
        return self._internal.local_south

    @property
    def local_up(self) -> np.ndarray:
        """
        Returns
        -------
        np.ndarray
            A unit vector pointing up (perpindicular to the ellipsoidal surface)
        """
        return self._internal.local_up

    @property
    def local_west(self) -> np.ndarray:
        """
        Returns
        -------
        np.ndarray
            A unit vector pointing in the local west direction
        """
        return self._internal.local_west

    def altitude_intercepts(
        self, altitude: float, observer: np.ndarray, look_vector: np.ndarray
    ) -> (np.ndarray, np.ndarray):  # type: ignore  # noqa: PGH003
        """
        Calculate the two intersections of a line of sight and an altitude.

        Parameters
        ----------
        altitude : float
            Altitude in meters.
        observer : np.ndarray
            Three element array containing the obsever position in geocentric coordinates.
        look_vector : np.ndarray
            Three element array containing a normalized look vector.

        Returns
        -------
        np.ndarray
            Three element array containing the first (entering) intercept in geocentric coordinates.
        np.ndarray
            Three element array containing the second (exiting) intercept in geocentric coordinates.

        Examples
        --------
        >>> import sasktran2 as sk
        >>> import numpy as np
        >>> geodetic = sk.WGS84()
        >>> look = geodetic.from_tangent_altitude(15322, [3.676013154788849600e+005, 1.009976313640051500e+006, \
                                                    -6.871601202127538600e+006], [0, 0, 1])
        >>> obs = geodetic.location
        >>> intercept1, intercept2 = geodetic.altitude_intercepts(16000, obs, look)
        >>> print(np.array_str(intercept1, precision=3))
        [ 1147302.059  3152186.5   -5425360.027]
        >>> print(np.array_str(intercept2, precision=3))
        [ 1201098.489  3299990.978 -5325574.803]
        """
        return self._internal.altitude_intercepts(
            altitude,
            np.atleast_1d(observer).astype(np.float64),
            np.atleast_1d(look_vector).astype(np.float64),
        )

    def osculating_spheroid(self) -> tuple[float, np.ndarray]:
        """
        Returns the osculating spheroid at the current location.

        The osculating spheroid is the best fit sphere to the geoid at the current location.

        Returns
        -------
        float
            The radius of the osculating spheroid in meters.
        np.ndarray
            A three element array containing the offset of the center of the geoid to the center of the osculating
            spheroid in geocentric coordinates.
        """
        result = self._internal.osculating_spheroid()

        return result[0], np.atleast_1d(result[1]).astype(np.float64)


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
