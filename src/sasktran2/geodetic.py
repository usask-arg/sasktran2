from __future__ import annotations

from sasktran2._core import Geodetic


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
