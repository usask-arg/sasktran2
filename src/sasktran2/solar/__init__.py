from __future__ import annotations

import abc

import pandas as pd

from .model import SolarModel

__all__ = ["SolarModel"]


class SolarGeometryHandlerBase:
    def __init__(self):
        """
        Interface class for the solar handler.  A solar handler is used to calculate the
        solar angles for a give time and location.
        """

    @abc.abstractmethod
    def target_solar_angles(
        self, latitude: float, longitude: float, altitude: float, time: pd.Timestamp
    ) -> tuple[float, float]:
        """
        Calculates the solar zenith and solar azimuth angles for a given location and time.

        Parameters
        ----------
        latitude : float
            Latitude in [degrees north]
        longitude : float
            Longitude in [degrees east]
        altitude : float
            Altius in [m]
        time : pd.Timestamp
            Time of interest

        Returns
        -------
        tuple[float, float]
            Solar zenith and solar azimuth angles in [degrees]. Solar azimuth is relative to true north, not
            to the observer. And is defined clockwise from north, such that the east is 90 degrees, south is 180 degrees, etc.
            The azimuthal direction points towards the sun.
        """


class SolarGeometryHandlerForced(SolarGeometryHandlerBase):
    def __init__(self, solar_zenith: float, solar_azimuth: float):
        """
        A solar handler where the solar angles are forced to be the same for all locations and times.

        Parameters
        ----------
        solar_zenith : float
            Solar zenith angle in [degrees]
        solar_azimuth : float
            Solar azimuth angle in [degrees]. Here 0 degrees corresponds to pointing true north, 90 degrees to east, etc.
            Note that this is NOT relative to the observer.  Note this is also pointing towards the sun, not away from the sun.
        """
        self.solar_zenith = solar_zenith
        self.solar_azimuth = solar_azimuth

    def target_solar_angles(
        self,
        latitude: float,  # noqa: ARG002
        longitude: float,  # noqa: ARG002
        altitude: float,  # noqa: ARG002
        time: pd.Timestamp,  # noqa: ARG002
    ) -> tuple[float, float]:
        return self.solar_zenith, self.solar_azimuth


class SolarGeometryHandlerAstropy(SolarGeometryHandlerBase):
    """
    Solar handler where the astropy package is used to calculate the solar angles.
    Must have astropy installed in the Python environment in order to use this solar handler.
    """

    def target_solar_angles(
        self, latitude: float, longitude: float, altitude: float, time: pd.Timestamp
    ) -> tuple[float, float]:
        try:
            from astropy.coordinates import AltAz, EarthLocation, get_sun
            from astropy.time import Time
        except ImportError:
            msg = "Astropy is required to use the astropy solar hanlder"
            raise ImportError(msg)  # noqa: B904

        # Location object
        loc = EarthLocation(lat=latitude, lon=longitude, height=altitude)
        sun = get_sun(Time(time))

        altaz = sun.transform_to(AltAz(obstime=time, location=loc))

        return 90.0 - altaz.alt.deg, altaz.az.deg
