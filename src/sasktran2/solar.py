import abc

import pandas as pd


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
            to the observer.
        """


class SolarHandlerForced(SolarGeometryHandlerBase):
    def __init__(self, solar_zenith: float, solar_azimuth: float):
        """
        A solar handler where the solar angles are forced to be the same for all locations and times.

        Parameters
        ----------
        solar_zenith : float
            Solar zenith angle in [degrees]
        solar_azimuth : float
            Solar azimuth angle in [degrees]. Here 0 degrees corresponds to pointing true north.
            Note that this is NOT relative to the observer.
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
    def target_solar_angles(
        self, latitude: float, longitude: float, altitude: float, time: pd.Timestamp
    ) -> tuple[float, float]:
        pass
