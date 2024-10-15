import abc

import pandas as pd


class SolarGeometryHandlerBase:
    def __init__(self):
        pass

    @abc.abstractmethod
    def target_solar_angles(
        self, latitude: float, longitude: float, altitude: float, time: pd.Timestamp
    ) -> tuple[float, float]:
        pass


class SolarHandlerForced(SolarGeometryHandlerBase):
    def __init__(self, solar_zenith: float, solar_azimuth: float):
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
