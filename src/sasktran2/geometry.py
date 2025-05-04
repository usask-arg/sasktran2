from __future__ import annotations

import numpy as np

from sasktran2._core_rust import GeometryType, InterpolationMethod, PyGeometry1D


class Geometry1D:
    _geometry: PyGeometry1D

    def __init__(
        self,
        cos_sza: float,
        solar_azimuth: float,
        earth_radius_m: float,
        altitude_grid_m: np.ndarray,
        interpolation_method: InterpolationMethod,
        geometry_type: GeometryType,
    ):
        self._geometry = PyGeometry1D(
            cos_sza,
            solar_azimuth,
            earth_radius_m,
            np.atleast_1d(altitude_grid_m).astype(np.float64),
            interpolation_method,
            geometry_type,
        )

    def altitudes(self) -> np.ndarray:
        return self._geometry.altitudes()

    @property
    def refractive_index(self) -> np.ndarray:
        return self._geometry.refractive_index

    @refractive_index.setter
    def refractive_index(self, value: np.ndarray) -> None:
        self._geometry.refractive_index = value
