from sasktran2._core_rust import PyGeometry1D
from sasktran2._core_rust import GeometryType, InterpolationMethod
import numpy as np

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