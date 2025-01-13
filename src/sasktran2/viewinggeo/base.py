from __future__ import annotations

import numpy as np
import xarray as xr

from sasktran2._core import Geometry1D


class ViewingGeometryContainer:
    def __init__(self, geometry_ds: xr.Dataset | None = None):
        self._geometry_ds = geometry_ds

    @property
    def geometry_ds(self) -> xr.Dataset:
        return self._geometry_ds

    def add_geometry_to_radiance(self, radiance: xr.Dataset) -> xr.Dataset:
        return xr.merge([self.geometry_ds, radiance])

    def recommended_earth_radius(self) -> float:
        return 6371000.0

    def recommended_cos_sza(self) -> float:
        raise NotImplementedError()

    def model_geometry(self, altitude_grid_m: np.ndarray) -> Geometry1D:
        raise NotImplementedError()
