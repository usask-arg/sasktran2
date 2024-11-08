import xarray as xr

from .ecef import ecef_to_sasktran2_ray

__all__ = ["ecef_to_sasktran2_ray"]


class ViewingGeometryContainer:
    def __init__(self, geometry_ds: xr.Dataset | None = None):
        self.geometry_ds = geometry_ds

    def add_geometry_to_radiance(self, radiance: xr.Dataset) -> xr.Dataset:
        return radiance
