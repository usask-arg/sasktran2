import xarray as xr

from sasktran2 import ViewingGeometry


class ViewingGeometryContainer(ViewingGeometry):
    def __init__(self, geometry_ds: xr.Dataset | None = None):
        super().__init__()
        self.geometry_ds = geometry_ds

    def add_geometry_to_radiance(self, radiance: xr.Dataset) -> xr.Dataset:
        return radiance
