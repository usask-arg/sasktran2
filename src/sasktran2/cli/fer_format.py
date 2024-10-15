from __future__ import annotations

import xarray as xr


class FrontEndRadiance:
    def __init__(self, sk2_radiance: xr.Dataset, geometry_ds: xr.Dataset):
        """
        A base container for the front-end radiance data format.

        Parameters
        ----------
        ds : xr.Dataset
            _description_
        """

        self._data = xr.merge([sk2_radiance, geometry_ds])

    @property
    def data(self) -> xr.Dataset:
        return self._data
