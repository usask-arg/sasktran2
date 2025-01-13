from __future__ import annotations

from pathlib import Path

import xarray as xr

from sasktran2.util import get_hapi

from .base import CachedDatabase


class HITRANLineDatabase(CachedDatabase):
    def __init__(
        self, db_root: Path | None = None, rel_path: Path | None = "hitran_lines"
    ):
        """
        Downloads the HITRAN Line databases from hitran using the hitran-api package.  Both the raw
        HITRAN data is stored for use in the hitran-api, and a netcdf file version for use in the
        SASKTRAN2 broadening routines.
        """
        CachedDatabase.__init__(self, db_root, rel_path=rel_path)

    def path(self, key: str, **kwargs) -> Path:
        data_file = self._db_root.joinpath(f"{key}.data")
        header_file = self._db_root.joinpath(f"{key}.header")
        if not (data_file.exists() and header_file.exists()):
            # download lines for this molecule
            self._download_line_db(key)

    def load_ds(self, key: str, **kwargs) -> xr.Dataset:
        nc_file = self._db_root.joinpath(f"{key}.nc")
        if nc_file.exists():
            return xr.open_dataset(nc_file)

        hapi = get_hapi()

        self.path(key)
        self.initialize_hapi(key)

        data_vars = {}
        for param in hapi.PARLIST_ALL:
            try:
                param_column = hapi.getColumn(key, param)
                data_vars[param] = ("line", param_column)
            except KeyError:  # param does not exist for this molecule
                continue

        ds = xr.Dataset(data_vars=data_vars)
        ds.to_netcdf(nc_file)
        return ds

    def clear(self):
        # deletes all .data and .header files in database folder
        for header_file in self._db_root.glob("*.header"):
            header_file.unlink()
        for data_file in self._db_root.glob("*.data"):
            data_file.unlink()

    def _download_line_db(self, molecule):
        hapi = get_hapi()

        self.initialize_hapi()

        # Get global isotopologue ids for this molecule
        mol_index = hapi.ISO_ID_INDEX["mol_name"]
        iso_ids = []
        for id in hapi.ISO_ID:
            if hapi.ISO_ID[id][mol_index] == molecule:
                iso_ids.append(id)

        # Download files
        hapi.fetch_by_ids(
            molecule, iso_ids, 0.0, 1.0e6, ParameterGroups=["par_line", "sdvoigt", "ht"]
        )

        return

    def initialize_hapi(self, table_name: str | None = None):
        """
        Initialize the HAPI database with one or no tables. This is to reduce the time and memory that hapi.db_begin
        uses when it is given a folder with many data files.
        """
        hapi = get_hapi()

        hapi.VARIABLES["BACKEND_DATABASE_NAME"] = str(self._db_root)
        hapi.LOCAL_TABLE_CACHE = {}
        if table_name is not None:
            hapi.storage2cache(table_name)

        return
