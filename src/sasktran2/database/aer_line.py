from __future__ import annotations

from os import stat
from pathlib import Path

import xarray as xr

from sasktran2.util import get_hapi

from .base import CachedDatabase

import pandas as pd
from pathlib import Path
import xarray as xr
import numpy as np
import shlex

def _read_line_file(file_path: Path) -> xr.Dataset:
    """
    Reads an AER linefile and puts it in the same format as the HITRAN optical reader

    Parameters
    ----------
    file_path : Path
        Path to the individual molecule datafile

    Returns
    -------
    xr.Dataset
    """
    # HITRAN F100 format is
    widths = np.array([2, 1, 12, 10, 10, 5, 5, 10, 4, 8, 3, 3, 9, 9, 3, 6])
    names = [
        "molec_id",
        "local_iso_id",
        "nu",
        "sw",
        "rsq",
        "gamma_air",
        "gamma_self",
        "elower",
        "n_air",
        "delta_air",
        "upper_quanta",
        "lower_quanta",
        "upper_local_q",
        "lower_local_q",
        "error_codes",
        "reference"
        ]

    data = pd.read_fwf(file_path,
                       widths=widths,
                       names=names,
                       #converters=converters,
                       comment=">",
                       on_bad_lines="skip",
                       error_bad_lines=False
                       )
    # Drop rows where mol_id is %% 
    data = data[data["molec_id"] != "%%"]

    # Rows where iso_id is NaN are extra information
    data = data.dropna(subset=["local_iso_id"])
    data = data[data["local_iso_id"] != "-"]

    # In Line intensity we have to replace the D with an E, if it was loaded in as a float
    if data["sw"].dtype != float:
        data.loc[:, "sw"] = data["sw"].str.replace("D", "E")

    # Now perform the conversions
    data["molec_id"] = data["molec_id"].astype(int)
    data["local_iso_id"] = data["local_iso_id"].astype(int)
    data["nu"] = data["nu"].astype(float)
    data["sw"] = data["sw"].astype(float)
    data["rsq"] = data["rsq"].astype(float)
    data["gamma_air"] = data["gamma_air"].astype(float)
    data["gamma_self"] = data["gamma_self"].astype(float)
    data["elower"] = data["elower"].astype(float)
    data["n_air"] = data["n_air"].astype(float)
    data["delta_air"] = data["delta_air"].astype(float)

    # Drop all other columns
    data = data[["molec_id", "local_iso_id", "nu", "sw", "rsq", "gamma_air", "gamma_self", "elower", "n_air", "delta_air"]]

    # Reset the index column
    data = data.reset_index(drop=True)

    # TODO: Speed dependent parameters, line coupling, and line mixing
    return xr.Dataset.from_dataframe(data)


class AERLineDatabase(CachedDatabase):
    def __init__(
        self, db_root: Path | None = None, rel_path: Path | None = "aer_lines", version: str = "3.8.1"
    ):
        """
        Loads in the AER line database from the AER line files.
        """
        CachedDatabase.__init__(self, db_root, rel_path=rel_path)

        self._hapi = get_hapi()
        self._version = version

        self._version_map = {"3.8.1": "5120012"}

    def _get_molecule_id(self, name: str):
        for i in range(1, 48):
            mol_name = self._hapi.moleculeName(i)

            if mol_name.lower() == name.lower():
                return i
        msg = f"Could not find molecule {name}"
        raise ValueError(msg)

    def path(self, key: str, **kwargs) -> Path:
        dir = self._db_root
        version_dir = dir.joinpath(f"aer_v_{self._version}")
        mol_id = self._get_molecule_id(key)
        mol_dir = version_dir.joinpath(f"line_files_By_Molecule/{mol_id:02d}_{key.upper()}")

        nc_file = mol_dir.joinpath(f"{mol_id:02d}_{key.upper()}.nc")

        if nc_file.exists():
            return nc_file

        # Else, check if the AER database is downloaded, if not download it

        if not version_dir.exists():
            from zenodo_get import zenodo_get
            try:
                zenodo_get(shlex.split(f'--record {self._version_map[self._version]} -o "{dir.as_posix()}"'))
            except ImportError as e:
                msg = "zenodo_get is required to download the AER line database"
                raise ImportError(msg) from e

            file = dir.joinpath(f"aer_v_{self._version}.tar.gz")
            # Extract the tar file
            import tarfile
            with tarfile.open(file, "r:gz") as tar:
                tar.extractall(dir)
            # Delete the tar file
            file.unlink()

        # Now get the raw data file
        mol_data_file = mol_dir.joinpath(f"{mol_id:02d}_{key.upper()}")
        if not mol_data_file.exists():
            msg = f"Could not find line file for molecule {key}, tried {mol_data_file.as_posix()}"
            raise FileNotFoundError(msg)

        # And convert to netcdf
        ds = _read_line_file(mol_data_file)
        ds.to_netcdf(nc_file)

        return nc_file

    def load_ds(self, key: str, **kwargs) -> xr.Dataset:
        return xr.open_dataset(self.path(key, **kwargs))


    def clear(self):
        # Delete the entire database folder
        if self._db_root.exists():
            for file in self._db_root.glob("*"):
                if file.is_dir():
                    for subfile in file.glob("*"):
                        subfile.unlink()
                    file.rmdir()
                else:
                    file.unlink()

if __name__ == "__main__":
    line_file = Path("~/Documents/data/aer_v_3.8.1/line_files_By_Molecule/02_CO2/02_CO2/").expanduser()

    _read_line_file(line_file)