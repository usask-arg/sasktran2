from __future__ import annotations

import shlex
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from sasktran2.database.base import CachedDatabase
from sasktran2.util import get_hapi


def _read_line_file_pd(file_path: Path) -> xr.Dataset:
    """
    Reads an AER linefile and puts it in the same format as the HITRAN optical reader.
    This function uses the pandas backend, it is faster but does not support all features

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
        "reference",
    ]

    data = pd.read_fwf(
        file_path,
        widths=widths,
        names=names,
        # converters=converters,
        comment=">",
        on_bad_lines="skip",
        error_bad_lines=False,
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
    data = data[
        [
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
        ]
    ]

    # Reset the index column
    data = data.reset_index(drop=True)

    # TODO: Speed dependent parameters, line coupling, and line mixing
    return xr.Dataset.from_dataframe(data)


def _read_line_file_py(file_path: Path) -> xr.Dataset:
    """
    Reads an AER linefile and puts it in the same format as the HITRAN optical reader.
    This function uses the python backend, it is slower but supports all features such as
    line coupling

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
    cols = np.concatenate(([0], np.cumsum(widths)))  # e.g. [ 0,  2,  3, 15, 25, ...]

    names = {
        "molec_id": int,
        "local_iso_id": int,
        "nu": float,
        "sw": float,
        "rsq": float,
        "gamma_air": float,
        "gamma_self": float,
        "elower": float,
        "n_air": float,
        "delta_air": float,
        "upper_quanta": str,
        "lower_quanta": str,
        "upper_local_q": str,
        "lower_local_q": str,
        "error_codes": str,
        "reference": str,
    }
    name_items = list(names.items())  # for indexing

    all_data = []

    with file_path.open() as f:
        for line in f:
            # Skip blank lines or comment lines
            if (not line.strip()) or line[0] in (">", "%"):
                continue

            # Slice out all fields at once
            fields = [
                line[cols[i] : cols[i + 1]].replace("D", "E")
                for i in range(len(widths))
            ]

            # Parse them according to the dictionary
            record = {}
            for (k, cast), raw in zip(name_items, fields, strict=True):
                record[k] = cast(raw)

            # Check if we need to read the line-coupling parameters
            if record["reference"][-2:-1] == "-":
                lc = next(f)  # read the next line
                record["Y_200"] = float(lc[3:15])
                record["G_200"] = float(lc[15:26])
                record["Y_250"] = float(lc[26:39])
                record["G_250"] = float(lc[39:50])
                record["Y_296"] = float(lc[50:63])
                record["G_296"] = float(lc[63:74])
                record["Y_340"] = float(lc[74:87])
                record["G_340"] = float(lc[87:98])
            else:
                record["Y_200"] = np.nan
                record["G_200"] = np.nan
                record["Y_250"] = np.nan
                record["G_250"] = np.nan
                record["Y_296"] = np.nan
                record["G_296"] = np.nan
                record["Y_340"] = np.nan
                record["G_340"] = np.nan

            all_data.append(record)

    # Convert to dataframe
    data = pd.DataFrame(all_data)

    # Reset the index column
    data = data.reset_index(drop=True)

    ds = xr.Dataset.from_dataframe(data)

    ds["Y"] = ds[["Y_200", "Y_250", "Y_296", "Y_340"]].to_array(dim="temperature")
    ds["G"] = ds[["G_200", "G_250", "G_296", "G_340"]].to_array(dim="temperature")

    ds = ds.drop_vars(
        ["Y_200", "G_200", "Y_250", "G_250", "Y_296", "G_296", "Y_340", "G_340"]
    )

    ds = ds.assign_coords(temperature=[200.0, 250.0, 296.0, 340.0])

    return ds.drop(
        [
            "error_codes",
            "reference",
            "upper_quanta",
            "lower_quanta",
            "upper_local_q",
            "lower_local_q",
        ]
    )


class AERLineDatabase(CachedDatabase):
    """
    Database for AER line files.  Currently only supports the 3.8.1 version of the AER line database
    """

    def __init__(
        self,
        db_root: Path | None = None,
        rel_path: Path | None = "aer_lines",
        version: str = "3.8.1",
    ):
        CachedDatabase.__init__(self, db_root, rel_path=rel_path)

        self._hapi = get_hapi()
        self._version = version

        self._version_map = {"3.8.1": "5120012"}

    def path(self, _key: str, **kwargs) -> Path:
        dir = self._db_root
        version_dir = dir.joinpath(f"aer_v_{self._version}")
        # check if the AER database is downloaded, if not download it

        if not version_dir.exists():
            from zenodo_get import zenodo_get

            try:
                zenodo_get(
                    shlex.split(
                        f'--record {self._version_map[self._version]} -o "{dir.as_posix()}"'
                    )
                )
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

        return version_dir

    def load_ds(self, _key: str, **kwargs) -> xr.Dataset:
        return None

    def clear(self):
        # Delete the entire database folder
        folder = self._db_root.joinpath(f"aer_v_{self._version}")
        if folder.exists():
            for file in folder.glob("*"):
                if file.is_dir():
                    for subfile in file.glob("*"):
                        subfile.unlink()
                    file.rmdir()
                else:
                    file.unlink()
