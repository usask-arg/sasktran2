from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate

from sasktran2.appconfig import database_root
from sasktran2.database.web import StandardDatabase


class RefractiveIndex:
    def __init__(
        self,
        refractive_index_fn: Callable[[float], complex],
        identifier: str,
        args: list[str] | None = None,
    ) -> None:
        """
        A generic implementation of a refractive index function. This class is a light wrapper
        on top of a callable function that returns the complex refractive index for a given
        wavelength, as well as providing a unique identifier for the refractive index.

        Parameters
        ----------
        refractive_index_fn : Callable[[float], complex]
            Function that takes in a wavelength in nm (and any keyword arguments) and returns the complex refractive index
        identifier : str
            A unique identifier for the refractive index
        args : list[str]
            Names of keyword arguments required by refractive_index_fn (default no arguments)
        """
        self._fn = refractive_index_fn
        self._identifier = identifier
        self._args = [] if args is None else args

    def refractive_index(self, wavelength_nm: np.ndarray, **kwargs) -> np.ndarray:
        """
        Returns the complex refractive index for a given wavelength

        Parameters
        ----------
        wavelength_nm : np.ndarray

        Returns
        -------
        np.ndarray
        """
        return self._fn(wavelength_nm, **kwargs)

    @property
    def refractive_index_fn(self):
        """
        Get the function that returns the complex refractive index
        """
        return self._fn

    @property
    def identifier(self):
        """
        Get the unique identifier for this refractive index
        """
        return self._identifier

    @property
    def args(self):
        """
        Get the list of any keyword arguments required by refractive_index_fn
        """
        return self._args


def _from_osiris_file(path: Path):
    """
    Reads in data from the old OSIRIS refractive index files and returns a callable function
    """
    data = pd.read_csv(path.as_posix(), header=None)

    return interpolate.interp1d(
        data.values[:, 0], data.values[:, 1] - 1j * data.values[:, 2]
    )


def _from_aria_file(key: str):
    """
    Reads in data from ARIA .ri file and returns a callable function
    Original function accessed from https://eodg.atm.ox.ac.uk/ARIA/media/files/read_ri.py on 2025-12-12
    """

    data_dir = database_root() / "refractive_index" / "aria"

    zip_file = data_dir / "ARIA.zip"
    url = "https://eodg.atm.ox.ac.uk/ARIA/data_files/ARIA.zip"
    if not zip_file.is_file():
        import requests

        data_dir.mkdir(parents=True, exist_ok=True)
        with requests.get(url) as response:
            response.raise_for_status()
            with zip_file.open("wb") as f:
                f.write(response.content)

    if not (data_dir / "data_files").is_dir():
        import tarfile

        with tarfile.open(zip_file, "r:gz") as tar:
            tar.extractall(data_dir)

    matches = sorted(data_dir.glob(f"**/{key}.ri"))
    if len(matches) == 0:
        msg = f"No ARIA data file found with key {key}"
        raise ValueError(msg)
    if len(matches) > 1:
        msg = f"Multiple ARIA data file found with key {key}"
        raise ValueError(msg)
    path = matches[0]

    class ReadError(Exception):
        def __init__(self, value):
            self.parameter = value

        def __str__(self):
            return repr(self.parameter)

    expected_header_names = [
        "FORMAT",
        "DESCRIPTION",
        "DISTRIBUTEDBY",
        "SUBSTANCE",
        "SAMPLEFORM",
        "TEMPERATURE",
        "CONCENTRATION",
        "REFERENCE",
        "DOI",
        "SOURCE",
        "CONTACT",
        "COMMENT",
    ]
    expected_column_names = ["wavl", "wavn", "n", "dn", "k", "dk"]
    with path.open() as f:
        t = f.readlines()
        t = [x.strip() for x in t]  # strips whitespace from beginning and end of lines

        out = {"header": {}, "data": {}}

        header_lines = 0
        data_lines = 0
        for line in t:
            if line[0] == "#":
                header_lines += 1
                if data_lines > 0:
                    msg = f"Incorrectly formatted file ({path}): Header not contiguous."
                    raise ReadError(msg)
            else:
                data_lines += 1
        for i in range(1, data_lines):  # ignore blank lines at end of file
            if any(char.isdigit() for char in t[-i]):
                break
            data_lines -= 1

        if header_lines == 0:
            msg = f"Incorrectly formatted file ({path}): No header."
            raise ReadError(msg)
        if data_lines == 0:
            msg = f"Incorrectly formatted file ({path}): No data."
            raise ReadError(msg)

        for i in range(header_lines):
            line = t[i][1:]  # strip leading '#'
            if line[0] != "#":
                tag_name = line.split("=", 1)[0].strip().upper()
                if tag_name not in expected_header_names:
                    logging.debug(
                        f'Unknown header tag "{tag_name}", so ignored (file: {path})'
                    )
                    continue
                try:  # ensure tag content is encoded consistently
                    tag_content = (
                        line.split("=", 1)[1].strip()
                        # .decode("utf8")
                        # .encode("utf8", "xmlcharrefreplace")
                    )
                except UnicodeDecodeError:
                    tag_content = (
                        line.split("=", 1)[1].strip()
                        # .decode("iso-8859-1")
                        # .encode("utf8", "xmlcharrefreplace")
                    )
                if tag_name in out["header"]:
                    tag_content = out["header"][tag_name] + " " + tag_content
                out["header"][tag_name] = tag_content
            elif tag_name in expected_header_names:
                tag_content = line[1:].strip()
                out["header"][tag_name] = out["header"][tag_name] + " " + tag_content

        if "FORMAT" not in out["header"]:
            msg = f"Incorrectly formatted file ({path}): No FORMAT tag in header."
            raise ReadError(msg)

        column_labels = out["header"]["FORMAT"].split()
        column_labels = [x.strip().lower() for x in column_labels]
        for cl in column_labels:
            if cl not in expected_column_names:
                logging.debug(f'Unknown column name "{cl}", so ignored (file: {path})')
                continue
            out["data"][cl] = []
        for ln in range(header_lines, data_lines):
            line = t[ln].split()
            line = [x.strip() for x in line]
            for c in range(len(column_labels)):
                if column_labels[c] in expected_column_names:
                    out["data"][column_labels[c]].append(float(line[c]))

        # add wavl & wavn columns if needed (wavl in micro-m, wavn in cm-1)
        if "wavn" not in out["data"]:
            out["data"]["wavn"] = [
                float(10000) / x if x != 0 else float("nan")
                for x in out["data"]["wavl"]
            ]
        if "wavl" not in out["data"]:
            out["data"]["wavl"] = [
                float(10000) / x if x != 0 else float("nan")
                for x in out["data"]["wavn"]
            ]

        col_lengths = []
        for col in out["data"]:
            col_lengths.append(len(out["data"][col]))
        if len(set(col_lengths)) > 1:
            msg = f"Incorrectly formatted file ({path}): Data columns have different lengths."
            raise ReadError(msg)

    wavl = np.array(out["data"]["wavl"])
    n = np.array(out["data"]["n"])
    k = np.array(out["data"]["k"])
    isort = np.argsort(wavl)

    return interpolate.interp1d(1e3 * wavl[isort], n[isort] - 1j * k[isort])


class H2SO4(RefractiveIndex):
    def __init__(self, source: str = "osiris") -> None:
        """
        A refractive index for H2SO4. The default source is the old OSIRIS data.
        """
        if source.lower() == "osiris":
            super().__init__(
                _from_osiris_file(
                    StandardDatabase().path("refractive_index/refrac_h2so4_osiris.txt")
                ),
                "H2SO4_osiris",
            )
        else:
            msg = "Only osiris source is supported for H2SO4 refractive index data"
            raise ValueError(msg)


class Dust(RefractiveIndex):
    def __init__(self, source: str = "osiris") -> None:
        """
        A refractive index for Dust. The default source is the old OSIRIS data.
        """
        if source.lower() == "osiris":
            super().__init__(
                _from_osiris_file(
                    StandardDatabase().path("refractive_index/refrac_dust_osiris.txt")
                ),
                "dust_osiris",
            )
        else:
            msg = "Only osiris source is supported for dust refractive index data"
            raise ValueError(msg)


class Ice(RefractiveIndex):
    def __init__(self, source: str = "osiris") -> None:
        """
        A refractive index for ice. The default source is the old OSIRIS data.
        """
        if source.lower() == "osiris":
            super().__init__(
                _from_osiris_file(
                    StandardDatabase().path("refractive_index/refrac_ice_osiris.txt")
                ),
                "ice_osiris",
            )
        else:
            msg = "Only osiris source is supported for ice refractive index data"
            raise ValueError(msg)


class Water(RefractiveIndex):
    def __init__(self, source: str = "osiris") -> None:
        """
        A refractive index for Water. The default source is the old OSIRIS data.
        """
        if source.lower() == "osiris":
            super().__init__(
                _from_osiris_file(
                    StandardDatabase().path("refractive_index/refrac_water_osiris.txt")
                ),
                "water_osiris",
            )
        elif source.lower() == "aria":
            super().__init__(_from_aria_file("H2O_Mcgarragh_2018"), "water_aria")
        else:
            msg = "Only osiris and aria sources are supported for water refractive index data"
            raise ValueError(msg)
