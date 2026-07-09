"""Local builders for Zenodo-hosted ozone absorption cross-section databases.

The public optical properties in :mod:`sasktran2.optical` expect the same compact
NetCDF layout used by the built-in DBM database: an ``xs`` variable in
``m^2/molecule`` with ``temperature_k`` and ``wavelength_nm`` coordinates.  The
Zenodo sources are not distributed in that layout, so these classes are small
cache managers that download the original files once and convert them into the
standard SASKTRAN2 absorber database format.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np
import xarray as xr

from .base import CachedDatabase
from .zenodo import download_zenodo_record


class _ZenodoO3Database(CachedDatabase):
    """Base class for one-file Zenodo O3 datasets with a derived NetCDF cache.

    ``_raw_filename`` is kept next to ``_processed_filename`` under the user's
    configured database root.  On first use we download the Zenodo record, locate
    the expected raw file, and then write the converted NetCDF file.  Later calls
    only open the converted cache.
    """

    _record_id: str
    _raw_filename: str
    _processed_filename: str

    def __init__(
        self,
        db_root: Path | None = None,
        rel_path: Path | None = Path("cross_sections/o3"),
    ) -> None:
        super().__init__(db_root, rel_path=rel_path)

    def path(self, _key: str = "", **kwargs) -> Path:
        output_file = self._db_root.joinpath(self._processed_filename)
        if not output_file.exists():
            raw_file = self._raw_file()
            self._convert(raw_file, output_file)

        return output_file

    def load_ds(self, key: str = "", **kwargs) -> xr.Dataset:
        return xr.open_dataset(self.path(key, **kwargs))

    def clear(self):
        for filename in (self._processed_filename, self._raw_filename):
            path = self._db_root.joinpath(filename)
            if path.exists():
                path.unlink()

    def _raw_file(self) -> Path:
        raw_file = self._db_root.joinpath(self._raw_filename)
        if raw_file.exists():
            return raw_file

        download_zenodo_record(self._record_id, self._db_root)

        if raw_file.exists():
            return raw_file

        matches = list(self._db_root.rglob(self._raw_filename))
        if matches:
            return matches[0]

        msg = (
            f"Could not find {self._raw_filename} after downloading Zenodo "
            f"record {self._record_id}"
        )
        raise OSError(msg)

    def _convert(self, raw_file: Path, output_file: Path) -> None:
        raise NotImplementedError


class O3BirkWagnerDatabase(_ZenodoO3Database):
    """
    Cached Birk and Wagner O3 absorption cross-section database.

    The source archive contains one ASCII table for each measured temperature.
    Each table has wavelength, absorption cross section, and uncertainty columns.
    Only the measured cross sections are stored in the SASKTRAN2 cache; the raw
    uncertainty columns and polynomial coefficient file remain available in the
    downloaded Zenodo archive if users need them for their own analysis. The
    wavelength grid is preserved as a vacuum wavelength grid; applying an
    air-to-vacuum conversion shifts the Huggins-band structure out of alignment
    with the independently vacuum-labelled Serdyuchenko dataset.
    """

    _record_id = "1485588"
    _raw_filename = "O3_UV_region_absorption_cross_section_database_13112018.zip"
    _processed_filename = "birk_wagner.nc"

    def _convert(self, raw_file: Path, output_file: Path) -> None:
        # The measured ACS files are named O3_ACS_193K.asc, O3_ACS_213K.asc,
        # etc.  The archive also contains a polynomial-coefficient file, which
        # is intentionally ignored here so the optical property matches the
        # discrete measured-temperature database style used by O3DBM.
        #
        # No air-to-vacuum conversion is applied to the wavelength column. In the
        # 315-345 nm overlap region, the measured structures align with the
        # Serdyuchenko vacuum wavelength grid as-is and become visibly shifted if
        # the Birk/Wagner grid is treated as air wavelengths.
        temperature_re = re.compile(r"O3_ACS_(\d+)K\.asc$")
        temperatures = []
        cross_sections = []
        wavelength_nm = None

        with zipfile.ZipFile(raw_file) as zf:
            source_files = [
                name
                for name in zf.namelist()
                if temperature_re.search(Path(name).name) is not None
            ]

            if not source_files:
                msg = f"No Birk/Wagner O3_ACS temperature files found in {raw_file}"
                raise OSError(msg)

            for source_file in source_files:
                match = temperature_re.search(Path(source_file).name)
                temperatures.append(float(match.group(1)))

                with zf.open(source_file) as f:
                    data = np.loadtxt(f, skiprows=2)

                if wavelength_nm is None:
                    wavelength_nm = data[:, 0]
                elif not np.allclose(wavelength_nm, data[:, 0], rtol=0, atol=1e-8):
                    msg = "Birk/Wagner O3 files do not share a common wavelength grid"
                    raise ValueError(msg)

                cross_sections.append(data[:, 1])

        sort_idx = np.argsort(temperatures)
        temperature_k = np.asarray(temperatures, dtype=float)[sort_idx]
        xs = np.asarray(cross_sections, dtype=float)[sort_idx] / 1e4

        ds = xr.Dataset(
            {
                "xs": (
                    ["temperature_k", "wavelength_nm"],
                    xs,
                    {"units": "m^2 molecule^-1"},
                )
            },
            coords={
                "temperature_k": temperature_k,
                "wavelength_nm": wavelength_nm.astype(float),
            },
            attrs={
                "source": "Birk and Wagner O3 UV absorption cross sections",
                "doi": "10.5281/zenodo.1485588",
            },
        )

        ds.to_netcdf(output_file)


class O3SerdyuchenkoDatabase(_ZenodoO3Database):
    """
    Cached Serdyuchenko-Gorshelev O3 absorption cross-section database.

    The source file is a single wide ASCII table. Column 1 is vacuum wavelength
    in nm, followed by cross sections for temperatures from 293 K down to 193 K.
    The converted cache sorts those temperatures into ascending order so the
    generic absorber database can interpolate consistently.
    """

    _record_id = "5793207"
    _raw_filename = "SerdyuchenkoGorshelev5digits_latest.dat"
    _processed_filename = "serdyuchenko.nc"

    def _convert(self, raw_file: Path, output_file: Path) -> None:
        # The file header states "NUMBER OF HEADER LINES: 45".  np.loadtxt then
        # reads a rectangular array with wavelength in column 0 and temperature
        # cross sections in columns 1:.
        data = np.loadtxt(raw_file, skiprows=45)

        temperature_k = np.array(
            [
                293.0,
                283.0,
                273.0,
                263.0,
                253.0,
                243.0,
                233.0,
                223.0,
                213.0,
                203.0,
                193.0,
            ]
        )
        sort_idx = np.argsort(temperature_k)

        ds = xr.Dataset(
            {
                "xs": (
                    ["temperature_k", "wavelength_nm"],
                    data[:, 1:].T[sort_idx] / 1e4,
                    {"units": "m^2 molecule^-1"},
                )
            },
            coords={
                "temperature_k": temperature_k[sort_idx],
                "wavelength_nm": data[:, 0].astype(float),
            },
            attrs={
                "source": "Serdyuchenko-Gorshelev O3 UV/VIS/NIR absorption cross sections",
                "doi": "10.5281/zenodo.5793207",
            },
        )

        ds.to_netcdf(output_file)
