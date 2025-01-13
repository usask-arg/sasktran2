from __future__ import annotations

import logging
import urllib.request
import zipfile
from pathlib import Path

import xarray as xr

from .base import CachedDatabase


class WebDatabase(CachedDatabase):
    def __init__(
        self, url: str, db_root: Path | None = None, rel_path: Path | None = None
    ) -> None:
        """
        A web database is a database that is downloaded from a URL and stored.  It must
        consist of a single file.

        Parameters
        ----------
        url : str
            URL Location for the file
        db_root : Path | None, optional
            Database root directory. Optional, by default None which will set it to the config file database_root
        rel_path : Path | None, optional
            Relative path to place the file inside the db_root, by default None which will place it in the root.
        """
        super().__init__(db_root)
        self._url = url

        # Get the directory to place the data files
        if rel_path is not None:
            self._data_directory = self._db_root.joinpath(rel_path)
        else:
            self._data_directory = self._db_root

        self._data_directory.mkdir(parents=True, exist_ok=True)

        # Get the filename from the url
        self._filename = Path(url).name

    def _post_process(self):
        """
        Called after the file is downloaded. Used by derived classes to post process the file
        """

    def load(self, **kwargs):
        """
        Downloads the database from the URL if it doesn't exist
        """
        output_file = self._data_directory.joinpath(self._filename)

        if output_file.exists():
            # Already have the database downloaded
            return

        try:
            urllib.request.urlretrieve(
                self._url,
                filename=output_file.as_posix(),
            )
        except Exception as e:
            logging.exception(e)

        self._post_process()

    def clear(self):
        output_file = self._data_directory.joinpath(self._filename)

        if output_file.exists():
            output_file.unlink()

    def output_file(self):
        return self._data_directory.joinpath(self._filename)

    def path(self, key: str, **kwargs) -> Path | None:
        return self._data_directory.joinpath(key)

    def load_ds(self, key: str, **kwargs):
        return xr.open_dataset(self.path(key, **kwargs))


class ZipWebDatabase(WebDatabase):
    def __init__(
        self, url: str, db_root: Path | None = None, rel_path: Path | None = None
    ) -> None:
        """
        A web database that is a zip file. It is downloaded from a URL and then extracted.

        Parameters
        ----------
        url : str
            url of the zip file
        db_root : Path | None, optional
            Database root directory. Optional, by default None which will set it to the config file database_root
        rel_path : Path | None, optional
            Relative path to place the file inside the db_root, by default None which will place it in the root.
        """
        super().__init__(url, db_root, rel_path)

    def _post_process(self):
        output_file = self._data_directory.joinpath(self._filename)

        try:
            with zipfile.ZipFile(output_file.as_posix(), "r") as zip_ref:
                zip_ref.extractall(self._db_root)
        except Exception as e:
            logging.exception(e)


class StandardDatabase(CachedDatabase):
    def __init__(self, version: str = "latest") -> None:
        """
        Standard databases that are downloaded from the sasktran website.

        Parameters
        ----------
        version : str, optional
            Specific version to use, by default "latest"
        """
        self._url = f"https://arg.usask.ca/sasktranfiles/sasktran2_db/v_{version}/"

    def path(self, key: str, **kwargs) -> Path | None:
        db = WebDatabase(self._url + key, rel_path=Path(key).parent)
        db.load()

        return db.output_file()

    def load_ds(self, key: str, **kwargs):
        return xr.open_dataset(self.path(key, **kwargs))

    def clear(self):
        msg = "The standard database can only be cleared manually"
        raise NotImplementedError(msg)
