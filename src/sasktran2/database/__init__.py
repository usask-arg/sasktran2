import abc
from pathlib import Path

import xarray as xr

from sasktran2 import appconfig


class AbstractDatabase(abc.ABC):
    """
    Defines the interface for database classes within sasktran2
    """

    @abc.abstractmethod
    def path(self, key: str, **kwargs) -> Path | None:
        pass

    @abc.abstractmethod
    def load_ds(self, key: str, **kwargs) -> xr.Dataset:
        pass


class CachedDatabase(AbstractDatabase):
    def __init__(self, db_root: Path | None = None) -> None:
        """
        A CachedDatabase is a database that caches data on the local file system.
        The root directly can optionally be specified, otherwise the default
        folder is used.

        Parameters
        ----------
        db_root : Path, optional
            The root directory to store the database, by default None
        """
        super().__init__()

        if db_root is None:
            self._db_root = appconfig.database_root()
        else:
            self._db_root = db_root

        if not self._db_root.exists():
            self._db_root.mkdir(parents=True)

    @abc.abstractmethod
    def clear(self):
        """
        Deletes the entire database
        """
