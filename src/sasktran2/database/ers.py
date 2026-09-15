"""Download and cache the versioned CAIRT Extended Reference Scenarios."""

from __future__ import annotations

import hashlib
import shutil
import tempfile
import urllib.request
from pathlib import Path

import xarray as xr

from .base import CachedDatabase

_RECORD_ID = "10022129"
_FILES = {"v07": ("CAIRT_ERS_v07.nc", "71a5ed74d7056538cfd6a99d20ca3599")}


def _md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class ERSDatabase(CachedDatabase):
    """Cache only the ERS NetCDF, verifying the published file checksum.

    Parameters
    ----------
    version : str, optional
        Pinned dataset version. Currently only ``"v07"`` is supported.
    db_root : Path or None, optional
        Database root; defaults to the configured SASKTRAN2 database directory.
        Files are stored under ``climatology/ers/<version>``.
    """

    def __init__(self, version: str = "v07", db_root: Path | None = None):
        if version not in _FILES:
            msg = f"Unsupported ERS version {version!r}; available versions: {list(_FILES)}"
            raise ValueError(msg)
        self.version = version
        self._filename, self._checksum = _FILES[version]
        super().__init__(
            None if db_root is None else Path(db_root),
            rel_path=Path("climatology/ers") / version,
        )

    def _verify(self, path: Path):
        if _md5(path) != self._checksum:
            msg = (
                f"ERS checksum mismatch for {path}; expected MD5 {self._checksum}. "
                "Remove the damaged cache file to download it again."
            )
            raise OSError(msg)

    def path(self, key: str = "", **kwargs) -> Path:
        """Return the verified NetCDF path, downloading it on first use."""
        if key not in ("", self._filename):
            msg = f"ERS contains only {self._filename}, not {key!r}"
            raise ValueError(msg)
        destination = self._db_root / self._filename
        if destination.exists():
            self._verify(destination)
            return destination

        url = f"https://zenodo.org/api/records/{_RECORD_ID}/files/{self._filename}/content"
        # A unique temporary file also makes concurrent first-use downloads safe.
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=self._db_root, delete=False) as stream:
                temporary = Path(stream.name)
                with urllib.request.urlopen(url, timeout=60) as response:
                    shutil.copyfileobj(response, stream)
            self._verify(temporary)
            temporary.replace(destination)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return destination

    def load_ds(self, key: str = "", **kwargs) -> xr.Dataset:
        """Load the raw arrays into memory and close the underlying file."""
        with xr.open_dataset(self.path(key), decode_timedelta=False) as dataset:
            return dataset.load()

    def clear(self):
        """Remove this version's cached NetCDF."""
        (self._db_root / self._filename).unlink(missing_ok=True)
