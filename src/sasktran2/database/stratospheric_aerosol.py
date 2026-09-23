"""Versioned, checksum-verified SAGE-derived aerosol reference data."""

from __future__ import annotations

import hashlib
import shutil
import tempfile
from importlib.resources import files
from pathlib import Path

import xarray as xr

from .base import CachedDatabase

_FILES = {
    "v1": (
        "stratospheric_aerosol_v1.nc",
        "3c949e2eeaff85318de9dc197cf3bf717d87a1f2bf748a7667cc633a55e6557f",
    )
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class StratosphericAerosolDatabase(CachedDatabase):
    """Cache the compact catalogue shipped with SASKTRAN2.

    No source archive or network download is needed. Both the packaged file
    and existing cache entries are checked against the pinned SHA-256 digest.
    ``db_root`` overrides the usual SASKTRAN2 database root. Versions are
    independent of the SAGE and USask source-product versions.
    """

    def __init__(self, version: str = "v1", db_root: str | Path | None = None):
        if version not in _FILES:
            msg = f"Unsupported aerosol catalogue version {version!r}; available: {list(_FILES)}"
            raise ValueError(msg)
        self.version = version
        self._filename, self._checksum = _FILES[version]
        super().__init__(
            None if db_root is None else Path(db_root),
            rel_path=Path("climatology/stratospheric_aerosol") / version,
        )

    def _verify(self, path: Path):
        if _sha256(path) != self._checksum:
            msg = f"Aerosol catalogue checksum mismatch for {path}; remove the damaged cache entry to restore it"
            raise OSError(msg)

    def path(self, key: str = "", **kwargs) -> Path:
        """Return the verified catalogue, populating the cache on first use."""
        if key not in ("", self._filename):
            msg = f"Aerosol catalogue contains only {self._filename}, not {key!r}"
            raise ValueError(msg)
        destination = self._db_root / self._filename
        if destination.exists():
            self._verify(destination)
            return destination
        resource = files("sasktran2").joinpath(
            "_data", "stratospheric_aerosol", self._filename
        )
        temporary = None
        try:
            with (
                resource.open("rb") as source,
                tempfile.NamedTemporaryFile(dir=self._db_root, delete=False) as stream,
            ):
                temporary = Path(stream.name)
                shutil.copyfileobj(source, stream)
            self._verify(temporary)
            temporary.replace(destination)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        return destination

    def load_ds(self, key: str = "", **kwargs) -> xr.Dataset:
        """Load raw catalogue arrays and close the file."""
        with xr.open_dataset(self.path(key)) as data:
            return data.load()

    def clear(self):
        """Remove only this version's cached file."""
        (self._db_root / self._filename).unlink(missing_ok=True)
