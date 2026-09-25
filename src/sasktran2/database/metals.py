from __future__ import annotations

from pathlib import Path

import xarray as xr

from sasktran2 import appconfig
from sasktran2.database.base import AbstractDatabase
from sasktran2.database.web import StandardDatabase


class MetalSpectroscopyDatabase(AbstractDatabase):
    """Cached spectroscopy for metal resonance and volume-emission properties.

    Files live below ``spectroscopy/metals/{atomic,molecular,emission}`` in the normal
    SASKTRAN2 database. Missing files are obtained through :class:`StandardDatabase`.
    An explicit ``db_root`` instead selects a local-only database, useful for
    offline work and independently prepared spectroscopy.
    See ``tools/spectroscopy`` for reproducible data preparation and provenance.
    """

    def __init__(self, db_root: str | Path | None = None):
        self._root = Path(db_root) if db_root is not None else appconfig.database_root()
        self._download_missing = db_root is None

    def path(self, key: str, *, kind: str = "atomic", **kwargs) -> Path:
        if kind not in {"atomic", "molecular", "emission"}:
            msg = "kind must be 'atomic', 'molecular', or 'emission'"
            raise ValueError(msg)
        if not key or any(
            c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"
            for c in key
        ):
            msg = "A species key must contain only letters, digits, and underscores"
            raise ValueError(msg)
        relative = Path("spectroscopy") / "metals" / kind / f"{key}.nc"
        path = self._root / relative
        if not path.is_file() and self._download_missing:
            path = StandardDatabase().path(relative.as_posix())
        if not path.is_file():
            msg = (
                f"Metal spectroscopy data are not installed at {path}. "
                "Prepare them with tools/spectroscopy, or provide db_filepath. "
                "Not every screened species has a usable spectroscopy file."
            )
            raise FileNotFoundError(msg)
        return path

    def load_ds(self, key: str, **kwargs) -> xr.Dataset:
        with xr.open_dataset(self.path(key, **kwargs)) as dataset:
            return dataset.load()

    def available_species(self, *, kind: str = "atomic") -> list[str]:
        """List locally installed species keys; availability is not detectability."""
        if kind not in {"atomic", "molecular", "emission"}:
            msg = "kind must be 'atomic', 'molecular', or 'emission'"
            raise ValueError(msg)
        folder = self._root / "spectroscopy" / "metals" / kind
        return sorted(path.stem for path in folder.glob("*.nc"))
