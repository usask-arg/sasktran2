from __future__ import annotations

import importlib
import shlex
from pathlib import Path


def download_zenodo_record(record_id: str, output_dir: Path) -> None:
    """
    Download a Zenodo record into a local directory using zenodo-get.

    zenodo-get changed its public API between major versions. This helper supports
    both the old ``zenodo_get(argv)`` entry point and the newer ``download`` helper.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        zenodo_get = importlib.import_module("zenodo_get")
    except ImportError as e:
        msg = (
            "zenodo_get is required to download this database from Zenodo. "
            "Install the optional dependency with `pip install sasktran2[zenodo]` "
            "or install `zenodo-get` directly."
        )
        raise ImportError(msg) from e

    if hasattr(zenodo_get, "download"):
        zenodo_get.download(str(record_id), output_dir=output_dir.as_posix())
    elif hasattr(zenodo_get, "zenodo_get"):
        zenodo_get.zenodo_get(
            shlex.split(f'--record {record_id} -o "{output_dir.as_posix()}"')
        )
    else:
        msg = "Installed zenodo_get package does not expose a supported download API"
        raise ImportError(msg)
