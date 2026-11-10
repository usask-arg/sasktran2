from __future__ import annotations

import importlib
import inspect
import shlex
from pathlib import Path

#: Retries after a failed file download, and the pause between them [s].
#: zenodo-get retries nothing by default, and a dropped connection mid-file is
#: otherwise fatal.
DOWNLOAD_RETRIES = 3
DOWNLOAD_RETRY_PAUSE_S = 10.0


def download_zenodo_record(record_id: str, output_dir: Path) -> None:
    """
    Download a Zenodo record into a local directory using zenodo-get.

    zenodo-get changed its public API between major versions. This helper supports
    both the old ``zenodo_get(argv)`` entry point and the newer ``download`` helper.
    Failed file downloads are retried :data:`DOWNLOAD_RETRIES` times.

    Raises
    ------
    ImportError
        If zenodo-get is not installed or has no supported download API.
    ConnectionError
        If Zenodo cannot be reached or a file still fails after the retries.
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

    try:
        if hasattr(zenodo_get, "download"):
            zenodo_get.download(
                str(record_id),
                output_dir=output_dir.as_posix(),
                **_supported(
                    zenodo_get.download,
                    retry_attempts=DOWNLOAD_RETRIES,
                    retry_pause=DOWNLOAD_RETRY_PAUSE_S,
                ),
            )
        elif hasattr(zenodo_get, "zenodo_get"):
            zenodo_get.zenodo_get(
                shlex.split(
                    f'--record {record_id} -o "{output_dir.as_posix()}" '
                    f"--retry {DOWNLOAD_RETRIES} --pause {DOWNLOAD_RETRY_PAUSE_S}"
                )
            )
        else:
            msg = (
                "Installed zenodo_get package does not expose a supported download API"
            )
            raise ImportError(msg)
    except (OSError, RuntimeError) as e:
        # zenodo-get 3 raises ConnectionError for unreachable metadata and
        # RuntimeError for a file that failed every attempt.
        msg = f"Could not download Zenodo record {record_id} to {output_dir}: {e}"
        raise ConnectionError(msg) from e


def _supported(function, **kwargs) -> dict:
    """The keyword arguments that ``function`` accepts."""
    try:
        parameters = inspect.signature(function).parameters
    except (TypeError, ValueError):
        return {}
    if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in parameters.values()):
        return kwargs
    return {k: v for k, v in kwargs.items() if k in parameters}
