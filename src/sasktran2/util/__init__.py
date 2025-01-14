# ruff: noqa: F401
from __future__ import annotations

import contextlib
import importlib.util
import io
import logging

from .._core import WignerD
from . import interpolation

_cached_hapi = None  # Cache variable to store the result of the import


def get_hapi():
    global _cached_hapi  # Use the cached result if available  # noqa: PLW0603
    if _cached_hapi is not None:
        return _cached_hapi

    try:
        # Suppress print messages during import
        with contextlib.redirect_stdout(io.StringIO()):
            import hapi  # Attempt to import the package
        _cached_hapi = hapi  # Cache the imported module
    except ImportError:
        _cached_hapi = None  # Cache the failure
        logging.error(
            "Optional dependency 'hapi' is not installed. Install it to use this functionality."
        )
        raise  # Re-raise the exception

    return _cached_hapi
