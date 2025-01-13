from __future__ import annotations

import itertools
import os

_counter = itertools.count()


def parameterized(names, params):
    def decorator(func):
        func.param_names = names
        func.params = params
        return func

    return decorator


def requires_dask():
    try:
        import dask  # noqa: F401
    except ImportError:
        raise NotImplementedError()  # noqa: B904


def _skip_slow():
    """
    Use this function to skip slow or highly demanding tests.

    Use it as a `Class.setup` method or a `function.setup` attribute.

    Examples
    --------
    >>> from . import _skip_slow
    >>> def time_something_slow():
    ...     pass
    ...
    >>> time_something.setup = _skip_slow
    """
    if os.environ.get("ASV_SKIP_SLOW", "0") == "1":
        msg = "Skipping this test..."
        raise NotImplementedError(msg)
