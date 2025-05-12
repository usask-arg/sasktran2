from __future__ import annotations

from sasktran2._core_rust import openmp_support_enabled as py_openmp_support_enabled


def openmp_support_enabled() -> bool:
    """
    Check if OpenMP support is enabled in the current build of sasktran2.

    Returns:
        bool: True if OpenMP support is enabled, False otherwise.
    """
    return py_openmp_support_enabled()
