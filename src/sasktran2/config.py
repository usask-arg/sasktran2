from __future__ import annotations

from sasktran2._core_rust import PyConfig


class Config:
    _config: PyConfig

    def __init__(self):
        """
        Object which stores all of the configuration settings for the radiative transfer calculation.
        """
        self._config = PyConfig()
