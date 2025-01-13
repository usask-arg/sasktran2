from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate

from sasktran2.database.web import StandardDatabase


class RefractiveIndex:
    def __init__(
        self, refractive_index_fn: Callable[[float], complex], identifier: str
    ) -> None:
        """
        A generic implementation of a refractive index function. This class is a light wrapper
        on top of a callable function that returns the complex refractive index for a given
        wavelength, as well as providing a unique identifier for the refractive index.

        Parameters
        ----------
        refractive_index_fn : Callable[[float], complex]
            Function that takes in a wavelength in nm and returns the complex refractive index
        identifier : str
            A unique identifier for the refractive index
        """
        self._fn = refractive_index_fn
        self._identifier = identifier

    def refractive_index(self, wavelength_nm: np.ndarray) -> np.ndarray:
        """
        Returns the complex refractive index for a given wavelength

        Parameters
        ----------
        wavelength_nm : np.ndarray

        Returns
        -------
        np.ndarray
        """
        return self._fn(wavelength_nm)

    @property
    def refractive_index_fn(self):
        """
        Get the function that returns the complex refractive index
        """
        return self._fn

    @property
    def identifier(self):
        """
        Get the unique identifier for this refractive index
        """
        return self._identifier


def _from_osiris_file(path: Path):
    """
    Reads in data from the old OSIRIS refractive index files and returns a callable function
    """
    data = pd.read_csv(path.as_posix(), header=None)

    return interpolate.interp1d(
        data.values[:, 0], data.values[:, 1] - 1j * data.values[:, 2]
    )


class H2SO4(RefractiveIndex):
    def __init__(self, source: str = "osiris") -> None:
        """
        A refractive index for H2SO4. The default source is the old OSIRIS data.
        """
        if source.lower() == "osiris":
            super().__init__(
                _from_osiris_file(
                    StandardDatabase().path("refractive_index/refrac_h2so4_osiris.txt")
                ),
                "H2SO4_osiris",
            )
        else:
            msg = "Only osiris source is supported for H2SO4 refractive index data"
            raise ValueError(msg)


class Dust(RefractiveIndex):
    def __init__(self, source: str = "osiris") -> None:
        """
        A refractive index for Dust. The default source is the old OSIRIS data.
        """
        if source.lower() == "osiris":
            super().__init__(
                _from_osiris_file(
                    StandardDatabase().path("refractive_index/refrac_dust_osiris.txt")
                ),
                "dust_osiris",
            )
        else:
            msg = "Only osiris source is supported for dust refractive index data"
            raise ValueError(msg)


class Ice(RefractiveIndex):
    def __init__(self, source: str = "osiris") -> None:
        """
        A refractive index for ice. The default source is the old OSIRIS data.
        """
        if source.lower() == "osiris":
            super().__init__(
                _from_osiris_file(
                    StandardDatabase().path("refractive_index/refrac_ice_osiris.txt")
                ),
                "ice_osiris",
            )
        else:
            msg = "Only osiris source is supported for ice refractive index data"
            raise ValueError(msg)


class Water(RefractiveIndex):
    def __init__(self, source: str = "osiris") -> None:
        """
        A refractive index for Water. The default source is the old OSIRIS data.
        """
        if source.lower() == "osiris":
            super().__init__(
                _from_osiris_file(
                    StandardDatabase().path("refractive_index/refrac_water_osiris.txt")
                ),
                "water_osiris",
            )
        else:
            msg = "Only osiris source is supported for water refractive index data"
            raise ValueError(msg)
