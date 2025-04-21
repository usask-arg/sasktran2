from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import PyRayleigh

from .base import Constituent


class Rayleigh(Constituent):
    _rayleigh: PyRayleigh

    def __init__(
        self,
        method: str = "bates",
        n2_percentage: float = 78.084,
        o2_percentage: float = 20.946,
        ar_percentage: float = 0.934,
        co2_percentage: float = 0.036,
        wavelengths_nm: np.ndarray = None,
        xs: np.ndarray = None,
        king_factor: np.ndarray = None,
    ) -> None:
        """
        An implementation of Rayleigh scattering.  Cross sections (and depolarization factors) can be
        calculated multiple ways, with the default method being that of 'bates'.

        Rayleigh scattering number density is estimated through the ideal gas law.

        This Constituent requires that the atmosphere object have `temperature_k`, `pressure_pa`, and
        `wavelength_nm` are all defined inside the :py:class:`sasktran2.Atmosphere` object.

        Parameters
        ----------
        method : str, default='bates'
            Method to use to calculate the cross section.  Supported methods are
            ['bates', 'manual'], by default 'bates'
        n2_percentage : float, optional
            Percentage of N2 in the atmosphere, by default 78.084
        o2_percentage : float, optional
            Percentage of O2 in the atmosphere, by default 20.946
        ar_percentage : float, optional
            Percentage of Ar in the atmosphere, by default 0.934
        co2_percentage : float, optional
            Percentage of CO2 in the atmosphere, by default 0.036
        wavelengths_nm : numpy.ndarray
            Wavelengths in nm to use for the cross section
        xs : numpy.ndarray
            Cross section in m2/molecule to use for the cross section
        king_factor : numpy.ndarray
            King factor to use for the cross section

        Raises
        ------
        ValueError
            If input method is not supported
        """
        self._rayleigh = PyRayleigh(
            method=method.lower(),
            n2_percentage=n2_percentage,
            o2_percentage=o2_percentage,
            ar_percentage=ar_percentage,
            co2_percentage=co2_percentage,
            wavelengths_nm=(
                wavelengths_nm.astype(np.float64)
                if wavelengths_nm is not None
                else None
            ),
            xs=xs.astype(np.float64) if xs is not None else None,
            king_factor=(
                king_factor.astype(np.float64) if king_factor is not None else None
            ),
        )

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere


        :meta private:
        """
        self._rayleigh.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        self._rayleigh.register_derivative(atmo, name)
