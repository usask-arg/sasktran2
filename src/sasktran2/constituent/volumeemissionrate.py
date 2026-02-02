from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import PyMonochromaticVolumeEmissionRate

from .base import Constituent


class MonochromaticVolumeEmissionRate(Constituent):
    _ver: PyMonochromaticVolumeEmissionRate

    def __init__(self, altitudes_m: np.ndarray, ver: np.ndarray, wavelength_nm: float):
        """
        A monochromatic volume emission rate (assumed to be a delta function at the specified wavelength).

        Because the volume emission rate is monochromatic, it only contributes to radiance at the specified wavelength.  Internally
        SASKTRAN2 will scale the volume emission rate by the wavelength bin width when calculating radiance so that integrals over the
        spectral radiance band produce the correct total radiance.

        Parameters
        ----------
        altitudes_m : np.ndarray
            Altitudes in [m] that correspond to the volume emission rate profile. Shape [N].
        ver : np.ndarray
            Volume emission rate profile. Units should be in (power / area / m).  Common units are ph cm^-2 / m s^-1 or
            W / m^3.  Note that when mixing photochemical calculations with solar sources, the units must match.  The default
            units for the solar spectrum are W / m^2 / nm, so VER in W / m^3 is recommended.  Also note that the VER is NOT
            per steradian.  A 4pi factor is applied internally when calculating radiance. Shape [N].
        wavelength_nm : float
            Wavelength in [nm] at which the volume emission rate is emitted.
        """
        self._ver = PyMonochromaticVolumeEmissionRate(
            np.atleast_1d(altitudes_m).astype(np.float64),
            np.atleast_1d(ver).astype(np.float64),
            np.float64(wavelength_nm),
        )

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere


        :meta private:
        """
        self._ver.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        self._ver.register_derivative(atmo, name)

    @property
    def ver(self):
        return self._ver.ver

    @ver.setter
    def ver(self, ver: np.array):
        self._ver.ver = ver

    @property
    def altitudes_m(self):
        return self._ver.altitudes_m
