from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import PyLineListVolumeEmissionRate

from .base import Constituent


class LineListVolumeEmissionRate(Constituent):
    _ver: PyLineListVolumeEmissionRate

    def __init__(
        self,
        altitudes_m: np.ndarray,
        photon_ver: np.ndarray,
        wavelengths_nm: np.ndarray,
        weights: np.ndarray,
        molecular_mass_g_per_mol: float | None = None,
    ):
        """
        A spectrally resolved line-list volume emission rate.

        Lines are Doppler broadened using the atmosphere temperature on the model altitude grid.
        The Doppler width uses ``molecular_mass_g_per_mol`` (O2 by default).

        Parameters
        ----------
        altitudes_m : np.ndarray
            Altitudes in [m] that correspond to the volume emission rate profile.
        photon_ver : np.ndarray
            Total photon volume emission rate profile in photons m^-3 s^-1.
        wavelengths_nm : np.ndarray
            Emission line wavelengths in [nm].
        weights : np.ndarray
            Relative line weights, ``(line,)`` or ``(altitude, line)``. Values are
            normalized internally.
        molecular_mass_g_per_mol : float, optional
            Emitter molar mass for Doppler broadening; O2 (31.9988) by default.
        """
        self._ver = PyLineListVolumeEmissionRate(
            np.atleast_1d(altitudes_m).astype(np.float64),
            np.atleast_1d(photon_ver).astype(np.float64),
            np.atleast_1d(wavelengths_nm).astype(np.float64),
            np.atleast_1d(weights).astype(np.float64),
            molecular_mass_g_per_mol=molecular_mass_g_per_mol,
        )

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        self._ver.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        self._ver.register_derivative(atmo, name)

    @property
    def photon_ver(self):
        return self._ver.photon_ver

    @photon_ver.setter
    def photon_ver(self, photon_ver: np.ndarray):
        self._ver.photon_ver = photon_ver

    @property
    def altitudes_m(self):
        return self._ver.altitudes_m

    @property
    def wavelengths_nm(self):
        return self._ver.wavelengths_nm

    @property
    def weights(self):
        return self._ver.weights
