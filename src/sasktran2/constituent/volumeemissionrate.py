from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import PyMonochromaticVolumeEmissionRate

from .base import Constituent


class MonochromaticVolumeEmissionRate(Constituent):
    _ver: PyMonochromaticVolumeEmissionRate

    def __init__(
        self,
        altitudes_m: np.ndarray,
        ver: np.ndarray,
        wavelength_nm: float,
        *,
        out_of_bounds_mode: str = "zero",
        line_shape: str = "delta",
        emitter_molecular_weight_g_per_mol: float | None = None,
    ):
        """
        A monochromatic volume emission rate.

        By default this is represented as a delta function at the specified wavelength.  Internally
        SASKTRAN2 will scale the volume emission rate by the wavelength bin width when calculating radiance so that integrals over the
        spectral radiance band produce the correct total radiance. If ``line_shape="doppler"``, the emitted line is instead broadened
        with a normalized thermal Doppler profile using the atmosphere temperature and the supplied emitter molecular weight.

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
        out_of_bounds_mode : str
            Altitude interpolation mode, either ``"zero"`` or ``"extend"``.
        line_shape : str
            Spectral line shape, either ``"delta"`` or ``"doppler"``.
        emitter_molecular_weight_g_per_mol : float | None
            Emitter molecular weight in [g mol^-1]. Required when ``line_shape="doppler"``.
        """
        self._ver = PyMonochromaticVolumeEmissionRate(
            np.atleast_1d(altitudes_m).astype(np.float64),
            np.atleast_1d(ver).astype(np.float64),
            np.float64(wavelength_nm),
            out_of_bounds_mode,
            line_shape,
            emitter_molecular_weight_g_per_mol,
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
