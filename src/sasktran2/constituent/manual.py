from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import PyManual

from .base import Constituent


class Manual(Constituent):
    _inner: PyManual

    def __init__(
        self,
        extinction: np.ndarray,
        ssa: np.ndarray,
        legendre_moments: np.ndarray | None = None,
        delta_scale: bool = False,
    ) -> None:
        """
        An implementation of a manual constituent where the user provides the extinction, single
        scattering albedo, and optionally the Legendre moments.

        Note that this is manual in the sense that the user provides all necessary atmospheric
        quantities directly, on the model grid.  No interpolation is done between levels, or between
        "wavelength" calculations.

        The legendre_moments must also be provided with the same number of moments as the model,
        including the batching of spherical legendre moments if using multiple stokes parameters.

        Parameters
        ----------
        extinction : numpy.ndarray
            Extinction cross section in m^-1. Shape [num_altitudes, num_wavelengths]
        ssa : numpy.ndarray
            Single scattering albedo (unitless). Shape [num_altitudes, num_wavelengths]
        legendre_moments : numpy.ndarray | None, optional
            Legendre moments (unitless), by default None. Shape [num_moments, num_altitudes, num_wavelengths]
        delta_scale : bool, optional
            Whether to apply delta-scaling to the scattering properties, by default False.

        """
        if delta_scale:
            f = legendre_moments[::4][-1] / (2.0 * legendre_moments.shape[0] / 4 + 1)
            extinction *= 1.0 - ssa * f
            ssa *= (1.0 - f) / (1.0 - ssa * f)

            for i in range(int(legendre_moments.shape[0] / 4)):
                # a1
                legendre_moments[4 * i] -= f * (2.0 * i + 1)
                legendre_moments[4 * i] /= 1.0 - f

                # a2
                legendre_moments[4 * i + 1] -= f * (2.0 * i + 1)
                legendre_moments[4 * i + 1] /= 1.0 - f

                # a3
                legendre_moments[4 * i + 2] -= f * (2.0 * i + 1)
                legendre_moments[4 * i + 2] /= 1.0 - f

                # b1
                legendre_moments[4 * i + 3] /= 1.0 - f

        self._inner = PyManual(
            extinction=extinction,
            ssa=ssa,
            legendre_moments=legendre_moments,
        )

    @property
    def extinction(self) -> np.ndarray:
        return self._inner.extinction

    @property
    def ssa(self) -> np.ndarray:
        return self._inner.ssa

    @property
    def leg_coeff(self) -> np.ndarray | None:
        return self._inner.leg_coeff

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere


        :meta private:
        """
        self._inner.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        self._inner.register_derivative(atmo, name)
