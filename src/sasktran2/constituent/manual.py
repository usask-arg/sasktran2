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
            Extinction cross section in m^-1. Shape ``(altitude, wavelength)``
            or ``(horizontal, altitude, wavelength)``. Altitude-only input is
            broadcast horizontally in a 2D atmosphere.
        ssa : numpy.ndarray
            Single scattering albedo with the same shape as ``extinction``.
        legendre_moments : numpy.ndarray | None, optional
            Legendre moments, by default None. Shape
            ``(moment, *extinction.shape)``.
        delta_scale : bool, optional
            Whether to apply delta-scaling to the scattering properties, by default False.

        """
        extinction = np.asarray(extinction)
        ssa = np.asarray(ssa)
        if extinction.shape != ssa.shape:
            msg = "extinction and ssa must have the same shape"
            raise ValueError(msg)
        if extinction.ndim not in (2, 3):
            msg = (
                "extinction and ssa must have shape (altitude, wavelength) or "
                "(horizontal, altitude, wavelength)"
            )
            raise ValueError(msg)

        self._input_volume_shape = (
            extinction.shape[:-1] if extinction.ndim == 3 else None
        )

        if legendre_moments is not None:
            legendre_moments = np.asarray(legendre_moments)
            if (
                legendre_moments.ndim != extinction.ndim + 1
                or legendre_moments.shape[1:] != extinction.shape
            ):
                msg = (
                    "legendre_moments must have shape (moment, *extinction.shape); "
                    f"got {legendre_moments.shape} for extinction shape "
                    f"{extinction.shape}"
                )
                raise ValueError(msg)

        if delta_scale and legendre_moments is None:
            msg = "legendre_moments must be provided when delta_scale is enabled"
            raise ValueError(msg)

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

        flat_extinction = extinction.reshape(-1, extinction.shape[-1])
        flat_ssa = ssa.reshape(-1, ssa.shape[-1])
        flat_legendre = (
            None
            if legendre_moments is None
            else legendre_moments.reshape(
                legendre_moments.shape[0], -1, legendre_moments.shape[-1]
            )
        )

        self._inner = PyManual(
            extinction=flat_extinction,
            ssa=flat_ssa,
            legendre_moments=flat_legendre,
        )

    @property
    def volume_spatial_mode(self) -> str:
        return (
            "native_2d" if self._input_volume_shape is not None else "altitude_profile"
        )

    @property
    def extinction(self) -> np.ndarray:
        extinction = self._inner.extinction
        if self._input_volume_shape is None:
            return extinction
        return extinction.reshape((*self._input_volume_shape, extinction.shape[-1]))

    @property
    def ssa(self) -> np.ndarray:
        ssa = self._inner.ssa
        if self._input_volume_shape is None:
            return ssa
        return ssa.reshape((*self._input_volume_shape, ssa.shape[-1]))

    @property
    def leg_coeff(self) -> np.ndarray | None:
        leg_coeff = self._inner.leg_coeff
        if leg_coeff is None or self._input_volume_shape is None:
            return leg_coeff
        return leg_coeff.reshape(
            (leg_coeff.shape[0], *self._input_volume_shape, leg_coeff.shape[-1])
        )

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere


        :meta private:
        """
        extinction = self.extinction
        if extinction.shape[-1] != atmo.num_wavel:
            msg = (
                "Manual constituent wavelength dimension does not match the atmosphere: "
                f"{extinction.shape[-1]} != {atmo.num_wavel}"
            )
            raise ValueError(msg)

        leg_coeff = self.leg_coeff
        if (
            leg_coeff is not None
            and leg_coeff.shape[0] != atmo.storage.leg_coeff.shape[0]
        ):
            msg = (
                "Manual constituent Legendre moment dimension does not match the "
                f"atmosphere: {leg_coeff.shape[0]} != "
                f"{atmo.storage.leg_coeff.shape[0]}"
            )
            raise ValueError(msg)

        if self._input_volume_shape is not None:
            if tuple(atmo.volume_shape) != tuple(self._input_volume_shape):
                msg = (
                    "Native 2D Manual constituent shape does not match the atmosphere: "
                    f"{self._input_volume_shape} != {atmo.volume_shape}"
                )
                raise ValueError(msg)
        elif extinction.shape[0] != atmo.volume_shape[-1]:
            msg = (
                "Altitude-profile Manual constituent does not match the atmosphere: "
                f"{extinction.shape[0]} altitudes != {atmo.volume_shape[-1]}"
            )
            raise ValueError(msg)

        self._inner.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        self._inner.register_derivative(atmo, name)
