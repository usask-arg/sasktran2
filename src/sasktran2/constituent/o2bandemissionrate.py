from __future__ import annotations

import numpy as np

from sasktran2._core_rust import PyO2BandEmissionRate
from sasktran2.database.hitran_line import HITRANLineDatabase

from .base import Constituent


class O2BandEmissionRate(Constituent):
    """O2 emission with an independently specified band photon VER profile.

    The VER is interpolated to the atmospheric grid before distributing it
    among rotational lines at the current ``atmo.temperature_k``. Rotational
    populations follow that temperature; vibrational band VERs are independent.
    Lines are Doppler broadened using the O2 molecular mass.

    The constituent provides ``wf_<name>_photon_ver`` on its VER altitude grid.
    When temperature derivatives are enabled, rotational redistribution and
    Doppler broadening contribute to ``wf_temperature_k`` at fixed band VER.
    Add an O2 absorber separately to include self-absorption and its derivatives.

    Parameters
    ----------
    altitudes_m : np.ndarray
        Increasing altitudes in meters for the VER profile.
    photon_ver : np.ndarray
        Total photon VER of this band in photons m^-3 s^-1, integrated over all
        directions and all its lines. The source includes the 1 / (4 pi) factor.
    band : str
        ``"0-0"`` or ``"1-1"`` for the A-band, or ``"1-0"`` for the B-band.
    line_weight_model : str
        ``"einstein_a_branching"`` or ``"hitran_line_strength"``.
    db : HITRANLineDatabase, optional
        Database containing the O2 emission lines.
    out_of_bounds_mode : str
        ``"zero"`` (default) or ``"extend"`` outside the VER altitude grid.
    """

    def __init__(
        self,
        altitudes_m: np.ndarray,
        photon_ver: np.ndarray,
        band: str = "0-0",
        line_weight_model: str = "einstein_a_branching",
        db: HITRANLineDatabase | None = None,
        out_of_bounds_mode: str = "zero",
    ):
        db = db or HITRANLineDatabase()
        self._emission = PyO2BandEmissionRate(
            np.asarray(altitudes_m, dtype=np.float64),
            np.asarray(photon_ver, dtype=np.float64),
            db.path("O2").as_posix(),
            band,
            line_weight_model,
            out_of_bounds_mode,
        )

    @classmethod
    def _from_native(cls, emission):
        obj = cls.__new__(cls)
        obj._emission = emission
        return obj

    def add_to_atmosphere(self, atmo):
        self._emission.add_to_atmosphere(atmo)

    def register_derivative(self, atmo, name):
        self._emission.register_derivative(atmo, name)

    @property
    def photon_ver(self):
        """Mutable photon VER profile in photons m^-3 s^-1."""
        return self._emission.photon_ver

    @photon_ver.setter
    def photon_ver(self, values):
        self._emission.photon_ver = np.asarray(values, dtype=np.float64)

    @property
    def altitudes_m(self):
        """Copy of the fixed VER altitude grid in meters."""
        return self._emission.altitudes_m

    @property
    def wavelengths_nm(self):
        return self._emission.wavelengths_nm

    @property
    def band(self):
        return self._emission.band

    def line_weights(self, temperature_k):
        """Normalized line weights for each supplied temperature (K)."""
        return self._emission.line_weights(
            np.atleast_1d(temperature_k).astype(np.float64)
        )
