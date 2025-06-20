from __future__ import annotations

import numpy as np

from sasktran2._core_rust import LineDatabaseType, PyLineAbsorber
from sasktran2.database.base import CachedDatabase
from sasktran2.database.hitran_line import HITRANLineDatabase
from sasktran2.optical.base import OpticalProperty
from sasktran2.util import get_hapi


class LineAbsorber(OpticalProperty):
    _internal: PyLineAbsorber

    def __init__(
        self,
        db_type: LineDatabaseType,
        db: CachedDatabase,
        molecule: str,
        cull_factor: float = 0.0,
        line_coupling=False,
        line_contribution_width: float = 25.0,
    ):
        hapi = get_hapi()

        # Ensure the db is downloaded
        db_path = db.path(molecule)

        self._internal = PyLineAbsorber(
            db_type,
            molecule,
            db_path.as_posix(),
            cull_factor,
            line_coupling,
            hapi.partitionSum,
            hapi.molecularMass,
            line_contribution_width,
        )

    def atmosphere_quantities(self, atmo, **kwargs):
        return self._internal.atmosphere_quantities(atmo)

    def _into_rust_object(self):
        return self._internal

    def cross_sections(self, wavelengths_nm, _altitudes_m, **kwargs):
        pressure_pa = kwargs.get("pressure_pa")
        if pressure_pa is None:
            msg = "Pressure must be provided in the kwargs"
            raise ValueError(msg)

        temperature_k = kwargs.get("temperature_k")
        if temperature_k is None:
            msg = "Temperature must be provided in the kwargs"
            raise ValueError(msg)

        p_self = kwargs.get("p_self")
        if p_self is None:
            p_self = np.zeros_like(pressure_pa)

        return self._internal.cross_section(
            1e7 / np.atleast_1d(wavelengths_nm).astype(float),
            np.atleast_1d(temperature_k).astype(float),
            np.atleast_1d(pressure_pa).astype(float),
            np.atleast_1d(p_self).astype(float),
        )


class HITRANAbsorber(LineAbsorber):
    def __init__(self, molecule: str, **kwargs):
        super().__init__(
            LineDatabaseType.HITRAN, HITRANLineDatabase(), molecule, **kwargs
        )
