from __future__ import annotations

import sasktran2 as sk
from sasktran2._core_rust import PyPopulationEmissionRate
from sasktran2.database.hitran_line import HITRANLineDatabase

from .base import Constituent


class PopulationEmissionRate(Constituent):
    _emission: PyPopulationEmissionRate

    def __init__(
        self,
        populations,
        species=("O2",),
        line_weight_model="einstein_a_branching",
        db: HITRANLineDatabase | None = None,
    ):
        """
        Photochemical population-to-emission constituent.

        The current implementation supports O2 A-band and B-band emission from population
        profiles such as the output of :class:`sasktran2.photchem.Yankovsky`.
        The population dataset must contain ``altitude``, ``temperature``, and
        ``O2(b)``. ``O2(b, v=1)`` and ``O2(b, v=2)`` are optional.

        Parameters
        ----------
        populations
            xarray Dataset containing population profiles on an altitude grid.
        species
            Species to include. Only ``"O2"`` is currently supported.
        line_weight_model
            ``"einstein_a_branching"`` or ``"hitran_line_strength"``.
        db
            HITRAN database used to load the O2 A-band line list.
        """
        db = db or HITRANLineDatabase()
        db_path = db.path("O2")
        if isinstance(species, str):
            species = [species]

        self._emission = PyPopulationEmissionRate(
            populations,
            db_path.as_posix(),
            list(species),
            line_weight_model,
        )

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        self._emission.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        self._emission.register_derivative(atmo, name)

    @property
    def photon_ver(self):
        return self._emission.photon_ver

    @property
    def altitudes_m(self):
        return self._emission.altitudes_m

    @property
    def wavelengths_nm(self):
        return self._emission.wavelengths_nm

    @property
    def weights(self):
        return self._emission.weights

    @property
    def num_line_list_emissions(self):
        return self._emission.num_line_list_emissions

    def line_list_photon_ver(self, index=0):
        return self._emission.line_list_photon_ver(index)

    def line_list_wavelengths_nm(self, index=0):
        return self._emission.line_list_wavelengths_nm(index)

    def line_list_weights(self, index=0):
        return self._emission.line_list_weights(index)
