from __future__ import annotations

import sasktran2 as sk
from sasktran2._core_rust import PyPopulationEmissionRate
from sasktran2.database.hitran_line import HITRANLineDatabase

from .base import Constituent
from .o2bandemissionrate import O2BandEmissionRate


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

        The current implementation supports O2 A-band and B-band emission from
        population profiles such as the output of
        :class:`sasktran2.photchem.Yankovsky`. The population dataset must
        contain ``altitude`` in meters, ``temperature`` in kelvin, and ``O2(b)``
        number density in m^-3. ``O2(b, v=1)`` and ``O2(b, v=2)`` are optional
        m^-3 number-density profiles.

        Internally, populations are multiplied by Einstein-A coefficients to
        obtain photon volume emission rates in photons m^-3 s^-1. The emitted
        source is assumed isotropic; the constituent applies the 1 / 4pi factor
        when adding the spectral source to the atmosphere.

        Each band's VER is interpolated onto the atmospheric grid, then its
        rotational line weights and Doppler widths are evaluated using the
        current atmospheric temperature. Temperature derivatives hold the
        supplied populations and band Einstein-A coefficients fixed.

        The inspection arrays (``photon_ver``, ``altitudes_m``,
        ``wavelengths_nm``, ``weights``, and the ``line_list_*`` methods) are
        read-only views of the construction-time A/B-band spectra using the
        dataset's temperature. Attempts to modify them raise ``ValueError``.
        Use :meth:`to_band_emissions` to obtain independent, mutable band VER
        constituents for retrievals.

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

    def to_band_emissions(self) -> dict[str, O2BandEmissionRate]:
        """Return independent band constituents initialized from populations.

        Keys are ``"0-0"``, ``"1-1"``, and ``"1-0"`` when their lines are
        available. Add these constituents instead of this population constituent
        to retrieve their VER profiles independently. Mutating the returned
        constituents does not change this object or the input population dataset.
        """
        return {
            emission.band: O2BandEmissionRate._from_native(emission)
            for emission in self._emission.to_band_emissions()
        }

    @property
    def photon_ver(self):
        """Read-only combined A-band photon VER at construction time."""
        return self._emission.photon_ver

    @property
    def altitudes_m(self):
        """Read-only altitude grid in meters."""
        return self._emission.altitudes_m

    @property
    def wavelengths_nm(self):
        """Read-only A-band line wavelengths in nanometers."""
        return self._emission.wavelengths_nm

    @property
    def weights(self):
        """Read-only A-band line weights at the input dataset's temperature."""
        return self._emission.weights

    @property
    def num_line_list_emissions(self):
        return self._emission.num_line_list_emissions

    def line_list_photon_ver(self, index=0):
        """Read-only photon VER of a construction-time line-list component."""
        return self._emission.line_list_photon_ver(index)

    def line_list_wavelengths_nm(self, index=0):
        """Read-only wavelengths of a construction-time line-list component."""
        return self._emission.line_list_wavelengths_nm(index)

    def line_list_weights(self, index=0):
        """Read-only weights of a construction-time line-list component."""
        return self._emission.line_list_weights(index)
