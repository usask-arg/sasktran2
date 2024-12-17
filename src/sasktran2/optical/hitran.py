import logging

import numpy as np

from sasktran2.atmosphere import Atmosphere
from sasktran2.database.hitran import HITRANLineDatabase
from sasktran2.optical.base import OpticalProperty, OpticalQuantities
from sasktran2.spectroscopy import voigt_broaden


class HITRANAbsorber(OpticalProperty):
    def __init__(self, molecule: str, **kwargs):
        self._line_db = HITRANLineDatabase().load_ds(molecule)
        self._molecule = molecule
        self._kwargs = kwargs

        if self._molecule.lower() == "co2":
            # CO2 is bugged
            self._line_db.local_iso_id.values[
                self._line_db.local_iso_id.values == 0
            ] = 10

    def atmosphere_quantities(self, atmo: Atmosphere, **kwargs) -> OpticalQuantities:
        try:
            from hapi import molecularMass, partitionSum
        except ImportError:
            msg = "hapi is required to calculate HITRAN cross sections"
            raise ImportError(msg)  # noqa: B904

        # Prepare inputs

        iso_keys = np.unique(self._line_db.local_iso_id.to_numpy())

        partition_ratio = np.zeros((len(atmo.temperature_k), len(iso_keys)))
        molecular_mass = np.zeros(len(iso_keys))

        for key in iso_keys:
            vals = partitionSum(
                self._line_db["molec_id"].to_numpy()[0],
                key,
                [*list(atmo.temperature_k), 296],
            )
            partition_ratio[:, key - 1] = vals[:-1] / vals[-1]
            molecular_mass[key - 1] = molecularMass(
                self._line_db["molec_id"].to_numpy()[0], key
            )

        sidx = np.argsort(atmo.wavenumbers_cminv)

        result = OpticalQuantities(
            extinction=np.zeros(
                (len(atmo.model_geometry.altitudes()), len(atmo.wavelengths_nm))
            ),
            ssa=np.zeros(
                (len(atmo.model_geometry.altitudes()), len(atmo.wavelengths_nm))
            ),
        )

        result2 = np.zeros(
            (len(atmo.wavelengths_nm), (len(atmo.model_geometry.altitudes()))),
            order="f",
        )

        logging.debug(f"Starting Broadening for {self._molecule}")  # noqa: G004
        voigt_broaden(
            self._line_db.nu.to_numpy(),
            self._line_db.sw.to_numpy(),
            self._line_db.elower.to_numpy(),
            self._line_db.gamma_air.to_numpy(),
            self._line_db.gamma_self.to_numpy(),
            self._line_db.delta_air.to_numpy(),
            self._line_db.n_air.to_numpy(),
            self._line_db.local_iso_id.to_numpy(),
            partition_ratio,
            molecular_mass,
            atmo.pressure_pa / 101325.0,
            np.ones(len(atmo.pressure_pa)) * 0,
            atmo.temperature_k,
            atmo.wavenumbers_cminv[sidx],
            result2,
            num_threads=atmo._config.num_threads,
            **self._kwargs,
        )

        result.extinction[:, sidx] = result2.T / 1e4

        logging.debug(f"Finished Broadening for {self._molecule}")  # noqa: G004

        return result

    def optical_derivatives(self, atmo: Atmosphere, **kwargs) -> dict:  # noqa: ARG002
        return {}

    def cross_sections(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs  # noqa: ARG002
    ) -> OpticalQuantities:
        msg = "Not Supported"
        raise NotImplementedError(msg)

    def cross_section_derivatives(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs  # noqa: ARG002
    ) -> dict:
        return {}
