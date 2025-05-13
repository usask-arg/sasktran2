from __future__ import annotations

import logging

import numpy as np

from sasktran2._core_rust import (
    voigt_broaden,
    voigt_broaden_uniform,
    voigt_broaden_with_line_coupling,
)
from sasktran2.atmosphere import Atmosphere
from sasktran2.database.hitran_line import HITRANLineDatabase
from sasktran2.optical.base import OpticalProperty, OpticalQuantities
from sasktran2.util import get_hapi


class LineAbsorber(OpticalProperty):
    def __init__(
        self, molecule: str, line_db, line_fn="voigt", line_coupling=False, **kwargs
    ):
        """
        Absorption cross sections calculated using discrete absorption lines from a HITRAN-like database
        and broadened using the internal SASKTRAN2 Voigt line shape function.

        Notes:

        A line database is not included with SASKTRAN2 and must be downloaded separately. The hitran-api package
        is required to use HITRAN lines, and the zenodo_get package is required to download the AER line database.

        Broadening is calculated "online" rather than using pre-computed tables. This is slower but allows for
        more flexibility in the calculation.

        Derivatives with respect to pressure/temperature are currently not supported

        Parameters
        ----------
        molecule : str
            HITRAN Molecule identifier, e.g. "CO2", "H2O", "O3", "CH4"
        line_db : LineDatabase
            Line database to use, one of HITRANLineDatabase or AERLineDatabase
        line_fn : str
            Line shape function to use, currently only "voigt" is supported
        kwargs : dict
            Additional keyword arguments to pass to the line broadening function
        """
        defaults = {
            "line_contribution_width": 25.0,
            "cull_factor": 0.0,
            "subtract_pedastal": False,
        }

        self._line_db = line_db.load_ds(molecule)
        self._molecule = molecule
        self._kwargs = defaults
        self._kwargs.update(kwargs)
        self._line_coupling = line_coupling

        if line_fn != "voigt":
            msg = "Only Voigt line shape is supported"
            raise ValueError(msg)

        if self._molecule.lower() == "co2":
            # CO2 is bugged
            self._line_db.local_iso_id.values[
                self._line_db.local_iso_id.values == 0
            ] = 10

    def atmosphere_quantities(self, atmo: Atmosphere, **kwargs) -> OpticalQuantities:
        result = OpticalQuantities(
            extinction=np.zeros(
                (len(atmo.model_geometry.altitudes()), len(atmo.wavelengths_nm)),
                order="F",
            ),
            ssa=np.zeros(
                (len(atmo.model_geometry.altitudes()), len(atmo.wavelengths_nm)),
                order="F",
            ),
        )

        result.extinction[:] = self.cross_sections(
            atmo.wavelengths_nm,
            atmo.model_geometry.altitudes(),
            temperature_k=atmo.temperature_k,
            pressure_pa=atmo.pressure_pa,
            num_threads=atmo._config.num_threads,
            **kwargs,
        )
        return result

    def optical_derivatives(self, atmo: Atmosphere, **kwargs) -> dict:  # noqa: ARG002
        return {}

    def cross_sections(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs  # noqa: ARG002
    ) -> OpticalQuantities:
        if "temperature_k" not in kwargs:
            msg = "Temperature must be provided to calculate the cross sections"
            raise ValueError(msg)
        if "pressure_pa" not in kwargs:
            msg = "Pressure must be provided to calculate the cross sections"
            raise ValueError(msg)

        num_threads = kwargs.get("num_threads", 1)

        hapi = get_hapi()

        wavenumbers_cminv = 1e7 / wavelengths_nm
        temperature_k = kwargs["temperature_k"]
        pressure_pa = kwargs["pressure_pa"]

        iso_keys = np.unique(self._line_db.local_iso_id.to_numpy())

        partition_ratio = np.zeros((len(temperature_k), len(iso_keys)))
        molecular_mass = np.zeros(len(iso_keys))

        for key in iso_keys:
            vals = hapi.partitionSum(
                self._line_db["molec_id"].to_numpy()[0],
                key,
                [*list(temperature_k), 296],
            )
            partition_ratio[:, key - 1] = vals[:-1] / vals[-1]
            molecular_mass[key - 1] = hapi.molecularMass(
                self._line_db["molec_id"].to_numpy()[0], key
            )

        sidx = np.argsort(wavenumbers_cminv)

        result = np.zeros(
            (len(pressure_pa), (len(wavenumbers_cminv))),
        )

        partial_pressures = kwargs.get("vmr", np.zeros(len(pressure_pa))) * pressure_pa
        logging.debug(f"Starting Broadening for {self._molecule}")

        if self._line_coupling:
            Y = (
                self._line_db["Y"]
                .fillna(0.0)
                .interp(temperature=temperature_k, kwargs={"fill_value": 0.0})
            )
            G = (
                self._line_db["G"]
                .fillna(0.0)
                .interp(temperature=temperature_k, kwargs={"fill_value": 0.0})
            )

            Y_arr = np.zeros((len(self._line_db.nu), len(temperature_k)))
            G_arr = np.zeros((len(self._line_db.nu), len(temperature_k)))
            Y_arr[:, :] = Y.to_numpy().T
            G_arr[:, :] = G.to_numpy().T

            voigt_broaden_with_line_coupling(
                self._line_db.nu.to_numpy(),
                self._line_db.sw.to_numpy(),
                self._line_db.elower.to_numpy(),
                self._line_db.gamma_air.to_numpy(),
                self._line_db.gamma_self.to_numpy(),
                self._line_db.delta_air.to_numpy(),
                self._line_db.n_air.to_numpy(),
                self._line_db.local_iso_id.to_numpy().astype(np.int32),
                partition_ratio,
                Y_arr,
                G_arr,
                molecular_mass,
                pressure_pa / 101325.0,
                np.ones(len(pressure_pa)) * 0,
                temperature_k,
                wavenumbers_cminv[sidx],
                result,
                num_threads=num_threads,
                **self._kwargs,
            )
            result[:, sidx] = result

            # Line coupling numerics can give very small negative values
            result[result < 0] = 0.0

            return result / 1e4

        # else
        # Check if the input is uniform wavenumber grid
        uniform_grid = np.allclose(
            np.diff(wavenumbers_cminv), np.diff(wavenumbers_cminv)[0], atol=1e-6
        )

        if uniform_grid:
            voigt_broaden_uniform(
                self._line_db.nu.to_numpy(),
                self._line_db.sw.to_numpy(),
                self._line_db.elower.to_numpy(),
                self._line_db.gamma_air.to_numpy(),
                self._line_db.gamma_self.to_numpy(),
                self._line_db.delta_air.to_numpy(),
                self._line_db.n_air.to_numpy(),
                self._line_db.local_iso_id.to_numpy().astype(np.int32),
                partition_ratio,
                molecular_mass,
                pressure_pa / 101325.0,
                partial_pressures / 101325.0,
                temperature_k.astype(np.float64),
                wavenumbers_cminv[sidx][0],
                wavenumbers_cminv[sidx][1] - wavenumbers_cminv[sidx][0],
                result,
                num_threads=num_threads,
                **self._kwargs,
            )

            result[:, sidx] = result

            return result / 1e4

        voigt_broaden(
            self._line_db.nu.to_numpy(),
            self._line_db.sw.to_numpy(),
            self._line_db.elower.to_numpy(),
            self._line_db.gamma_air.to_numpy(),
            self._line_db.gamma_self.to_numpy(),
            self._line_db.delta_air.to_numpy(),
            self._line_db.n_air.to_numpy(),
            self._line_db.local_iso_id.to_numpy().astype(np.int32),
            partition_ratio,
            molecular_mass,
            pressure_pa / 101325.0,
            np.ones(len(pressure_pa)) * 0,
            temperature_k,
            wavenumbers_cminv[sidx],
            result,
            num_threads=num_threads,
            **self._kwargs,
        )
        result[:, sidx] = result

        return result / 1e4

    def cross_section_derivatives(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs  # noqa: ARG002
    ) -> dict:
        return {}


class HITRANAbsorber(LineAbsorber):
    def __init__(self, molecule: str, **kwargs):
        super().__init__(molecule, HITRANLineDatabase(), **kwargs)
