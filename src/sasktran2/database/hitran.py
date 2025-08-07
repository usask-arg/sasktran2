from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import xarray as xr

from sasktran2.optical.aerline import AERLineAbsorber
from sasktran2.optical.database import (
    OpticalDatabaseGenericAbsorber,
)
from sasktran2.optical.hitran import HITRANAbsorber
from sasktran2.units import wavenumber_cminv_to_wavlength_nm
from sasktran2.util import get_hapi

from .base import CachedDatabase
from .hitran_line import HITRANLineDatabase

PRESSURE_GRID = (
    np.array(
        [
            1.03181655e03,
            9.99910226e02,
            9.72231997e02,
            9.44960640e02,
            9.18098202e02,
            8.91646636e02,
            8.65607893e02,
            8.39983721e02,
            8.14775620e02,
            7.89985094e02,
            7.65613488e02,
            7.41661804e02,
            7.18131094e02,
            6.95022257e02,
            6.72335891e02,
            6.50072549e02,
            6.28232578e02,
            6.06816233e02,
            5.85823558e02,
            5.65254408e02,
            5.45108580e02,
            5.25385676e02,
            5.06085044e02,
            4.87205987e02,
            4.68747652e02,
            4.50708890e02,
            4.33088502e02,
            4.15885138e02,
            3.99097246e02,
            3.82723078e02,
            3.66760783e02,
            3.51208311e02,
            3.36063513e02,
            3.21323986e02,
            3.06987184e02,
            2.93050454e02,
            2.79510897e02,
            2.66365512e02,
            2.53611100e02,
            2.41244361e02,
            2.29261844e02,
            2.17659797e02,
            2.06434425e02,
            1.95581772e02,
            1.85097641e02,
            1.74977731e02,
            1.65217644e02,
            1.55812775e02,
            1.46758379e02,
            1.38049500e02,
            1.29681042e02,
            1.21647854e02,
            1.13944585e02,
            1.06565734e02,
            9.95056250e01,
            9.27584798e01,
            8.63183864e01,
            8.01792901e01,
            7.43350047e01,
            6.87792048e01,
            6.35054499e01,
            5.85071744e01,
            5.37776926e01,
            4.93102089e01,
            4.50978073e01,
            4.11334673e01,
            3.74100678e01,
            3.39203935e01,
            3.06571276e01,
            2.76128649e01,
            2.47801240e01,
            2.21513493e01,
            1.97189141e01,
            1.74751173e01,
            1.54122078e01,
            1.35223946e01,
            1.17978408e01,
            1.02306666e01,
            8.81297618e00,
            7.53686238e00,
            6.39440891e00,
            5.37770790e00,
            4.47887025e00,
            3.69003722e00,
            3.00339331e00,
            2.41118145e00,
            1.90571877e00,
            1.47941371e00,
            1.12478403e00,
            8.34477191e-01,
            6.01292929e-01,
            4.18208541e-01,
            2.78406802e-01,
            1.75308171e-01,
            1.02608084e-01,
            5.43204621e-02,
        ]
    )
    * 100
)
TEMP_GRID = np.arange(190, 311, 10)


class HITRANDatabase(CachedDatabase, OpticalDatabaseGenericAbsorber):
    def __init__(
        self,
        molecule: str,
        start_wavenumber: float,
        end_wavenumber: float,
        wavenumber_resolution: float,
        reduction_factor: int,
        num_threads: int = 1,
        db_root: Path | None = None,
        backend: str = "sasktran2",
        profile: str = "voigt",
        line_list: str = "hitran",
    ) -> None:
        """


        Parameters
        ----------
        molecule : str
            The HITRAN molecule identifier, e.g. "CO2", "H2O", "O3", "CH4"
        start_wavenumber : float
            The start wavenumber of the database in cm^-1
        end_wavenumber : float
            The end wavenumber of the database in cm^-1
        wavenumber_resolution : float
            The wavenumber resolution of the database in cm^-1
        num_threads: int, optional
            The number of threads to use in the cross section calculation. Default is 1.
        db_root : Path, optional
            The root directory to store the database, by default None
        backend : str, optional
            The backend to use, by default "sasktran2" which uses the internal broadening algorithms.
            options are ["sasktran2", "sasktran_legacy", "hapi"]. Default is "sasktran2".
            Note: "sasktran_legacy" requires the sasktran package to be installed.
            Both the "sasktran2" and "hapi" backends use the hapi package to access the HITRAN database.
        profile : str, optional
            The line shape profile to use in the cross section calculation. Currently, only the "hapi" backend supports options other
            than the default "voigt". Supported options are ["voigt", "sdvoigt", "ht", "priority", "lorentz", "doppler"].
            Default is "voigt".
        line_list : str, optional
            The line list database to use in the cross section calculation. Only the "sasktran2" backend supports options other
            than the default "hitran". Supported options are ["hitran", "aer"].
            Default is "hitran".
        """

        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        self._profile = profile.lower()
        if backend == "sasktran_legacy" and self._profile != "voigt":
            logging.warning(
                "%s line profile is not supported by sasktran_legacy backend. Switching to voigt",
                self._profile,
            )
            self._profile = "voigt"

        self._line_list = line_list.lower()
        if backend != "sasktran2" and line_list != "hitran":
            logging.warning(
                "%s line list is not supported by %s backend. Switching to hitran",
                line_list,
                backend,
            )
            self._line_list = "hitran"

        hasher = hashlib.sha1()
        encoded = json.dumps(
            {
                "start_wavenumber": start_wavenumber,
                "end_wavenumber": end_wavenumber,
                "wavenumber_resolution": wavenumber_resolution,
                "reduction_factor": reduction_factor,
                "profile": self._profile,
                "line_list": self._line_list,
            },
            sort_keys=True,
            cls=NumpyEncoder,
        ).encode()
        hasher.update(encoded)
        identifier = hasher.hexdigest()
        self._backend = backend
        self._start_wavenumber = start_wavenumber
        self._end_wavenumber = end_wavenumber
        self._wavenumber_resolution = wavenumber_resolution
        self._reduction_factor = reduction_factor
        self._molecule = molecule
        self._ltol = 1e-9
        self._wavenumber_wing = 50
        self._wavenumber_margin = 10
        self._num_threads = num_threads

        CachedDatabase.__init__(
            self,
            db_root=db_root,
            rel_path=Path("hitran").joinpath(self._molecule, self._backend),
        )

        self._data_file = self._db_root.joinpath(f"{identifier}.nc")

        OpticalDatabaseGenericAbsorber.__init__(self, self.path())

    def generate(self):
        if self._backend == "sasktran2":
            self._generate_sasktran2()
        elif self._backend == "sasktran_legacy":
            self._generate_sasktran_legacy()
        elif self._backend == "hapi":
            self._generate_hapi()
        else:
            msg = f"Invalid backend {self._backend}"
            raise ValueError(msg)

    def clear(self):
        if self._data_file.exists():
            self._data_file.unlink()

    def path(self, key: str | None = None, **kwargs) -> Path | None:  # noqa: ARG002
        if not self._data_file.exists():
            self.generate()
        return self._data_file

    def load_ds(self, key: str | None = None, **kwargs) -> xr.Dataset:
        return xr.open_dataset(self.path(key, **kwargs))

    def _generate_sasktran_legacy(self):
        try:
            from sasktran import ClimatologyUserDefined, HITRANChemical
        except ImportError as err:
            msg = "sasktran_legacy is required to generate Mie databases, try pip install sasktran"
            raise ImportError(msg) from err

        hires_wavenumber_grid = np.arange(
            self._start_wavenumber, self._end_wavenumber, self._wavenumber_resolution
        )

        xs = np.zeros((len(PRESSURE_GRID), len(TEMP_GRID), len(hires_wavenumber_grid)))

        for idx, pres in enumerate(PRESSURE_GRID):
            for idy, temp in enumerate(TEMP_GRID):
                xs[idx, idy] = (
                    HITRANChemical(self._molecule, line_tolerance=self._ltol)
                    .calculate_cross_sections(
                        ClimatologyUserDefined(
                            altitudes=[0, 100000],
                            values={
                                "SKCLIMATOLOGY_PRESSURE_PA": [pres, pres],
                                "SKCLIMATOLOGY_TEMPERATURE_K": [temp, temp],
                            },
                        ),
                        0,
                        0,
                        10000,
                        54372,
                        wavenumber_cminv_to_wavlength_nm(hires_wavenumber_grid),
                    )
                    .total
                    / 1e4
                )

        ds = xr.Dataset(
            {"xs": (["pressure", "temperature", "wavenumber_cminv"], xs)},
            coords={
                "pressure": PRESSURE_GRID,
                "temperature": TEMP_GRID,
                "wavenumber_cminv": hires_wavenumber_grid,
            },
        )

        ds.to_netcdf(self._data_file)

    def _generate_hapi(self):
        line_db = HITRANLineDatabase()
        line_db.path(self._molecule)

        hapi = get_hapi()

        abs_cal_fn = {
            "voigt": hapi.absorptionCoefficient_Voigt,
            "sdvoigt": hapi.absorptionCoefficient_SDVoigt,
            "ht": hapi.absorptionCoefficient_HT,
            "priority": hapi.absorptionCoefficient_Priority,
            "lorentz": hapi.absorptionCoefficient_Lorentz,
            "doppler": hapi.absorptionCoefficient_Doppler,
        }
        if self._profile not in abs_cal_fn:
            msg = f"Invalid line profile {self._profile}"
            raise ValueError(msg)

        line_db.initialize_hapi(self._molecule)

        hapi.select(
            self._molecule,
            DestinationTableName="spectral_window",
            Conditions=(
                "between",
                "nu",
                self._start_wavenumber - self._wavenumber_margin,
                self._end_wavenumber + self._wavenumber_margin,
            ),
        )

        hires_wavenumber_grid = np.arange(
            self._start_wavenumber, self._end_wavenumber, self._wavenumber_resolution
        )

        xs = np.zeros((len(PRESSURE_GRID), len(TEMP_GRID), len(hires_wavenumber_grid)))

        for idx, pres in enumerate(PRESSURE_GRID):
            for idy, temp in enumerate(TEMP_GRID):
                _, xs_hapi = abs_cal_fn[self._profile](
                    SourceTables="spectral_window",
                    Environment={"T": temp, "p": pres / 101325.0},
                    WavenumberGrid=hires_wavenumber_grid.tolist(),
                    WavenumberWing=self._wavenumber_wing,
                )
                xs[idx, idy] = xs_hapi / 1e4

        ds = xr.Dataset(
            {"xs": (["pressure", "temperature", "wavenumber_cminv"], xs)},
            coords={
                "pressure": PRESSURE_GRID,
                "temperature": TEMP_GRID,
                "wavenumber_cminv": hires_wavenumber_grid,
            },
        )

        ds.to_netcdf(self._data_file)

    def _generate_sasktran2(self):
        hires_wavenumber_grid = np.arange(
            self._start_wavenumber, self._end_wavenumber, self._wavenumber_resolution
        )

        if self._line_list == "hitran":
            calculator = HITRANAbsorber(
                self._molecule, line_contribution_width=self._wavenumber_wing / 2
            )
        elif self._line_list == "aer":
            calculator = AERLineAbsorber(
                self._molecule, line_contribution_width=self._wavenumber_wing / 2
            )
        else:
            msg = f"Invalid line_list {self._line_list}"
            raise ValueError(msg)

        xs = np.zeros((len(PRESSURE_GRID), len(TEMP_GRID), len(hires_wavenumber_grid)))

        for i, temperature in enumerate(TEMP_GRID):
            xs[:, i, :] = calculator.cross_sections(
                wavenumber_cminv_to_wavlength_nm(hires_wavenumber_grid),
                np.ones(len(PRESSURE_GRID)) * 0,
                temperature_k=np.ones(len(PRESSURE_GRID)) * temperature,
                pressure_pa=PRESSURE_GRID,
                num_threads=self._num_threads,
            )
        ds = xr.Dataset(
            {"xs": (["pressure", "temperature", "wavenumber_cminv"], xs)},
            coords={
                "pressure": PRESSURE_GRID,
                "temperature": TEMP_GRID,
                "wavenumber_cminv": hires_wavenumber_grid,
            },
        )

        ds.to_netcdf(self._data_file)
