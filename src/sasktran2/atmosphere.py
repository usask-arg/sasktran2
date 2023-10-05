import abc

import numpy as np

import sasktran2 as sk
from sasktran2.units import wavlength_nm_to_wavenumber_cminv, wavenumber_cminv_to_wavlength_nm


class Atmosphere:
    def __init__(self, numwavel: int, model_geometry: sk.Geometry1D, config: sk.Config, nstokes: int = 1):
        self._nstokes = nstokes

        if nstokes == 1:
            self._atmosphere = sk.AtmosphereStokes_1(numwavel, model_geometry, config, True)
        elif nstokes == 3:
            self._atmosphere = sk.AtmosphereStokes_3(numwavel, model_geometry, config, True)

        self._storage_needs_reset = False

        self._pressure_pa = None
        self._temperature_k = None

        self._wavelengths_nm = None
        self._wavenumbers_cminv = None

        self._model_geometry = model_geometry
        self._config = config

        self._constituents = dict()

    @property
    def model_geometry(self):
        return self._model_geometry

    @property
    def nstokes(self):
        return self._nstokes

    @property
    def storage(self):
        return self._atmosphere.storage

    @property
    def temperature_k(self):
        return self._temperature_k

    @temperature_k.setter
    def temperature_k(self, temp: np.array):
        self._temperature_k = temp

    @property
    def pressure_pa(self):
        return self._pressure_pa

    @pressure_pa.setter
    def pressure_pa(self, pres: np.array):
        self._pressure_pa = pres

    @property
    def wavelengths_nm(self):
        return self._wavelengths_nm

    @wavelengths_nm.setter
    def wavelengths_nm(self, wav: np.array):
        self._wavelengths_nm = wav
        self._wavenumbers_cminv = wavlength_nm_to_wavenumber_cminv(wav)

    @property
    def wavenumber_cminv(self):
        return self._wavenumber_cminv

    @wavenumber_cminv.setter
    def wavenumber_cminv(self, wav: np.array):
        self._wavenumber_cminv = wav
        self._wavelengths_nm = wavenumber_cminv_to_wavlength_nm(wav)

    def _zero_storage(self):
        self.storage.ssa[:] = 0
        self.storage.total_extinction[:] = 0
        for phase in self.storage.phase:
            phase.leg_coeff[:] = 0

    def internal_object(self):
        print('making atmo')
        if len(self._constituents) > 0:
            # Using the constituent interface
            if self._storage_needs_reset:
                self._zero_storage()
            for name, constituent in self._constituents.items():
                constituent.add_to_atmosphere(self)

            self.storage.leg_coeff /= self.storage.ssa[np.newaxis, :, :]

            self.storage.ssa /= self.storage.total_extinction

        self._storage_needs_reset = True
        print('atmo made')
        return self._atmosphere

    def __setitem__(self, item, value):
        self._constituents[item] = value

    def __getitem__(self, item):
        return self._constituents.get(item)

    def register_derivative(self):
        pass
