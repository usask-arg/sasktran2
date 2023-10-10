from copy import copy
from dataclasses import dataclass
from typing import Union

import numpy as np

import sasktran2 as sk
from sasktran2.units import (
    wavenumber_cminv_to_wavlength_nm,
    wavlength_nm_to_wavenumber_cminv,
)


@dataclass
class NativeGridDerivative:
    """
    Internal input object that defines the model input quantities necessary to calculate a derivative.
    This mapping is from the native model grid to the native grid, it does not change the gridding.
    """
    d_extinction: np.ndarray
    d_ssa: np.ndarray
    d_leg_coeff: np.ndarray = None


class DerivativeMapping:
    def __init__(self, native_grid_mapping: NativeGridDerivative, summable: bool = False):
        """
        A class which defines a mapping from the internal model derivative quantities (derivatives with respect to,
        extinction, single scatter albedo, legendre coefficients) and user input quantities (e.g. VMR at a specific level).

        Parameters
        ----------
        native_grid_mapping : NativeGridDerivative
            A mapping from the native atmospheric grid to the same grid.
        summable : bool, optional
            True if this quantity should be accumulated and summed across each constituent. For example, each constituent may have a derivative with respect
            to a quantity like atmospheric temperature or pressure, by default False
        """
        self._native_grid_mapping = native_grid_mapping
        self._summable = summable

    @property
    def native_grid_mapping(self) -> NativeGridDerivative:
        """
        The mapping on the native grid.

        Returns
        -------
        NativeGridDerivative
        """
        return self._native_grid_mapping

    @property
    def summable(self) -> bool:
        """
        True if the mapping should be summed and accumulated.

        Returns
        -------
        bool
        """
        return self._summable

    def map_derivative(self):
        pass



class Atmosphere:
    def __init__(self, numwavel: int, model_geometry: sk.Geometry1D, config: sk.Config, calculate_derivatives: bool = True):
        """
        The main specification for the atmospheric state.

        Parameters
        ----------
        numwavel : int
            Number of wavelengths to include in the calculation.  Note that wavelength is a dummy variable, this dimension can be
            anaything.
        model_geometry : sk.Geometry1D
            The geometry defining where the atmospheric quantities are specified on.
        config : sk.Config
            Main configuration object.
        calculate_derivatives : bool, optional
            Whether or not the model should calculate derivatives with respect to atmospheric quantities., by default True
        """
        self._nstokes = config.num_stokes
        self._calculate_derivatives = calculate_derivatives

        if self._nstokes == 1:
            self._atmosphere = sk.AtmosphereStokes_1(numwavel, model_geometry, config, calculate_derivatives)
        elif self._nstokes == 3:
            self._atmosphere = sk.AtmosphereStokes_3(numwavel, model_geometry, config, calculate_derivatives)

        self._storage_needs_reset = False

        self._pressure_pa = None
        self._temperature_k = None

        self._wavelengths_nm = None
        self._wavenumbers_cminv = None

        self._model_geometry = model_geometry
        self._config = config

        self._constituents = {}
        self._derivs = {}
        self._unscaled_ssa = None
        self._unscaled_extinction = None

    @property
    def model_geometry(self) -> sk.Geometry1D:
        """
        The model geometry object

        Returns
        -------
        sk.Geometry1D
        """        
        return self._model_geometry

    @property
    def nstokes(self) -> int:
        """
        The number of stokes parameters the atmosphere contains information to calculate

        Returns
        -------
        int
        """
        return self._nstokes

    @property
    def storage(self) -> Union[sk.AtmosphereStorageStokes_1, sk.AtmosphereStorageStokes_3]:
        """
        The internal object which contains the atmosphere extinction, single scatter albedo, and legendre coefficients.

        Returns
        -------
        Union[sk.AtmosphereGridStorageStokes_1, sk.AtmosphereGridStorageStokes_3]
        """
        return self._atmosphere.storage

    @property
    def temperature_k(self) -> np.array:
        """
        Atmospheric temperature in [K] on the same grid as :py:property:model_geometry

        Returns
        -------
        np.array
        """
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
        self.storage.leg_coeff[:] = 0

    @property
    def deriv_mappings(self):
        return self._derivs

    @property
    def unscaled_ssa(self):
        return self._unscaled_ssa

    @property
    def unscaled_extinction(self):
        return self._unscaled_extinction

    def internal_object(self):
        if len(self._constituents) > 0:
            # Using the constituent interface
            if self._storage_needs_reset:
                self._zero_storage()
            for _, constituent in self._constituents.items():
                constituent.add_to_atmosphere(self)

            self.storage.leg_coeff /= self.storage.ssa[np.newaxis, :, :]

            self.storage.ssa /= self.storage.total_extinction

            self._derivs = {}
            if self._calculate_derivatives:
                for name, constituent in self._constituents.items():
                    self._derivs[name] = constituent.register_derivative(self)
        else:
            # using the raw interface
            if self._calculate_derivatives and len(self._derivs) == 0:
                self._derivs['raw'] = {}
                self._derivs['raw']['extinction'] = DerivativeMapping(NativeGridDerivative(d_extinction=np.ones_like(self.storage.total_extinction),
                                                                                            d_ssa=np.zeros_like(self.storage.ssa)), summable=True)
                self._derivs['raw']['ssa'] = DerivativeMapping(NativeGridDerivative(d_extinction=np.zeros_like(self.storage.total_extinction),
                                                                                    d_ssa=np.ones_like(self.storage.ssa)), summable=True)

                # TODO: Also calculate dlegendre for raw? probably...

        # Store the unscaled optical properties for use in the derivative mappings
        self._unscaled_ssa = copy(self.storage.ssa)
        self._unscaled_extinction = copy(self.storage.total_extinction)

        self._storage_needs_reset = True
        return self._atmosphere

    def __setitem__(self, item, value):
        self._constituents[item] = value

    def __getitem__(self, item):
        return self._constituents.get(item)
