import logging
from copy import copy
from dataclasses import dataclass
from typing import List, Optional, Union

import numpy as np
import xarray as xr
from scipy import sparse

import sasktran2 as sk
from sasktran2.polarization import LegendreStorageView
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

    d_extinction: np.ndarray = None
    d_ssa: np.ndarray = None
    d_leg_coeff: np.ndarray = None
    scat_factor: np.ndarray = None
    scat_deriv_index: int = None
    d_albedo: np.ndarray = None


class DerivativeMapping:
    def __init__(
        self,
        native_grid_mapping: NativeGridDerivative,
        summable: bool = False,
        log_radiance_space: bool = False,
        name_prefix: str = "wf_",
    ):
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
        log_radiance_space : bool, optional
            True if the derivative should be returned in log_radiance space, i.e. if we should divide the absolute derivative by the radiance
        name_prefix: str, optional
            Prefix to put infront of the variable in the output dataset, default is "wf_".
        """
        self._native_grid_mapping = native_grid_mapping
        self._summable = summable
        self._log_radiance_space = log_radiance_space
        self._name_prefix = name_prefix

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
    def is_surface_derivative(self) -> bool:
        return self._native_grid_mapping.d_albedo is not None

    @property
    def summable(self) -> bool:
        """
        True if the mapping should be summed and accumulated.

        Returns
        -------
        bool
        """
        return self._summable

    @property
    def log_radiance_space(self) -> bool:
        """
        True if the derivative should be returned in log_radiance space, i.e. if we should divide the absolute derivative by the radiance

        Returns
        -------
        bool
        """
        return self._log_radiance_space

    @property
    def name_prefix(self) -> str:
        return self._name_prefix

    def map_derivative(self, data: np.ndarray, dimensions: List[str]):
        return xr.DataArray(
            data,
            dims=dimensions,
        )

    def name(self) -> Optional[str]:
        return None


class InterpolatedDerivativeMapping(DerivativeMapping):
    def __init__(
        self,
        native_grid_mapping: NativeGridDerivative,
        interpolating_matrix: np.ndarray,
        interp_dim="altitude",
        result_dim="interp_altitude",
        summable: bool = False,
        log_radiance_space: bool = False,
        name: Optional[str] = None,
    ):
        """
        A class which defines a mapping from internal model derivative quantities to user input quantities
        that are not on the native model grid

        Parameters
        ----------
        native_grid_mapping : NativeGridDerivative
            Mapping of the quantity in question from the native grid to the native grid
        interpolating_matrix : np.ndarray
            An interpolating matrix such that user quantity on the native grid can be calculated by multiplying
            the matrix by the user input quantity
        interp_dim : str, optional
            Dimension in the native mapping the interpolation is done over, by default "altitude"
        result_dim : str, optional
            string to name the resulting dimension, by default "interp_altitude"
        summable : bool, optional
            See :py:class:`DerivativeMapping`, by default False
        log_radiance_space : bool, optional
            See :py:class:`DerivativeMapping`, by default False
        name : Optional[str], optional
            A name for the derivative mapping, if set to None then the name is inferred from the key
            the mapping resides in.  This is only useful if you have two derivatives with the same key
            that should be summed together, by default None
        """
        super().__init__(native_grid_mapping, summable, log_radiance_space)
        self._name = name

        if interp_dim != result_dim:
            self._xr_interpolator = xr.DataArray(
                interpolating_matrix, dims=[interp_dim, result_dim]
            )
            self._rename_map = {}
        else:
            self._xr_interpolator = xr.DataArray(
                interpolating_matrix, dims=["tempDIM", result_dim]
            )
            self._rename_map = {interp_dim: "tempDIM"}

        self._sparse_mapping = sparse.csc_matrix(self._xr_interpolator.to_numpy())

    def map_derivative(self, data: np.ndarray, dimensions: List[str]):
        new_dims = copy(dimensions)
        new_dims[-1] = self._xr_interpolator.dims[-1]

        return xr.DataArray(
            (data.reshape(-1, data.shape[-1]) @ self._sparse_mapping).reshape(
                np.concatenate((data.shape[:-1], [-1]))
            ),
            dims=new_dims,
        ).transpose(*[new_dims[-1], *new_dims[:-1]])

    def name(self) -> Optional[str]:
        return self._name


class SurfaceDerivativeMapping(DerivativeMapping):
    def __init__(
        self,
        native_grid_mapping: NativeGridDerivative,
        interpolating_matrix: np.ndarray,
        interp_dim="wavelength",
        result_dim="interp_wavelength",
        summable: bool = False,
        log_radiance_space: bool = False,
    ):
        """
        A class which defines a mapping from internal model surface derivative quantities to user input quantities
        that are not on the native model grid

        Parameters
        ----------
        native_grid_mapping : NativeGridDerivative
            Mapping of the quantity in question from the native grid to the native grid
        interpolating_matrix : np.ndarray
            An interpolating matrix such that user quantity on the native grid can be calculated by multiplying
            the matrix by the user input quantity
        interp_dim : str, optional
            Dimension in the native mapping the interpolation is done over, by default "wavelength"
        result_dim : str, optional
            string to name the resulting dimension, by default "interp_wavelength"
        summable : bool, optional
            See :py:class:`DerivativeMapping`, by default False
        log_radiance_space : bool, optional
            See :py:class:`DerivativeMapping`, by default False
        """
        super().__init__(native_grid_mapping, summable, log_radiance_space)

        self._interp_dim = interp_dim
        self._result_dim = result_dim

        self._xr_interpolator = xr.DataArray(
            interpolating_matrix, dims=["tempDIM", result_dim]
        )

    def map_derivative(self, data: np.ndarray, dimensions: List[str]):
        return xr.DataArray(
            np.einsum(
                "ijk, il->lijk", data, self._xr_interpolator.to_numpy(), optimize=True
            ),
            dims=[self._result_dim, *dimensions],
        )


class Atmosphere:
    def __init__(
        self,
        model_geometry: sk.Geometry1D,
        config: sk.Config,
        wavelengths_nm: np.array = None,
        wavenumber_cminv: np.array = None,
        numwavel: Optional[int] = None,
        calculate_derivatives: bool = True,
        pressure_derivative: bool = True,
        temperature_derivative: bool = True,
    ):
        """
        The main specification for the atmospheric state.

        See :py:attr:`~storage` for details on how the atmospheric state parameters are stored.

        See :py:attr:`~surface` for details on how the surface parameters are stored.

        Parameters
        ----------
        model_geometry : sk.Geometry1D
            The geometry defining where the atmospheric quantities are specified on.
        config : sk.Config
            Main configuration object.
        wavelengths_nm : np.array, optional
            Array of wavelengths to specify the atmosphere at in [nm].  One of wavelengths_nm, wavenumber_cminv, numwavel must be set.
        wavenumber_cminv : np.array, optional
            Array of wavenumbers to specify the atmosphere at in [:math:`\\text{cm}^{-1}`]. One of wavelengths_nm, wavenumber_cminv, numwavel must be set.
        numwavel : int, optional
            Number of wavelengths to include in the calculation.  Note that wavelength is a dummy variable, this dimension can be
            anaything. One of wavelengths_nm, wavenumber_cminv, numwavel must be set.
        calculate_derivatives : bool, optional
            Whether or not the model should calculate derivatives with respect to atmospheric quantities., by default True
        pressure_derivative: bool, optional
            Whether or not the model should calculate derivatives with respect to pressure., by default True
        temperature_derivative: bool, optional
            Whether or not the model should calculate derivatives with respect to temperature., by default True
        """
        self._wavelengths_nm = None
        self._wavenumbers_cminv = None

        # Set wavelengths
        if (
            (wavelengths_nm is None)
            and (wavenumber_cminv is None)
            and (numwavel is None)
        ):
            msg = "One of wavelengths_nm, wavenumber_cminv, numwavel must be set when constructing the Atmosphere object"
            raise ValueError(msg)

        if wavelengths_nm is not None:
            self.wavelengths_nm = wavelengths_nm

        if wavenumber_cminv is not None:
            self.wavenumbers_cminv = wavenumber_cminv

        nwavel = len(self.wavelengths_nm) if numwavel is None else numwavel

        self._nstokes = config.num_stokes
        self._calculate_derivatives = calculate_derivatives
        self._pressure_derivative = pressure_derivative
        self._temperature_derivative = temperature_derivative

        if self._nstokes == 1:
            self._atmosphere = sk.AtmosphereStokes_1(
                nwavel, model_geometry, config, calculate_derivatives
            )
        elif self._nstokes == 3:
            self._atmosphere = sk.AtmosphereStokes_3(
                nwavel, model_geometry, config, calculate_derivatives
            )

        self._leg_coeff = LegendreStorageView(
            self._atmosphere.storage.leg_coeff, self._nstokes
        )

        self.surface.albedo = np.zeros(nwavel)

        self._storage_needs_reset = False

        self._pressure_pa = None
        self._temperature_k = None

        self._model_geometry = model_geometry
        self._config = config

        self._constituents = {}
        self._derivs = {}
        self._unscaled_ssa = None
        self._unscaled_extinction = None
        self._applied_delta_m_order = None

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
    def applied_delta_m_order(self) -> Optional[int]:
        """
        The order of the applied delta_m scaling. Can be None if no scaling has been applied
        """
        return self._applied_delta_m_order

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
    def calculate_temperature_derivative(self) -> bool:
        """
        True if we are calculating the derivative with respect to temperature

        Returns
        -------
        bool
        """
        return self._temperature_derivative

    @property
    def calculate_pressure_derivative(self) -> bool:
        """
        True if we are calculating the derivative with respect to pressure

        Returns
        -------
        bool
        """
        return self._pressure_derivative

    @property
    def storage(
        self,
    ) -> Union[sk.AtmosphereStorageStokes_1, sk.AtmosphereStorageStokes_3]:
        """
        The internal object which contains the atmosphere extinction, single scatter albedo, and legendre coefficients.

        Returns
        -------
        Union[sk.AtmosphereGridStorageStokes_1, sk.AtmosphereGridStorageStokes_3]
        """
        return self._atmosphere.storage

    @property
    def surface(self) -> sk.Surface:
        """
        The surface object

        Returns
        -------
        sk.Surface
        """
        return self._atmosphere.surface

    @property
    def temperature_k(self) -> np.array:
        """
        Atmospheric temperature in [K] on the same grid as :py:attr:`~model_geometry`

        Returns
        -------
        np.array
        """
        return self._temperature_k

    @temperature_k.setter
    def temperature_k(self, temp: np.array):
        self._temperature_k = temp

    @property
    def pressure_pa(self) -> np.array:
        """
        Pressure in [Pa] on the same grid as :py:attr:`~model_geometry`

        Returns
        -------
        np.array
        """
        return self._pressure_pa

    @pressure_pa.setter
    def pressure_pa(self, pres: np.array):
        self._pressure_pa = pres

    @property
    def wavelengths_nm(self) -> Optional[np.array]:
        """
        The wavelengths in [nm] the atmosphere is specified at.  This is an
        optional property, it may be None.

        Returns
        -------
        Optional[np.array]
        """
        return self._wavelengths_nm

    @wavelengths_nm.setter
    def wavelengths_nm(self, wav: np.array):
        self._wavelengths_nm = wav
        self._wavenumbers_cminv = wavlength_nm_to_wavenumber_cminv(wav)

    @property
    def wavenumbers_cminv(self) -> Optional[np.array]:
        """
        The wavenumbers in [:math:`\\text{cm}^{-1}`].  This is an optional property, it may be set to None

        Returns
        -------
        Optional[np.array]
        """
        return self._wavenumbers_cminv

    @wavenumbers_cminv.setter
    def wavenumbers_cminv(self, wav: np.array):
        self._wavenumbers_cminv = wav
        self._wavelengths_nm = wavenumber_cminv_to_wavlength_nm(wav)

    def _zero_storage(self):
        """
        Sets the internal storage object to 0.  Typically used in-between calculations to reset the
        atmosphere without reconstructing it.
        """
        self.storage.ssa[:] = 0
        self.storage.total_extinction[:] = 0
        self.storage.leg_coeff[:] = 0

    @property
    def deriv_mappings(self) -> dict:
        """
        A nested dictionary of :py:class:`sasktran2.atmosphere.DerivativeMapping` objects.

        Returns
        -------
        dict
        """
        return self._derivs

    @property
    def unscaled_ssa(self) -> np.ndarray:
        """
        The unscaled single scatter albedo.  After performing the calculation, `storage.ssa` may be
        modified because of the delta-m scaling.  This property stores the original unscaled single
        scatter albedo.

        Returns
        -------
        np.ndarray
        """
        return self._unscaled_ssa

    @property
    def unscaled_extinction(self) -> np.ndarray:
        """
        The unscaled extinction. After performing the calculation, `storage.total_extinction` may be
        modified because of the delta-m scaling. This property stores the original unscaled single
        scatter albedo.

        Returns
        -------
        np.ndarray
        """
        return self._unscaled_extinction

    @property
    def leg_coeff(self) -> LegendreStorageView:
        return self._leg_coeff

    def internal_object(self) -> Union[sk.AtmosphereStokes_1, sk.AtmosphereStokes_3]:
        """
        The internal `pybind11` object that can be used to perform the radiative transfer calculation.
        Calling this method will trigger a construction of the atmosphere object if necessary, and then
        return back the internal object.

        Returns
        -------
        Union[sk.AtmosphereStokes_1, sk.AtmosphereStokes_3]
        """
        num_scat_derivs = 0

        if len(self._constituents) > 0:
            logging.debug("Setting atmosphere from constituents")
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
                    self._derivs[name] = constituent.register_derivative(self, name)

                    if self._derivs[name] is not None:
                        for _, deriv in self._derivs[name].items():
                            if deriv.native_grid_mapping.d_leg_coeff is not None:
                                num_scat_derivs += 1
            logging.debug("Finished setting atmosphere from constituents")
        else:
            # using the raw interface
            if self._calculate_derivatives and len(self._derivs) == 0:
                self._derivs["raw"] = {}
                self._derivs["raw"]["extinction"] = DerivativeMapping(
                    NativeGridDerivative(
                        d_extinction=np.ones_like(self.storage.total_extinction),
                        d_ssa=np.zeros_like(self.storage.ssa),
                    ),
                    summable=True,
                )
                self._derivs["raw"]["ssa"] = DerivativeMapping(
                    NativeGridDerivative(
                        d_extinction=np.zeros_like(self.storage.total_extinction),
                        d_ssa=np.ones_like(self.storage.ssa),
                    ),
                    summable=True,
                )

                for i in range(self.storage.leg_coeff.shape[0]):
                    if i == 0:
                        # No derivative for the first leg_coeff
                        continue

                    d_leg_coeff = np.zeros_like(self.storage.leg_coeff)
                    d_leg_coeff[i, :] = 1
                    self._derivs["raw"][f"leg_coeff_{i}"] = DerivativeMapping(
                        NativeGridDerivative(
                            d_extinction=np.zeros_like(self.storage.total_extinction),
                            d_ssa=np.zeros_like(self.storage.ssa),
                            d_leg_coeff=d_leg_coeff,
                            scat_factor=np.ones_like(self.storage.ssa),
                        ),
                        summable=True,
                    )
                    num_scat_derivs += 1

        # Now we need to resize the phase derivative storage if necessary, and set the scattering derivatives
        if num_scat_derivs > 0:
            self.storage.resize_derivatives(num_scat_derivs)

            scat_index = 0
            for _, mappings in self._derivs.items():
                if mappings is not None:
                    for _, mapping in mappings.items():
                        if (
                            mapping is not None
                            and mapping.native_grid_mapping.d_leg_coeff is not None
                        ):
                            self.storage.d_leg_coeff[
                                :, :, :, scat_index
                            ] = mapping.native_grid_mapping.d_leg_coeff
                            mapping.native_grid_mapping.scat_deriv_index = scat_index
                            scat_index += 1

        # Store the unscaled optical properties for use in the derivative mappings
        self._unscaled_ssa = copy(self.storage.ssa)
        self._unscaled_extinction = copy(self.storage.total_extinction)

        # Apply delta scaling if required
        if self._config.delta_m_scaling:
            # Find the order of the scaling
            if self._config.num_streams == self._leg_coeff.a1.shape[0]:
                # Can't apply scaling
                logging.info(
                    "Delta-m Scaling NOT applied since the number of MS streams is equal to the number of user input legendre coefficients"
                )
            else:
                self._atmosphere.apply_delta_m_scaling(self._config.num_streams)
                self._applied_delta_m_order = self._config.num_streams

        self._storage_needs_reset = True
        return self._atmosphere

    def __setitem__(self, item, value):
        self._constituents[item] = value

    def __getitem__(self, item):
        return self._constituents.get(item)
