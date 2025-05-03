from __future__ import annotations

import logging
from copy import copy
from dataclasses import dataclass

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import EmissionSource, PyAtmosphere
from sasktran2.polarization import LegendreStorageView
from sasktran2.units import (
    wavenumber_cminv_to_wavlength_nm,
    wavlength_nm_to_wavenumber_cminv,
)
from sasktran2.util.state import EquationOfState


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
    d_brdf: np.ndarray = None


class Atmosphere:
    def __init__(
        self,
        model_geometry: sk.Geometry1D,
        config: sk.Config,
        wavelengths_nm: np.array = None,
        wavenumber_cminv: np.array = None,
        numwavel: int | None = None,
        calculate_derivatives: bool = True,
        pressure_derivative: bool = True,
        temperature_derivative: bool = True,
        specific_humidity_derivative: bool = True,
        legendre_derivative: bool = True,
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
        legendre_derivative: bool, optional
            Whether or not the model should calculate derivatives with respect to the legendre coefficients., by default True
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
            self.wavelengths_nm = wavelengths_nm.astype(np.float64)

        if wavenumber_cminv is not None:
            self.wavenumbers_cminv = wavenumber_cminv.astype(np.float64)

        nwavel = len(self.wavelengths_nm) if numwavel is None else numwavel

        self._nstokes = config.num_stokes
        self._calculate_derivatives = calculate_derivatives
        self._pressure_derivative = pressure_derivative
        self._temperature_derivative = temperature_derivative
        self._legendre_derivative = legendre_derivative
        self._specific_humidity_derivative = specific_humidity_derivative

        self._atmosphere = PyAtmosphere(
            nwavel,
            len(model_geometry.altitudes()),
            config.num_singlescatter_moments,
            calculate_derivatives,
            config.emission_source != EmissionSource.NoSource,
            config.num_stokes,
        )

        self._leg_coeff = LegendreStorageView(
            self._atmosphere.storage.leg_coeff, self._nstokes
        )

        self._storage_needs_reset = False

        self._equation_of_state = EquationOfState()

        self._model_geometry = model_geometry
        self._config = config

        self._constituents = {}
        self._derivs = {}
        self._unscaled_ssa = None
        self._unscaled_extinction = None
        self._applied_delta_m_order = None
        self._nwavel = nwavel

    @property
    def num_wavel(self) -> int:
        """
        The number of wavelengths the atmosphere is specified at

        Returns
        -------
        int
        """
        return self._nwavel

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
    def applied_delta_m_order(self) -> int | None:
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
    def calculate_specific_humidity_derivative(self) -> bool:
        """
        True if we are calculating the derivative with respect to specific humidity

        Returns
        -------
        bool
        """
        return self._specific_humidity_derivative

    @property
    def storage(
        self,
    ) -> sk.AtmosphereStorageStokes_1 | sk.AtmosphereStorageStokes_3:
        """
        The internal object which contains the atmosphere extinction, single scatter albedo, and legendre coefficients.

        Returns
        -------
        Union[sk.AtmosphereGridStorageStokes_1, sk.AtmosphereGridStorageStokes_3]
        """
        return self._atmosphere.storage

    @property
    def surface(self) -> sk.SurfaceStokes_1 | sk.SurfaceStokes_3:
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
        return self._equation_of_state.temperature_k

    @temperature_k.setter
    def temperature_k(self, temp: np.array):
        self._equation_of_state.temperature_k = temp

    @property
    def pressure_pa(self) -> np.array:
        """
        Pressure in [Pa] on the same grid as :py:attr:`~model_geometry`

        Returns
        -------
        np.array
        """
        return self._equation_of_state.pressure_pa

    @property
    def specific_humidity(self) -> np.array:
        """
        Specific humidity on the same grid as :py:attr:`~model_geometry`

        Returns
        -------
        np.array
        """
        return self._equation_of_state.specific_humidity

    @property
    def state_equation(self) -> EquationOfState:
        """
        The equation of state object which contains the temperature, pressure, and specific humidity

        Returns
        -------
        EquationOfState
        """
        return self._equation_of_state

    @specific_humidity.setter
    def specific_humidity(self, sh: np.array):
        self._equation_of_state.specific_humidity = sh

    @pressure_pa.setter
    def pressure_pa(self, pres: np.array):
        self._equation_of_state.pressure_pa = pres

    @property
    def wavelengths_nm(self) -> np.ndarray | None:
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
    def wavenumbers_cminv(self) -> np.ndarray | None:
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
        self.storage.set_zero()
        self.surface.set_zero()

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

    def internal_object(self) -> sk.AtmosphereStokes_1 | sk.AtmosphereStokes_3:
        """
        The internal `pybind11` object that can be used to perform the radiative transfer calculation.
        Calling this method will trigger a construction of the atmosphere object if necessary, and then
        return back the internal object.

        Returns
        -------
        Union[sk.AtmosphereStokes_1, sk.AtmosphereStokes_3]
        """

        if len(self._constituents) > 0:
            logging.debug("Setting atmosphere from constituents")
            # Using the constituent interface
            if self._storage_needs_reset:
                self._zero_storage()
            for _, constituent in self._constituents.items():
                constituent.add_to_atmosphere(self)

            self.storage.normalize_by_extinctions()

            self._derivs = {}
            if self._calculate_derivatives:
                for name, constituent in self._constituents.items():
                    constituent.register_derivative(self, name)

            logging.debug("Finished setting atmosphere from constituents")
        else:
            # using the raw interface
            if self._calculate_derivatives and len(self._derivs) == 0:
                deriv_mapping = self.storage.get_derivative_mapping("wf_extinction")
                deriv_mapping.d_extinction[:] = 1
                deriv_mapping.d_ssa[:] = 0.0
                deriv_mapping.interp_dim = "altitude"

                deriv_mapping = self.storage.get_derivative_mapping("wf_ssa")
                deriv_mapping.d_extinction[:] = 0.0
                deriv_mapping.d_ssa[:] = 1.0
                deriv_mapping.interp_dim = "altitude"

                deriv_mapping = self.surface.get_derivative_mapping("wf_albedo")
                deriv_mapping.d_brdf[:] = 1.0

                if self._legendre_derivative:
                    for i in range(self.storage.leg_coeff.shape[0]):
                        if i == 0:
                            # No derivative for the first leg_coeff
                            continue

                        deriv_mapping = self.storage.get_derivative_mapping(
                            f"wf_leg_coeff_{i}"
                        )
                        deriv_mapping.d_leg_coeff[i, :] = 1
                        deriv_mapping.d_extinction[:] = 0
                        deriv_mapping.d_ssa[:] = 0
                        deriv_mapping.scat_factor[:] = 1
                        deriv_mapping.interp_dim = "altitude"

        # Now we need to resize the phase derivative storage if necessary, and set the scattering derivatives
        self.storage.finalize_scattering_derivatives(0)

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
