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


@dataclass(frozen=True)
class _AtmosphereSpatialLayout:
    """Internal description of the volume-atmosphere storage layout."""

    num_horizontal: int
    num_altitudes: int
    geometry_is_2d: bool

    @classmethod
    def from_geometry(
        cls, geometry: sk.Geometry1D | sk.Geometry2D
    ) -> _AtmosphereSpatialLayout:
        if isinstance(geometry, sk.Geometry2D):
            num_horizontal, num_altitudes = geometry.shape
            return cls(num_horizontal, num_altitudes, True)

        return cls(1, len(geometry.altitudes()), False)

    @property
    def is_2d(self) -> bool:
        return self.geometry_is_2d

    @property
    def shape(self) -> tuple[int, ...]:
        if self.is_2d:
            return (self.num_horizontal, self.num_altitudes)
        return (self.num_altitudes,)

    @property
    def num_locations(self) -> int:
        return self.num_horizontal * self.num_altitudes

    def native_altitudes(self, geometry: sk.Geometry1D | sk.Geometry2D) -> np.ndarray:
        altitudes = np.asarray(geometry.altitudes(), dtype=np.float64)
        if self.is_2d:
            return np.tile(altitudes, self.num_horizontal)
        return altitudes

    def validate_state(self, value: np.ndarray, name: str) -> np.ndarray:
        array = np.asarray(value)
        valid_shapes = {(self.num_altitudes,)}
        if self.is_2d:
            valid_shapes.add((self.num_horizontal, self.num_altitudes))

        if array.shape not in valid_shapes:
            expected = " or ".join(str(shape) for shape in sorted(valid_shapes))
            msg = f"{name} must have shape {expected}; got {array.shape}"
            raise ValueError(msg)
        return array

    def native_state(self, value: np.ndarray | None, name: str) -> np.ndarray | None:
        if value is None:
            return None

        array = self.validate_state(value, name)
        if self.is_2d and array.ndim == 1:
            array = np.broadcast_to(array, (self.num_horizontal, self.num_altitudes))
        return np.ascontiguousarray(array.reshape(-1), dtype=np.float64)

    def reshape_native(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values)
        if values.shape[0] != self.num_locations:
            msg = (
                "The first dimension must be the flattened atmosphere location "
                f"dimension ({self.num_locations}); got {values.shape[0]}"
            )
            raise ValueError(msg)
        return values.reshape((*self.shape, *values.shape[1:]))


@dataclass(frozen=True)
class _NativeAtmosphereInputs:
    altitudes_m: np.ndarray
    state_equation: EquationOfState

    def state(self, name: str) -> np.ndarray | None:
        return getattr(self.state_equation, name)


class _ExpandedProfileGeometry:
    """Geometry facade used by altitude-only constituents in a 2D atmosphere."""

    def __init__(
        self,
        geometry: sk.Geometry2D,
        native_altitudes_m: np.ndarray,
    ) -> None:
        self._geometry = geometry
        self._native_altitudes = native_altitudes_m

    def altitudes(self) -> np.ndarray:
        return self._native_altitudes

    def __getattr__(self, name):
        return getattr(self._geometry, name)


class _ExpandedProfileAtmosphere:
    """Atmosphere facade that broadcasts altitude-only inputs onto 2D locations."""

    def __init__(
        self, atmosphere: Atmosphere, native_inputs: _NativeAtmosphereInputs
    ) -> None:
        self._atmosphere = atmosphere
        self._model_geometry = _ExpandedProfileGeometry(
            atmosphere.model_geometry, native_inputs.altitudes_m
        )
        self._native_inputs = native_inputs

    @property
    def model_geometry(self):
        return self._model_geometry

    @property
    def state_equation(self) -> EquationOfState:
        return self._native_inputs.state_equation

    @property
    def pressure_pa(self) -> np.ndarray | None:
        return self._native_inputs.state("pressure_pa")

    @property
    def temperature_k(self) -> np.ndarray | None:
        return self._native_inputs.state("temperature_k")

    @property
    def specific_humidity(self) -> np.ndarray | None:
        return self._native_inputs.state("specific_humidity")

    def _native_state(self, name: str) -> np.ndarray | None:
        return self._native_inputs.state(name)

    def _native_altitudes(self) -> np.ndarray:
        return self._native_inputs.altitudes_m

    def _native_state_equation(self) -> EquationOfState:
        return self._native_inputs.state_equation

    def __getattr__(self, name):
        return getattr(self._atmosphere, name)


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
        model_geometry: sk.Geometry1D | sk.Geometry2D,
        config: sk.Config,
        wavelengths_nm: np.array = None,
        wavenumber_cminv: np.array = None,
        numwavel: int | None = None,
        calculate_derivatives: bool = True,
        pressure_derivative: bool = True,
        temperature_derivative: bool = True,
        specific_humidity_derivative: bool = True,
        legendre_derivative: bool = True,
        spectral_grid: sk.basis.Grid | None = None,  # noqa: ARG002
    ):
        """
        The main specification for the atmospheric state.

        See :py:attr:`~storage` for details on how the atmospheric state parameters are stored.

        See :py:attr:`~surface` for details on how the surface parameters are stored.

        Parameters
        ----------
        model_geometry : sk.Geometry1D | sk.Geometry2D
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

        grid = None
        wavenumber_space = None

        if wavelengths_nm is not None:
            self.wavelengths_nm = wavelengths_nm.astype(np.float64)
            grid = sk.basis.Grid.from_triangles(self.wavelengths_nm)
            wavenumber_space = False

        if wavenumber_cminv is not None:
            self.wavenumbers_cminv = wavenumber_cminv.astype(np.float64)
            grid = sk.basis.Grid.from_triangles(self.wavenumbers_cminv)
            wavenumber_space = True

        nwavel = len(self.wavelengths_nm) if numwavel is None else numwavel

        self._nstokes = config.num_stokes
        self._spectral_integration_mode = config.spectral_grid_mode
        self._calculate_derivatives = calculate_derivatives
        self._pressure_derivative = pressure_derivative
        self._temperature_derivative = temperature_derivative
        self._legendre_derivative = legendre_derivative
        self._specific_humidity_derivative = specific_humidity_derivative

        self._spatial_layout = _AtmosphereSpatialLayout.from_geometry(model_geometry)

        self._atmosphere = PyAtmosphere(
            nwavel,
            self._spatial_layout.num_locations,
            config.num_singlescatter_moments,
            calculate_derivatives,
            config.emission_source != EmissionSource.NoSource,
            config.num_stokes,
            grid._internal_object() if grid is not None else None,
            grid._internal_object() if grid is not None else None,
            wavenumber_space,
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
        self._derivative_output_shapes = {}
        self._spatial_state_derivatives = set()
        self._native_input_cache = None

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
    def model_geometry(self) -> sk.Geometry1D | sk.Geometry2D:
        """
        The model geometry object

        Returns
        -------
        sk.Geometry1D | sk.Geometry2D
        """
        return self._model_geometry

    @property
    def volume_shape(self) -> tuple[int, ...]:
        """Shape of the volume atmosphere, excluding wavelength."""
        return self._spatial_layout.shape

    @property
    def num_locations(self) -> int:
        """Number of flattened volume-atmosphere locations."""
        return self._spatial_layout.num_locations

    def reshape_native(self, values: np.ndarray) -> np.ndarray:
        """Restore flattened native-location data to the geometry shape."""
        return self._spatial_layout.reshape_native(values)

    def derivative_output_shape(self, mapping_name: str) -> tuple[int, ...] | None:
        """Structured parameter shape for a native-location derivative mapping."""
        return self._derivative_output_shapes.get(mapping_name)

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
    def calculate_derivatives(self) -> bool:
        """Whether the atmosphere registers and calculates derivatives."""
        return self._calculate_derivatives

    @property
    def revision(self) -> int:
        """Monotonic revision of the built native atmosphere state."""
        return self._atmosphere.revision

    def mark_changed(self) -> None:
        """Indicate that directly mutable native atmosphere storage changed.

        Constituent-backed atmospheres call this automatically when their
        native representation is rebuilt. Users of the raw storage interface
        must call it after modifying a borrowed NumPy storage view.
        """
        self._atmosphere.mark_changed()

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
        self._equation_of_state.temperature_k = self._validate_state(
            temp, "temperature_k"
        )

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
        self._equation_of_state.specific_humidity = self._validate_state(
            sh, "specific_humidity"
        )

    @pressure_pa.setter
    def pressure_pa(self, pres: np.array):
        self._equation_of_state.pressure_pa = self._validate_state(pres, "pressure_pa")

    def _validate_state(self, value: np.ndarray | None, name: str) -> np.ndarray | None:
        if value is None:
            return None
        if isinstance(self._model_geometry, sk.Geometry2D):
            return self._spatial_layout.validate_state(value, name)
        return value

    def _native_state(self, name: str) -> np.ndarray | None:
        if self._native_input_cache is not None:
            return self._native_input_cache.state(name)
        if not self._spatial_layout.is_2d:
            return getattr(self._equation_of_state, name)
        return self._spatial_layout.native_state(
            getattr(self._equation_of_state, name), name
        )

    def _native_altitudes(self) -> np.ndarray:
        if self._native_input_cache is not None:
            return self._native_input_cache.altitudes_m
        if not self._spatial_layout.is_2d:
            return self._model_geometry.altitudes()
        return self._spatial_layout.native_altitudes(self._model_geometry)

    def _native_state_equation(self) -> EquationOfState:
        if self._native_input_cache is not None:
            return self._native_input_cache.state_equation
        if not self._spatial_layout.is_2d:
            return self._equation_of_state
        state = EquationOfState(self._equation_of_state.molar_mass_dry_air)
        state.pressure_pa = self._native_state("pressure_pa")
        state.temperature_k = self._native_state("temperature_k")
        state.specific_humidity = self._native_state("specific_humidity")
        return state

    def _make_native_input_cache(self) -> _NativeAtmosphereInputs:
        state = EquationOfState(self._equation_of_state.molar_mass_dry_air)
        for name in ("pressure_pa", "temperature_k", "specific_humidity"):
            value = self._spatial_layout.native_state(
                getattr(self._equation_of_state, name), name
            )
            setattr(state, name, value)
        return _NativeAtmosphereInputs(
            self._spatial_layout.native_altitudes(self._model_geometry), state
        )

    def _constituent_atmosphere(self, constituent, profile_atmosphere):
        spatial_mode = getattr(constituent, "volume_spatial_mode", "altitude_profile")
        if spatial_mode not in ("altitude_profile", "native_2d"):
            msg = (
                f"Unsupported volume_spatial_mode {spatial_mode!r}; expected "
                "'altitude_profile' or 'native_2d'"
            )
            raise ValueError(msg)

        if not isinstance(self._model_geometry, sk.Geometry2D):
            return self
        if spatial_mode == "native_2d":
            return self
        return profile_atmosphere

    def _finalize_spatial_derivatives(self) -> None:
        self._derivative_output_shapes = {}
        if not self._spatial_layout.is_2d:
            return

        state_assign_names = {
            "pressure_pa": "wf_pressure_pa",
            "temperature_k": "wf_temperature_k",
            "specific_humidity": "wf_specific_humidity",
        }
        for mapping_name in self.storage.derivative_mapping_names():
            mapping = self.storage.get_derivative_mapping(mapping_name)
            interpolator = np.asarray(mapping.interpolator)
            if mapping.interp_dim == "location" and (
                interpolator.size == 0
                or (
                    interpolator.ndim == 2
                    and interpolator.shape[1] == self._spatial_layout.num_locations
                )
            ):
                self._derivative_output_shapes[mapping_name] = self.volume_shape

            state_name = next(
                (
                    name
                    for name, assign_name in state_assign_names.items()
                    if mapping.assign_name == assign_name
                ),
                None,
            )
            if state_name is None:
                continue

            self._spatial_state_derivatives.add(mapping_name)

            state_value = getattr(self._equation_of_state, state_name)
            if state_value is None and state_name != "specific_humidity":
                continue

            if interpolator.size != 0 and (
                interpolator.ndim != 2
                or interpolator.shape[1] != self._spatial_layout.num_locations
            ):
                continue

            if state_value is None or np.asarray(state_value).ndim == 1:
                if interpolator.size == 0:
                    broadcast = np.zeros(
                        (
                            self._spatial_layout.num_locations,
                            self._spatial_layout.num_altitudes,
                        )
                    )
                    location = np.arange(self._spatial_layout.num_locations)
                    broadcast[
                        location, location % self._spatial_layout.num_altitudes
                    ] = 1
                    mapping.interpolator = broadcast
                else:
                    mapping.interpolator = interpolator.reshape(
                        interpolator.shape[0],
                        self._spatial_layout.num_horizontal,
                        self._spatial_layout.num_altitudes,
                    ).sum(axis=1)
                mapping.interp_dim = "altitude"
            else:
                if interpolator.size != 0:
                    mapping.interpolator = interpolator
                mapping.interp_dim = "location"
                self._derivative_output_shapes[mapping_name] = self.volume_shape

    def _prepare_spatial_derivatives(self) -> None:
        if not self._spatial_layout.is_2d:
            return
        for mapping_name in self._spatial_state_derivatives:
            mapping = self.storage.get_derivative_mapping(mapping_name)
            mapping.clear_interpolator()
            mapping.interp_dim = "location"

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
    def spectral_integration_mode(self) -> sk.SpectralGridMode:
        """
        The spectral integration mode for the atmosphere

        Returns
        -------
        sk.SpectralGridMode
        """
        return self._spectral_integration_mode

    @property
    def leg_coeff(self) -> LegendreStorageView:
        return self._leg_coeff

    def _into_rust_object(self) -> PyAtmosphere:
        return self._atmosphere

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
            # Rebuilding mutates shared native storage. Invalidate existing
            # linearization sessions before any mutation, including rebuilds
            # that subsequently fail.
            self.mark_changed()
            # Using the constituent interface
            if self._storage_needs_reset:
                self._zero_storage()

            # Mark the storage dirty before constituent work begins. If a
            # constituent fails, the next call must rebuild from zero rather
            # than accumulating on a partially populated atmosphere.
            self._storage_needs_reset = True

            profile_atmosphere = self
            try:
                if isinstance(self._model_geometry, sk.Geometry2D):
                    native_inputs = self._make_native_input_cache()
                    self._native_input_cache = native_inputs
                    profile_atmosphere = _ExpandedProfileAtmosphere(self, native_inputs)

                for _, constituent in self._constituents.items():
                    constituent.add_to_atmosphere(
                        self._constituent_atmosphere(constituent, profile_atmosphere)
                    )

                self.storage.normalize_by_extinctions()

                self._derivs = {}
                if self._calculate_derivatives:
                    self._prepare_spatial_derivatives()
                    for name, constituent in self._constituents.items():
                        constituent.register_derivative(
                            self._constituent_atmosphere(
                                constituent, profile_atmosphere
                            ),
                            name,
                        )
                    self._finalize_spatial_derivatives()
            except Exception:
                self._zero_storage()
                self._derivs = {}
                self._derivative_output_shapes = {}
                raise
            finally:
                self._native_input_cache = None

            logging.debug("Finished setting atmosphere from constituents")
        else:
            # using the raw interface
            if self._calculate_derivatives and len(self._derivs) == 0:
                self._derivative_output_shapes = {}
                deriv_mapping = self.storage.get_derivative_mapping("wf_extinction")
                deriv_mapping.d_extinction[:] = 1
                deriv_mapping.d_ssa[:] = 0.0
                deriv_mapping.interp_dim = (
                    "location" if self._spatial_layout.is_2d else "altitude"
                )
                if self._spatial_layout.is_2d:
                    self._derivative_output_shapes["wf_extinction"] = self.volume_shape

                deriv_mapping = self.storage.get_derivative_mapping("wf_ssa")
                deriv_mapping.d_extinction[:] = 0.0
                deriv_mapping.d_ssa[:] = 1.0
                deriv_mapping.interp_dim = (
                    "location" if self._spatial_layout.is_2d else "altitude"
                )
                if self._spatial_layout.is_2d:
                    self._derivative_output_shapes["wf_ssa"] = self.volume_shape

                deriv_mapping = self.surface.get_derivative_mapping("wf_albedo")
                deriv_mapping.d_brdf[:] = 1.0

                if self._config.emission_source != sk.EmissionSource.NoSource:
                    deriv_mapping = self.storage.get_derivative_mapping("wf_emission")
                    deriv_mapping.d_extinction[:] = 0.0
                    deriv_mapping.d_ssa[:] = 0.0
                    deriv_mapping.d_emission[:] = 1.0
                    deriv_mapping.interp_dim = (
                        "location" if self._spatial_layout.is_2d else "altitude"
                    )
                    if self._spatial_layout.is_2d:
                        self._derivative_output_shapes["wf_emission"] = (
                            self.volume_shape
                        )

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
                        deriv_mapping.interp_dim = (
                            "location" if self._spatial_layout.is_2d else "altitude"
                        )
                        if self._spatial_layout.is_2d:
                            self._derivative_output_shapes[f"wf_leg_coeff_{i}"] = (
                                self.volume_shape
                            )

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
