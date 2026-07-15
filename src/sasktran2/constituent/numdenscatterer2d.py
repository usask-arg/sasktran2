from __future__ import annotations

from typing import Any

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import PyNumberDensityScatterer2D
from sasktran2.atmosphere import Atmosphere
from sasktran2.optical.base import OpticalProperty

from .base import Constituent


class NumberDensityScatterer2D(Constituent):
    """A scatterer specified by number density on a Geometry2D grid.

    The number density and all spatially varying optical-property arguments are
    stored in ``(horizontal, altitude)`` order. They are used directly on the
    native :class:`~sasktran2.Geometry2D` grid without interpolation or an
    expanded identity mapping.

    Parameters
    ----------
    optical_property : OpticalProperty
        Property supplying extinction, scattering, and phase-function data.
    number_density : numpy.ndarray
        Particle number density in ``m^-3`` with shape
        ``(horizontal, altitude)``.
    **kwargs
        Additional spatial inputs for ``optical_property``. Each input must be
        either a scalar or an array with the same shape as ``number_density``.
        Scalars are applied uniformly to all locations.
    """

    _constituent: PyNumberDensityScatterer2D

    def __init__(
        self,
        optical_property: OpticalProperty,
        number_density: np.ndarray,
        **kwargs,
    ) -> None:
        super().__init__()
        number_density = self._validate_native_profile(number_density, "number_density")

        self._volume_shape = number_density.shape
        self._optical_property = optical_property
        self._kwargs = {
            name: self._validate_aux_input(value, name)
            for name, value in kwargs.items()
        }
        self._vertical_deriv_factor = np.ones(self._volume_shape, dtype=np.float64)
        self._d_vertical_deriv_factor: dict[str, np.ndarray] = {}
        self._wf_name = "number_density"
        self._constituent = PyNumberDensityScatterer2D(
            optical_property,
            np.ascontiguousarray(number_density).reshape(-1),
            **self._flat_kwargs(),
        )
        self._sync_constituent_state()

    @staticmethod
    def _validate_native_profile(value: np.ndarray, name: str) -> np.ndarray:
        value = np.asarray(value, dtype=np.float64)
        if value.ndim != 2:
            msg = f"{name} must have shape (horizontal, altitude); got {value.shape}"
            raise ValueError(msg)
        if 0 in value.shape:
            msg = f"{name} horizontal and altitude dimensions must both be non-empty"
            raise ValueError(msg)
        return np.ascontiguousarray(value)

    def _validate_aux_input(self, value: np.ndarray, name: str) -> np.ndarray:
        value = np.asarray(value, dtype=np.float64)
        if value.ndim == 0:
            return np.full(self._volume_shape, value.item(), dtype=np.float64)
        if value.shape != self._volume_shape:
            msg = (
                f"{name} must be scalar or have shape {self._volume_shape}; "
                f"got {value.shape}"
            )
            raise ValueError(msg)
        return np.ascontiguousarray(value)

    def _flat_kwargs(self) -> dict[str, np.ndarray]:
        return {
            name: np.ascontiguousarray(value).reshape(-1)
            for name, value in self._kwargs.items()
        }

    def _sync_constituent_state(self) -> None:
        self._constituent.aux_inputs = self._flat_kwargs()
        self._constituent.vertical_deriv_factor = np.ascontiguousarray(
            self._vertical_deriv_factor
        ).reshape(-1)
        self._constituent.d_vertical_deriv_factor = {
            name: np.ascontiguousarray(value).reshape(-1)
            for name, value in self._d_vertical_deriv_factor.items()
        }
        self._constituent.wf_name = self._wf_name

    def __getattr__(self, name: str) -> Any:
        if name in self.__dict__.get("_kwargs", {}):
            return self._kwargs[name]
        msg = f"{type(self).__name__!s} has no attribute {name!r}"
        raise AttributeError(msg)

    def __setattr__(self, name: str, value: Any) -> None:
        if name in self.__dict__.get("_kwargs", {}):
            self._kwargs[name] = self._validate_aux_input(value, name)
            self._sync_constituent_state()
        else:
            super().__setattr__(name, value)

    @property
    def volume_spatial_mode(self) -> str:
        return "native_2d"

    @property
    def number_density(self) -> np.ndarray:
        """Particle number density in ``(horizontal, altitude)`` order."""
        return self._constituent.number_density.reshape(self._volume_shape)

    @number_density.setter
    def number_density(self, number_density: np.ndarray) -> None:
        number_density = np.asarray(number_density, dtype=np.float64)
        if number_density.shape != self._volume_shape:
            msg = (
                f"number_density must retain shape {self._volume_shape}; "
                f"got {number_density.shape}"
            )
            raise ValueError(msg)
        self._constituent.number_density = np.ascontiguousarray(number_density).reshape(
            -1
        )

    def _validate_atmosphere(self, atmo: Atmosphere) -> None:
        if not isinstance(atmo.model_geometry, sk.Geometry2D):
            msg = f"{type(self).__name__} requires an atmosphere using Geometry2D"
            raise TypeError(msg)
        if tuple(atmo.volume_shape) != self._volume_shape:
            msg = (
                f"{type(self).__name__} shape does not match the atmosphere: "
                f"{self._volume_shape} != {atmo.volume_shape}"
            )
            raise ValueError(msg)

    def add_to_atmosphere(self, atmo: Atmosphere) -> None:
        self._validate_atmosphere(atmo)
        self._sync_constituent_state()
        self._constituent.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: Atmosphere, name: str) -> None:
        self._validate_atmosphere(atmo)
        self._sync_constituent_state()
        self._constituent.register_derivative(atmo, name)


class ExtinctionScatterer2D(NumberDensityScatterer2D):
    """A Geometry2D scatterer normalized to an extinction profile.

    This is the native-grid counterpart to
    :class:`~sasktran2.constituent.ExtinctionScatterer`. The supplied
    extinction is converted to number density using the optical property's
    cross section at ``extinction_wavelength_nm``. The conversion is evaluated
    independently at every Geometry2D location, including the supplied optical
    property arguments.

    Parameters
    ----------
    optical_property : OpticalProperty
        Property supplying extinction, scattering, and phase-function data.
    extinction_per_m : numpy.ndarray
        Extinction in ``m^-1`` with shape ``(horizontal, altitude)``.
    extinction_wavelength_nm : float
        Wavelength in nanometres at which ``extinction_per_m`` is specified.
    **kwargs
        Additional optical-property inputs. Shape rules are identical to
        :class:`NumberDensityScatterer2D`.
    """

    def __init__(
        self,
        optical_property: OpticalProperty,
        extinction_per_m: np.ndarray,
        extinction_wavelength_nm: float,
        **kwargs,
    ) -> None:
        extinction_per_m = self._validate_native_profile(
            extinction_per_m, "extinction_per_m"
        )
        extinction_wavelength_nm = float(extinction_wavelength_nm)
        if not np.isfinite(extinction_wavelength_nm) or extinction_wavelength_nm <= 0:
            msg = "extinction_wavelength_nm must be finite and positive"
            raise ValueError(msg)

        self._extinction_per_m = extinction_per_m.copy()
        self._extinction_wavelength_nm = extinction_wavelength_nm
        self._extinction_to_numden_factors: np.ndarray | None = None
        super().__init__(
            optical_property,
            np.zeros_like(extinction_per_m, dtype=np.float64),
            **kwargs,
        )
        self._wf_name = "extinction"
        self._sync_constituent_state()

    def _update_number_density(self, atmo: Atmosphere) -> None:
        native_altitudes = np.asarray(atmo._native_altitudes(), dtype=np.float64)
        factors = np.asarray(
            self._optical_property.cross_sections(
                np.array([self._extinction_wavelength_nm]),
                altitudes_m=native_altitudes,
                **self._flat_kwargs(),
            ).extinction,
            dtype=np.float64,
        ).reshape(-1)
        expected_size = int(np.prod(self._volume_shape))
        if factors.size != expected_size:
            msg = (
                "Optical-property reference cross section has spatial size "
                f"{factors.size}; expected {expected_size}"
            )
            raise ValueError(msg)
        if np.any(~np.isfinite(factors)) or np.any(factors <= 0):
            msg = "Reference extinction cross sections must be finite and positive"
            raise ValueError(msg)

        self._extinction_to_numden_factors = factors.reshape(self._volume_shape)
        self._vertical_deriv_factor = 1.0 / self._extinction_to_numden_factors
        self.number_density = (
            self._extinction_per_m / self._extinction_to_numden_factors
        )

        derivatives = self._optical_property.cross_section_derivatives(
            np.array([self._extinction_wavelength_nm]),
            altitudes_m=native_altitudes,
            **self._flat_kwargs(),
        )
        self._d_vertical_deriv_factor = {}
        for name, derivative in derivatives.items():
            derivative_array = np.asarray(derivative, dtype=np.float64)
            if derivative_array.size != expected_size:
                msg = (
                    f"Cross-section derivative {name!r} has spatial size "
                    f"{derivative_array.size}; expected {expected_size}"
                )
                raise ValueError(msg)
            derivative_array = derivative_array.reshape(self._volume_shape)
            self._d_vertical_deriv_factor[name] = (
                -derivative_array * self._vertical_deriv_factor**2
            )

    @property
    def extinction_per_m(self) -> np.ndarray:
        """Reference extinction in ``(horizontal, altitude)`` order."""
        return self._extinction_per_m

    @extinction_per_m.setter
    def extinction_per_m(self, extinction_per_m: np.ndarray) -> None:
        extinction_per_m = np.asarray(extinction_per_m, dtype=np.float64)
        if extinction_per_m.shape != self._volume_shape:
            msg = (
                f"extinction_per_m must retain shape {self._volume_shape}; "
                f"got {extinction_per_m.shape}"
            )
            raise ValueError(msg)
        self._extinction_per_m = np.ascontiguousarray(extinction_per_m)

    @property
    def extinction_wavelength_nm(self) -> float:
        """Wavelength in nanometres used to normalize the extinction."""
        return self._extinction_wavelength_nm

    def add_to_atmosphere(self, atmo: Atmosphere) -> None:
        self._validate_atmosphere(atmo)
        self._update_number_density(atmo)
        super().add_to_atmosphere(atmo)
