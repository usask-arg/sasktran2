from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import PyVMRAbsorber2D
from sasktran2.atmosphere import Atmosphere
from sasktran2.optical.base import OpticalProperty

from .base import Constituent


class VMRAbsorber2D(Constituent):
    """An absorber specified by VMR on the native Geometry2D grid.

    The VMR is stored in ``(horizontal, altitude)`` order and must match the
    :class:`~sasktran2.Geometry2D` used by the atmosphere. No spatial
    interpolation or horizontal broadcasting is performed.

    Parameters
    ----------
    optical_property : OpticalProperty
        Property that supplies absorption cross sections.
    vmr : numpy.ndarray
        Volume mixing ratio with shape ``(horizontal, altitude)``.
    """

    _constituent: PyVMRAbsorber2D

    def __init__(self, optical_property: OpticalProperty, vmr: np.ndarray) -> None:
        super().__init__()
        vmr = np.asarray(vmr, dtype=np.float64)
        if vmr.ndim != 2:
            msg = f"vmr must have shape (horizontal, altitude); got {vmr.shape}"
            raise ValueError(msg)
        if 0 in vmr.shape:
            msg = "vmr horizontal and altitude dimensions must both be non-empty"
            raise ValueError(msg)

        self._volume_shape = vmr.shape
        self._constituent = PyVMRAbsorber2D(
            optical_property, np.ascontiguousarray(vmr).reshape(-1)
        )
        self._optical_property = optical_property

    @property
    def volume_spatial_mode(self) -> str:
        return "native_2d"

    @property
    def vmr(self) -> np.ndarray:
        """Volume mixing ratio in ``(horizontal, altitude)`` order."""
        return self._constituent.vmr.reshape(self._volume_shape)

    @vmr.setter
    def vmr(self, vmr: np.ndarray) -> None:
        vmr = np.asarray(vmr, dtype=np.float64)
        if vmr.shape != self._volume_shape:
            msg = f"vmr must retain shape {self._volume_shape}; got {vmr.shape}"
            raise ValueError(msg)
        self._constituent.vmr = np.ascontiguousarray(vmr).reshape(-1)

    def _validate_atmosphere(self, atmo: Atmosphere) -> None:
        if not isinstance(atmo.model_geometry, sk.Geometry2D):
            msg = "VMRAbsorber2D requires an atmosphere using Geometry2D"
            raise TypeError(msg)
        if tuple(atmo.volume_shape) != self._volume_shape:
            msg = (
                "VMRAbsorber2D shape does not match the atmosphere: "
                f"{self._volume_shape} != {atmo.volume_shape}"
            )
            raise ValueError(msg)

    def add_to_atmosphere(self, atmo: Atmosphere) -> None:
        self._validate_atmosphere(atmo)
        self._constituent.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: Atmosphere, name: str) -> None:
        self._validate_atmosphere(atmo)
        self._constituent.register_derivative(atmo, name)
