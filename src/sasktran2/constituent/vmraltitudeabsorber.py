from __future__ import annotations

import numpy as np

from sasktran2._core_rust import PyVMRAltitudeAbsorber
from sasktran2.atmosphere import Atmosphere
from sasktran2.optical.base import OpticalProperty

from .base import Constituent


class VMRAltitudeAbsorber(Constituent):
    _constituent: PyVMRAltitudeAbsorber

    def __init__(
        self,
        optical_property: OpticalProperty,
        altitudes_m: np.array,
        vmr: np.array,
        out_of_bounds_mode: str = "zero",
    ) -> None:
        """
        An atmospheric constituent that is specified through volume mixing ratio (VMR) on an altitude grid.
        The altitude grid need not match the global atmospheric grid, in the case they are different the VMRs
        will be linearly interpolated to the atmospheric grid locations.

        The constituent has an associated OpticalProperty that defines the absorption cross section as a function
        of wavelength.

        Parameters
        ----------
        optical_property : OpticalProperty
            An object that provides absorption cross sections as a function of wavelength
        altitudes_m : np.array
            Altitudes in [m] that the VMR is specified on
        vmr : np.array
            Volume mixing ratio
        out_of_bounds_mode : str, optional
            One of ['zero', 'extend'], 'zero' will set the VMR to be 0 on any altitudes that are
            out of bounds.  'extend' will extend the last or first VMR value.  By default 'zero'
        """
        super().__init__()

        self._constituent = PyVMRAltitudeAbsorber(
            optical_property,
            altitudes_m.astype(np.float64),
            vmr.astype(np.float64),
            out_of_bounds_mode,
        )

        self._temp = optical_property

    @property
    def vmr(self):
        return self._constituent.vmr

    @vmr.setter
    def vmr(self, vmr: np.array):
        self._constituent.vmr = vmr

    @property
    def altitudes_m(self):
        return self._constituent.altitudes_m

    def add_to_atmosphere(self, atmo: Atmosphere):
        self._constituent.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: Atmosphere, name: str):
        self._constituent.register_derivative(atmo, name)
