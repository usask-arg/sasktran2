from __future__ import annotations

import abc
from dataclasses import dataclass

import numpy as np

from sasktran2.atmosphere import Atmosphere, NativeGridDerivative


@dataclass
class OpticalQuantities:
    extinction: np.ndarray = None
    ssa: np.ndarray = None
    a1: np.ndarray = None
    a2: np.ndarray = None
    a3: np.ndarray = None
    a4: np.ndarray = None
    b1: np.ndarray = None
    b2: np.ndarray = None


class OpticalProperty(abc.ABC):
    @abc.abstractmethod
    def atmosphere_quantities(self, atmo: Atmosphere, **kwargs) -> OpticalQuantities:
        pass

    def optical_derivatives(self, atmo: Atmosphere, **kwargs) -> dict:  # noqa: ARG002
        return {}

    def cross_sections(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs  # noqa: ARG002
    ) -> OpticalQuantities:
        msg = "Not Supported"
        raise NotImplementedError(msg)

    def cross_section_derivatives(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs  # noqa: ARG002
    ) -> dict:
        return {}

    def _into_rust_object(self):
        return None

    def __add__(self, other):
        return AdditiveOpticalProperty(self, other)


class AdditiveOpticalProperty(OpticalProperty):
    def __init__(self, property1: OpticalProperty, property2: OpticalProperty):
        self._property1 = property1
        self._property2 = property2

    def atmosphere_quantities(self, atmo: Atmosphere, **kwargs) -> OpticalQuantities:
        res1 = self._property1.atmosphere_quantities(atmo, **kwargs)
        res2 = self._property2.atmosphere_quantities(atmo, **kwargs)
        return OpticalQuantities(
            extinction=res1.extinction + res2.extinction, ssa=res1.ssa + res2.ssa
        )

    def optical_derivatives(self, atmo: Atmosphere, **kwargs) -> dict:
        res1 = self._property1.optical_derivatives(atmo, **kwargs)
        res2 = self._property2.optical_derivatives(atmo, **kwargs)

        all_keys = set(res1.keys()).union(set(res2.keys()))
        result = {}
        for k in all_keys:
            if k not in res1:
                # Just in res2
                result[k] = res2[k]
            elif k not in res2:
                # Just in res1
                result[k] = res1[k]
            else:
                result[k] = NativeGridDerivative(
                    d_extinction=res1[k].d_extinction + res2[k].d_extinction
                )

        return result
