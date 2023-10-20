import abc
from dataclasses import dataclass

import numpy as np

from sasktran2.atmosphere import Atmosphere


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

    def cross_sections(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs  # noqa: ARG002
    ) -> OpticalQuantities:
        msg = "Not Supported"
        raise NotImplementedError(msg)
