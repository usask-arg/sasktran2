import abc
import numpy as np
from dataclasses import dataclass


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
    def atmosphere_quantities(self, atmo: Atmosphere) -> OpticalQuantities:
        pass
