import abc
from sasktran2 import Atmosphere


class Constituent(abc.ABC):
    @abc.abstractmethod
    def name(self) -> str:
        pass

    @abc.abstractmethod
    def add_to_atmosphere(self, atmo: Atmosphere):
        pass

    @abc.abstractmethod
    def register_derivative(self, atmo: Atmosphere):
        pass