import abc

import numpy as np
from scipy.stats import lognorm, rv_continuous


class ParticleSizeDistribution(abc.ABC):
    def __init__(self, identifier: str) -> None:
        """
        Abstract class to define particle size distributions that Mie parameters can be
        integrated over.  This class is a light wrapper on top of scipy.stats.rv_continuous
        which adds some additional information.

        Parameters
        ----------
        identifier : str
            A unique identifier for the distribution
        """
        self._identifier = identifier

    @abc.abstractmethod
    def distribution(self, **kwargs) -> rv_continuous:
        """
        Returns back the scipy object representing this distribution

        Returns
        -------
        rv_continuous
        """
        return self._distribution

    @property
    def identifier(self) -> str:
        """
        Get the unique identifier for this distribution

        Returns
        -------
        str
        """
        return self._identifier


class LogNormalDistribution(ParticleSizeDistribution):
    def __init__(self) -> None:
        """
        A log normal particle size distribution, defined by two parameters, the median radius and mode width
        """
        super().__init__("lognormal")

    def distribution(self, **kwargs):
        return lognorm(np.log(kwargs["mode_width"]), scale=kwargs["median_radius"])

    @staticmethod
    def args():
        return ["median_radius", "mode_width"]
