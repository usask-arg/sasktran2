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

    @abc.abstractmethod
    def args(self):
        """
        A list of arguments that are required to define this distribution when calling distribution
        """

    def freeze(self, **kwargs):
        """
        Freeze some of the arguments of this distribution. E.g. if `y` is an argument of this distrubtion, calling
        `freeze(y=5)` will return a new distribution that is the same as this one, but with `y` frozen to 5.

        Returns
        -------
        ParticleSizeDistribution
            A new distribution with some of args frozen
        """
        return FrozenDistribution(self, kwargs)


class LogNormalDistribution(ParticleSizeDistribution):
    def __init__(self) -> None:
        """
        A log normal particle size distribution, defined by two parameters, the median radius and mode width
        """
        super().__init__("lognormal")

    def distribution(self, **kwargs):
        return lognorm(np.log(kwargs["mode_width"]), scale=kwargs["median_radius"])

    def args(self):
        return ["median_radius", "mode_width"]


class FrozenDistribution(ParticleSizeDistribution):
    def __init__(
        self, base_distribution: ParticleSizeDistribution, frozen_parameters: dict
    ) -> None:
        """
        A particle size distribution that is frozen in time, useful for testing

        Parameters
        ----------
        base_distribution : ParticleSizeDistribution
            The distribution to freeze
        """
        identifier = f"frozen_{base_distribution.identifier}"
        for key, value in frozen_parameters.items():
            identifier += f"_{key}_{value}"

            if key not in base_distribution.args():
                msg = f"Frozen key {key} not in base distribution args"
                raise ValueError(msg)

        super().__init__(identifier)
        self._distribution = base_distribution

        self._frozen_parameters = frozen_parameters
        self._args = [
            arg for arg in base_distribution.args() if arg not in frozen_parameters
        ]

    def distribution(self, **kwargs):
        return self._distribution.distribution(**{**self._frozen_parameters, **kwargs})

    def args(self):
        return self._args
