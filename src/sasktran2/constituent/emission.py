import numpy as np

import sasktran2 as sk
from sasktran2.atmosphere import Atmosphere
from sasktran2.constants import K_BOLTZMANN, PLANCK, SPEED_OF_LIGHT

from .base import Constituent


def planck_blackbody_radiance(
    temperature_k: np.ndarray, wavelengths_nm: np.ndarray
) -> np.ndarray:
    """
    Calculates the Planck function for a given set of temperatures and wavelengths.

    Parameters
    ----------
    temperature_k : np.ndarray
        Temperatures in [K]
    wavelengths_nm : np.ndarray
        Wavelengths in [nm] to calculate the radiance at

    Returns
    -------
    np.ndarray
        The blackbody radiance with units of W / (m^2 nm sr). The shape of the returned
        array is shape(len(temperature_k), len(wavelengths_nm)).
    """
    wavelengths_m = wavelengths_nm * 1e-9
    exponent = (
        PLANCK
        * SPEED_OF_LIGHT
        / (wavelengths_m * K_BOLTZMANN * temperature_k[:, np.newaxis])
    )
    return (
        (2 * PLANCK * SPEED_OF_LIGHT**2 / wavelengths_m**5)
        / (np.exp(exponent) - 1)
        * 1e-9
    )


class ThermalEmission(Constituent):
    def __init__(self):
        """
        An implementation of thermal emissions calculated from the Planck function. The emission is
        calculated with units of [W / (m^2 nm sr)].

        This Constituent requires that the atmosphere object have `temperature_k` and
        `wavelength_nm` defined inside the :py:class:`sasktran2.Atmosphere` object.
        """
        return

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere

        :meta private:
        """
        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the ThermalEmission constituent"
            raise ValueError(msg)

        if atmo.temperature_k is None:
            msg = "It is required to set the temperature_k property in the Atmosphere object to use the ThermalEmission Constituent"
            raise ValueError(msg)

        # calculate radiance in W / (m^2 nm sr)
        atmo.storage.emission_source += planck_blackbody_radiance(
            atmo.temperature_k, atmo.wavelengths_nm
        )

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        return super().register_derivative(atmo=atmo, name=name)


class SurfaceThermalEmission(Constituent):
    def __init__(
        self,
        temperature_k: float,
        emissivity: np.ndarray,
    ) -> None:
        """
        A thermally emissive surface that is defined by a temperature and emissivity. The
        emission is calculated as the product of emissivity and the Planck function.

        This can either operate in a "scalar" mode where the emissivity is constant in wavelength,
        or a "native" mode where the emissivity is defined on the same grid as the atmosphere.

        Parameters
        ----------
        temperature_k : float
            Surface temperature.
        emissivity : np.ndarray
            Surface emissivity. Can be a scalar to indicate it is constant in wavelength. If set to an
            array it must match the atmosphere wavelength grid dimenstion.
        """
        Constituent.__init__(self)
        self._emissivity = np.atleast_1d(emissivity)
        self._temperature_k = temperature_k

    @property
    def temperature_k(self) -> float:
        return self._temperature_k

    @temperature_k.setter
    def temperature_k(self, temperature_k: float):
        self._temperature_k = temperature_k

    @property
    def emissivity(self) -> np.ndarray:
        return self._emissivity

    @emissivity.setter
    def emissivity(self, emissivity: np.ndarray):
        self._emissivity = np.atleast_1d(emissivity)

    def add_to_atmosphere(self, atmo: Atmosphere):
        atmo.surface.emission[:] += (
            self.emissivity
            * planck_blackbody_radiance(
                np.atleast_1d(self.temperature_k), atmo.wavelengths_nm
            ).flatten()
        )

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        return super().register_derivative(atmo=atmo, name=name)
