from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import PyThermalEmission
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


def d_planck_blackbody_radiance_d_temperature(
    temperature_k: np.ndarray, wavelengths_nm: np.ndarray
) -> np.ndarray:
    """
    Calculates the derivative of the Planck function for a given set of temperatures and wavelengths
    with respect to the temperature parameter

    Parameters
    ----------
    temperature_k : np.ndarray
        Temperatures in [K]
    wavelengths_nm : np.ndarray
        Wavelengths in [nm] to calculate the radiance at

    Returns
    -------
    np.ndarray
        The derivative of the blackbody radiance with respect to temperature with units of W / (m^2 nm sr K). The shape of the returned
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
        / (np.exp(exponent) - 1) ** 2
        * (
            PLANCK
            * SPEED_OF_LIGHT
            / (wavelengths_m * K_BOLTZMANN * temperature_k[:, np.newaxis] ** 2)
        )
        * np.exp(exponent)
        * 1e-9
    )


class ThermalEmission(Constituent):
    _thermal_emission: PyThermalEmission

    def __init__(self):
        """
        An implementation of thermal emissions calculated from the Planck function. The emission is
        calculated with units of [W / (m^2 nm sr)].

        This Constituent requires that the atmosphere object have `temperature_k` and
        `wavelength_nm` defined inside the :py:class:`sasktran2.Atmosphere` object.
        """
        self._thermal_emission = PyThermalEmission()

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere


        :meta private:
        """
        self._thermal_emission.add_to_atmosphere(atmo)

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        self._thermal_emission.register_derivative(atmo, name)


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

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        atmo.surface.emission[:] += (
            self.emissivity
            * planck_blackbody_radiance(
                np.atleast_1d(self.temperature_k), atmo.wavelengths_nm
            ).flatten()
        )

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        planck = planck_blackbody_radiance(
            np.atleast_1d(self.temperature_k), atmo.wavelengths_nm
        )

        # Surface temperature derivative
        deriv_mapping = atmo.surface.get_derivative_mapping(f"wf_{name}_temperature_k")
        deriv_mapping.d_emission[:] = (
            self.emissivity
            * d_planck_blackbody_radiance_d_temperature(
                np.atleast_1d(self.temperature_k), atmo.wavelengths_nm
            )
        )
        deriv_mapping.interp_dim = "dummy"
        deriv_mapping.interpolator = np.ones((len(atmo.wavelengths_nm), 1))

        # Surface emissivity derivative
        deriv_mapping = atmo.surface.get_derivative_mapping(f"wf_{name}_emissivity")
        deriv_mapping.d_emission[:] = planck.flatten()
        deriv_mapping.interp_dim = "dummy"

        # If emissivity is scalear, we interpolate
        if len(self.emissivity) == 1:
            deriv_mapping.interpolator = np.ones((len(atmo.wavelengths_nm), 1))

        # Else don't need an interpolator
