from .base import Constituent
from sasktran2.optical.rayleigh import rayleigh_cross_section_bates
import sasktran2 as sk
from sasktran2.optical import pressure_temperature_to_numberdensity
from sasktran2.units import celsius_to_kelvin
import numpy as np


class Rayleigh(Constituent):
    def __init__(self, method: str ='bates', name='rayleigh'):
        """
        An implementation of Rayleigh scattering.  Cross sections (and depolarization factors) can be
        calculated multiple ways, with the default method being that of 'bates'.  

        Rayleigh scattering number density is estimated through the ideal gas law.

        This Constituent requires that the atmosphere object have `temperature_k`, `pressure_pa`, and
        `wavelength_nm` are all defined inside the :py:class:`sasktran2.Atmosphere` object.

        Parameters
        ----------
        method : str, optional
            Method to use to calculate the cross section.  Supported methods are
            ['bates'], by default 'bates'
        name : str, optional
            Name to give the constituent, by default 'rayleigh'

        Raises
        ------
        ValueError
            If input method is not supported
        """
        if method.lower() == 'bates':
            self._rayleigh_cross_fn = rayleigh_cross_section_bates
        else:
            msg = 'Method must be bates'
            raise ValueError(msg)

        self._name = name

    def name(self) -> str:
        return self._name

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere

 
        :meta private:
        """
        if atmo.wavelengths_nm is None:
            msg = 'It is required to give the Atmosphere object wavelengths to use the Rayleigh constituent'
            raise ValueError(msg)

        if atmo.pressure_pa is None:
            msg = 'It is required to set the pressure_pa property in the Atmosphere object to use the Rayleigh Constituent'
            raise ValueError(msg)

        if atmo.temperature_k is None:
            msg = 'It is required to set the temperature_k property in the Atmosphere object to use the Rayleigh Constituent'
            raise ValueError(msg)

        # Get the number density from the atmosphere object at the grid points
        num_dens = pressure_temperature_to_numberdensity(atmo.pressure_pa, atmo.temperature_k)

        scattering_xs, king_factor = self._rayleigh_cross_fn(atmo.wavelengths_nm / 1000)

        # King factor F = (6 + 3 d) / (6 - 7 d)
        # So then (6 - 7 d) * F = 6 + 3 d
        # d (3 + 7 F) = 6(F - 1)
        # d = 6 (F - 1) / (3 + 7 F)
        delta = 6 * (king_factor - 1) / (3 + 7 * king_factor)

        ray_ext = np.outer(num_dens, scattering_xs)
        # Start by adding in the extinction
        atmo.storage.total_extinction += ray_ext

        # SSA temporarily stores the scattering extinction
        atmo.storage.ssa += ray_ext

        if atmo.nstokes == 1:
            atmo.storage.leg_coeff[0] += ray_ext
            atmo.storage.leg_coeff[2] += ray_ext * ((1 - delta) / (2 + delta))[np.newaxis, :]
        else:
            msg = 'NSTOKES={} not currently implemented for Rayleigh constituent'.format(atmo.nstokes)
            raise ValueError(msg)


    def register_derivative(self, atmo: sk.Atmosphere):
        pass