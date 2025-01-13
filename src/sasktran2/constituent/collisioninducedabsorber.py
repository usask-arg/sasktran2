from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2.optical import pressure_temperature_to_numberdensity

from ..optical.base import OpticalProperty
from .base import Constituent


class CollisionInducedAbsorber(Constituent):
    def __init__(self, optical_property: OpticalProperty, name: str):
        """
        An implementation of collision-induced absorption for quantities with known constant mole
        fractions.

        Air number density is estimated through the ideal gas law. Currently the mole fraction is applied
        to the air number density, though in the future it should be applied to dry air only.

        This Constituent requires that the atmosphere object have `temperature_k`, `pressure_pa`, and
        `wavelength_nm` all defined inside the :py:class:`sasktran2.Atmosphere` object.

        Parameters
        ----------
        optical_property : OpticalProperty
            Optical property that computes the collision-induced cross section in [m^5].
        name : str
            Molecule names to load mole fractions.  Currently supported: O2O2
        """
        self._optical_property = optical_property
        if name.lower() == "o2o2":
            self._fraction_product = 0.20964**2
        else:
            msg = f"Unknown name '{name}'"
            raise ValueError(msg)

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere

        :meta private:
        """
        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the O2O2 constituent"
            raise ValueError(msg)

        if atmo.pressure_pa is None:
            msg = "It is required to set the pressure_pa property in the Atmosphere object to use the O2O2 Constituent"
            raise ValueError(msg)

        if atmo.temperature_k is None:
            msg = "It is required to set the temperature_k property in the Atmosphere object to use the O2O2 Constituent"
            raise ValueError(msg)

        # Get the number density from the atmosphere object at the grid points (m^-3)
        num_dens = pressure_temperature_to_numberdensity(
            atmo.pressure_pa, atmo.temperature_k
        )

        self._optical_quants = self._optical_property.atmosphere_quantities(atmo)

        atmo.storage.total_extinction += (  # m^-1
            self._optical_quants.extinction  # m^5 (cross section)
            * (self._fraction_product * num_dens**2)[:, np.newaxis]
        )

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        number_density, dN_dP, dN_dT = pressure_temperature_to_numberdensity(
            atmo.pressure_pa, atmo.temperature_k, include_derivatives=True
        )
        derivs = {}

        deriv_names = []
        d_vals = []
        if atmo.calculate_pressure_derivative:
            deriv_names.append("pressure_pa")
            d_vals.append(dN_dP)
        if atmo.calculate_temperature_derivative:
            deriv_names.append("temperature_k")
            d_vals.append(dN_dT)

        # Contributions from the change in number density via pressure/temperature
        for deriv_name, dN_dX in zip(deriv_names, d_vals, strict=False):
            dk_dX = (
                2
                * self._fraction_product
                * (number_density * dN_dX)[:, np.newaxis]
                * self._optical_quants.extinction
            )
            deriv_mapping = atmo.storage.get_derivative_mapping(
                f"wf_{name}_{deriv_name}"
            )
            deriv_mapping.d_extinction[:] += dk_dX
            deriv_mapping.d_ssa[:] += (
                dk_dX
                * (self._optical_quants.ssa - atmo.storage.ssa)
                / atmo.storage.total_extinction
            )
            deriv_mapping.interpolator = np.eye(len(number_density))
            deriv_mapping.interp_dim = "altitude"
            deriv_mapping.assign_name = f"wf_{deriv_name}"

        if len(deriv_names) > 0:
            optical_derivs = self._optical_property.optical_derivatives(atmo=atmo)

            for key, val in optical_derivs.items():
                deriv_mapping = atmo.storage.get_derivative_mapping(
                    f"wf_{name}_{deriv_name}_{key}_xs"
                )
                deriv_mapping.d_extinction[:] += val.d_extinction
                deriv_mapping.d_ssa[:] += (
                    val.d_extinction
                    * (self._optical_quants.ssa - atmo.storage.ssa)
                    / atmo.storage.total_extinction
                )
                deriv_mapping.interpolator = (
                    np.eye(len(number_density))
                    * (self._fraction_product * number_density**2)[np.newaxis, :]
                )
                deriv_mapping.interp_dim = "altitude"
                deriv_mapping.assign_name = f"wf_{deriv_name}"

        return derivs
