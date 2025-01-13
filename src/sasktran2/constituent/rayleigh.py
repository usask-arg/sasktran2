from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2.optical.rayleigh import rayleigh_cross_section_bates
from sasktran2.polarization import LegendreStorageView

from .base import Constituent


class Rayleigh(Constituent):
    def __init__(self, method: str = "bates", method_kwargs: dict | None = None):
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

        method_kwargs: dict, optional
            kwargs that can be passed to the method

        Raises
        ------
        ValueError
            If input method is not supported
        """
        self._method_kwargs = method_kwargs
        if method.lower() == "bates":
            self._rayleigh_cross_fn = rayleigh_cross_section_bates
            if self._method_kwargs is not None:
                self._fn_kwargs = self._method_kwargs
            else:
                self._fn_kwargs = {}
        elif method.lower() == "manual":

            def temp_fn(wv_micron: np.array):
                xs = np.interp(
                    wv_micron * 1000,
                    self._method_kwargs["wavelength_nm"],
                    self._method_kwargs["xs"],
                )
                king = np.interp(
                    wv_micron * 1000,
                    self._method_kwargs["wavelength_nm"],
                    self._method_kwargs["king_factor"],
                )

                return xs, king

            self._rayleigh_cross_fn = temp_fn
            self._fn_kwargs = {}
        else:
            msg = "Method must be bates or manual"
            raise ValueError(msg)

        self._ray_ext = None
        self._delta = None

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere


        :meta private:
        """
        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the Rayleigh constituent"
            raise ValueError(msg)

        if atmo.pressure_pa is None:
            msg = "It is required to set the pressure_pa property in the Atmosphere object to use the Rayleigh Constituent"
            raise ValueError(msg)

        if atmo.temperature_k is None:
            msg = "It is required to set the temperature_k property in the Atmosphere object to use the Rayleigh Constituent"
            raise ValueError(msg)

        # Get the number density from the atmosphere object at the grid points
        num_dens = atmo.state_equation.air_numberdensity["N"]

        scattering_xs, king_factor = self._rayleigh_cross_fn(
            atmo.wavelengths_nm / 1000, **self._fn_kwargs
        )

        # King factor F = (6 + 3 d) / (6 - 7 d)
        # So then (6 - 7 d) * F = 6 + 3 d
        # d (3 + 7 F) = 6(F - 1)
        # d = 6 (F - 1) / (3 + 7 F)
        self._delta = 6 * (king_factor - 1) / (3 + 7 * king_factor)

        self._ray_ext = np.outer(num_dens, scattering_xs)
        # Start by adding in the extinction
        atmo.storage.total_extinction += self._ray_ext

        # SSA temporarily stores the scattering extinction
        atmo.storage.ssa += self._ray_ext

        atmo.leg_coeff.a1[0] += self._ray_ext
        atmo.leg_coeff.a1[2] += (
            self._ray_ext * ((1 - self._delta) / (2 + self._delta))[np.newaxis, :]
        )

        if atmo.nstokes == 3:
            atmo.leg_coeff.a2[2] += (
                self._ray_ext
                * 6
                * ((1 - self._delta) / (2 + self._delta))[np.newaxis, :]
            )

            atmo.leg_coeff.b1[2] += (
                self._ray_ext
                * np.sqrt(6.0)
                * ((1 - self._delta) / (2 + self._delta))[np.newaxis, :]
            )
        elif atmo.nstokes == 4:
            msg = (
                "NSTOKES={} not currently implemented for Rayleigh constituent".format(
                    atmo.nstokes
                )
            )
            raise ValueError(msg)

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        num_dens = atmo.state_equation.air_numberdensity
        N = num_dens["N"]

        derivs = {}

        # dI/dP = dI/dN * dN/dP, and dI/dN is the extinction
        # It's easier to just treat this as a number density derivative using the Interpolated Derivative Mapping
        # and an Identity matrix as the mapping matrix

        xs = self._ray_ext / N[:, np.newaxis]

        deriv_names = []
        d_vals = []
        if atmo.calculate_pressure_derivative:
            deriv_names.append("pressure_pa")
            d_vals.append(num_dens["dN_dP"])
        if atmo.calculate_temperature_derivative:
            deriv_names.append("temperature_k")
            d_vals.append(num_dens["dN_dT"])

        for deriv_name, vert_factor in zip(deriv_names, d_vals, strict=True):
            mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_{deriv_name}")
            mapping.d_extinction[:] += xs
            mapping.d_ssa[:] += (
                xs * (1 - atmo.storage.ssa) / atmo.storage.total_extinction
            )
            mapping.d_leg_coeff[:] -= atmo.storage.leg_coeff
            mapping.scat_factor[:] += xs / (
                atmo.storage.ssa * atmo.storage.total_extinction
            )

            mapping.interpolator = np.eye(len(N)) * vert_factor[np.newaxis, :]
            mapping.interp_dim = "altitude"
            mapping.assign_name = f"wf_{deriv_name}"

            d_leg_coeff = LegendreStorageView(mapping.d_leg_coeff, atmo.nstokes)

            d_leg_coeff.a1[0][:] += 1
            d_leg_coeff.a1[2][:] += ((1 - self._delta) / (2 + self._delta))[
                np.newaxis, :
            ]

            if atmo.nstokes >= 3:
                d_leg_coeff.a2[2][:] += (
                    6 * ((1 - self._delta) / (2 + self._delta))[np.newaxis, :]
                )

                d_leg_coeff.b1[2][:] += (
                    np.sqrt(6.0)
                    * ((1 - self._delta) / (2 + self._delta))[np.newaxis, :]
                )

        return derivs
