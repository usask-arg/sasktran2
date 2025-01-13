from __future__ import annotations

import logging

import numpy as np

from sasktran2._core import LinearizedMie
from sasktran2.atmosphere import Atmosphere
from sasktran2.mie.distribution import ParticleSizeDistribution, integrate_mie
from sasktran2.mie.refractive import RefractiveIndex
from sasktran2.optical.base import OpticalProperty, OpticalQuantities
from sasktran2.polarization import LegendreStorageView


class Mie(OpticalProperty):
    def __init__(
        self,
        psize_distribution: ParticleSizeDistribution,
        refractive_index: RefractiveIndex,
    ):
        """
        Mie scattering optical property where the Mie calculations are done on the fly.

        Note, this is much slower than using a precomputed Mie table, and in most cases it is
        recommended to use the MieDatabase class instead.

        Parameters
        ----------
        psize_distribution : ParticleSizeDistribution
            Particle size distribution to use for the Mie calculations
        refractive_index : RefractiveIndex
            Refraction index to use for the Mie calculations
        """
        self._psize_distribution = psize_distribution

        self._geometry_dependent = len(psize_distribution.args()) > 0

        self._refrac_index = refractive_index

        self._calculation_ds = {}
        self._xs_ds = {}

    def _internal_generate(self, wavelengths_nm, result_ds, **kwargs):
        if len(wavelengths_nm) > 20:
            logging.warning(
                "Calculating Mie scattering parameters for a large number of wavelengths. You may want to use the MieDatabase class instead."
            )

        dist_args = self._psize_distribution.args()

        compute_coeff = "num_legendre" in kwargs

        for arg in dist_args:
            if arg not in kwargs:
                msg = f"Missing argument {arg} for particle size distribution"
                raise ValueError(msg)

        if self._geometry_dependent:
            # Find the unique values of the distribution arguments
            arg_tuples = np.array([kwargs[arg] for arg in dist_args]).T
            unique_args = np.unique(arg_tuples, axis=0)

            if len(unique_args) > 20:
                logging.warning(
                    "Calculating Mie scattering parameters for a large number of particle size distribution arguments. You may want to use the MieDatabase class instead."
                )

            for args in unique_args:
                # Generate the Mie scattering cross sections for each unique set of distribution arguments
                ds = integrate_mie(
                    LinearizedMie(kwargs.get("num_threads", 1)),
                    self._psize_distribution.distribution(
                        **dict(zip(dist_args, args, strict=True))
                    ),
                    self._refrac_index.refractive_index_fn,
                    wavelengths_nm,
                    compute_coeffs=compute_coeff,
                    num_coeffs=kwargs.get("num_legendre", 10),
                )

                result_ds[tuple(args)] = ds
        else:
            ds = integrate_mie(
                LinearizedMie(kwargs.get("num_threads", 1)),
                self._psize_distribution.distribution(),
                self._refrac_index.refractive_index_fn,
                wavelengths_nm,
                compute_coeffs=compute_coeff,
                num_coeffs=kwargs.get("num_legendre", 10),
            )

            result_ds[()] = ds

    def atmosphere_quantities(self, atmo: Atmosphere, **kwargs) -> OpticalQuantities:
        self._internal_generate(
            atmo.wavelengths_nm,
            self._calculation_ds,
            **{
                **kwargs,
                "num_legendre": atmo.leg_coeff.a1.shape[0],
                "num_threads": atmo._config.num_threads,
            },
        )

        # Then calculate the atmosphere quantities
        quants = OpticalQuantities(
            extinction=np.zeros(
                (len(atmo.model_geometry.altitudes()), len(atmo.wavelengths_nm))
            ),
            ssa=np.zeros(
                (len(atmo.model_geometry.altitudes()), len(atmo.wavelengths_nm))
            ),
        )

        quants.leg_coeff = np.zeros_like(atmo.storage.leg_coeff)

        leg_coeff = LegendreStorageView(quants.leg_coeff, atmo.nstokes)

        if self._geometry_dependent:
            arg_tuples = np.array(
                [kwargs[arg] for arg in self._psize_distribution.args()]
            ).T
            for i, arg in enumerate(arg_tuples):
                ds = self._calculation_ds[tuple(arg)].sel(
                    wavelength=atmo.wavelengths_nm
                )

                # Convert nm^2 to m^2
                quants.extinction[i, :] = ds["xs_total"].to_numpy() / 1e18
                quants.ssa[i, :] = ds["xs_scattering"].to_numpy() / 1e18

                leg_coeff.a1[:, i, :] = ds["lm_a1"].to_numpy().T
                if atmo.nstokes == 3:
                    leg_coeff.a2[:, i, :] = ds["lm_a2"].to_numpy().T
                    leg_coeff.b1[:, i, :] = ds["lm_b1"].to_numpy().T
                    leg_coeff.a3[:, i, :] = ds["lm_a3"].to_numpy().T
        else:
            ds = self._calculation_ds[()].sel(wavelength=atmo.wavelengths_nm)

            # Convert nm^2 to m^2
            quants.extinction[:] = ds["xs_total"].to_numpy()[np.newaxis, :] / 1e18
            quants.ssa[:] = ds["xs_scattering"].to_numpy()[np.newaxis, :] / 1e18

            leg_coeff.a1[:] = ds["lm_a1"].to_numpy().T[:, np.newaxis, :]
            if atmo.nstokes == 3:
                leg_coeff.a2[:] = ds["lm_a2"].to_numpy().T[:, np.newaxis, :]
                leg_coeff.b1[:] = ds["lm_b1"].to_numpy().T[:, np.newaxis, :]
                leg_coeff.a3[:] = ds["lm_a3"].to_numpy().T[:, np.newaxis, :]

        quants.extinction[np.isnan(quants.extinction)] = 0
        quants.ssa[np.isnan(quants.ssa)] = 0

        return quants

    def cross_sections(self, wavelengths_nm, altitudes_m, **kwargs):
        self._internal_generate(wavelengths_nm, self._xs_ds, **kwargs)

        # Then calculate the atmosphere quantities
        quants = OpticalQuantities(
            extinction=np.zeros((len(altitudes_m), len(wavelengths_nm))),
            ssa=np.zeros((len(altitudes_m), len(wavelengths_nm))),
        )

        if self._geometry_dependent:
            arg_tuples = np.array(
                [kwargs[arg] for arg in self._psize_distribution.args()]
            ).T
            for i, arg in enumerate(arg_tuples):
                ds = self._xs_ds[tuple(arg)].sel(wavelength=wavelengths_nm)

                # Convert nm^2 to m^2
                quants.extinction[i, :] = ds["xs_total"].to_numpy() / 1e18
                quants.ssa[i, :] = (
                    ds["xs_scattering"].to_numpy() / 1e18
                ) / quants.extinction[i, :]

            quants.extinction[np.isnan(quants.extinction)] = 0
            quants.ssa[np.isnan(quants.ssa)] = 0
        else:
            ds = self._xs_ds[()].sel(wavelength=wavelengths_nm)

            # Convert nm^2 to m^2
            quants.extinction[:] = ds["xs_total"].to_numpy() / 1e18
            quants.ssa[:] = (ds["xs_scattering"].to_numpy() / 1e18) / quants.extinction

            quants.extinction[np.isnan(quants.extinction)] = 0
            quants.ssa[np.isnan(quants.ssa)] = 0

        return quants
