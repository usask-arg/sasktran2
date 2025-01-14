from __future__ import annotations

import abc
import logging
import time

import numpy as np
import scipy.integrate as integrate
import xarray as xr
from scipy.stats import lognorm, rv_continuous

from sasktran2 import mie
from sasktran2.legendre import compute_greek_coefficients
from sasktran2.mie import MieOutput


def _post_process(total: dict, wavelengths):
    """
    Internal method to post-process the Mie solution.  This involves transforming the quantities that we integrate over
    to the relevant bulk properties and applying the necessary normalizations.
        Parameters
    ----------
    total : dict
        Dictionary containing calculated and integrated mie parameters
    wavelengths : np.array
        Wavelengths of interest

    """
    k = 2 * np.pi / wavelengths
    c = 4 * np.pi / (2 * k**2 * total["xs_scattering"])

    total["p11"] *= c
    total["p12"] *= c
    total["p33"] *= c
    total["p34"] *= c

    return


def _integrable_parameters(output: MieOutput, pdf, x):
    """
    Helper function to set up the parameters to be integrated

    Parameters
    ----------
    output : MieOutput
        MieOutput structure containing angles, siza parameters, refractive index, and calculated mie parameters
    pdf : np.array
        Probability distribution evaluated for radii x
    x : np.array
        Radii of the distribution
    Returns
    -------
    Dict
        Dictionary containing the keys 'p11', 'p12', 'p33', 'p34', 'Cext', 'Csca', and 'Cabs'
    """

    out = {}
    s1 = output.values.S1
    s2 = output.values.S2

    # Phase function elements
    # Note that P22 = P11 and P44 = p33 so we have 6 independent elements
    out["p11"] = np.abs(s1) ** 2 + np.abs(s2) ** 2
    out["p12"] = np.abs(s1) ** 2 - np.abs(s2) ** 2
    out["p33"] = np.real(s1 * np.conj(s2) + s2 * np.conj(s1))
    out["p34"] = np.real(-1j * (s1 * np.conj(s2) - s2 * np.conj(s1)))

    # Cross sections in units of input radius**2
    out["Cext"] = output.values.Qext * np.pi * np.square(x)
    out["Csca"] = output.values.Qsca * np.pi * np.square(x)
    out["Cabs"] = (output.values.Qext - output.values.Qsca) * np.pi * np.square(x)

    for key in ["Cext", "Csca", "Cabs"]:
        out[key] *= pdf

    for key in ["p11", "p12", "p33", "p34"]:
        out[key] *= np.transpose([pdf])

    return out


def integrate_mie(
    my_mie: mie,
    prob_dist: rv_continuous,
    refrac_index_fn,
    wavelengths,
    num_angles=1801,
    num_quad=1024,
    maxintquantile=0.99999,
    compute_coeffs=False,
    num_coeffs=64,
):
    """
    Integrates the Mie parameters over an arbitrary particle size distribution, returning cross sections and phase matrices
      at a set of wavelengths.

    Note that the units of the input parameters are left arbitrary, but they must be consistent between each other,
    i.e., if wavelengths is specified in nm it is expected that the particle size distribution will be as a function
    of nm.

    Parameters
    ----------
    mie : sk.Mie
        Internal Mie algorithms to use
    prob_dist : rv_continuous
        Instance of a probability distribution to integrate over.  Distribution must be specified in the same units
        as wavelengths
    refrac_index_fn : fn
        Function taking in one argument (wavelength) and returning a complex number for refractive index at that
        wavelength.  Function argument should be in the same units as wavelengths.
    wavelengths : np.array
        Wavelengths to calculate for.  Units should be consistent with prob_dist
    num_angles: int, Optional
        Number of angles at which to evaluate p11, p12, p33, p34
    num_quad : int, Optional
        Number of quadrature points to use when integrating over particle size.  Default is 1024
    maxintquantile : float
        Used to determine the maximum radius to use in integration.  A value of 0.99 means that at least 99% of the
        probability distribution area multipied by radius**2 will be included within the integration bounds.  Default
        0.99999
    compute_coeffs : bool
        Option to also compute and return the greek coefficients. Default False.
    num_coeffes : int
        Optional parameter, maximum number of coefficients to return in the expansion. Default 64.
    Returns
    -------
    xr.Dataset
        Dataset containing the keys 'p11', 'p12', 'p33', 'p34', 'xs_total', 'xs_absorption', and 'xs_scattering'
        All scattering cross sections will be in units of wavelength**2. If compute_coeffs is True, will also
        contain keys 'lm_a1', 'lm_a2', 'lm_a3', 'lm_a4', 'lm_b1', and 'lm_b2'.
    """
    t_full = time.time()

    angles = np.linspace(0, 180, num_angles)

    # TODO FROM HERE
    norm = integrate.quad(
        lambda x: prob_dist.pdf(x) * x**2, 0, 1e25, points=(prob_dist.mean())
    )[0]

    def pdf(x):
        return prob_dist.pdf(x) * x**2 / norm

    # We have to determine the maximum R to integrate to, this is apparently very challenging, what we do is
    # calculate the CDF of repeatadly large values until it is greater than our threshold
    max_r = prob_dist.mean()
    while (
        integrate.quad(pdf, 0, max_r * 2, points=(prob_dist.mean()))[0]
        - integrate.quad(pdf, 0, max_r, points=(prob_dist.mean()))[0]
    ) > (1 - maxintquantile):
        max_r *= 2

    from scipy.special import roots_legendre

    x, w = roots_legendre(num_quad)
    x = 0.5 * (x + 1) * max_r
    w *= max_r / 2
    # TODO TO HERE, need to fix

    all_output = xr.Dataset(
        {
            "p11": (["wavelength", "angle"], np.zeros((len(wavelengths), len(angles)))),
            "p12": (["wavelength", "angle"], np.zeros((len(wavelengths), len(angles)))),
            "p33": (["wavelength", "angle"], np.zeros((len(wavelengths), len(angles)))),
            "p34": (["wavelength", "angle"], np.zeros((len(wavelengths), len(angles)))),
            "xs_total": (["wavelength"], np.zeros(len(wavelengths))),
            "xs_scattering": (["wavelength"], np.zeros(len(wavelengths))),
            "xs_absorption": (["wavelength"], np.zeros(len(wavelengths))),
        },
        coords={"wavelength": wavelengths, "angle": angles},
    )

    for idx, wavelength in enumerate(wavelengths):
        n = np.cdouble(refrac_index_fn(wavelength))

        size_param = 2 * np.pi * x / wavelength
        t = time.time()
        output = my_mie.calculate(size_param, n, np.cos(angles * np.pi / 180), False)
        logging.debug(f"Internal Mie Calculation {time.time() - t}")
        params = _integrable_parameters(output, prob_dist.pdf(x), x)

        # performing the integral
        all_output["p11"].values[idx, :] = np.sum(
            params["p11"] * np.transpose([w]), axis=0
        )
        all_output["p12"].values[idx, :] = np.sum(
            params["p12"] * np.transpose([w]), axis=0
        )
        all_output["p33"].values[idx, :] = np.sum(
            params["p33"] * np.transpose([w]), axis=0
        )
        all_output["p34"].values[idx, :] = np.sum(
            params["p34"] * np.transpose([w]), axis=0
        )
        all_output["xs_total"].values[idx] = np.sum(params["Cext"] * w)
        all_output["xs_scattering"].values[idx] = np.sum(params["Csca"] * w)
        all_output["xs_absorption"].values[idx] = np.sum(params["Cabs"] * w)
    _post_process(all_output, wavelengths)

    if compute_coeffs:
        # Note that P22 = P11 and P44 = p33
        lm_a1, lm_a2, lm_a3, lm_a4, lm_b1, lm_b2 = compute_greek_coefficients(
            p11=all_output["p11"].data,
            p12=all_output["p12"].data,
            p22=all_output["p11"].data,
            p33=all_output["p33"].data,
            p34=all_output["p34"].data,
            p44=all_output["p33"].data,
            angle_grid=angles,
            num_coeff=num_coeffs,
        )
        coeffs_output = xr.Dataset(
            {
                "lm_a1": (["wavelength", "legendre"], lm_a1),
                "lm_a2": (["wavelength", "legendre"], lm_a2),
                "lm_a3": (["wavelength", "legendre"], lm_a3),
                "lm_a4": (["wavelength", "legendre"], lm_a4),
                "lm_b1": (["wavelength", "legendre"], lm_b1),
                "lm_b2": (["wavelength", "legendre"], lm_b2),
            },
            coords={"wavelength": wavelengths, "legendre": np.arange(num_coeffs)},
        )
        all_output = all_output.merge(coeffs_output)

    logging.debug(f"Total Mie Calculation {time.time() - t_full}")
    return all_output


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
