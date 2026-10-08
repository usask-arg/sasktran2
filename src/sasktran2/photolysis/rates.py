"""Photolysis and photoexcitation rates from an actinic flux."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass

import numpy as np
import xarray as xr

QuantumYield = float | Callable[[np.ndarray, np.ndarray], np.ndarray]


@dataclass(frozen=True)
class Photolysis:
    """Rate of a continuum (or resolved-line) absorption process, per molecule.

    ``J(z) = integral of F(lambda, z) sigma(lambda, z) phi(lambda, T(z)) d lambda``
    over ``wavelength_range_nm``, with ``F`` the actinic flux and ``sigma``
    the absorber's cross section from the flux calculation.

    Parameters
    ----------
    name
        Rate name, e.g. a mechanism rate input such as ``"J_O3_O1D"``.
    absorber
        Species whose cross section is used; must be in the flux dataset.
    quantum_yield
        Constant, or a function of ``(wavelength_nm, temperature_k)`` arrays
        that broadcast to ``(wavelength, altitude)``.
    wavelength_range_nm
        Inclusive integration range; the whole grid when None.
    max_grid_spacing_nm
        If set, the wavelength grid must be at least this fine over the
        range, e.g. so that absorption lines are resolved.
    reference
        Source of the cross section and quantum yield, for provenance.
    """

    name: str
    absorber: str
    quantum_yield: QuantumYield = 1.0
    wavelength_range_nm: tuple[float, float] | None = None
    max_grid_spacing_nm: float | None = None
    reference: str = ""


@dataclass(frozen=True)
class LinePhotolysis:
    """Rate of absorption in a single unresolved solar line, per molecule.

    ``J(z) = phi * sigma_eff * F_line * F(lambda_0, z) / F_toa(lambda_0)``:
    the top-of-atmosphere line-integrated photon flux ``F_line`` scaled by
    the actinic-flux transmission at the line centre.
    """

    name: str
    line_center_nm: float
    cross_section_m2: float
    toa_line_flux_photons_m2_s: float
    quantum_yield: float = 1.0
    reference: str = ""


def _quantum_yield(reaction: Photolysis, wavelength_nm, temperature_k) -> np.ndarray:
    if callable(reaction.quantum_yield):
        return reaction.quantum_yield(
            wavelength_nm[:, np.newaxis], temperature_k[np.newaxis, :]
        )
    return np.asarray(reaction.quantum_yield, dtype=float)


def _continuum_rate(flux: xr.Dataset, reaction: Photolysis) -> np.ndarray:
    wavelength = flux["wavelength"].to_numpy()
    if reaction.absorber not in flux["species"]:
        msg = (
            f"{reaction.name}: no cross section for {reaction.absorber} in the flux "
            f"dataset; it has {list(flux['species'].values)}"
        )
        raise ValueError(msg)

    lo, hi = reaction.wavelength_range_nm or (wavelength[0], wavelength[-1])
    selected = (wavelength >= lo) & (wavelength <= hi)
    if selected.sum() < 2:
        msg = f"{reaction.name}: fewer than two wavelengths in {lo}-{hi} nm"
        raise ValueError(msg)
    w = wavelength[selected]
    if reaction.max_grid_spacing_nm is not None and np.diff(
        w
    ).max() > reaction.max_grid_spacing_nm * (1 + 1e-9):
        msg = (
            f"{reaction.name}: needs a wavelength grid of at most "
            f"{reaction.max_grid_spacing_nm} nm in {lo}-{hi} nm, got up to "
            f"{np.diff(w).max():.4g} nm"
        )
        raise ValueError(msg)

    # Discrete-ordinates fluxes can be slightly negative where they are
    # vanishingly small; negative flux or cross section is numerical noise.
    actinic = np.maximum(flux["actinic_flux"].to_numpy()[selected, :], 0.0)
    sigma = np.maximum(
        flux["cross_section"].sel(species=reaction.absorber).to_numpy()[selected, :],
        0.0,
    )
    phi = _quantum_yield(reaction, w, flux["temperature_k"].to_numpy())
    return np.trapezoid(actinic * sigma * phi, w, axis=0)


def _line_rate(flux: xr.Dataset, reaction: LinePhotolysis) -> np.ndarray:
    wavelength = flux["wavelength"].to_numpy()
    if not wavelength[0] <= reaction.line_center_nm <= wavelength[-1]:
        msg = (
            f"{reaction.name}: line at {reaction.line_center_nm} nm is outside the "
            f"wavelength grid ({wavelength[0]}-{wavelength[-1]} nm)"
        )
        raise ValueError(msg)
    transmission = np.maximum(
        (
            flux["actinic_flux"].interp(wavelength=reaction.line_center_nm)
            / flux["solar_flux"].interp(wavelength=reaction.line_center_nm)
        ).to_numpy(),
        0.0,
    )
    return (
        reaction.quantum_yield
        * reaction.cross_section_m2
        * reaction.toa_line_flux_photons_m2_s
        * transmission
    )


def photolysis_rates(
    flux: xr.Dataset, reactions: Iterable[Photolysis | LinePhotolysis]
) -> xr.Dataset:
    """Per-molecule rates [s^-1] on the flux dataset's altitude grid.

    ``flux`` is the output of :meth:`ActinicFlux.calculate`. Each reaction
    gives one variable, named by the reaction's ``name``.
    """
    rates = xr.Dataset(coords={"altitude": flux["altitude"]})
    for reaction in reactions:
        if reaction.name in rates:
            msg = f"Rate '{reaction.name}' is defined twice"
            raise ValueError(msg)
        if isinstance(reaction, LinePhotolysis):
            values = _line_rate(flux, reaction)
        else:
            values = _continuum_rate(flux, reaction)
        rates[reaction.name] = ("altitude", values)
        rates[reaction.name].attrs = {"units": "s^-1", "reference": reaction.reference}
    return rates
