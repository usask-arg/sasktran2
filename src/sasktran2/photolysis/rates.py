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
    the line-integrated photon flux ``F_line`` at 1 AU, scaled by the
    Earth-Sun distance and the actinic-flux transmission at the line centre.
    """

    name: str
    line_center_nm: float
    cross_section_m2: float
    toa_line_flux_photons_m2_s: float
    quantum_yield: float = 1.0
    reference: str = ""


# Chabrillat and Kockarts (1997), Table 1: the solar Lyman-alpha line reduction
# factor R(N) = sum b_i exp(-c_i N) and the O2 effective photolysis cross
# section R(N) sigma(N) = sum d_i exp(-e_i N), for slant O2 column N [cm^-2].
_CK_B = (6.8431e-01, 2.29841e-01, 8.65412e-02)
_CK_C_CM2 = (8.22114e-21, 1.77556e-20, 8.22112e-21)
_CK_D_CM2 = (6.0073e-21, 4.28569e-21, 1.28059e-20)
_CK_E_CM2 = (8.21666e-21, 1.63296e-20, 4.85121e-17)


@dataclass(frozen=True)
class LymanAlphaPhotolysis:
    """O2 photolysis by the solar Lyman-alpha line, per molecule.

    ``J(z) = phi * F_line * sum_i d_i exp(-e_i N(z))``, the parameterisation
    of Chabrillat and Kockarts (1997), with ``N`` the slant O2 column from the
    flux dataset and ``F_line`` the line-integrated photon flux at 1 AU,
    scaled by the Earth-Sun distance. The effective cross section falls with
    depth because the surviving photons are in the line wings, where O2
    absorbs weakly.

    Chabrillat, S., and G. Kockarts (1997), Simple parameterization of the
    absorption of the solar Lyman-alpha line, Geophys. Res. Lett., 24(21),
    2659-2662.
    """

    name: str
    toa_line_flux_photons_m2_s: float
    quantum_yield: float = 1.0
    reference: str = "Chabrillat and Kockarts (1997), Geophys. Res. Lett., 24, 2659"


def lyman_alpha_reduction_factor(o2_column_m2) -> np.ndarray:
    """Fraction of the solar Lyman-alpha line remaining after a slant O2 column [m^-2]."""
    n_cm2 = np.asarray(o2_column_m2, dtype=float) * 1.0e-4
    return sum(b * np.exp(-c * n_cm2) for b, c in zip(_CK_B, _CK_C_CM2, strict=True))


def lyman_alpha_o2_rate_per_photon_m2(o2_column_m2) -> np.ndarray:
    """O2 photolysis rate per unit line flux, ``R(N) sigma(N)`` [m^2], at slant column [m^-2]."""
    n_cm2 = np.asarray(o2_column_m2, dtype=float) * 1.0e-4
    return 1.0e-4 * sum(
        d * np.exp(-e * n_cm2) for d, e in zip(_CK_D_CM2, _CK_E_CM2, strict=True)
    )


def _distance_factor(flux: xr.Dataset) -> float:
    return 1.0 / float(flux.attrs.get("earth_sun_distance_au", 1.0)) ** 2


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
        * _distance_factor(flux)
        * transmission
    )


def _lyman_alpha_rate(flux: xr.Dataset, reaction: LymanAlphaPhotolysis) -> np.ndarray:
    if "O2" not in flux["species"]:
        msg = f"{reaction.name}: the flux dataset has no O2 slant column"
        raise ValueError(msg)
    column = flux["slant_column"].sel(species="O2").to_numpy()
    return (
        reaction.quantum_yield
        * reaction.toa_line_flux_photons_m2_s
        * _distance_factor(flux)
        * lyman_alpha_o2_rate_per_photon_m2(column)
    )


def photolysis_rates(
    flux: xr.Dataset,
    reactions: Iterable[Photolysis | LinePhotolysis | LymanAlphaPhotolysis],
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
        elif isinstance(reaction, LymanAlphaPhotolysis):
            values = _lyman_alpha_rate(flux, reaction)
        else:
            values = _continuum_rate(flux, reaction)
        rates[reaction.name] = ("altitude", values)
        rates[reaction.name].attrs = {"units": "s^-1", "reference": reaction.reference}
    return rates
