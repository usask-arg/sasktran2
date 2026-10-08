"""Actinic flux at the resolution of TUV, with its O2 band parameterisations.

:class:`TUVActinicFlux` reproduces the spectral and vertical treatment of
TUV-x in its v5.4 configuration, with SASKTRAN2's discrete-ordinates solver
in place of TUV's two-stream solver:

- the 156 wavelength bins of TUV-x v5.4 (120-735 nm);
- homogeneous layers between the altitude levels, with layer columns
  integrated assuming exponential variation;
- O2 in the Lyman-alpha bin from Chabrillat and Kockarts (1997), and in the
  17 Schumann-Runge band bins from Koppers and Murtagh (1996). Both give
  effective cross sections that depend on the slant O2 column, converted to
  layer optical depths as TUV does;
- rates are sums over bins of flux x cross section x quantum yield.

Cross sections and the solar spectrum elsewhere come either from
SASKTRAN2's own optical properties and TSIS-1 HSRS, averaged over each bin
(``data="sasktran2"``), or from TUV-x v5.4 itself (``data="tuv-x"``), which
with ``num_streams=2`` is the closest replica of TUV-x.

The TUV-x tables, ``photolysis/tuvx_v54.nc`` in the sasktran2 standard
database, are built by ``tools/nlte/build_tuvx_v54.py``. TUV-x and its data
are Copyright (C) 2020 National Center for Atmospheric Research, Apache
License 2.0.
"""

from __future__ import annotations

import functools
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import xarray as xr
from numpy.polynomial import chebyshev
from scipy import constants
from scipy.integrate import cumulative_trapezoid

import sasktran2 as sk
from sasktran2.constants import K_BOLTZMANN
from sasktran2.database.web import StandardDatabase
from sasktran2.optical.rayleigh import rayleigh_cross_section_bates

from .flux import EARTH_RADIUS_M, slant_columns
from .rates import _CK_B, _CK_C_CM2, _CK_D_CM2, _CK_E_CM2

TUVX_V54_DATABASE_KEY = "photolysis/tuvx_v54.nc"
HSRS_DATABASE_KEY = "solar/solar_irradiance_hsrs_2022_11_30_extended.nc"

#: Absorbers used with ``data="sasktran2"``, keyed by species id. O2 is the
#: ultraviolet continuum only; its Lyman-alpha and Schumann-Runge band bins
#: are parameterised, and TUV includes no O2 absorption beyond 243 nm.
TUV_OPTICAL_PROPERTIES = {
    "O3": sk.optical.O3DBM,
    "O2": sk.optical.O2UV,
    "NO2": sk.optical.NO2Vandaele,
}

_CM2 = 1.0e-4
# TUV's guards in the Lyman-alpha parameterisation.
_LA_TINY = 1.0e-100
_LA_MIN_CROSS_SECTION_M2 = 1.0e-20 * _CM2
_LA_LARGE_OPTICAL_DEPTH = 1000.0
# Relative column change below which TUV treats two levels as one.
_SRB_PRECISION = 1.0e-7


@functools.cache
def tuvx_v54_tables() -> xr.Dataset:
    """The TUV-x v5.4 tables from the sasktran2 standard database.

    Raises
    ------
    OSError
        If the tables cannot be found
    """
    path = StandardDatabase().path(TUVX_V54_DATABASE_KEY)
    if path is None or not path.exists():
        msg = f"Could not find the TUV-x v5.4 tables at {path}"
        raise OSError(msg)
    return xr.load_dataset(path)


def tuvx_v54_wavelength_edges() -> np.ndarray:
    """Edges [nm] of the 156 TUV-x v5.4 wavelength bins."""
    return tuvx_v54_tables()["wavelength_edge"].to_numpy()


def chabrillat_kockarts_cross_section(o2_column_m2) -> np.ndarray:
    """O2 effective Lyman-alpha cross section [m^2] at slant O2 columns [m^-2].

    ``R(N) sigma(N) / R(N)`` of Chabrillat and Kockarts (1997), with TUV's
    floor of 1e-20 cm^2 where the line is extinguished.
    """
    n_cm2 = np.asarray(o2_column_m2, dtype=float) * _CM2
    reduction = _ck_sum(_CK_B, _CK_C_CM2, n_cm2)
    rate = _ck_sum(_CK_D_CM2, _CK_E_CM2, n_cm2)
    valid = (reduction > _LA_TINY) & (rate > _LA_TINY)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(valid, rate / reduction * _CM2, _LA_MIN_CROSS_SECTION_M2)


def _ck_sum(weights, coefficients, n_cm2) -> np.ndarray:
    with np.errstate(over="ignore"):
        return sum(
            w * np.exp(-c * n_cm2) for w, c in zip(weights, coefficients, strict=True)
        )


def koppers_murtagh_cross_section(
    o2_column_m2, temperature_k, tables: xr.Dataset | None = None
) -> np.ndarray:
    """O2 effective cross sections [m^2] in the 17 Schumann-Runge band bins.

    Koppers and Murtagh (1996) Chebyshev fits in the log slant O2 column and
    temperature, with TUV's handling of columns outside the fitted range.
    Columns and temperatures are on altitude levels, ordered upward; the
    result is (level, band).
    """
    tables = tuvx_v54_tables() if tables is None else tables
    n_cm2 = np.atleast_1d(np.asarray(o2_column_m2, dtype=float)) * _CM2
    temperature_k = np.broadcast_to(temperature_k, n_cm2.shape)
    lower, upper = (float(v) for v in tables.attrs["srb_log_column_limits"])
    t0 = float(tables.attrs["srb_reference_temperature_k"])
    a = tables["srb_chebyshev_a"].to_numpy().copy()
    b = tables["srb_chebyshev_b"].to_numpy().copy()
    # TUV's Chebyshev series halves the first coefficient.
    a[0] *= 0.5
    b[0] *= 0.5

    num_levels = n_cm2.size
    sigma = np.tile(tables["srb_default_cross_section"].to_numpy(), (num_levels, 1))
    column = np.maximum(n_cm2, np.exp(lower))
    with np.errstate(divide="ignore"):
        x = np.log(column)
    small = n_cm2 < np.exp(lower)
    fitted = ~small & (x <= upper)
    y = (2.0 * x[fitted] - (lower + upper)) / (upper - lower)
    sigma[fitted] = (
        np.exp(
            chebyshev.chebval(y, a).T * (temperature_k[fitted, np.newaxis] - t0)
            + chebyshev.chebval(y, b).T
        )
        * _CM2
    )
    # Below the fitted range (from the bottom up), take the first fitted
    # level; above it, the last level before the column becomes too small.
    deep = np.flatnonzero(~small & (x > upper))
    num_deep = min(deep[-1] + 1, num_levels - 1) if deep.size else 0
    if num_deep:
        sigma[:num_deep] = sigma[num_deep]
    if small.any():
        top = np.flatnonzero(small)[0] - 1
        if top >= 0:
            sigma[top + 1 :] = sigma[top]
    return sigma


def _layer_columns(altitude_m, density_m3) -> np.ndarray:
    """Columns [m^-2] of each layer, assuming exponential variation within it."""
    dz = np.diff(altitude_m)
    lo, hi = density_m3[..., :-1], density_m3[..., 1:]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log(lo / hi)
        return np.where(
            (lo > 0) & (hi > 0) & (np.abs(ratio) > 1e-8),
            dz * (lo - hi) / ratio,
            dz * 0.5 * (lo + hi),
        )


def _secant(slant_air_m2, layer_air_m2) -> np.ndarray:
    """TUV's ratio of slant to vertical air column at each level."""
    vertical = np.concatenate([np.cumsum(layer_air_m2[::-1])[::-1], [0.0]])
    secant = np.empty_like(vertical)
    with np.errstate(divide="ignore", invalid="ignore"):
        secant[:-1] = slant_air_m2[:-1] / vertical[:-1]
    # Isotropic value where the sun is blocked; the top level copies the one
    # below.
    secant[:-1] = np.where(np.isfinite(slant_air_m2[:-1]), secant[:-1], 2.0)
    secant[-1] = secant[-2]
    return secant


def _lyman_alpha_optical_depth(o2_column_m2, secant) -> np.ndarray:
    """Vertical O2 optical depth of each layer in the Lyman-alpha bin."""
    reduction = _ck_sum(_CK_B, _CK_C_CM2, np.asarray(o2_column_m2) * _CM2)
    lower, upper = reduction[:-1], reduction[1:]
    with np.errstate(divide="ignore", invalid="ignore"):
        tau = np.log(upper) / secant[1:] - np.log(lower) / secant[:-1]
    return np.where((lower > _LA_TINY) & (upper > 0.0), tau, _LA_LARGE_OPTICAL_DEPTH)


def _schumann_runge_optical_depth(o2_column_m2, sigma_m2, secant, tables) -> np.ndarray:
    """Vertical O2 optical depth (layer, band) in the Schumann-Runge band bins.

    The slant optical depth of a layer integrates sigma(N) dN with sigma a
    power law in N across the layer, as in TUV.
    """
    lower_limit = float(tables.attrs["srb_log_column_limits"][0])
    column = np.maximum(np.asarray(o2_column_m2) * _CM2, np.exp(lower_limit))[
        :, np.newaxis
    ]
    sigma = sigma_m2 / _CM2
    absorbed = sigma * column
    num_levels = column.shape[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        power = np.log(sigma[1:] / sigma[:-1]) / np.log(column[1:] / column[:-1])
        tau = np.abs(absorbed[1:] - absorbed[:-1]) / (1.0 + power)
    tau = 2.0 * tau / (secant[:-1, np.newaxis] + secant[1:, np.newaxis])
    # TUV's value where the column barely changes, not divided by the secant.
    same = np.abs(1.0 - column[1:] / column[:-1]) <= 2.0 * _SRB_PRECISION
    tau = np.where(same, absorbed[1:] / (num_levels - 1), tau)
    # Near the lower end of the fitted columns sigma can fall faster than 1/N,
    # where TUV's layer optical depth turns slightly negative (~1e-6).
    return np.maximum(tau, 0.0)


def _rayleigh_second_moment(wavelength_nm) -> np.ndarray:
    """Second Legendre moment of the Rayleigh phase function, from the King factor."""
    _, king = rayleigh_cross_section_bates(np.asarray(wavelength_nm) / 1.0e3)
    depolarization = 6.0 * (king - 1.0) / (3.0 + 7.0 * king)
    gamma = depolarization / (2.0 - depolarization)
    return (1.0 - gamma) / (2.0 * (1.0 + 2.0 * gamma))


@functools.cache
def _hsrs_photon_cumulative() -> tuple[np.ndarray, np.ndarray]:
    ds = xr.load_dataset(StandardDatabase().path(HSRS_DATABASE_KEY))
    wavelength = ds["wavelength"].to_numpy()
    photons = ds["irradiance"].to_numpy() / (
        constants.h * constants.c / (wavelength * 1.0e-9)
    )
    return wavelength, cumulative_trapezoid(photons, wavelength, initial=0.0)


def hsrs_photons_per_bin(edges_nm) -> np.ndarray:
    """TSIS-1 HSRS photon flux [photons m^-2 s^-1] integrated over each bin, at 1 AU."""
    wavelength, cumulative = _hsrs_photon_cumulative()
    return np.diff(np.interp(edges_nm, wavelength, cumulative))


@dataclass(frozen=True)
class TUVXQuantumYield:
    """A TUV-x v5.4 quantum yield, for :class:`Photolysis` on TUV-x bins.

    Looks up ``table`` (``"o3_o1d"`` or ``"o3_o3p"``) by bin and interpolates
    linearly in temperature, as TUV-x does. Only valid at bin centres.
    """

    table: str

    def __call__(self, wavelength_nm, temperature_k) -> np.ndarray:
        tables = tuvx_v54_tables()
        values = tables[f"{self.table}_quantum_yield"].to_numpy()
        t_grid = tables["temperature_k"].to_numpy()
        edges = tables["wavelength_edge"].to_numpy()
        wavelength_nm, temperature_k = np.broadcast_arrays(wavelength_nm, temperature_k)
        index = np.clip(
            np.searchsorted(edges, wavelength_nm, side="right") - 1,
            0,
            edges.size - 2,
        )
        t = np.clip(temperature_k, t_grid[0], t_grid[-1])
        j = np.clip(np.searchsorted(t_grid, t, side="right") - 1, 0, t_grid.size - 2)
        weight = (t - t_grid[j]) / (t_grid[j + 1] - t_grid[j])
        return (1.0 - weight) * values[j, index] + weight * values[j + 1, index]


class TUVActinicFlux:
    """Actinic flux on the TUV-x v5.4 wavelength bins, with TUV's O2 parameterisations.

    Parameters
    ----------
    altitudes_m
        Altitude levels [m]; the layers between them are homogeneous.
    data
        ``"sasktran2"`` for SASKTRAN2's cross sections (averaged over each
        bin) and the TSIS-1 HSRS solar spectrum, or ``"tuv-x"`` for TUV-x
        v5.4's O2, O3 and Rayleigh cross sections, extraterrestrial flux and
        O3 quantum yields (see :class:`TUVXQuantumYield`). Absorbers other
        than O2 and O3 always use SASKTRAN2's cross sections.
    num_streams
        Discrete-ordinates streams. TUV-x v5.4 uses a two-stream
        delta-Eddington solver; ``num_streams=2`` is the closest match.
    optical_properties
        SASKTRAN2 cross sections by species id, replacing entries of
        :data:`TUV_OPTICAL_PROPERTIES`.
    cross_section_resolution_nm
        Sampling used to average SASKTRAN2 cross sections over each bin.
    num_threads
        Engine threads.
    """

    def __init__(
        self,
        altitudes_m,
        data: str = "sasktran2",
        num_streams: int = 4,
        optical_properties: Mapping | None = None,
        cross_section_resolution_nm: float = 0.05,
        num_threads: int = 8,
    ):
        if data not in ("sasktran2", "tuv-x"):
            msg = f"data must be 'sasktran2' or 'tuv-x', not {data!r}"
            raise ValueError(msg)
        self.altitudes_m = np.asarray(altitudes_m, dtype=float)
        self.data = data
        self.num_streams = num_streams
        self.optical_properties = dict(optical_properties or {})
        self.cross_section_resolution_nm = cross_section_resolution_nm
        self.num_threads = num_threads

    def _optical_property(self, species: str):
        if species not in self.optical_properties:
            self.optical_properties[species] = TUV_OPTICAL_PROPERTIES[species]()
        return self.optical_properties[species]

    def _binned_cross_sections(self, species, edges, altitude, temperature, pressure):
        """Bin-averaged SASKTRAN2 cross sections [m^2], (bin, altitude)."""
        widths = np.diff(edges)
        counts = np.maximum(
            np.ceil(widths / self.cross_section_resolution_nm).astype(int), 1
        )
        fine = np.concatenate(
            [
                lo + width * (np.arange(n) + 0.5) / n
                for lo, width, n in zip(edges[:-1], widths, counts, strict=True)
            ]
        )
        config = sk.Config()
        geometry = sk.Geometry1D(
            1.0,
            0.0,
            EARTH_RADIUS_M,
            altitude,
            sk.InterpolationMethod.LinearInterpolation,
            sk.GeometryType.PlaneParallel,
        )
        atmo = sk.Atmosphere(geometry, config, fine, calculate_derivatives=False)
        atmo.temperature_k = temperature
        atmo.pressure_pa = pressure
        xs = self._optical_property(species).atmosphere_quantities(atmo).extinction
        starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
        return np.add.reduceat(xs, starts, axis=1).T / counts[:, np.newaxis]

    def calculate(
        self,
        atmosphere: xr.Dataset,
        cos_sza: float,
        albedo: float = 0.0,
        earth_sun_distance_au: float = 1.0,
    ) -> xr.Dataset:
        """Actinic flux for one solar zenith angle.

        Takes the same atmosphere as :meth:`ActinicFlux.calculate` and
        returns the same variables, on bin centres with a ``wavelength_edge``
        coordinate; ``actinic_flux`` and ``solar_flux`` are bin averages
        [photons m^-2 s^-1 nm^-1]. In the parameterised O2 bins,
        ``cross_section`` is the effective cross section at each level.
        """
        tables = tuvx_v54_tables()
        edges = tables["wavelength_edge"].to_numpy()
        centres = tables["wavelength"].to_numpy()
        widths = np.diff(edges)
        lyman_alpha = tables["lyman_alpha_bin"].to_numpy().astype(bool)
        srb = tables["schumann_runge_bin"].to_numpy().astype(bool)

        z = self.altitudes_m
        z_mid = 0.5 * (z[1:] + z[:-1])
        source_altitude = atmosphere["altitude"].to_numpy()

        def levels(name, at, log=False):
            values = atmosphere[name].to_numpy()
            if log:
                return np.exp(
                    np.interp(at, source_altitude, np.log(np.maximum(values, 1e-300)))
                )
            return np.interp(at, source_altitude, values)

        temperature = levels("temperature_k", z)
        temperature_mid = levels("temperature_k", z_mid)
        pressure = levels("pressure_pa", z, log=True)
        pressure_mid = levels("pressure_pa", z_mid, log=True)
        air = pressure / (K_BOLTZMANN * temperature)
        air_layer = _layer_columns(z, air)

        known = [
            *TUV_OPTICAL_PROPERTIES,
            *(s for s in self.optical_properties if s not in TUV_OPTICAL_PROPERTIES),
        ]
        absorbers = [s for s in known if s in atmosphere]
        densities = {s: levels(s, z, log=True) for s in absorbers}
        layer_columns = {s: _layer_columns(z, densities[s]) for s in absorbers}

        slant = slant_columns(
            z, np.array([air, densities.get("O2", np.zeros_like(z))]), cos_sza
        )
        secant = _secant(slant[0], air_layer)

        # Cross sections at levels (for rates) and layers (for optical depth).
        combined = np.empty(2 * z.size - 1)
        combined[0::2], combined[1::2] = z, z_mid
        combined_t = np.empty_like(combined)
        combined_t[0::2], combined_t[1::2] = temperature, temperature_mid
        combined_p = np.empty_like(combined)
        combined_p[0::2], combined_p[1::2] = pressure, pressure_mid
        sigma_level, sigma_layer = {}, {}
        for species in absorbers:
            if self.data == "tuv-x" and species == "O3":
                t_grid = tables["temperature_k"].to_numpy()
                table = tables["o3_cross_section"].to_numpy()
                both = np.array(
                    [
                        np.interp(
                            np.clip(combined_t, t_grid[0], t_grid[-1]), t_grid, col
                        )
                        for col in table.T
                    ]
                )
            elif self.data == "tuv-x" and species == "O2":
                values = np.nan_to_num(tables["o2_cross_section"].to_numpy())
                both = np.repeat(values[:, np.newaxis], combined.size, axis=1)
            else:
                both = self._binned_cross_sections(
                    species, edges, combined, combined_t, combined_p
                )
            sigma_level[species], sigma_layer[species] = both[:, 0::2], both[:, 1::2]

        tau = np.zeros((centres.size, z.size - 1))
        for species in absorbers:
            tau += sigma_layer[species] * layer_columns[species]
        if "O2" in absorbers:
            o2_slant = slant[1]
            km = koppers_murtagh_cross_section(o2_slant, temperature, tables)
            sigma_level["O2"][lyman_alpha] = chabrillat_kockarts_cross_section(o2_slant)
            sigma_level["O2"][srb] = km.T
            o2_tau = {
                "la": _lyman_alpha_optical_depth(o2_slant, secant),
                "srb": _schumann_runge_optical_depth(o2_slant, km, secant, tables).T,
            }
            parameterised = sigma_layer["O2"] * layer_columns["O2"]
            tau[lyman_alpha] += o2_tau["la"] - parameterised[lyman_alpha]
            tau[srb] += o2_tau["srb"] - parameterised[srb]

        if self.data == "tuv-x":
            rayleigh = tables["rayleigh_cross_section"].to_numpy()
            solar_per_bin = tables["solar_flux"].to_numpy()
        else:
            fine = np.linspace(0.0, 1.0, 21)[1:-1]
            rayleigh = np.mean(
                [
                    rayleigh_cross_section_bates((edges[:-1] + widths * f) / 1.0e3)[0]
                    for f in fine
                ],
                axis=0,
            )
            solar_per_bin = hsrs_photons_per_bin(edges)
        tau_rayleigh = rayleigh[:, np.newaxis] * air_layer[np.newaxis, :]
        total = tau + tau_rayleigh

        config = sk.Config()
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        config.flux_types = [sk.FluxType.Actinic]
        config.num_streams = self.num_streams
        config.num_forced_azimuth = 1
        config.num_threads = self.num_threads
        geometry = sk.Geometry1D(
            cos_sza,
            0.0,
            EARTH_RADIUS_M,
            z,
            sk.InterpolationMethod.LowerInterpolation,
            sk.GeometryType.PseudoSpherical,
        )
        viewing = sk.ViewingGeometry()
        for altitude in z:
            viewing.add_flux_observer(sk.FluxObserverSolar(cos_sza, altitude))
        atmo = sk.Atmosphere(geometry, config, centres, calculate_derivatives=False)
        atmo.temperature_k = temperature
        atmo.pressure_pa = pressure

        # Each level holds the layer above it; the top level is never used.
        extinction = np.empty((z.size, centres.size))
        extinction[:-1] = (total / np.diff(z)[np.newaxis, :]).T
        extinction[-1] = extinction[-2]
        with np.errstate(divide="ignore", invalid="ignore"):
            ssa = np.where(total > 0, tau_rayleigh / total, 0.0).T
        ssa = np.vstack([ssa, ssa[-1:]])
        moments = np.zeros((atmo.storage.leg_coeff.shape[0], *extinction.shape))
        moments[0] = 1.0
        moments[2] = _rayleigh_second_moment(centres)[np.newaxis, :]
        atmo["tuv"] = sk.constituent.Manual(extinction, ssa, moments)
        atmo["surface"] = sk.constituent.LambertianSurface(albedo)

        engine = sk.Engine(config, geometry, viewing)
        radiance = engine.calculate_radiance(atmo)
        solar = solar_per_bin / widths / earth_sun_distance_au**2

        ds = xr.Dataset(
            {
                "actinic_flux": (
                    ("wavelength", "altitude"),
                    radiance["actinic_flux"].to_numpy() * solar[:, np.newaxis],
                ),
                "solar_flux": (("wavelength",), solar),
                "cross_section": (
                    ("species", "wavelength", "altitude"),
                    (
                        np.stack([sigma_level[s] for s in absorbers])
                        if absorbers
                        else np.zeros((0, centres.size, z.size))
                    ),
                ),
                "slant_column": (
                    ("species", "altitude"),
                    (
                        slant_columns(
                            z, np.array([densities[s] for s in absorbers]), cos_sza
                        )
                        if absorbers
                        else np.zeros((0, z.size))
                    ),
                ),
                "temperature_k": (("altitude",), temperature),
            },
            coords={
                "altitude": z,
                "wavelength": centres,
                "wavelength_edge": ("wavelength_edge", edges),
                "species": absorbers,
            },
            attrs={
                "cos_sza": cos_sza,
                "albedo": albedo,
                "earth_sun_distance_au": earth_sun_distance_au,
                "num_streams": self.num_streams,
                "spectral_sampling": "bins",
                "data": self.data,
            },
        )
        ds["actinic_flux"].attrs["units"] = "photons m^-2 s^-1 nm^-1"
        ds["solar_flux"].attrs["units"] = "photons m^-2 s^-1 nm^-1"
        ds["cross_section"].attrs["units"] = "m^2"
        ds["slant_column"].attrs["units"] = "m^-2"
        ds["wavelength_edge"].attrs["units"] = "nm"
        return ds
