from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.constants import c, h, k, physical_constants
from scipy.special import wofz

from sasktran2.atmosphere import NativeGridDerivative
from sasktran2.database.metals import MetalSpectroscopyDatabase
from sasktran2.optical.base import OpticalProperty
from sasktran2.polarization import LegendreStorageView

_AMU = physical_constants["atomic mass constant"][0]
_INTEGRATED_CROSS_SECTION = (
    np.pi * physical_constants["classical electron radius"][0] * c
)
_C2 = h * c * 100 / k  # kelvin cm


def resonance_polarizability(lower_j, upper_j) -> np.ndarray:
    """Electronic/rotational W2 for an isolated electric-dipole return transition.

    This is ``3 (2 Ju + 1) {1 1 2; Ju Ju Jl}**2``. It assumes an
    unpolarized lower level, no magnetic field, and no hyperfine or J-state
    interference. Nuclear hyperfine quantum numbers must not be substituted
    without also supplying their component strengths and frequencies.
    """
    lower, upper = np.broadcast_arrays(
        np.asarray(lower_j, dtype=float), np.asarray(upper_j, dtype=float)
    )
    delta = upper - lower
    valid = (
        np.isfinite(lower)
        & np.isfinite(upper)
        & (lower >= 0)
        & (upper >= 0)
        & (lower * 2 == np.round(lower * 2))
        & (upper * 2 == np.round(upper * 2))
        & np.isin(delta, [-1, 0, 1])
        & ((lower + upper) > 0)
    )
    if not np.all(valid):
        msg = "J values must be nonnegative half-integers for an allowed E1 transition"
        raise ValueError(msg)
    result = np.zeros_like(lower)
    plus = delta == 1
    j = lower[plus]
    result[plus] = (2 * j + 5) * (j + 2) / (10 * (j + 1) * (2 * j + 1))
    same = delta == 0
    j = lower[same]
    result[same] = (2 * j - 1) * (2 * j + 3) / (10 * j * (j + 1))
    minus = delta == -1
    j = lower[minus]
    result[minus] = (2 * j - 3) * (j - 1) / (10 * j * (2 * j + 1))
    return result


def resonance_phase_function(cos_angle, polarizability):
    """Scalar phase function with unit spherical mean: ``1 + W2/2 P2(mu)``."""
    mu = np.asarray(cos_angle, dtype=float)
    w2 = np.asarray(polarizability, dtype=float)
    if not np.all(np.isfinite(mu) & (np.abs(mu) <= 1)):
        msg = "cos_angle must be finite and between -1 and 1"
        raise ValueError(msg)
    if not np.all(np.isfinite(w2) & (w2 >= 0) & (w2 <= 1)):
        msg = "polarizability must be finite and between 0 and 1"
        raise ValueError(msg)
    return 1 + w2 * (3 * mu**2 - 1) / 4


@dataclass
class ResonanceQuantities:
    """Cross sections in m² and phase coefficients.

    ``cross_sections`` returns dimensionless albedo in ``ssa``. Following the
    engine's optical-property interface, ``atmosphere_quantities`` instead puts
    the scattering cross section in that field.
    """

    extinction: np.ndarray
    ssa: np.ndarray
    leg_coeff: np.ndarray
    polarizability: np.ndarray


def _phase_coefficients(w2, num_stokes, num_moments, *, derivative=False):
    if num_stokes not in (1, 3):
        msg = "Resonance phase matrices support num_stokes=1 or 3"
        raise ValueError(msg)
    if not isinstance(num_moments, int | np.integer) or num_moments < 3:
        msg = "At least three moments are required for resonance scattering"
        raise ValueError(msg)
    result = np.zeros(((1 if num_stokes == 1 else 4) * num_moments, *w2.shape))
    view = LegendreStorageView(result, num_stokes)
    if not derivative:
        view.a1[0] = 1
    view.a1[2] = w2 / 2
    if num_stokes == 3:
        view.a2[2] = 3 * w2
        view.b1[2] = np.sqrt(6) * w2 / 2
    return result


class LineResonance(OpticalProperty):
    """LTE line extinction and same-transition resonance scattering.

    Vacuum wavelengths are used throughout. Each Voigt profile is normalized
    in frequency, with Doppler standard deviation ``nu/c sqrt(k T / mass)``
    and natural Lorentz HWHM ``(Gamma_upper + Gamma_lower)/(4 pi)``.
    The integrated absorption cross section is ``pi r_e c f`` times the
    LTE population of the lower state. No pressure broadening, collisions,
    stimulated emission, hyperfine/isotope components, magnetic effects, or
    frequency redistribution is added to the source data.

    The scattering albedo includes only radiative return to the *same* lower
    state, ``A_ul / Gamma_upper``. Other radiative branches are extinction
    losses in this monochromatic approximation; fluorescence feeding other
    wavelengths requires a separate source calculation. Incomplete NIST
    decay sums make the estimated elastic return an upper bound.

    Parameters
    ----------
    db_filepath, db
        Exactly one NetCDF file or in-memory line dataset. See the metal
        resonance user guide for the schema and source assumptions.
    line_wing_cutoff_nm
        Symmetric hard truncation about each vacuum line center, default
        0.1 nm, without renormalization. Set to None for complete Voigt wings.
        This is a computational cutoff, not instrumental broadening.
    """

    def __init__(
        self,
        db_filepath: str | Path | None = None,
        *,
        db: xr.Dataset | None = None,
        line_wing_cutoff_nm: float | None = 0.1,
    ):
        if (db_filepath is None) == (db is None):
            msg = "Provide exactly one of db_filepath or db"
            raise ValueError(msg)
        if db is None:
            with xr.open_dataset(db_filepath) as source:
                db = source.load()
        self._database = db.copy(deep=True)
        self._cutoff = line_wing_cutoff_nm
        if self._cutoff is not None and (
            not np.isfinite(self._cutoff) or self._cutoff <= 0
        ):
            msg = "line_wing_cutoff_nm must be positive or None"
            raise ValueError(msg)
        self._validate_database()

    @property
    def database(self) -> xr.Dataset:
        """A copy of the source lines, populations, and provenance metadata."""
        return self._database.copy(deep=True)

    def _validate_database(self):
        db = self._database
        required = (
            "wavelength_nm",
            "oscillator_strength",
            "lower_energy_cminv",
            "lower_j",
            "upper_j",
            "einstein_a_s",
            "upper_total_a_s",
        )
        for name in required:
            if name not in db or db[name].dims != ("line",):
                msg = f"{name} must be present with dimension 'line'"
                raise ValueError(msg)
            values = db[name].values
            if not np.all(np.isfinite(values)) or np.any(values < 0):
                msg = f"{name} must contain finite nonnegative values"
                raise ValueError(msg)
        if db.sizes["line"] == 0:
            msg = "A resonance database must contain at least one line"
            raise ValueError(msg)
        for name in (
            "wavelength_nm",
            "oscillator_strength",
            "einstein_a_s",
            "upper_total_a_s",
        ):
            if np.any(db[name].values <= 0):
                msg = f"{name} must be strictly positive"
                raise ValueError(msg)
        self._mass = float(db.attrs.get("mass_amu", np.nan))
        if not np.isfinite(self._mass) or self._mass <= 0:
            msg = "The database requires a positive mass_amu attribute"
            raise ValueError(msg)
        if (
            db.attrs.get("wavelength_medium", db.wavelength_nm.attrs.get("medium"))
            != "vacuum"
        ):
            msg = "Resonance line wavelengths must be in vacuum"
            raise ValueError(msg)
        if np.any(db.einstein_a_s.values > db.upper_total_a_s.values * (1 + 1e-10)):
            msg = "upper_total_a_s cannot be smaller than einstein_a_s"
            raise ValueError(msg)
        for name in ("lower_total_a_s", "lower_statistical_weight"):
            if name in db:
                values = db[name].values
                if (
                    db[name].dims != ("line",)
                    or not np.all(np.isfinite(values))
                    or np.any(values < 0)
                ):
                    msg = f"{name} must be a finite nonnegative line array"
                    raise ValueError(msg)
        self._wavelength = db.wavelength_nm.values.astype(float)
        self._frequency = c * 1e9 / self._wavelength
        self._f = db.oscillator_strength.values.astype(float)
        self._energy = db.lower_energy_cminv.values.astype(float)
        self._weight = (
            db.lower_statistical_weight.values.astype(float)
            if "lower_statistical_weight" in db
            else 2 * db.lower_j.values + 1
        )
        if np.any(self._weight <= 0):
            msg = "Lower statistical weights must be positive"
            raise ValueError(msg)
        self._w2 = resonance_polarizability(db.lower_j.values, db.upper_j.values)
        lower_decay = db.lower_total_a_s.values if "lower_total_a_s" in db else 0
        self._gamma = (db.upper_total_a_s.values + lower_decay) / (4 * np.pi)
        self._branch = np.minimum(db.einstein_a_s.values / db.upper_total_a_s.values, 1)
        self._has_states = "energy_cminv" in db and "statistical_weight" in db
        if self._has_states:
            energies = db.energy_cminv.values
            weights = db.statistical_weight.values
            if (
                energies.ndim != 1
                or weights.shape != energies.shape
                or not energies.size
                or not np.all(np.isfinite(energies))
                or np.any(energies < 0)
                or not np.all(np.isfinite(weights))
                or np.any(weights <= 0)
                or not np.any(energies == 0)
            ):
                msg = "Partition states require nonnegative energies including zero and positive weights"
                raise ValueError(msg)
        else:
            if "partition_temperature_k" not in db or "partition_function" not in db:
                msg = "Provide partition states or a partition_temperature_k/partition_function table"
                raise ValueError(msg)
            grid = db.partition_temperature_k.values
            values = db.partition_function.values
            if (
                grid.ndim != 1
                or values.shape != grid.shape
                or len(grid) < 2
                or not np.all(np.isfinite(grid))
                or not np.all(np.isfinite(values))
                or np.any(grid <= 0)
                or np.any(values <= 0)
                or np.any(np.diff(grid) <= 0)
            ):
                msg = "Partition table must have increasing positive temperatures and positive values"
                raise ValueError(msg)
        if "upper_decay_data_complete" in db and not np.all(
            db.upper_decay_data_complete.values
        ):
            warnings.warn(
                "Some upper-state decay sums are incomplete: natural widths are lower "
                "bounds and elastic-return albedos are upper bounds. See dataset provenance.",
                UserWarning,
                stacklevel=3,
            )

    def _partition(self, temperature):
        db = self._database
        if self._has_states:
            q = np.zeros_like(temperature)
            dq = np.zeros_like(temperature)
            energies = db.energy_cminv.values
            weights = db.statistical_weight.values
            for start in range(0, len(energies), 4096):
                e = _C2 * energies[start : start + 4096, np.newaxis]
                terms = weights[start : start + 4096, np.newaxis] * np.exp(
                    -e / temperature
                )
                q += terms.sum(axis=0)
                dq += (terms * e / temperature**2).sum(axis=0)
            return q, dq / q
        grid = db.partition_temperature_k.values
        values = db.partition_function.values
        if np.any((temperature < grid[0]) | (temperature > grid[-1])):
            msg = "Temperature is outside the supplied partition-function table"
            raise ValueError(msg)
        index = np.clip(
            np.searchsorted(grid, temperature, side="right") - 1, 0, len(grid) - 2
        )
        slope = np.log(values[index + 1] / values[index]) / np.log(
            grid[index + 1] / grid[index]
        )
        return values[index] * (temperature / grid[index]) ** slope, slope / temperature

    def _temperature(self, temperature_k, size):
        if temperature_k is None:
            msg = "temperature_k is required for metal resonance cross sections"
            raise ValueError(msg)
        temperature = np.asarray(temperature_k, dtype=float)
        if temperature.ndim > 1 or temperature.size not in (1, size):
            msg = "temperature_k must be scalar or have one value per altitude"
            raise ValueError(msg)
        temperature = np.broadcast_to(temperature.reshape(-1), (size,))
        if not np.all(np.isfinite(temperature) & (temperature > 0)):
            msg = "temperature_k must be finite and positive"
            raise ValueError(msg)
        for attr, comparison in (
            ("temperature_min_k", np.less),
            ("temperature_max_k", np.greater),
        ):
            if attr in self._database.attrs and np.any(
                comparison(temperature, self._database.attrs[attr])
            ):
                msg = f"Temperature violates the database {attr}={self._database.attrs[attr]} limit"
                raise ValueError(msg)
        return temperature

    def line_strengths(self, temperature_k):
        """Frequency-integrated LTE cross sections [m² Hz], shape (T, line)."""
        temperature = self._temperature(temperature_k, np.size(temperature_k))
        q, _ = self._partition(temperature)
        population = (
            self._weight
            * np.exp(-_C2 * self._energy / temperature[:, None])
            / q[:, None]
        )
        return _INTEGRATED_CROSS_SECTION * self._f * population

    def _evaluate(
        self,
        wavelengths_nm,
        altitudes_m,
        temperature_k,
        *,
        derivatives=False,
        native=False,
    ):
        wavelengths = np.atleast_1d(np.asarray(wavelengths_nm, dtype=float))
        altitudes = np.atleast_1d(np.asarray(altitudes_m, dtype=float))
        if (
            wavelengths.ndim != 1
            or not wavelengths.size
            or not np.all(np.isfinite(wavelengths) & (wavelengths > 0))
        ):
            msg = "wavelengths_nm must be a nonempty vector of positive vacuum wavelengths"
            raise ValueError(msg)
        for attr, comparison in (
            ("wavelength_min_nm", np.less),
            ("wavelength_max_nm", np.greater),
        ):
            if attr in self._database.attrs and np.any(
                comparison(wavelengths, self._database.attrs[attr])
            ):
                msg = (
                    f"Wavelength is outside database coverage ({attr}="
                    f"{self._database.attrs[attr]}). Missing spectroscopy "
                    "cannot be assumed to have zero cross section."
                )
                raise ValueError(msg)
        if (
            altitudes.ndim != 1
            or not altitudes.size
            or not np.all(np.isfinite(altitudes))
        ):
            msg = "altitudes_m must be a nonempty finite vector"
            raise ValueError(msg)
        temperature = self._temperature(temperature_k, len(altitudes))
        q, dlogq = self._partition(temperature)
        order = np.argsort(wavelengths)
        sorted_wavelengths = wavelengths[order]
        frequencies = c * 1e9 / sorted_wavelengths
        shape = (len(altitudes), len(wavelengths))
        total, scatter, weighted = (np.zeros(shape) for _ in range(3))
        if derivatives:
            dtotal, dscatter, dweighted = (np.zeros(shape) for _ in range(3))
        thermal = np.sqrt(k * temperature / (self._mass * _AMU)) / c
        if self._cutoff is None:
            left_indices = np.zeros(len(self._wavelength), dtype=int)
            right_indices = np.full(len(self._wavelength), len(wavelengths), dtype=int)
        else:
            left_indices = np.searchsorted(
                sorted_wavelengths, self._wavelength - self._cutoff
            )
            right_indices = np.searchsorted(
                sorted_wavelengths, self._wavelength + self._cutoff, side="right"
            )
        for i in np.flatnonzero(right_indices > left_indices):
            left, right = left_indices[i], right_indices[i]
            sigma = self._frequency[i] * thermal[:, None]
            z = (
                frequencies[None, left:right] - self._frequency[i] + 1j * self._gamma[i]
            ) / (np.sqrt(2) * sigma)
            w = wofz(z)
            profile = w.real / (np.sqrt(2 * np.pi) * sigma)
            population = (
                self._weight[i] * np.exp(-_C2 * self._energy[i] / temperature) / q
            )
            strength = _INTEGRATED_CROSS_SECTION * self._f[i] * population[:, None]
            cross = strength * profile
            total[:, left:right] += cross
            scatter[:, left:right] += cross * self._branch[i]
            weighted[:, left:right] += cross * self._branch[i] * self._w2[i]
            if derivatives:
                # w + z w' loses precision far from line center. Its asymptotic
                # expansion starts at z^-3 because the z^-1 terms cancel.
                combination = w + z * (-2 * z * w + 2j / np.sqrt(np.pi))
                far = np.abs(z) > 20
                if np.any(far):
                    invz = 1 / z[far]
                    power = invz**3
                    series = np.zeros_like(power)
                    coefficient = 1.0
                    for n in range(1, 8):
                        coefficient *= (2 * n - 1) / 2
                        series += -2 * n * coefficient * power
                        power *= invz**2
                    combination[far] = 1j / np.sqrt(np.pi) * series
                dprofile = -combination.real / (
                    2 * temperature[:, None] * np.sqrt(2 * np.pi) * sigma
                )
                dlog_population = _C2 * self._energy[i] / temperature**2 - dlogq
                dcross = strength * (dprofile + profile * dlog_population[:, None])
                dtotal[:, left:right] += dcross
                dscatter[:, left:right] += dcross * self._branch[i]
                dweighted[:, left:right] += dcross * self._branch[i] * self._w2[i]
        ssa = np.divide(scatter, total, out=np.zeros(shape), where=total > 0)
        w2 = np.divide(weighted, scatter, out=np.zeros(shape), where=scatter > 0)
        inverse = np.argsort(order)
        if not derivatives:
            return (
                total[:, inverse],
                (scatter if native else ssa)[:, inverse],
                w2[:, inverse],
            )
        dssa = np.divide(
            dscatter - ssa * dtotal, total, out=np.zeros(shape), where=total > 0
        )
        dw2 = np.divide(
            dweighted - w2 * dscatter, scatter, out=np.zeros(shape), where=scatter > 0
        )
        return (
            dtotal[:, inverse],
            (dscatter if native else dssa)[:, inverse],
            dw2[:, inverse],
        )

    def cross_sections(
        self,
        wavelengths_nm,
        altitudes_m,
        *,
        temperature_k=None,
        num_stokes=1,
        num_moments=3,
        **kwargs,
    ) -> ResonanceQuantities:
        extinction, ssa, w2 = self._evaluate(wavelengths_nm, altitudes_m, temperature_k)
        return ResonanceQuantities(
            extinction, ssa, _phase_coefficients(w2, num_stokes, num_moments), w2
        )

    def cross_section_derivatives(
        self, wavelengths_nm, altitudes_m, *, temperature_k=None, **kwargs
    ):
        extinction, _, _ = self._evaluate(
            wavelengths_nm, altitudes_m, temperature_k, derivatives=True
        )
        return {"temperature_k": extinction.flatten()}

    @staticmethod
    def _atmosphere_inputs(atmo, kwargs):
        temperature = kwargs.get("temperature_k", atmo._native_state("temperature_k"))
        return {
            "wavelengths_nm": atmo.wavelengths_nm,
            "altitudes_m": atmo._native_altitudes(),
            "temperature_k": temperature,
        }

    def atmosphere_quantities(self, atmo, **kwargs):
        """Engine quantities: ``ssa`` stores scattering cross section in m²."""
        extinction, scattering, w2 = self._evaluate(
            **self._atmosphere_inputs(atmo, kwargs), native=True
        )
        return ResonanceQuantities(
            extinction,
            scattering,
            _phase_coefficients(w2, atmo.nstokes, atmo.leg_coeff.a1.shape[0]),
            w2,
        )

    def optical_derivatives(self, atmo, **kwargs):
        extinction, ssa, w2 = self._evaluate(
            **self._atmosphere_inputs(atmo, kwargs), derivatives=True, native=True
        )
        return {
            "temperature_k": NativeGridDerivative(
                d_extinction=extinction,
                d_ssa=ssa,
                d_leg_coeff=_phase_coefficients(
                    w2, atmo.nstokes, atmo.leg_coeff.a1.shape[0], derivative=True
                ),
            )
        }


class AtomicResonance(LineResonance):
    """A neutral or ionic species from the cached NIST-derived line database.

    Use NIST species keys such as ``Na_I``, ``Mg_II``, or ``Fe_I``. Species
    number density is the total of the represented charge state, with LTE
    populations among the supplied atomic levels. No ionization equilibrium
    is imposed. Fine-structure lines are unresolved in isotope/hyperfine space.
    """

    _species = None
    _kind = "atomic"

    def __init__(self, species=None, db_filepath=None, *, db=None, **kwargs):
        species = species or self._species
        if db_filepath is None and db is None:
            if species is None:
                msg = "An atomic or molecular species key is required"
                raise ValueError(msg)
            db_filepath = MetalSpectroscopyDatabase().path(
                species.replace(" ", "_"), kind=self._kind
            )
        super().__init__(db_filepath, db=db, **kwargs)


class MolecularResonance(AtomicResonance):
    """ExoMol LTE absorption and elastic return for one molecular isotopologue.

    This property does not represent total solar-pumped molecular fluorescence
    or chemiluminescence. Each dataset specifies its isotopologue and license;
    no terrestrial isotope-abundance scaling is applied automatically.
    """

    _kind = "molecular"


class Sodium(AtomicResonance):
    """Na I resonance lines; unresolved isotope/hyperfine structure."""

    _species = "Na_I"


class Potassium(AtomicResonance):
    """K I resonance lines; unresolved isotope/hyperfine structure."""

    _species = "K_I"


class Lithium(AtomicResonance):
    """Li I resonance lines; unresolved isotope/hyperfine structure."""

    _species = "Li_I"


class Magnesium(AtomicResonance):
    """Mg I resonance lines."""

    _species = "Mg_I"


class MagnesiumIon(AtomicResonance):
    """Mg II resonance lines."""

    _species = "Mg_II"


class Calcium(AtomicResonance):
    """Ca I resonance lines."""

    _species = "Ca_I"


class CalciumIon(AtomicResonance):
    """Ca II resonance lines, with loss to other radiative branches."""

    _species = "Ca_II"


class Iron(AtomicResonance):
    """Fe I low-state resonance lines."""

    _species = "Fe_I"


class IronIon(AtomicResonance):
    """Fe II low-state resonance lines."""

    _species = "Fe_II"


class Aluminium(AtomicResonance):
    """Al I low-state resonance lines."""

    _species = "Al_I"


class Nickel(AtomicResonance):
    """Ni I low-state resonance lines."""

    _species = "Ni_I"


class Chromium(AtomicResonance):
    """Cr I low-state resonance lines."""

    _species = "Cr_I"


class Manganese(AtomicResonance):
    """Mn I low-state resonance lines."""

    _species = "Mn_I"


class Titanium(AtomicResonance):
    """Ti I low-state resonance lines."""

    _species = "Ti_I"


class Cobalt(AtomicResonance):
    """Co I low-state resonance lines."""

    _species = "Co_I"


class Copper(AtomicResonance):
    """Cu I low-state resonance lines."""

    _species = "Cu_I"


class Zinc(AtomicResonance):
    """Zn I intercombination line; strength reconciled against the NIST A value."""

    _species = "Zn_I"


class Silicon(AtomicResonance):
    """Si I low-state resonance lines; silicon is a meteoric metalloid."""

    _species = "Si_I"


class Rubidium(AtomicResonance):
    """Rb I lines, a speculative trace meteoric candidate."""

    _species = "Rb_I"


class Strontium(AtomicResonance):
    """Sr I lines, a speculative trace meteoric candidate."""

    _species = "Sr_I"


class Barium(AtomicResonance):
    """Ba I lines, a speculative trace meteoric candidate."""

    _species = "Ba_I"


class AluminiumOxide(MolecularResonance):
    """27Al16O LTE extinction and elastic return from ExoMol ATP."""

    _species = "AlO"


class MagnesiumOxide(MolecularResonance):
    """24Mg16O LTE extinction and elastic return from ExoMol LiTY."""

    _species = "MgO"


class CalciumOxide(MolecularResonance):
    """40Ca16O LTE extinction and elastic return from ExoMol VBATHY."""

    _species = "CaO"


class TitaniumOxide(MolecularResonance):
    """48Ti16O LTE extinction and elastic return from ExoMol TOTO."""

    _species = "TiO"
