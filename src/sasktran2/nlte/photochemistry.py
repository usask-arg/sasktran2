"""Photochemically excited emitters for a radiative-transfer atmosphere.

:func:`add_photochemical_species` connects the pieces: the atmosphere's
state and absorbers drive :class:`sasktran2.photolysis.ActinicFlux`, the
resulting photolysis and photoexcitation rates drive a bundled kinetic
mechanism through :func:`solve`, and the excited-state populations become
emission constituents of the atmosphere.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
import xarray as xr

import sasktran2 as sk

from . import Mechanism, solve

#: Volume mixing ratios used for background species that are neither in
#: ``background`` nor absorbers of the atmosphere.
DEFAULT_VMR = {"N2": 0.7808, "CO2": 4.2e-4}

#: Alternative names accepted in ``background``.
BACKGROUND_ALIASES = {"O": "O(3P)"}

#: Absorbers passed to the actinic-flux calculation when available.
PHOTOLYSIS_ABSORBERS = ("O2", "O3", "NO2", "N2", "O(3P)")


@dataclass(frozen=True)
class _Emitter:
    description: str
    #: Bundled kinetic mechanism, for emitters driven by state populations.
    mechanism: str | None = None
    #: Band constituent factory by mechanism transition id.
    bands: dict[str, Callable[[np.ndarray, np.ndarray], object]] | None = None
    #: Background species for emitters driven directly by the actinic flux.
    absorber: str | None = None
    #: Atomic lines by mechanism transition id: (vacuum wavelength [nm],
    #: emitter molar mass [g mol^-1]).
    lines: dict[str, tuple[float, float]] | None = None


class _Combined(sk.constituent.base.Constituent):
    """Several emission constituents added as one."""

    def __init__(self, parts: dict[str, object]):
        self.parts = parts

    def add_to_atmosphere(self, atmo):
        for part in self.parts.values():
            part.add_to_atmosphere(atmo)

    def register_derivative(self, atmo, name: str):
        for key, part in self.parts.items():
            part.register_derivative(atmo, f"{name}_{key}")


def _o2_band(band: str):
    def make(altitude_m, photon_ver):
        return sk.constituent.O2BandEmissionRate(altitude_m, photon_ver, band=band)

    return make


#: Emitters :func:`add_photochemical_species` can add, by state id.
PHOTOCHEMICAL_SPECIES = {
    "O2(b)": _Emitter(
        mechanism="oxygen",
        bands={
            "o2b0_a_band": _o2_band("0-0"),
            "o2b1_x1": _o2_band("1-1"),
            "o2b1_b_band": _o2_band("1-0"),
        },
        description=(
            "O2 b-X emission: the A band (0-0, 762 nm) with the 1-1 hot band, "
            "and the B band (1-0, 688 nm)"
        ),
    ),
    "O(1S)": _Emitter(
        mechanism="oxygen_green",
        lines={"o1s_green_line": (557.8888, 15.999), "o1s_297": (297.3159, 15.999)},
        description=(
            "The O(1S) green line (557.7 nm) and 297.2 nm line, from O2 "
            "photodissociation, the Barth mechanism and a supplied ionospheric "
            "production"
        ),
    ),
    "OH(A)": _Emitter(
        absorber="OH",
        description=(
            "OH A-X solar resonance fluorescence (bands near 282-350 nm, mostly "
            "0-0 at 308 nm) from thermal OH(X, v=0); needs OH in `background`"
        ),
    ),
}

_NOT_YET_AVAILABLE = {
    "O2(a)": "the O2 a-X (1.27 um) band emission is not implemented yet",
}


def _oxygen_rates(flux: xr.Dataset) -> xr.Dataset:
    from sasktran2 import photolysis

    return photolysis.photolysis_rates(flux, photolysis.presets.oxygen_photolysis())


def _green_rates(flux: xr.Dataset) -> xr.Dataset:
    from sasktran2 import photolysis

    return photolysis.photolysis_rates(flux, photolysis.presets.green_line_photolysis())


_RATE_FUNCTIONS = {"oxygen": _oxygen_rates, "oxygen_green": _green_rates}
#: Mechanisms that need another's solution, and the states they take from it.
_DEPENDS_ON = {"oxygen_green": ("oxygen", ("O2(a)",))}


def _canonical(name: str) -> str:
    known = {key.lower().replace(" ", ""): key for key in PHOTOCHEMICAL_SPECIES}
    key = name.lower().replace(" ", "")
    if key in known:
        return known[key]
    pending = {k.lower(): (k, why) for k, why in _NOT_YET_AVAILABLE.items()}
    if key in pending:
        state, why = pending[key]
        msg = f"{state} cannot be added yet: {why}"
        raise NotImplementedError(msg)
    msg = (
        f"Unknown photochemical species {name!r}; available: "
        f"{', '.join(PHOTOCHEMICAL_SPECIES)}"
    )
    raise ValueError(msg)


def _background_profiles(background: xr.Dataset | None) -> dict[str, xr.DataArray]:
    if background is None:
        return {}
    if "altitude" not in background.coords:
        msg = "background must have an 'altitude' coordinate [m]"
        raise ValueError(msg)
    profiles = {}
    for name, values in background.data_vars.items():
        profiles[BACKGROUND_ALIASES.get(str(name), str(name))] = values
    return profiles


def add_photochemical_species(
    atmosphere: sk.Atmosphere,
    species: Sequence[str],
    *,
    cos_sza: float,
    background: xr.Dataset | None = None,
    albedo: float = 0.3,
    earth_sun_distance_au: float = 1.0,
    actinic_flux=None,
    rates: xr.Dataset | None = None,
    ionospheric_o1s_production=None,
) -> xr.Dataset:
    """Solve the photochemistry of ``species`` and add their emission to ``atmosphere``.

    The atmosphere must have its temperature and pressure set. Background
    number densities for the kinetics and the actinic flux are taken, in
    order, from ``background``, from the atmosphere's
    :class:`~sasktran2.constituent.VMRAltitudeAbsorber` constituents with the
    same name (e.g. ``atmosphere["O3"]``), and for N2 and CO2 from
    :data:`DEFAULT_VMR`. Atomic oxygen, ``"O(3P)"`` (or ``"O"``), always comes
    from ``background``.

    The photochemistry is solved on the atmosphere's altitude grid, or on the
    grid of ``actinic_flux`` if one is given; the emission constituents
    interpolate onto the model grid. Emission only reaches the radiance when
    the engine's ``config.emission_source`` is
    ``sk.EmissionSource.VolumeEmissionRate``.

    Parameters
    ----------
    atmosphere
        A one-dimensional :class:`sasktran2.Atmosphere`.
    species
        Emitting states, from :data:`PHOTOCHEMICAL_SPECIES` (case
        insensitive). Each is added as the constituent ``"<state> emission"``.
    cos_sza
        Cosine of the solar zenith angle for the photochemistry, normally at
        the tangent point.
    background
        Number densities [m^-3] by species id on an ``altitude`` coordinate
        [m]. Interpolated in log space.
    albedo
        Lambertian surface albedo for the actinic flux.
    earth_sun_distance_au
        Earth-Sun distance for the actinic flux.
    actinic_flux
        A :class:`sasktran2.photolysis.ActinicFlux` to use instead of the
        default one on the atmosphere's altitude grid, e.g. to set the grid
        or the wavelength sampling.
    rates
        Precomputed rate inputs of the mechanisms [s^-1] on the photochemistry
        grid; skips the actinic-flux calculation.
    ionospheric_o1s_production
        For ``"O(1S)"``: O(1S) volume production [m^-3 s^-1] from ionospheric
        processes (N2(A) + O, photoelectron impact, O2+ + e), on the
        photochemistry grid. These dominate above about 105 km; without them
        the green line is underestimated there.

    Returns
    -------
    xr.Dataset
        The solution of :func:`solve` on the photochemistry grid, with the
        rate inputs that drove it.

    Examples
    --------
    A limb A-band calculation with daytime O2(b) emission::

        config = sk.Config()
        config.emission_source = sk.EmissionSource.VolumeEmissionRate
        atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths)
        atmosphere.temperature_k, atmosphere.pressure_pa = temperature, pressure
        atmosphere["O2"] = sk.constituent.VMRAltitudeAbsorber(o2_optics, z, o2_vmr)
        atmosphere["O3"] = sk.constituent.VMRAltitudeAbsorber(o3_optics, z, o3_vmr)
        atmosphere["rayleigh"] = sk.constituent.Rayleigh()
        atmosphere["solar"] = sk.constituent.SolarIrradiance(photon_units=True)

        background = xr.Dataset({"O": ("altitude", o_density)}, coords={"altitude": z})
        solution = sk.nlte.add_photochemical_species(
            atmosphere, ["O2(b)"], cos_sza=0.6, background=background
        )
        radiance = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    The emission is in photon units, so the solar spectrum must be too.
    """
    from sasktran2 import photolysis

    emitters = {_canonical(name): None for name in species}
    requested = {
        PHOTOCHEMICAL_SPECIES[name].mechanism
        for name in emitters
        if PHOTOCHEMICAL_SPECIES[name].mechanism is not None
    }
    order = []
    for name in sorted(requested, key=lambda m: m in _DEPENDS_ON):
        if name in _DEPENDS_ON and _DEPENDS_ON[name][0] not in order:
            order.append(_DEPENDS_ON[name][0])
        if name not in order:
            order.append(name)
    mechanisms = {name: Mechanism.bundled(name) for name in order}
    fluorescent = [n for n in emitters if PHOTOCHEMICAL_SPECIES[n].absorber]

    if not isinstance(atmosphere.model_geometry, sk.Geometry1D):
        msg = "add_photochemical_species supports one-dimensional atmospheres only"
        raise NotImplementedError(msg)
    if atmosphere.temperature_k is None or atmosphere.pressure_pa is None:
        msg = "set the atmosphere's temperature_k and pressure_pa first"
        raise ValueError(msg)

    model_altitude = np.asarray(atmosphere.model_geometry.altitudes(), dtype=float)
    altitude = (
        model_altitude
        if actinic_flux is None
        else np.asarray(actinic_flux.altitudes_m, dtype=float)
    )

    def log_interp(at, source_altitude, values):
        return np.exp(
            np.interp(at, source_altitude, np.log(np.maximum(values, 1e-300)))
        )

    temperature = np.interp(altitude, model_altitude, atmosphere.temperature_k)
    pressure = log_interp(altitude, model_altitude, atmosphere.pressure_pa)
    air = log_interp(
        altitude,
        model_altitude,
        atmosphere.state_equation.dry_air_numberdensity["N"],
    )

    profiles = _background_profiles(background)

    def density(name: str) -> np.ndarray | None:
        if name in profiles:
            values = profiles[name]
            return log_interp(
                altitude, values["altitude"].to_numpy(), values.to_numpy()
            )
        constituent = atmosphere[name]
        if isinstance(constituent, sk.constituent.VMRAltitudeAbsorber):
            vmr = np.interp(altitude, constituent.altitudes_m, constituent.vmr)
            return vmr * air
        if name in DEFAULT_VMR:
            return DEFAULT_VMR[name] * air
        return None

    from_solutions = {
        state for name in order if name in _DEPENDS_ON for state in _DEPENDS_ON[name][1]
    }
    needed = [
        *(
            b
            for m in mechanisms.values()
            for b in m.background
            if b not in from_solutions and b != "M"
        ),
        *(PHOTOCHEMICAL_SPECIES[n].absorber for n in fluorescent),
    ]
    densities = {name: density(name) for name in dict.fromkeys(needed)}
    missing = [name for name, values in densities.items() if values is None]
    if missing:
        msg = (
            f"No number densities for {', '.join(missing)}: give them in "
            "`background`, or add them to the atmosphere as VMRAltitudeAbsorber "
            "constituents"
        )
        raise ValueError(msg)

    state = {
        "temperature_k": ("altitude", temperature),
        "pressure_pa": ("altitude", pressure),
    }
    chemistry = xr.Dataset(
        {**state, **{name: ("altitude", values) for name, values in densities.items()}},
        coords={"altitude": altitude},
    )

    flux = None
    if fluorescent or (mechanisms and rates is None):
        absorbers = {
            name: values
            for name in PHOTOLYSIS_ABSORBERS
            if (values := density(name)) is not None
        }
        optics = xr.Dataset(
            {**state, **{name: ("altitude", v) for name, v in absorbers.items()}},
            coords={"altitude": altitude},
        )
        calculator = actinic_flux or photolysis.ActinicFlux(altitude)
        flux = calculator.calculate(
            optics,
            cos_sza=cos_sza,
            albedo=albedo,
            earth_sun_distance_au=earth_sun_distance_au,
        )
    if mechanisms and rates is None:
        rates = xr.merge([_RATE_FUNCTIONS[name](flux) for name in order])
    if "oxygen_green" in mechanisms:
        if ionospheric_o1s_production is None:
            warnings.warn(
                "No ionospheric O(1S) production given; the green line is "
                "underestimated above about 100 km",
                stacklevel=2,
            )
            production = np.zeros_like(altitude)
        else:
            production = np.asarray(ionospheric_o1s_production, dtype=float)
        rates = rates.assign(
            P_O1S_ION=("altitude", production / np.maximum(densities["O(3P)"], 1.0))
        )

    results = []
    solutions = {}
    for name, mechanism in mechanisms.items():
        inputs = chemistry.assign(M=("altitude", air))
        if name in _DEPENDS_ON:
            parent, states = _DEPENDS_ON[name]
            for state in states:
                inputs[state] = solutions[parent]["density"].sel(state=state)
        solutions[name] = solve(mechanism, inputs, rates)
        solution = solutions[name]
        if len(mechanisms) > 1:
            solution = solution.rename(relative_residual=f"relative_residual_{name}")
        results += [
            solution,
            xr.Dataset(
                {
                    rate: ("altitude", np.asarray(rates[rate]))
                    for rate in mechanism.rate_inputs
                    if rate in rates
                },
                coords={"altitude": altitude},
            ),
        ]

    config = getattr(atmosphere, "_config", None)
    if (
        config is not None
        and config.emission_source != sk.EmissionSource.VolumeEmissionRate
    ):
        warnings.warn(
            "Photochemical emission reaches the radiance only with "
            "config.emission_source = sk.EmissionSource.VolumeEmissionRate",
            stacklevel=2,
        )
    for name in emitters:
        emitter = PHOTOCHEMICAL_SPECIES[name]
        if emitter.absorber is not None:
            ver, constituent = _fluorescence(
                name, flux, temperature, densities[emitter.absorber], altitude
            )
            atmosphere[f"{name} emission"] = constituent
            results.append(
                xr.Dataset(
                    {f"{name} photon_ver": ("altitude", ver)},
                    coords={"altitude": altitude},
                )
            )
            continue
        solution = solutions[emitter.mechanism]
        parts = {
            transition: make(
                altitude,
                solution["photon_ver"].sel(transition=transition).to_numpy(),
            )
            for transition, make in (emitter.bands or {}).items()
        }
        for transition, (wavelength, mass) in (emitter.lines or {}).items():
            parts[transition] = sk.constituent.MonochromaticVolumeEmissionRate(
                altitude,
                solution["photon_ver"].sel(transition=transition).to_numpy(),
                wavelength,
                line_shape="doppler",
                emitter_molecular_weight_g_per_mol=mass,
            )
        atmosphere[f"{name} emission"] = _Combined(parts)

    return xr.merge(results, compat="no_conflicts", join="outer")


def _fluorescence(name, flux, temperature, density, altitude):
    from . import fluorescence

    if name != "OH(A)":
        msg = f"no fluorescence model for {name}"
        raise NotImplementedError(msg)
    lines = fluorescence.oh_ax_lines()
    at_lines = fluorescence.actinic_flux_at_lines(flux, lines.wavelength_nm)
    ver, weights = fluorescence.oh_ax_fluorescence(
        lines, temperature, density, at_lines
    )
    constituent = sk.constituent.LineListVolumeEmissionRate(
        altitude,
        ver,
        lines.wavelength_nm,
        weights,
        molecular_mass_g_per_mol=fluorescence.OH_MOLAR_MASS,
    )
    return ver, constituent
