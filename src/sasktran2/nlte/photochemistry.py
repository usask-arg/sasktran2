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
PHOTOLYSIS_ABSORBERS = ("O2", "O3", "NO2")


@dataclass(frozen=True)
class _Emitter:
    mechanism: str
    #: Band constituent factory by mechanism transition id.
    bands: dict[str, Callable[[np.ndarray, np.ndarray], object]]
    description: str


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
}

_NOT_YET_AVAILABLE = {
    "O2(a)": "the O2 a-X (1.27 um) band emission is not implemented yet",
    "O(1S)": "the oxygen mechanism has no O(1S) source yet",
}


def _oxygen_rates(flux: xr.Dataset) -> xr.Dataset:
    from sasktran2 import photolysis

    return photolysis.photolysis_rates(flux, photolysis.presets.oxygen_photolysis())


_RATE_FUNCTIONS = {"oxygen": _oxygen_rates}


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
        Precomputed rate inputs of the mechanism [s^-1] on the photochemistry
        grid; skips the actinic-flux calculation.

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
    mechanisms = {PHOTOCHEMICAL_SPECIES[name].mechanism for name in emitters}
    if len(mechanisms) != 1:
        msg = f"species from several mechanisms cannot be combined yet: {mechanisms}"
        raise NotImplementedError(msg)
    (mechanism_name,) = mechanisms
    mechanism = Mechanism.bundled(mechanism_name)

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

    densities = {name: density(name) for name in mechanism.background}
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

    if rates is None:
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
        rates = _RATE_FUNCTIONS[mechanism_name](flux)

    solution = solve(mechanism, chemistry, rates)

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
        atmosphere[f"{name} emission"] = _Combined(
            {
                transition: make(
                    altitude,
                    solution["photon_ver"].sel(transition=transition).to_numpy(),
                )
                for transition, make in emitter.bands.items()
            }
        )

    rate_inputs = xr.Dataset(
        {
            name: ("altitude", np.asarray(rates[name]))
            for name in mechanism.rate_inputs
            if name in rates
        },
        coords={"altitude": altitude},
    )
    return xr.merge([solution, rate_inputs])
