"""CAIRT Extended Reference Scenarios (ERS) v7.

Select latitude-band and seasonal reference atmospheres from
https://doi.org/10.5281/zenodo.10022129 (Quentin Errera, CC BY 4.0).
Several long-lived gases are scaled to expected 2030 abundances. These are
multiannual reference scenarios, not atmospheres for an arbitrary calendar year.

Selection is exact: latitude centers are -80, -45, 0, 45, 80 degrees; months
are 1, 4, 7, 10; local times are 9.5 and 21.5 hours. Raw and selected profiles
preserve source values, including known suspect data. Selected profiles warn
about suspect fields, and the constituent helpers validate gas VMRs before use.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import xarray as xr

import sasktran2 as sk
from sasktran2.constituent import VMRAltitudeAbsorber
from sasktran2.database.ers import ERSDatabase, _md5
from sasktran2.optical.base import OpticalProperty

_GASES = (
    "C2H2",
    "C2H6",
    "CH4",
    "CO",
    "H2O",
    "HCN",
    "HO2NO2",
    "N2O",
    "N2O5",
    "NH3",
    "NO",
    "NO2",
    "O3",
    "PAN",
    "HNO3",
    "ClO",
    "ClONO2",
    "BrONO2",
    "CCl4",
    "CFC11",
    "CFC12",
    "CO2",
    "HCFC22",
    "SF6",
    "OCS",
    "SO2",
    "O1D",
    "O",
    "O2",
    "CF4",
    "HDO",
)
_SPECIES = {name.upper(): name.lower() for name in _GASES}
_SPECIES.update(
    {
        "HNO4": "ho2no2",
        "F11": "cfc11",
        "CFCL3": "cfc11",
        "F12": "cfc12",
        "CF2CL2": "cfc12",
        "F22": "hcfc22",
        "CHCLF2": "hcfc22",
        "F14": "cf4",
    }
)
_STATE_NAMES = {
    "temperature_mean": "temperature_k",
    "temperature_std": "temperature_std_k",
    "pressure_mean": "pressure_pa",
    "surface_pressure_mean": "surface_pressure_pa",
}


def _species_name(species: str) -> str:
    try:
        return _SPECIES[species.upper()]
    except KeyError:
        msg = f"Unsupported ERS gas {species!r}; available gases: {', '.join(_GASES)}"
        raise ValueError(msg) from None


def load_dataset(
    version: str = "v07",
    *,
    path: str | Path | None = None,
    db_root: str | Path | None = None,
) -> xr.Dataset:
    """Return raw ERS arrays, loaded into memory with the file closed.

    Parameters
    ----------
    version : str, optional
        Dataset schema/version, currently ``"v07"`` only.
    path : str or Path or None, optional
        Explicit local NetCDF. This bypasses downloading and the pinned checksum,
        allowing offline subsets and modified files. The actual checksum is
        recorded in the returned metadata.
    db_root : str or Path or None, optional
        Override the configured download cache root. Cannot be combined with path.

    Returns
    -------
    xarray.Dataset
        Original coordinates, variables and units, with source provenance added.
        Month and hour remain numeric. No negative values or zeros are corrected.
    """
    if version != "v07":
        msg = f"Unsupported ERS version {version!r}; available versions: ['v07']"
        raise ValueError(msg)
    if path is not None and db_root is not None:
        msg = "Specify either path or db_root, not both"
        raise ValueError(msg)
    source = Path(path) if path is not None else ERSDatabase(version, db_root).path()
    with xr.open_dataset(source, decode_timedelta=False) as data:
        dataset = data.load()
    dataset.attrs.update(
        ers_version=version,
        source_file=source.name,
        source_md5=_md5(source),
        source_doi="10.5281/zenodo.10022129",
        source_mode="local" if path else "verified_cache",
        source_attribution="Quentin Errera; CAIRT Extended Reference Scenarios v7; CC BY 4.0",
        scaled_to_2030="CH4, N2O, CO2, SF6, CCl4, HCFC22, CFC11, CFC12",
        vmr_basis="Mole fraction in air; wet/dry convention not explicitly specified by source",
    )
    return dataset


def _altitudes(data: xr.Dataset) -> np.ndarray:
    altitude = data["altitude_m"].to_numpy()
    if (
        altitude.ndim != 1
        or len(altitude) < 2
        or not np.isfinite(altitude).all()
        or np.any(np.diff(altitude) <= 0)
    ):
        msg = "ERS altitude must be finite and strictly increasing with at least two levels"
        raise ValueError(msg)
    return altitude


def _units(variable: xr.DataArray, expected: str):
    if variable.attrs.get("units") != expected:
        msg = f"ERS {variable.name} must have units {expected!r}, got {variable.attrs.get('units')!r}"
        raise ValueError(msg)


def profile(
    *,
    month: int,
    latitude_degrees: float,
    local_time_hours: float,
    solar_activity: str = "minimum",
    volcanic_activity: str = "background",
    species: str | Sequence[str] | None = None,
    version: str = "v07",
    path: str | Path | None = None,
    db_root: str | Path | None = None,
) -> xr.Dataset:
    """Return an exact ERS scenario, preserving means and standard deviations.

    Parameters
    ----------
    month : int
        One of 1, 4, 7, 10. No year dependence is implied.
    latitude_degrees : float
        One of -80, -45, 0, 45, 80: centers of broad latitude bands.
    local_time_hours : float
        Either 9.5 or 21.5. Variables without an hour dimension retain their
        monthly mean. These times do not always correspond to day and night.
    solar_activity : str, optional
        ``"minimum"`` (default) or ``"maximum"``.
    volcanic_activity : str, optional
        ``"background"`` (default) or ``"enhanced"``.
    species : str or sequence of str or None, optional
        Gas names, case insensitive. None includes all available gases plus
        condensed sulfuric acid and air molar mass. State fields and latitude
        bounds are included regardless of species selection.
    version, path, db_root
        See :func:`load_dataset`.

    Returns
    -------
    xarray.Dataset
        Profiles on ``altitude_m``. State variables are ``temperature_k``,
        ``temperature_std_k``, ``pressure_pa``, and ``surface_pressure_pa``.
        Gases retain lower-case ``<species>_mean`` / ``<species>_std`` names and
        mol/mol units. Each variable's ``quality_flags`` attribute describes
        detected issues; values are not repaired. Standard deviations represent
        climatological spread, not uncertainty in the mean or covariance.

    Warns
    -----
    UserWarning
        If a returned field has suspect values. The warning names the affected
        fields; selecting specific species avoids warnings about unused gases
        and air molar mass.
    """
    result = _select_profile(
        month=month,
        latitude_degrees=latitude_degrees,
        local_time_hours=local_time_hours,
        solar_activity=solar_activity,
        volcanic_activity=volcanic_activity,
        species=species,
        version=version,
        path=path,
        db_root=db_root,
    )
    descriptions = {
        "nonfinite_values": "contains nonfinite values",
        "negative_values": "contains negative values",
        "all_zero_profile": "is entirely zero and requires verification before use",
        "source_molar_mass_requires_verification": "has suspect altitude dependence in ERS v7; verify before using it to calculate density",
    }
    issues = [
        f"{name} {descriptions[flag]}"
        for name, variable in result.data_vars.items()
        for flag in variable.attrs["quality_flags"].split(",")
        if flag
    ]
    if issues:
        warnings.warn(
            "Suspect ERS profile data: "
            + "; ".join(issues)
            + ". Source values are unchanged; inspect these fields before use.",
            UserWarning,
            stacklevel=2,
        )
    return result


def _select_profile(
    *,
    month: int,
    latitude_degrees: float,
    local_time_hours: float,
    solar_activity: str = "minimum",
    volcanic_activity: str = "background",
    species: str | Sequence[str] | None = None,
    version: str = "v07",
    path: str | Path | None = None,
    db_root: str | Path | None = None,
) -> xr.Dataset:
    """Select source values; callers warn or validate the fields they use."""
    requested = [species] if isinstance(species, str) else species
    stems = (
        None
        if requested is None
        else list(dict.fromkeys(_species_name(s) for s in requested))
    )
    categories = {
        "solmin_solmax": {"minimum": 0, "maximum": 1},
        "volcanism": {"background": 0, "enhanced": 1},
    }
    for dimension, choice in [
        ("solmin_solmax", solar_activity),
        ("volcanism", volcanic_activity),
    ]:
        if choice not in categories[dimension]:
            msg = f"Invalid ERS {dimension} choice {choice!r}; choose from {list(categories[dimension])}"
            raise ValueError(msg)
    selectors = {
        "month": month,
        "lat": latitude_degrees,
        "hour": local_time_hours,
        "solmin_solmax": categories["solmin_solmax"][solar_activity],
        "volcanism": categories["volcanism"][volcanic_activity],
    }
    data = load_dataset(version, path=path, db_root=db_root)
    for dimension, value in selectors.items():
        if dimension not in data.coords or data[dimension].dims != (dimension,):
            msg = f"ERS requires a one-dimensional {dimension} coordinate"
            raise ValueError(msg)
        values = data[dimension].to_numpy()
        if np.count_nonzero(values == value) != 1:
            msg = f"ERS requires exact {dimension} selection; got {value!r}, available: {values.tolist()}"
            raise ValueError(msg)
    _units(data["lev"], "km")
    if stems is None:
        stems = [name.lower() for name in _GASES if f"{name.lower()}_mean" in data]
    names = [name for name in _STATE_NAMES if name in data]
    for stem in stems:
        for suffix in ("mean", "std"):
            name = f"{stem}_{suffix}"
            if name not in data:
                msg = f"ERS file does not contain {name}"
                raise ValueError(msg)
            _units(data[name], "mol mol-1")
            names.append(name)
    for name in names:
        if name in _STATE_NAMES:
            _units(data[name], "K" if name.startswith("temperature") else "Pa")
    if "lat_bnds" in data:
        names.append("lat_bnds")
    if requested is None:
        names.extend(
            name
            for name in ("airmolmass", "h2so4m_c_mean", "h2so4m_c_std")
            if name in data
        )
    # Dataset.sel applies each selector only to variables with that dimension.
    result = (
        data[names]
        .sel({k: v for k, v in selectors.items() if k in data[names].dims})
        .copy(deep=True)
    )
    result = result.rename(
        {"lev": "altitude_m", **{k: v for k, v in _STATE_NAMES.items() if k in result}}
    )
    result = result.assign_coords(altitude_m=result.altitude_m.astype(float) * 1000)
    result.altitude_m.attrs = {"units": "m", "standard_name": "altitude"}
    _altitudes(result)
    for name, variable in result.data_vars.items():
        flags = []
        if not np.isfinite(variable).all():
            flags.append("nonfinite_values")
        if name != "lat_bnds" and bool((variable < 0).any()):
            flags.append("negative_values")
        if name in ("o_mean", "o1d_mean", "h2so4m_c_mean") and bool(
            (variable == 0).all()
        ):
            flags.append("all_zero_profile")
        if name == "airmolmass":
            flags.append("source_molar_mass_requires_verification")
        variable.attrs["quality_flags"] = ",".join(flags)
    result.attrs.update(
        selection="exact",
        selected_month=float(month),
        selected_latitude_degrees=float(latitude_degrees),
        selected_local_time_hours=float(local_time_hours),
        solar_activity=solar_activity,
        volcanic_activity=volcanic_activity,
    )
    return result


def _validated_vmr(data: xr.Dataset, stem: str, negative_vmr: str) -> np.ndarray:
    if negative_vmr not in ("raise", "clip"):
        msg = "negative_vmr must be 'raise' or 'clip'"
        raise ValueError(msg)
    variable = data[f"{stem}_mean"]
    if variable.dims != ("altitude_m",):
        msg = f"ERS {stem} must be a single altitude profile after selection"
        raise ValueError(msg)
    vmr = variable.to_numpy().copy()
    if not np.isfinite(vmr).all() or np.any(vmr > 1):
        msg = f"ERS {stem} contains nonfinite VMRs or VMRs greater than one"
        raise ValueError(msg)
    if stem in ("o", "o1d") and np.all(vmr == 0):
        msg = f"ERS {stem} is an all-zero atomic-oxygen scenario of uncertain validity; inspect profile() before use"
        raise ValueError(msg)
    if np.any(vmr < 0):
        if negative_vmr == "raise":
            msg = f"ERS {stem} contains negative VMRs; inspect profile() or explicitly use negative_vmr='clip'"
            raise ValueError(msg)
        warnings.warn(
            f"Clipped {np.count_nonzero(vmr < 0)} negative ERS {stem} VMRs to zero",
            UserWarning,
            stacklevel=3,
        )
        vmr = np.maximum(vmr, 0)
    return vmr


def constituent(
    species: str,
    optical_property: OpticalProperty,
    *,
    month: int,
    latitude_degrees: float,
    local_time_hours: float,
    solar_activity: str = "minimum",
    volcanic_activity: str = "background",
    negative_vmr: str = "raise",
    out_of_bounds_mode: str = "zero",
    version: str = "v07",
    path: str | Path | None = None,
    db_root: str | Path | None = None,
) -> VMRAltitudeAbsorber:
    """Create a gas constituent on the native ERS altitude grid.

    Scenario and source arguments are described in :func:`profile`.
    ``negative_vmr`` is ``"raise"`` (default) or ``"clip"`` (warn and set negative
    values to zero). All-zero O/O1D scenarios are rejected. ``out_of_bounds_mode``
    is ``"zero"`` or ``"extend"``, as in :class:`~sasktran2.constituent.VMRAltitudeAbsorber`.

    Gas VMR is passed through without a wet/dry conversion. Use with unset or
    zero specific humidity; the source convention requires further verification.
    HDO requires suitable isotope-specific optical properties supplied by the
    caller. Condensed H2SO4 is available through :func:`profile`, not this helper.
    """
    if out_of_bounds_mode not in ("zero", "extend"):
        msg = "out_of_bounds_mode must be 'zero' or 'extend'"
        raise ValueError(msg)
    data = _select_profile(
        month=month,
        latitude_degrees=latitude_degrees,
        local_time_hours=local_time_hours,
        solar_activity=solar_activity,
        volcanic_activity=volcanic_activity,
        species=[species],
        version=version,
        path=path,
        db_root=db_root,
    )
    return VMRAltitudeAbsorber(
        optical_property,
        _altitudes(data),
        _validated_vmr(data, _species_name(species), negative_vmr),
        out_of_bounds_mode,
    )


def add_to_atmosphere(
    atmosphere: sk.Atmosphere,
    species: Mapping[str, OpticalProperty],
    *,
    month: int,
    latitude_degrees: float,
    local_time_hours: float,
    solar_activity: str = "minimum",
    volcanic_activity: str = "background",
    set_pressure_temperature: bool = True,
    negative_vmr: str = "raise",
    out_of_bounds_mode: str = "raise",
    version: str = "v07",
    path: str | Path | None = None,
    db_root: str | Path | None = None,
) -> None:
    """Add ERS gas constituents and optionally pressure/temperature to an atmosphere.

    ``species`` maps gas names to caller-supplied optical properties. Scenario
    and source arguments are described in :func:`profile`. The file is loaded
    once and all requested inputs are validated before updating the atmosphere.

    Temperature and VMR are interpolated linearly; pressure is interpolated in
    log space. Source levels bracketing the atmosphere altitude range are kept,
    so unrelated negative VMRs above that range do not block a calculation.
    ``negative_vmr`` is ``"raise"`` or ``"clip"`` (with a warning).

    ``out_of_bounds_mode`` is ``"raise"`` (default) for model altitudes outside
    the source grid, or ``"extend"`` to hold boundary VMR, temperature, and
    pressure constant. No physical extrapolation is inferred.

    Specific humidity is not populated from H2O. Nonzero existing humidity is
    rejected when adding gases, to avoid an unverified wet/dry VMR correction.
    ``set_pressure_temperature=False`` leaves the existing state unchanged.
    The supplied surface pressure and air molar mass are not used to reconstruct
    the vertical pressure profile. An empty species mapping can set state only.
    """
    if out_of_bounds_mode not in ("raise", "extend"):
        msg = "out_of_bounds_mode must be 'raise' or 'extend'"
        raise ValueError(msg)
    if negative_vmr not in ("raise", "clip"):
        msg = "negative_vmr must be 'raise' or 'clip'"
        raise ValueError(msg)
    if (
        species
        and atmosphere.specific_humidity is not None
        and np.any(atmosphere.specific_humidity != 0)
    ):
        msg = "ERS gas helpers require unset or zero specific_humidity until the source wet/dry convention is verified"
        raise ValueError(msg)
    data = _select_profile(
        month=month,
        latitude_degrees=latitude_degrees,
        local_time_hours=local_time_hours,
        solar_activity=solar_activity,
        volcanic_activity=volcanic_activity,
        species=list(species),
        version=version,
        path=path,
        db_root=db_root,
    )
    source_altitudes = _altitudes(data)
    target = np.asarray(atmosphere.model_geometry.altitudes())
    if (
        target.ndim != 1
        or target.size == 0
        or not np.isfinite(target).all()
        or np.any(np.diff(target) <= 0)
    ):
        msg = "Atmosphere altitudes must be finite and strictly increasing"
        raise ValueError(msg)
    if out_of_bounds_mode == "raise" and (
        target[0] < source_altitudes[0] or target[-1] > source_altitudes[-1]
    ):
        msg = f"Atmosphere altitudes are outside ERS range {source_altitudes[0]} to {source_altitudes[-1]} m"
        raise ValueError(msg)
    lower = np.clip(
        np.searchsorted(source_altitudes, target[0], side="right") - 1,
        0,
        len(source_altitudes) - 2,
    )
    upper = np.clip(
        np.searchsorted(source_altitudes, target[-1]), 1, len(source_altitudes) - 1
    )
    data = data.isel(altitude_m=slice(lower, upper + 1))
    altitudes = _altitudes(data)
    prepared = {
        name: VMRAltitudeAbsorber(
            optical,
            altitudes,
            _validated_vmr(data, _species_name(name), negative_vmr),
            "extend",
        )
        for name, optical in species.items()
    }
    if set_pressure_temperature:
        for name in ("temperature_k", "pressure_pa"):
            if name not in data or data[name].dims != ("altitude_m",):
                msg = f"ERS requires a {name} altitude profile"
                raise ValueError(msg)
            if not np.isfinite(data[name]).all() or bool((data[name] <= 0).any()):
                msg = f"ERS {name} must be finite and positive"
                raise ValueError(msg)
        if np.any(np.diff(data.pressure_pa) >= 0):
            msg = "ERS pressure must decrease with altitude"
            raise ValueError(msg)
        temperature = np.interp(target, altitudes, data.temperature_k)
        pressure = np.exp(np.interp(target, altitudes, np.log(data.pressure_pa)))
    if set_pressure_temperature:
        atmosphere.temperature_k = temperature
        atmosphere.pressure_pa = pressure
    for name, absorber in prepared.items():
        atmosphere[name] = absorber
