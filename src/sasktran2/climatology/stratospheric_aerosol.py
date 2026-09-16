"""Add stratospheric sulfate aerosol using SAGE III-ISS reference scenarios.

Use :func:`constituent` to add aerosol to an atmosphere, :func:`scenarios` to
list the twelve available cases, and :func:`profile` to inspect or customize
extinction and particle size. Choose low, typical, elevated or extreme loading
for southern midlatitudes, the tropics or northern midlatitudes.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr
from scipy.integrate import trapezoid
from scipy.ndimage import gaussian_filter1d

import sasktran2 as sk
from sasktran2.constituent import ExtinctionScatterer
from sasktran2.database.stratospheric_aerosol import (
    StratosphericAerosolDatabase,
    _sha256,
)
from sasktran2.optical.base import OpticalProperty


def load_dataset(
    version: str = "v1",
    *,
    path: str | Path | None = None,
    db_root: str | Path | None = None,
) -> xr.Dataset:
    """Load unsmoothed observations, formal errors, QA flags and provenance.

    The catalogue is included with SASKTRAN2; no source archive is needed.
    Use ``path`` to load a custom NetCDF catalogue or ``db_root`` to choose
    the cache directory. These two options are mutually exclusive.

    ``observed_valid`` marks each scenario's contiguous reliable interval on
    the shared ``altitude_m`` coordinate; padding outside it is missing data.
    All nine source extinction channels are retained. Source radius units are
    normalized to nm and extinction to m^-1. Formal source errors are neither
    propagated through smoothing nor assigned to the modeled extensions.
    """
    if version != "v1":
        msg = f"Unsupported aerosol catalogue version {version!r}; available: ['v1']"
        raise ValueError(msg)
    if path is not None and db_root is not None:
        msg = "Specify either path or db_root, not both"
        raise ValueError(msg)
    source = (
        Path(path)
        if path is not None
        else StratosphericAerosolDatabase(version, db_root).path()
    )
    with xr.open_dataset(source) as data:
        result = data.load()
    if (
        result.attrs.get("catalogue_version") != version
        or result.attrs.get("schema_version") != 1
    ):
        msg = "Incompatible stratospheric aerosol catalogue schema/version"
        raise ValueError(msg)
    result.attrs.update(
        source_file=source.name,
        source_sha256=_sha256(source),
        source_mode="local" if path is not None else "verified_cache",
    )
    return result


def scenarios(
    version: str = "v1",
    *,
    path: str | Path | None = None,
    db_root: str | Path | None = None,
) -> xr.Dataset:
    """List case names and metadata, including loading percentiles and bounds.

    Names are ``sh_midlat_*``, ``tropical_*``, and ``nh_midlat_*``, with suffixes
    ``low``, ``typical``, ``elevated``, and ``extreme``. Loading labels refer to
    the 18-30 km 756 nm optical-depth distribution of the screened archive.
    Selection is exact; no nearest latitude or date selection is performed.
    See :func:`load_dataset` for version and local-file options.
    """
    data = load_dataset(version, path=path, db_root=db_root)
    return data[
        [name for name, value in data.data_vars.items() if value.dims == ("scenario",)]
    ]


def _grid(values, name):
    result = np.asarray(values, dtype=float)
    if (
        result.ndim != 1
        or len(result) < 2
        or not np.isfinite(result).all()
        or np.any(np.diff(result) <= 0)
    ):
        msg = f"{name} must contain at least two finite, strictly increasing altitudes"
        raise ValueError(msg)
    return result


def _positive(value, name, *, allow_zero=False):
    try:
        result = float(value)
    except (TypeError, ValueError):
        result = np.nan
    if not np.isfinite(result) or result < 0 or (result == 0 and not allow_zero):
        msg = f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}"
        raise ValueError(msg)
    return result


def _observations(data, scenario):
    if "scenario" not in data.coords or data.scenario.dims != ("scenario",):
        msg = "Catalogue requires a one-dimensional scenario coordinate"
        raise ValueError(msg)
    names = data.scenario.values.tolist()
    if names.count(scenario) != 1:
        msg = f"Select an exact aerosol scenario; got {scenario!r}, available: {names}"
        raise ValueError(msg)
    selected = data.sel(scenario=scenario)
    required = {
        "observed_valid",
        "observed_bottom_m",
        "observed_top_m",
        "source_aerosol_flag",
    }
    if missing := required - set(selected.variables):
        msg = f"Aerosol catalogue is missing required fields: {sorted(missing)}"
        raise ValueError(msg)
    for name, units in [
        ("altitude_m", "m"),
        ("wavelength_nm", "nm"),
        ("raw_extinction_per_m", "m-1"),
        ("raw_median_radius_nm", "nm"),
        ("reference_upper_scale_height_m", "m"),
    ]:
        if name not in selected or selected[name].attrs.get("units") != units:
            msg = f"Aerosol catalogue {name} must have units {units!r}"
            raise ValueError(msg)
    z = _grid(selected.altitude_m, "Source altitude_m")
    if selected.observed_valid.dims != ("altitude_m",):
        msg = "observed_valid must be an altitude mask"
        raise ValueError(msg)
    valid = selected.observed_valid.values
    indices = np.flatnonzero(valid == 1)
    if (
        not np.isin(valid, [0, 1]).all()
        or len(indices) < 2
        or np.any(np.diff(indices) != 1)
    ):
        msg = "Observed interval must be contiguous; internal gaps cannot be extended or smoothed"
        raise ValueError(msg)
    selected = selected.isel(altitude_m=indices)
    z = z[indices]
    if not np.allclose(np.diff(z), np.diff(z)[0], rtol=1e-8, atol=1e-6):
        msg = "Catalogue observations must use a regular native altitude grid"
        raise ValueError(msg)
    if (
        float(selected.observed_bottom_m) != z[0]
        or float(selected.observed_top_m) != z[-1]
    ):
        msg = "Catalogue bounds do not match its observed interval"
        raise ValueError(msg)
    if (
        data.attrs.get("reference_wavelength_nm") != 756
        or data.attrs.get("mode_width") != 1.6
    ):
        msg = "Expected 756 nm extinction and fixed lognormal width 1.6"
        raise ValueError(msg)
    if (
        selected.raw_extinction_per_m.dims
        != (
            "altitude_m",
            "wavelength_nm",
        )
        or selected.raw_median_radius_nm.dims != ("altitude_m",)
        or selected.source_aerosol_flag.dims != ("altitude_m", "wavelength_nm")
    ):
        msg = "Catalogue extinction/radius dimensions are invalid"
        raise ValueError(msg)
    extinction = selected.raw_extinction_per_m.sel(wavelength_nm=756).values
    radius = selected.raw_median_radius_nm.values
    if (
        not np.isfinite(extinction).all()
        or np.any(extinction <= 0)
        or not np.isfinite(radius).all()
        or np.any((radius <= 10) | (radius >= 590))
    ):
        msg = "Observed extinction must be positive and finite; median radii must lie strictly within 10-590 nm"
        raise ValueError(msg)
    flags = selected.source_aerosol_flag.sel(wavelength_nm=[756, 869, 1021, 1543])
    if not np.isin(flags, [2, 3]).all():
        msg = "Observed interval contains cloud/invalid source aerosol flags"
        raise ValueError(msg)
    return selected, z, extinction, radius


def profile(
    scenario: str,
    *,
    altitudes_m: np.ndarray | None = None,
    smoothing_fwhm_m: float = 1500.0,
    lower_extension: str = "exponential_to_zero",
    upper_extension: str = "exponential",
    lower_scale_height_m: float = 2000.0,
    upper_scale_height_m: float | str = "reference",
    ground_altitude_m: float = 0.0,
    version: str = "v1",
    path: str | Path | None = None,
    db_root: str | Path | None = None,
) -> xr.Dataset:
    """Prepare one paired extinction/radius scenario on a model altitude grid.

    Parameters
    ----------
    scenario : str
        Exact case name returned by :func:`scenarios`.
    altitudes_m : array or None
        Strictly increasing finite altitudes in metres. By default, use 0-100
        km at 500 m spacing, including the ground and native observed levels.
    smoothing_fwhm_m : float
        Smoothing width (Gaussian FWHM) in metres, default 1500 m. Smooths
        log extinction while preserving optical depth over the observed
        interval. Zero disables smoothing. Radius is unchanged.
    lower_extension : str
        ``"exponential_to_zero"`` (default) or ``"zero"`` outside the core.
    upper_extension : str
        ``"exponential"`` (default) or ``"zero"`` outside the core.
    lower_scale_height_m : float
        Positive downward-taper scale height, default 2000 m.
    upper_scale_height_m : float or str
        ``"reference"`` reads the fixed regular-condition scale height for the
        latitude band, shared by all loading tiers. A positive number overrides
        it. Changing smoothing never refits this reference.
    ground_altitude_m : float
        Finite ground height below the observed interval; extinction is exactly
        zero at/below this height. The lower taper models this stratospheric
        component only, not an independent tropospheric aerosol population.
    version, path, db_root
        See :func:`load_dataset`.

    Returns
    -------
    xarray.Dataset
        ``extinction_per_m``, ``median_radius_nm`` and ``region`` on
        ``altitude_m``. Raw observations/errors/flags and smoothed native
        extinction use the separate ``observed_altitude_m`` coordinate.
        Radius is held at the respective endpoint in both extensions.
        ``core_aod``, ``lower_extension_aod`` and
        ``upper_extension_aod_to_infinity`` describe the continuous prepared
        profile, independently of the requested output grid. Linear resampling
        onto a coarser grid need not preserve its numerical integral.

    Notes
    -----
    All extensions are modeled, including the numerical continuation to 100 km.
    Their extinction/particle size are not observations. Boundaries are
    continuous for exponential extensions, but derivatives need not match.
    Zero extensions may be discontinuous. No internal gap filling is performed.
    """
    if lower_extension not in (
        "exponential_to_zero",
        "zero",
    ) or upper_extension not in ("exponential", "zero"):
        msg = "Invalid extension: lower must be 'exponential_to_zero'/'zero', upper 'exponential'/'zero'"
        raise ValueError(msg)
    width = _positive(smoothing_fwhm_m, "smoothing_fwhm_m", allow_zero=True)
    lower_h = _positive(lower_scale_height_m, "lower_scale_height_m")
    ground = float(ground_altitude_m)
    if not np.isfinite(ground):
        msg = "ground_altitude_m must be finite"
        raise ValueError(msg)
    data = load_dataset(version, path=path, db_root=db_root)
    raw, z, extinction, radius = _observations(data, scenario)
    if ground >= z[0]:
        msg = "ground_altitude_m must be below the observed interval"
        raise ValueError(msg)
    reference = (
        isinstance(upper_scale_height_m, str) and upper_scale_height_m == "reference"
    )
    upper_h = _positive(
        (
            float(raw.reference_upper_scale_height_m)
            if reference
            else upper_scale_height_m
        ),
        "upper_scale_height_m",
    )
    target = (
        np.unique(np.concatenate((np.arange(0.0, 100001.0, 500.0), z, [ground])))
        if altitudes_m is None
        else _grid(altitudes_m, "altitudes_m")
    )
    core_aod = float(trapezoid(extinction, z))
    smooth = extinction.copy()
    correction = 1.0
    if width > 0:
        smooth = np.exp(
            gaussian_filter1d(
                np.log(extinction),
                width / np.sqrt(8 * np.log(2)) / (z[1] - z[0]),
                mode="reflect",
            )
        )
        correction = core_aod / float(trapezoid(smooth, z))
        smooth *= correction
    prepared = np.interp(target, z, smooth)
    below, above = target < z[0], target > z[-1]
    prepared[below | above] = 0
    lower_aod = upper_aod = 0.0
    if lower_extension == "exponential_to_zero":
        region = below & (target > ground)
        span = z[0] - ground
        # Equivalent to expm1((z-ground)/H)/expm1(span/H), without overflow.
        denominator = -np.expm1(-span / lower_h)
        prepared[region] = (
            smooth[0]
            * np.exp((target[region] - z[0]) / lower_h)
            * (-np.expm1(-(target[region] - ground) / lower_h))
            / denominator
        )
        x = span / lower_h
        # Series avoids subtracting almost equal terms when H >> span.
        integral = (
            span * (0.5 - x / 12 + x**3 / 720)
            if x < 1e-3
            else lower_h - span * np.exp(-x) / denominator
        )
        lower_aod = float(smooth[0] * integral)
    if upper_extension == "exponential":
        prepared[above] = smooth[-1] * np.exp(-(target[above] - z[-1]) / upper_h)
        upper_aod = float(smooth[-1] * upper_h)
    prepared[target <= ground] = 0
    result = raw.rename(altitude_m="observed_altitude_m").copy(deep=True)
    result["smoothed_observed_extinction_per_m"] = (
        "observed_altitude_m",
        smooth,
        {"units": "m-1"},
    )
    result = result.assign_coords(altitude_m=("altitude_m", target, {"units": "m"}))
    result["extinction_per_m"] = ("altitude_m", prepared, {"units": "m-1"})
    result["median_radius_nm"] = (
        "altitude_m",
        np.interp(target, z, radius),
        {"units": "nm"},
    )
    result["region"] = (
        "altitude_m",
        np.where(below, -1, np.where(above, 1, 0)).astype(np.int8),
        {
            "flag_values": [-1, 0, 1],
            "flag_meanings": "lower_extension observed_core upper_extension",
        },
    )
    for name, value in [
        ("core_aod", core_aod),
        ("lower_extension_aod", lower_aod),
        ("upper_extension_aod_to_infinity", upper_aod),
        ("total_aod_to_infinity", core_aod + lower_aod + upper_aod),
    ]:
        result[name] = xr.DataArray(value, attrs={"units": "1", "wavelength_nm": 756.0})
    result.attrs.update(data.attrs)
    result.attrs.update(
        smoothing_fwhm_m=width,
        smoothing_aod_correction=correction,
        lower_extension=lower_extension,
        upper_extension=upper_extension,
        lower_scale_height_m=lower_h,
        upper_scale_height_m=upper_h,
        upper_scale_height_source="reference" if reference else "override",
        ground_altitude_m=ground,
        radius_extension="constant endpoint",
        uncertainty_treatment="Raw formal uncertainties only; not propagated through smoothing or extensions",
    )
    return result


def constituent(
    scenario: str,
    optical_property: OpticalProperty | None = None,
    **profile_options,
) -> ExtinctionScatterer:
    """Create a sulfate aerosol constituent to add to an atmosphere.

    For example, add a typical tropical scenario on the atmosphere's grid::

        atmosphere["aerosol"] = sk.climatology.stratospheric_aerosol.constituent(
            "tropical_typical", altitudes_m=atmosphere.model_geometry.altitudes()
        )

    ``profile_options`` accepts every keyword of :func:`profile`, including
    ``altitudes_m`` (pass the atmosphere grid for direct evaluation there).
    The default optical property calculates sulfate Mie scattering on demand,
    with the OSIRIS H2SO4 refractive index and fixed lognormal width 1.6.
    Its refractive-index file uses the standard SASKTRAN2 download cache.

    For repeated/many-wavelength calculations pass a :class:`~sasktran2.database.MieDatabase`
    with the same composition and width and a ``median_radius`` argument in
    nm. The caller is responsible for custom optical-property compatibility.
    The constituent is zero outside the prepared grid. Its vertical derivative
    is with respect to prepared extinction, not the original observations.
    """
    data = profile(scenario, **profile_options)
    if optical_property is None:
        optical_property = sk.optical.Mie(
            sk.mie.LogNormalDistribution().freeze(mode_width=1.6),
            sk.mie.refractive.H2SO4(),
        )
    return ExtinctionScatterer(
        optical_property,
        data.altitude_m.values,
        data.extinction_per_m.values,
        756.0,
        "zero",
        median_radius=data.median_radius_nm.values,
    )
