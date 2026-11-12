"""Run TUV-x (v5.4 configuration) on a given atmosphere, for comparison with sasktran2.

Needs the ``musica`` package, which is not a sasktran2 dependency; run this
script in its own environment::

    python -m venv musica-venv && musica-venv/bin/pip install musica xarray netcdf4
    musica-venv/bin/python tools/nlte/tuvx_reference.py inputs.nc outputs.nc

``inputs.nc`` holds, on ``altitude`` [m]: ``temperature_k`` and number
densities ``air``, ``O2``, ``O3`` [m^-3]; scalar attributes ``sza_deg``,
``albedo`` and ``earth_sun_distance_au``; and optionally
``solar_flux_per_bin`` [photons m^-2 s^-1] on TUV-x's wavelength bins
(``wavelength_edge`` [nm]), which replaces TUV-x's extraterrestrial flux.
Aerosols are switched off.

``outputs.nc`` holds the O2 and O3 photolysis rates [s^-1] on the TUV-x
height edges and the actinic flux [photons m^-2 s^-1 nm^-1] by bin, height
and component (direct, upwelling, downwelling).
"""

from __future__ import annotations

import sys

import musica.tuvx.v54 as v54
import numpy as np
import xarray as xr
from musica.tuvx.grid_map import GridMap
from musica.tuvx.profile import Profile
from musica.tuvx.profile_map import ProfileMap
from musica.tuvx.radiator import Radiator
from musica.tuvx.radiator_map import RadiatorMap
from musica.tuvx.tuvx import TUVX

REACTIONS = {
    "J_O3_O1D": "O3+hv->O2+O(1D)",
    "J_O3_O3P": "O3+hv->O2+O(3P)",
    "J_O2": "O2+hv->O+O",
}


def layer_columns_cm2(edges_km, density_cm3) -> np.ndarray:
    """Columns [cm^-2] of each layer, assuming exponential variation within it."""
    dz_cm = np.diff(edges_km) * 1.0e5
    lo, hi = density_cm3[:-1], density_cm3[1:]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log(lo / hi)
        column = np.where(
            (lo > 0) & (hi > 0) & (np.abs(ratio) > 1e-8),
            dz_cm * (lo - hi) / ratio,
            dz_cm * 0.5 * (lo + hi),
        )
    return column


def build_calculator(inputs: xr.Dataset) -> TUVX:
    """A TUV-x v5.4 calculator whose profiles hold ``inputs``.

    Profiles must be complete when the calculator is constructed. Their
    arrays are zero-copy views, and writing them afterwards bypasses TUV-x's
    profile update: the radiators see the new layer densities, but the
    Lyman-alpha and Schumann-Runge band parameterisations keep the O2 column
    of the default v5.4 profile.
    """
    grids = GridMap()
    grids["height", "km"] = v54.height_grid()
    grids["wavelength", "nm"] = v54.wavelength_grid()
    heights = grids["height", "km"]
    wavelengths = grids["wavelength", "nm"]
    edges_km = np.array(heights.edges, copy=True)
    mid_km = 0.5 * (edges_km[:-1] + edges_km[1:])
    altitude_km = inputs["altitude"].to_numpy() / 1.0e3

    def interp_log(values, at_km):
        return np.exp(np.interp(at_km, altitude_km, np.log(np.maximum(values, 1e-30))))

    profiles = ProfileMap()
    for species in ("air", "O2", "O3"):
        density = inputs[species].to_numpy() * 1.0e-6
        edge_values = interp_log(density, edges_km)
        above = altitude_km >= edges_km[-1]
        column_above = np.trapezoid(density[above], altitude_km[above] * 1.0e5)
        # Beyond the input grid, extend with the scale height of its top layer.
        scale_height_cm = (
            (altitude_km[-1] - altitude_km[-2])
            * 1.0e5
            / max(np.log(density[-2] / density[-1]), 1e-6)
            if density[-1] > 0
            else 0.0
        )
        profiles[species, "molecule cm-3"] = Profile(
            name=species,
            units="molecule cm-3",
            grid=heights,
            edge_values=edge_values,
            midpoint_values=interp_log(density, mid_km),
            layer_densities=layer_columns_cm2(edges_km, edge_values),
            exo_layer_density=float(column_above + density[-1] * scale_height_cm),
        )

    temperature = inputs["temperature_k"].to_numpy()
    profiles["temperature", "K"] = Profile(
        name="temperature",
        units="K",
        grid=heights,
        edge_values=np.interp(edges_km, altitude_km, temperature),
        midpoint_values=np.interp(mid_km, altitude_km, temperature),
    )

    num_bins = len(wavelengths.edges) - 1
    albedo = float(inputs.attrs["albedo"])
    profiles["surface albedo", "none"] = Profile(
        name="surface albedo",
        units="none",
        grid=wavelengths,
        edge_values=np.full(num_bins + 1, albedo),
        midpoint_values=np.full(num_bins, albedo),
    )
    if "solar_flux_per_bin" in inputs:
        if not np.allclose(inputs["wavelength_edge"].to_numpy(), wavelengths.edges):
            raise ValueError("solar_flux_per_bin is not on the TUV-x wavelength bins")
        per_bin_cm2 = inputs["solar_flux_per_bin"].to_numpy() * 1.0e-4
        edge_values = np.empty(num_bins + 1)
        edge_values[1:-1] = 0.5 * (per_bin_cm2[:-1] + per_bin_cm2[1:])
        edge_values[0], edge_values[-1] = per_bin_cm2[0], per_bin_cm2[-1]
        profiles["extraterrestrial flux", "photon cm-2 s-1"] = Profile(
            name="extraterrestrial flux",
            units="photon cm-2 s-1",
            grid=wavelengths,
            edge_values=edge_values,
            midpoint_values=per_bin_cm2,
        )
    else:
        profiles["extraterrestrial flux", "photon cm-2 s-1"] = v54.profile(
            "extraterrestrial flux", wavelengths
        )

    shape = (num_bins, len(heights.edges) - 1)
    radiators = RadiatorMap()
    radiators["aerosol"] = Radiator(
        name="aerosol",
        height_grid=heights,
        wavelength_grid=wavelengths,
        optical_depths=np.zeros(shape),
        single_scattering_albedos=np.zeros(shape),
        asymmetry_factors=np.zeros(shape),
    )
    return TUVX(
        grid_map=grids,
        profile_map=profiles,
        radiator_map=radiators,
        config_path=v54.config_file_path(),
    )


def main(inputs_path: str, outputs_path: str) -> None:
    inputs = xr.open_dataset(inputs_path)
    tuvx = build_calculator(inputs)
    # MUSICA arrays are views into memory owned by the map objects: keep the
    # maps alive and copy.
    grid_map, profiles = tuvx.get_grid_map(), tuvx.get_profile_map()
    edges_km = np.array(grid_map["height", "km"].edges, copy=True)
    wavelength_edges = np.array(grid_map["wavelength", "nm"].edges, copy=True)

    result = tuvx.run(
        sza=np.radians(float(inputs.attrs["sza_deg"])),
        earth_sun_distance=float(inputs.attrs["earth_sun_distance_au"]),
    )
    rates = result["photolysis_rate_constants"]
    out = xr.Dataset(
        {
            name: (("altitude_km",), rates.sel(reaction=reaction).to_numpy())
            for name, reaction in REACTIONS.items()
        },
        coords={"altitude_km": edges_km},
    )
    # TUV-x returns the actinic flux normalised to the top-of-atmosphere direct
    # beam; scale by the extraterrestrial flux per nm [photons cm^-2 s^-1 nm^-1].
    extraterrestrial = np.array(
        profiles["extraterrestrial flux", "photon cm-2 s-1"].midpoint_values, copy=True
    ) / np.diff(wavelength_edges)
    distance = float(inputs.attrs["earth_sun_distance_au"])
    out["actinic_flux"] = (
        ("wavelength", "altitude_km", "component"),
        result["actinic_flux"].to_numpy()
        * (extraterrestrial / distance**2)[:, np.newaxis, np.newaxis]
        * 1.0e4,
    )
    out = out.assign_coords(
        wavelength=result["wavelength_midpoint"].to_numpy(),
        wavelength_edge=("wavelength_edge", wavelength_edges),
        component=["direct", "upwelling", "downwelling"],
    )
    out["actinic_flux"].attrs["units"] = "photons m^-2 s^-1 nm^-1"
    out.attrs.update(inputs.attrs)
    out.to_netcdf(outputs_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
