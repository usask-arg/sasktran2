# ruff: noqa: EM101, EM102, T201
"""Run one tangent-policy Geometry1D calculation per OMPS limb image."""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import sasktran2 as sk
import xarray as xr
from omps_orbital_plane_2d import (
    DEFAULT_ANC,
    DEFAULT_L1G,
    OmpsInputs,
    SelectedViewing,
    load_omps_inputs,
    select_viewing_geometry,
)


@dataclass
class ImageRayPolicies:
    viewing: sk.ViewingGeometry
    earth_radius_m: float
    reference_cos_sza: float
    tangent_altitude_m: np.ndarray
    tangent_altitude_residual_m: np.ndarray
    solar_azimuth_deg: np.ndarray
    viewing_azimuth_deg: np.ndarray
    relative_azimuth_deg: np.ndarray
    scattering_angle_deg: np.ndarray
    cos_sza: np.ndarray


def _wrap_degrees(angle_deg: float | np.ndarray) -> float | np.ndarray:
    return (np.asarray(angle_deg) + 180.0) % 360.0 - 180.0


def _sun_vector_ecef(
    solar_handler: sk.solar.SolarGeometryHandlerBase,
    latitude_deg: float,
    longitude_deg: float,
    altitude_m: float,
    timestamp: np.datetime64,
) -> tuple[np.ndarray, float, float]:
    geoid = sk.WGS84()
    geoid.from_lat_lon_alt(latitude_deg, longitude_deg, altitude_m)
    solar_zenith_deg, solar_azimuth_deg = solar_handler.target_solar_angles(
        latitude_deg,
        longitude_deg,
        altitude_m,
        pd.Timestamp(timestamp),
    )
    zenith = np.deg2rad(solar_zenith_deg)
    azimuth = np.deg2rad(solar_azimuth_deg)
    north = -geoid.local_south
    east = -geoid.local_west
    sun = np.cos(zenith) * geoid.local_up + np.sin(zenith) * (
        np.cos(azimuth) * north + np.sin(azimuth) * east
    )
    sun /= np.linalg.norm(sun)
    return sun, float(solar_zenith_deg), float(solar_azimuth_deg % 360.0)


def _azimuth_degrees(vector: np.ndarray, geoid: sk.Geodetic) -> float:
    north = -geoid.local_south
    east = -geoid.local_west
    return float(np.rad2deg(np.arctan2(vector.dot(east), vector.dot(north))) % 360.0)


def image_ray_policies(
    observer_positions_ecef_m: np.ndarray,
    look_directions_ecef: np.ndarray,
    nominal_tangent_altitude_m: np.ndarray,
    sun_vector_ecef: np.ndarray,
) -> ImageRayPolicies:
    """Convert one ECEF image to tangent policies for a spherical 1D engine."""
    viewing = sk.ViewingGeometry()
    size = len(nominal_tangent_altitude_m)
    tangent_altitude_m = np.empty(size)
    tangent_altitude_residual_m = np.empty(size)
    solar_azimuth_deg = np.empty(size)
    viewing_azimuth_deg = np.empty(size)
    relative_azimuth_deg = np.empty(size)
    scattering_angle_deg = np.empty(size)
    cos_sza = np.empty(size)
    surface_radii_m = np.empty(size)

    observer_geoid = sk.WGS84()
    observer_geoid.from_xyz(observer_positions_ecef_m[0])
    observer_altitude_m = float(observer_geoid.altitude)

    for index, (observer, look, nominal_altitude) in enumerate(
        zip(
            observer_positions_ecef_m,
            look_directions_ecef,
            nominal_tangent_altitude_m,
            strict=True,
        )
    ):
        tangent_geoid = sk.WGS84()
        tangent_geoid.from_tangent_point(observer, look)
        tangent_altitude_m[index] = tangent_geoid.altitude
        tangent_altitude_residual_m[index] = tangent_geoid.altitude - nominal_altitude
        cos_sza[index] = np.clip(tangent_geoid.local_up.dot(sun_vector_ecef), -1.0, 1.0)
        solar_azimuth_deg[index] = _azimuth_degrees(sun_vector_ecef, tangent_geoid)
        viewing_azimuth_deg[index] = _azimuth_degrees(look, tangent_geoid)
        relative_azimuth_deg[index] = float(
            _wrap_degrees(solar_azimuth_deg[index] - viewing_azimuth_deg[index])
        )
        # This is the same propagation-direction angle used by PhaseHandler:
        # acos((-sun) dot (-observer-to-scene look)) == acos(sun dot look).
        scattering_angle_deg[index] = np.rad2deg(
            np.arccos(np.clip(sun_vector_ecef.dot(look), -1.0, 1.0))
        )

        surface_geoid = sk.WGS84()
        surface_geoid.from_lat_lon_alt(
            tangent_geoid.latitude, tangent_geoid.longitude, 0.0
        )
        surface_radii_m[index] = np.linalg.norm(surface_geoid.location)
        viewing.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m=float(tangent_altitude_m[index]),
                relative_azimuth=np.deg2rad(relative_azimuth_deg[index]),
                observer_altitude_m=observer_altitude_m,
                cos_sza=float(cos_sza[index]),
            )
        )

    reference_index = int(np.argmin(np.abs(tangent_altitude_m - 35_000.0)))
    return ImageRayPolicies(
        viewing=viewing,
        earth_radius_m=float(np.mean(surface_radii_m)),
        reference_cos_sza=float(cos_sza[reference_index]),
        tangent_altitude_m=tangent_altitude_m,
        tangent_altitude_residual_m=tangent_altitude_residual_m,
        solar_azimuth_deg=solar_azimuth_deg,
        viewing_azimuth_deg=viewing_azimuth_deg,
        relative_azimuth_deg=relative_azimuth_deg,
        scattering_angle_deg=scattering_angle_deg,
        cos_sza=cos_sza,
    )


def calculate_image(
    data: OmpsInputs,
    selected: SelectedViewing,
    scan_index: int,
    wavelengths_nm: np.ndarray,
    config: sk.Config,
    solar_handler: sk.solar.SolarGeometryHandlerBase,
    ozone_optical: sk.optical.OpticalProperty,
) -> tuple[np.ndarray, ImageRayPolicies, float, float]:
    los_indices = np.flatnonzero(selected.scan_indices == scan_index)
    if len(los_indices) == 0:
        raise ValueError(f"Scan {scan_index} has no selected lines of sight")

    sun, ephemeris_sza_deg, ephemeris_saa_deg = _sun_vector_ecef(
        solar_handler,
        float(data.tangent_latitude_35km_deg[scan_index]),
        float(data.tangent_longitude_35km_deg[scan_index]),
        35_000.0,
        data.times[scan_index],
    )
    policies = image_ray_policies(
        selected.viewing.observer_positions_ecef_m[los_indices],
        selected.viewing.look_directions_ecef[los_indices],
        selected.observed_tangent_altitude_m[los_indices],
        sun,
    )
    geometry = sk.Geometry1D(
        cos_sza=policies.reference_cos_sza,
        solar_azimuth=0.0,
        earth_radius_m=policies.earth_radius_m,
        altitude_grid_m=data.atmosphere_altitude_m,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    pressure_pa = np.ascontiguousarray(data.pressure_pa[scan_index], dtype=np.float64)
    temperature_k = np.ascontiguousarray(
        data.temperature_k[scan_index], dtype=np.float64
    )
    ozone_vmr = np.ascontiguousarray(
        data.ozone_density_cm3[scan_index] / data.air_density_cm3[scan_index],
        dtype=np.float64,
    )
    if config.los_refraction:
        geometry.refractive_index = np.ascontiguousarray(
            sk.optical.refraction.ciddor_index_of_refraction(
                temperature_k,
                pressure_pa,
                np.zeros_like(pressure_pa),
                400.0,
                600.0,
            ),
            dtype=np.float64,
        )

    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=wavelengths_nm,
        calculate_derivatives=False,
    )
    atmosphere.pressure_pa = pressure_pa
    atmosphere.temperature_k = temperature_k
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        ozone_optical,
        data.atmosphere_altitude_m,
        ozone_vmr,
    )
    # Keep the native default solar irradiance of one. The resulting engine
    # output is directly comparable to OMPS SunNormalizedRadiance (L / E0).
    result = sk.Engine(config, geometry, policies.viewing).calculate_radiance(
        atmosphere
    )
    sun_normalized_radiance = result["radiance"].sel(stokes="I", drop=True).values
    return sun_normalized_radiance, policies, ephemeris_sza_deg, ephemeris_saa_deg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l1g", type=Path, default=DEFAULT_L1G)
    parser.add_argument("--anc", type=Path, default=DEFAULT_ANC)
    parser.add_argument("--slit", type=int, default=1)
    parser.add_argument("--scan-start", type=int, default=0)
    parser.add_argument("--num-scans", type=int)
    parser.add_argument("--scan-step", type=int, default=1)
    parser.add_argument("--tangent-altitude-min-km", type=float, default=10.0)
    parser.add_argument("--tangent-altitude-max-km", type=float, default=60.0)
    parser.add_argument("--wavelengths-nm", type=float, nargs="+", default=[350.0])
    parser.add_argument("--no-refraction", action="store_true")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.perf_counter()
    data = load_omps_inputs(args.l1g, args.anc, args.slit)
    if args.num_scans is None:
        scan_indices = np.arange(args.scan_start, len(data.times), args.scan_step)
    else:
        scan_indices = args.scan_start + args.scan_step * np.arange(args.num_scans)
    if np.any(scan_indices < 0) or np.any(scan_indices >= len(data.times)):
        raise ValueError("Requested scan indices lie outside the OMPS orbit")

    requested_wavelengths_nm = np.asarray(args.wavelengths_nm)
    wavelength_indices = np.asarray(
        [
            int(np.argmin(np.abs(data.wavelengths_nm - value)))
            for value in requested_wavelengths_nm
        ]
    )
    if len(np.unique(wavelength_indices)) != len(wavelength_indices):
        raise ValueError("Requested wavelengths map to duplicate OMPS channels")
    wavelengths_nm = data.wavelengths_nm[wavelength_indices]
    selected = select_viewing_geometry(
        data,
        scan_indices,
        args.tangent_altitude_min_km * 1e3,
        args.tangent_altitude_max_km * 1e3,
        wavelength_indices,
        sk.WGS84(),
    )

    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.occultation_source = sk.OccultationSource.NoSource
    config.los_refraction = not args.no_refraction
    solar_handler = sk.solar.SolarGeometryHandlerAstropy()
    ozone_optical = sk.optical.O3DBM()

    num_los = len(selected.scan_indices)
    modeled = np.full((len(wavelengths_nm), num_los), np.nan)
    policy_fields = {
        name: np.full(num_los, np.nan)
        for name in (
            "policy_tangent_altitude_m",
            "tangent_altitude_residual_m",
            "solar_azimuth_deg",
            "viewing_azimuth_deg",
            "relative_azimuth_deg",
            "scattering_angle_deg",
            "cos_sza",
            "image_earth_radius_m",
        )
    }
    ephemeris_sza_deg = np.full(num_los, np.nan)
    ephemeris_saa_deg = np.full(num_los, np.nan)

    for image_number, scan_index in enumerate(scan_indices):
        los_indices = np.flatnonzero(selected.scan_indices == scan_index)
        sun_normalized_radiance, policies, image_sza_deg, image_saa_deg = (
            calculate_image(
                data,
                selected,
                int(scan_index),
                wavelengths_nm,
                config,
                solar_handler,
                ozone_optical,
            )
        )
        modeled[:, los_indices] = sun_normalized_radiance
        policy_fields["policy_tangent_altitude_m"][
            los_indices
        ] = policies.tangent_altitude_m
        policy_fields["tangent_altitude_residual_m"][
            los_indices
        ] = policies.tangent_altitude_residual_m
        policy_fields["solar_azimuth_deg"][los_indices] = policies.solar_azimuth_deg
        policy_fields["viewing_azimuth_deg"][los_indices] = policies.viewing_azimuth_deg
        policy_fields["relative_azimuth_deg"][
            los_indices
        ] = policies.relative_azimuth_deg
        policy_fields["scattering_angle_deg"][
            los_indices
        ] = policies.scattering_angle_deg
        policy_fields["cos_sza"][los_indices] = policies.cos_sza
        policy_fields["image_earth_radius_m"][los_indices] = policies.earth_radius_m
        ephemeris_sza_deg[los_indices] = image_sza_deg
        ephemeris_saa_deg[los_indices] = image_saa_deg
        if image_number % 20 == 0 or image_number + 1 == len(scan_indices):
            print(
                f"Completed image {image_number + 1}/{len(scan_indices)} "
                f"(scan {scan_index})"
            )

    solar_irradiance = sk.solar.SolarModel().irradiance(wavelengths_nm)
    comparison = xr.Dataset(
        {
            "modeled_sun_normalized_radiance_1d": (
                ("wavelength", "los"),
                modeled,
            ),
            "observed_sun_normalized_radiance": (
                ("wavelength", "los"),
                selected.observed_sun_normalized_radiance,
            ),
            "modeled_radiance_1d": (
                ("wavelength", "los"),
                modeled * solar_irradiance[:, np.newaxis],
            ),
            "observed_radiance": (
                ("wavelength", "los"),
                selected.observed_radiance,
            ),
            "generic_solar_irradiance": ("wavelength", solar_irradiance),
            **{name: ("los", value) for name, value in policy_fields.items()},
            "ephemeris_solar_zenith_35km_deg": ("los", ephemeris_sza_deg),
            "ephemeris_solar_azimuth_35km_deg": ("los", ephemeris_saa_deg),
            "omps_solar_zenith_35km_deg": (
                "los",
                data.solar_zenith_35km_deg[selected.scan_indices],
            ),
            "omps_solar_azimuth_35km_deg": (
                "los",
                data.solar_azimuth_35km_deg[selected.scan_indices],
            ),
        },
        coords={
            "wavelength": wavelengths_nm,
            "los": np.arange(num_los),
            "time": ("los", selected.viewing.times),
            "scan_index": ("los", selected.scan_indices),
            "observed_tangent_altitude_m": (
                "los",
                selected.observed_tangent_altitude_m,
            ),
        },
    )
    comparison["sun_normalized_model_minus_observation_1d"] = (
        comparison.modeled_sun_normalized_radiance_1d
        - comparison.observed_sun_normalized_radiance
    )
    comparison["sun_normalized_model_to_observation_ratio_1d"] = (
        comparison.modeled_sun_normalized_radiance_1d
        / comparison.observed_sun_normalized_radiance
    )
    comparison["radiance_model_to_observation_ratio_1d"] = (
        comparison.modeled_radiance_1d / comparison.observed_radiance
    )
    comparison.attrs.update(
        l1g_file=str(args.l1g),
        ancillary_file=str(args.anc),
        modeled_measurement=(
            "Sun-normalized radiance calculated with unit incident solar irradiance"
        ),
        observed_normalized_dataset="GRIDDED_DATA/SunNormalizedRadiance",
        solar_geometry_source="Astropy ephemeris",
        relative_azimuth_convention=(
            "wrap(total solar azimuth - observer-to-scene viewing azimuth)"
        ),
        scattering_angle_convention=(
            "acos(sun_vector_ecef dot observer_to_scene_look_ecef), matching "
            "SASKTRAN2 PhaseHandler propagation directions"
        ),
    )

    center_los = np.asarray(
        [
            np.flatnonzero(selected.scan_indices == scan)[
                np.argmin(
                    np.abs(
                        selected.observed_tangent_altitude_m[
                            selected.scan_indices == scan
                        ]
                        - 35_000.0
                    )
                )
            ]
            for scan in scan_indices
        ]
    )
    ratio = comparison.sun_normalized_model_to_observation_ratio_1d.values
    valid_omps_azimuth = np.abs(comparison.omps_solar_azimuth_35km_deg.values) > 0
    azimuth_delta = _wrap_degrees(
        comparison.ephemeris_solar_azimuth_35km_deg.values
        - comparison.omps_solar_azimuth_35km_deg.values
    )
    print(f"Wavelengths [nm]: {wavelengths_nm.tolist()}")
    print(f"Images: {len(scan_indices)}; LOS: {num_los}")
    print(f"Elapsed: {time.perf_counter() - started:.3f} s")
    print(
        "Sun-normalized model/observation ratio percentiles: "
        f"{np.nanpercentile(ratio, [5, 25, 50, 75, 95]).tolist()}"
    )
    print(
        "Center-ray scattering angle: "
        f"{comparison.scattering_angle_deg.values[center_los[0]]:.3f} deg at "
        f"scan {scan_indices[0]}, "
        f"{comparison.scattering_angle_deg.values[center_los[-1]]:.3f} deg at "
        f"scan {scan_indices[-1]}"
    )
    if np.any(valid_omps_azimuth):
        print(
            "Maximum Astropy-OMPS azimuth difference where OMPS is non-zero: "
            f"{np.nanmax(np.abs(azimuth_delta[valid_omps_azimuth])):.6f} deg"
        )
    else:
        print("No populated OMPS solar azimuths in the selected images")
    print(
        "Maximum tangent-policy altitude residual: "
        f"{np.nanmax(np.abs(policy_fields['tangent_altitude_residual_m'])):.6f} m"
    )
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        comparison.to_netcdf(args.output)
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
