# ruff: noqa: EM101, EM102, T201
"""Run an orbital-plane 2D calculation using a local OMPS LP orbit.

The defaults use a locally materialized 2020 orbit from::

    ~/Documents/data/omps/l1g
    ~/Documents/data/omps/anc

The script intentionally keeps the data adaptation visible. It is an example
and an interface probe, not a retrieval-quality OMPS L1 reader.
"""

from __future__ import annotations

import argparse
import threading
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import psutil
import sasktran2 as sk
import xarray as xr

DEFAULT_L1G = (
    Path.home()
    / "Documents/data/omps/l1g"
    / "OMPS-NPP_LP-L1G-EV_v2.6_2020m0901t000735_o45835_2022m1008t105415.h5"
)
DEFAULT_ANC = (
    Path.home()
    / "Documents/data/omps/anc"
    / "OMPS-NPP_LP-L1-ANC_v2.0_2020m0901t000735_o45835_2022m1008t070135.h5"
)


@dataclass
class OmpsInputs:
    times: np.ndarray
    spacecraft_latitude_deg: np.ndarray
    spacecraft_longitude_deg: np.ndarray
    spacecraft_altitude_m: np.ndarray
    tangent_latitude_35km_deg: np.ndarray
    tangent_longitude_35km_deg: np.ndarray
    solar_zenith_35km_deg: np.ndarray
    solar_azimuth_35km_deg: np.ndarray
    tangent_height_m: np.ndarray
    wavelengths_nm: np.ndarray
    observed_radiance: np.ndarray
    observed_sun_normalized_radiance: np.ndarray
    atmosphere_altitude_m: np.ndarray
    pressure_pa: np.ndarray
    temperature_k: np.ndarray
    air_density_cm3: np.ndarray
    ozone_density_cm3: np.ndarray


@dataclass
class SelectedViewing:
    viewing: sk.OrbitalPlaneViewingGeometry
    scan_indices: np.ndarray
    observed_tangent_altitude_m: np.ndarray
    observed_radiance: np.ndarray
    observed_sun_normalized_radiance: np.ndarray


@dataclass(frozen=True)
class PhaseMeasurement:
    name: str
    elapsed_s: float
    rss_before_bytes: int
    rss_after_bytes: int
    peak_rss_bytes: int


class ResourceProfiler:
    """Measure wall time and native process RSS for coarse script phases."""

    def __init__(self) -> None:
        self._process = psutil.Process()
        self.initial_rss_bytes = self._rss_bytes()
        self.measurements: list[PhaseMeasurement] = []

    def _rss_bytes(self) -> int:
        return int(self._process.memory_info().rss)

    @contextmanager
    def phase(self, name: str) -> Iterator[None]:
        rss_before = self._rss_bytes()
        peak_rss = [rss_before]
        finished = threading.Event()

        def sample_memory() -> None:
            while not finished.wait(0.01):
                peak_rss[0] = max(peak_rss[0], self._rss_bytes())

        sampler = threading.Thread(target=sample_memory, daemon=True)
        sampler.start()
        start = time.perf_counter()
        try:
            yield
        finally:
            elapsed_s = time.perf_counter() - start
            finished.set()
            sampler.join()
            rss_after = self._rss_bytes()
            peak_rss[0] = max(peak_rss[0], rss_after)
            self.measurements.append(
                PhaseMeasurement(
                    name=name,
                    elapsed_s=elapsed_s,
                    rss_before_bytes=rss_before,
                    rss_after_bytes=rss_after,
                    peak_rss_bytes=peak_rss[0],
                )
            )

    def print_summary(self, setup_phase_count: int) -> None:
        mib = 1024**2
        setup = self.measurements[:setup_phase_count]
        calculations = self.measurements[setup_phase_count:]
        setup_elapsed_s = sum(item.elapsed_s for item in setup)
        setup_rss_bytes = setup[-1].rss_after_bytes

        print("Timing and process memory:")
        for item in self.measurements:
            print(
                f"  {item.name:<32} {item.elapsed_s:8.3f} s; "
                f"RSS end={item.rss_after_bytes / mib:8.1f} MiB; "
                f"phase delta={(item.rss_after_bytes - item.rss_before_bytes) / mib:+7.1f} MiB; "
                f"phase peak={item.peak_rss_bytes / mib:8.1f} MiB"
            )
        print(
            f"  {'setup total':<32} {setup_elapsed_s:8.3f} s; "
            f"RSS growth={(setup_rss_bytes - self.initial_rss_bytes) / mib:+.1f} MiB"
        )
        repeated = [
            item for item in calculations if item.name.startswith("repeat calculation")
        ]
        if repeated:
            print(
                f"  {'repeated calculation mean':<32} "
                f"{np.mean([item.elapsed_s for item in repeated]):8.3f} s"
            )


class OmpsSolarHandler:
    """Return the OMPS 35-km solar angles nearest to a group time."""

    def __init__(
        self,
        times: np.ndarray,
        solar_zenith_deg: np.ndarray,
        solar_azimuth_deg: np.ndarray,
    ) -> None:
        self._times_ns = times.astype("datetime64[ns]").astype(np.int64)
        self._solar_zenith_deg = solar_zenith_deg
        self._solar_azimuth_deg = solar_azimuth_deg

    def target_solar_angles(
        self,
        _latitude: float,
        _longitude: float,
        _altitude: float,
        time: pd.Timestamp,
    ) -> tuple[float, float]:
        time_ns = np.datetime64(time.to_datetime64(), "ns").astype(np.int64)
        index = int(np.argmin(np.abs(self._times_ns - time_ns)))
        return (
            float(self._solar_zenith_deg[index]),
            float(self._solar_azimuth_deg[index] % 360.0),
        )


def _decode_utc(values: np.ndarray) -> np.ndarray:
    timestamps = []
    for value in values:
        text = value.decode().rstrip("\x00").removesuffix("Z")
        timestamps.append(np.datetime64(text, "ns"))
    return np.asarray(timestamps)


def load_omps_inputs(l1g_path: Path, anc_path: Path, slit: int) -> OmpsInputs:
    """Load center-slit geometry, measurements, and ancillary profiles."""
    if not l1g_path.exists():
        raise FileNotFoundError(l1g_path)
    if not anc_path.exists():
        raise FileNotFoundError(anc_path)
    if slit not in (0, 1, 2):
        raise ValueError("OMPS slit must be 0, 1, or 2")

    with h5py.File(l1g_path) as l1g:
        geo = l1g["GEOLOCATION_DATA"]
        gridded = l1g["GRIDDED_DATA"]
        times = _decode_utc(gridded["DateTimeUTC"][:])
        spacecraft_latitude_deg = np.asarray(geo["SpacecraftLatitude"][:])
        spacecraft_longitude_deg = np.asarray(geo["SpacecraftLongitude"][:])
        spacecraft_altitude_m = np.asarray(geo["SpacecraftAltitude"][:, slit]) * 1e3
        tangent_latitude_35km_deg = np.asarray(geo["Latitude_35km"][:, slit])
        tangent_longitude_35km_deg = np.asarray(geo["Longitude_35km"][:, slit])
        solar_zenith_35km_deg = np.asarray(geo["SolarZenithAngle_35km"][:, slit])
        solar_azimuth_35km_deg = np.asarray(geo["SolarAzimuth_35km"][:, slit])
        tangent_height_m = np.asarray(gridded["TangentHeight"][:, slit, :]) * 1e3
        wavelengths_nm = np.asarray(gridded["WavelengthGrid"][:]) * 1e3
        observed_radiance = np.asarray(gridded["Radiance"][:, slit, :, :])
        normalized_name = (
            "SunNormalizedRadiance"
            if "SunNormalizedRadiance" in gridded
            else "Reflectance"
        )
        observed_sun_normalized_radiance = np.asarray(
            gridded[normalized_name][:, slit, :, :]
        )

    with h5py.File(anc_path) as anc:
        atmosphere_altitude_m = np.asarray(anc["GEOLOCATION_DATA/Altitude"][:]) * 1e3
        pressure_pa = np.asarray(anc["GRIDDED_DATA/Pressure"][:, slit, :]) * 100.0
        temperature_k = np.asarray(anc["GRIDDED_DATA/Temperature"][:, slit, :])
        air_density_cm3 = np.asarray(anc["GRIDDED_DATA/AirDensity"][:, slit, :])
        ozone_density_cm3 = np.asarray(anc["GRIDDED_DATA/O3Density"][:, slit, :])

    expected_profiles = len(times)
    for name, values in {
        "pressure": pressure_pa,
        "temperature": temperature_k,
        "air density": air_density_cm3,
        "ozone density": ozone_density_cm3,
    }.items():
        if values.shape != (expected_profiles, len(atmosphere_altitude_m)):
            raise ValueError(
                f"ANC {name} shape {values.shape} does not match the L1G orbit"
            )
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError(f"ANC {name} contains invalid or fill values")

    return OmpsInputs(
        times=times,
        spacecraft_latitude_deg=spacecraft_latitude_deg,
        spacecraft_longitude_deg=spacecraft_longitude_deg,
        spacecraft_altitude_m=spacecraft_altitude_m,
        tangent_latitude_35km_deg=tangent_latitude_35km_deg,
        tangent_longitude_35km_deg=tangent_longitude_35km_deg,
        solar_zenith_35km_deg=solar_zenith_35km_deg,
        solar_azimuth_35km_deg=solar_azimuth_35km_deg,
        tangent_height_m=tangent_height_m,
        wavelengths_nm=wavelengths_nm,
        observed_radiance=observed_radiance,
        observed_sun_normalized_radiance=observed_sun_normalized_radiance,
        atmosphere_altitude_m=atmosphere_altitude_m,
        pressure_pa=pressure_pa,
        temperature_k=temperature_k,
        air_density_cm3=air_density_cm3,
        ozone_density_cm3=ozone_density_cm3,
    )


def _ecef_location(
    geoid: sk.Geodetic, latitude_deg: float, longitude_deg: float, altitude_m: float
) -> np.ndarray:
    geoid.from_lat_lon_alt(latitude_deg, longitude_deg, altitude_m)
    return geoid.location.copy()


def select_viewing_geometry(
    data: OmpsInputs,
    scan_indices: np.ndarray,
    tangent_altitude_min_m: float,
    tangent_altitude_max_m: float,
    wavelength_indices: np.ndarray,
    geoid: sk.Geodetic,
) -> SelectedViewing:
    """Reconstruct ECEF rays and sample corresponding measurements."""
    observers: list[np.ndarray] = []
    look_directions: list[np.ndarray] = []
    times: list[np.datetime64] = []
    los_scan_indices: list[int] = []
    observed_altitudes: list[float] = []
    radiances: list[np.ndarray] = []
    sun_normalized_radiances: list[np.ndarray] = []

    ray_geoid = sk.WGS84()
    for scan_index in scan_indices:
        observer = _ecef_location(
            geoid,
            data.spacecraft_latitude_deg[scan_index],
            data.spacecraft_longitude_deg[scan_index],
            data.spacecraft_altitude_m[scan_index],
        )
        tangent_35km = _ecef_location(
            geoid,
            data.tangent_latitude_35km_deg[scan_index],
            data.tangent_longitude_35km_deg[scan_index],
            35_000.0,
        )
        boresight = tangent_35km - observer
        boresight /= np.linalg.norm(boresight)

        vertical_indices = np.flatnonzero(
            (data.tangent_height_m[scan_index] >= tangent_altitude_min_m)
            & (data.tangent_height_m[scan_index] <= tangent_altitude_max_m)
        )
        if len(vertical_indices) == 0:
            raise ValueError(
                f"Scan {scan_index} has no tangent samples between "
                f"{tangent_altitude_min_m / 1e3:g} and "
                f"{tangent_altitude_max_m / 1e3:g} km"
            )

        for vertical_index in vertical_indices:
            tangent_altitude_m = data.tangent_height_m[scan_index, vertical_index]
            look = ray_geoid.from_tangent_altitude(
                float(tangent_altitude_m), observer, boresight
            )
            profile_radiance = data.observed_radiance[scan_index, vertical_index, :]
            measured = profile_radiance[wavelength_indices].astype(np.float64)
            measured[measured <= 0] = np.nan
            profile_sun_normalized = data.observed_sun_normalized_radiance[
                scan_index, vertical_index, :
            ]
            measured_sun_normalized = profile_sun_normalized[wavelength_indices].astype(
                np.float64
            )
            measured_sun_normalized[measured_sun_normalized <= 0] = np.nan

            observers.append(observer)
            look_directions.append(look)
            times.append(data.times[scan_index])
            los_scan_indices.append(int(scan_index))
            observed_altitudes.append(float(tangent_altitude_m))
            radiances.append(measured)
            sun_normalized_radiances.append(measured_sun_normalized)

    viewing = sk.OrbitalPlaneViewingGeometry(
        times,
        np.asarray(observers),
        np.asarray(look_directions),
        vertical_slice=np.asarray(los_scan_indices),
        geoid=geoid,
    )
    return SelectedViewing(
        viewing=viewing,
        scan_indices=np.asarray(los_scan_indices),
        observed_tangent_altitude_m=np.asarray(observed_altitudes),
        observed_radiance=np.asarray(radiances).T,
        observed_sun_normalized_radiance=np.asarray(sun_normalized_radiances).T,
    )


def _surface_directions(data: OmpsInputs, geoid: sk.Geodetic) -> np.ndarray:
    locations = [
        _ecef_location(
            geoid,
            latitude,
            longitude,
            0.0,
        )
        for latitude, longitude in zip(
            data.tangent_latitude_35km_deg,
            data.tangent_longitude_35km_deg,
            strict=True,
        )
    ]
    directions = np.asarray(locations)
    return directions / np.linalg.norm(directions, axis=1, keepdims=True)


def map_ancillary_state(
    data: OmpsInputs,
    geometry: sk.OrbitalPlaneGeometry,
    geoid: sk.Geodetic,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Map ANC profiles along the ordered track without crossing orbit branches."""
    source_directions = _surface_directions(data, geoid)
    target_directions = geometry.ground_track_ecef_m
    target_directions /= np.linalg.norm(target_directions, axis=1, keepdims=True)
    similarities = target_directions @ source_directions.T
    source_rows = np.empty(len(source_directions), dtype=int)
    minimum_row = 0
    for source_index in range(len(source_directions)):
        source_rows[source_index] = minimum_row + int(
            np.argmax(similarities[minimum_row:, source_index])
        )
        minimum_row = source_rows[source_index]
    nearest_scan = np.argmin(
        np.abs(np.arange(len(target_directions))[:, np.newaxis] - source_rows), axis=1
    )

    pressure_pa = data.pressure_pa[nearest_scan]
    temperature_k = data.temperature_k[nearest_scan]
    air_density_cm3 = data.air_density_cm3[nearest_scan]
    ozone_vmr = data.ozone_density_cm3[nearest_scan] / air_density_cm3
    return pressure_pa, temperature_k, ozone_vmr, nearest_scan


def make_atmosphere(
    geometry: sk.OrbitalPlaneGeometry,
    config: sk.Config,
    wavelengths_nm: np.ndarray,
    pressure_pa: np.ndarray,
    temperature_k: np.ndarray,
    ozone_vmr: np.ndarray,
    *,
    calculate_derivatives: bool,
) -> tuple[sk.Atmosphere, sk.constituent.VMRAbsorber2D]:
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=wavelengths_nm,
        calculate_derivatives=calculate_derivatives,
    )
    atmosphere.pressure_pa = pressure_pa
    atmosphere.temperature_k = temperature_k
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    ozone = sk.constituent.VMRAbsorber2D(sk.optical.O3DBM(), ozone_vmr)
    atmosphere["ozone"] = ozone
    # Keep the native default solar irradiance of one. The resulting engine
    # output is directly comparable to OMPS SunNormalizedRadiance (L / E0).
    return atmosphere, ozone


def make_comparison(
    result: xr.Dataset,
    selected: SelectedViewing,
    nearest_anc_scan: np.ndarray,
    engine: sk.OrbitalPlaneEngine,
    l1g_path: Path,
    anc_path: Path,
) -> xr.Dataset:
    modeled_sun_normalized = result["radiance"].sel(stokes="I", drop=True)
    solar_irradiance = xr.DataArray(
        sk.solar.SolarModel().irradiance(modeled_sun_normalized.wavelength.values),
        dims=("wavelength",),
        coords={"wavelength": modeled_sun_normalized.wavelength},
    )
    comparison = xr.Dataset(
        {
            "modeled_sun_normalized_radiance": modeled_sun_normalized,
            "observed_sun_normalized_radiance": (
                ("wavelength", "los"),
                selected.observed_sun_normalized_radiance,
            ),
            "modeled_radiance": modeled_sun_normalized * solar_irradiance,
            "observed_radiance": (
                ("wavelength", "los"),
                selected.observed_radiance,
            ),
            "generic_solar_irradiance": solar_irradiance,
        }
    )
    comparison["sun_normalized_model_minus_observation"] = (
        comparison.modeled_sun_normalized_radiance
        - comparison.observed_sun_normalized_radiance
    )
    comparison["sun_normalized_model_to_observation_ratio"] = (
        comparison.modeled_sun_normalized_radiance
        / comparison.observed_sun_normalized_radiance
    )
    comparison["radiance_model_to_observation_ratio"] = (
        comparison.modeled_radiance / comparison.observed_radiance
    )
    comparison = comparison.assign_coords(
        time=("los", selected.viewing.times),
        scan_index=("los", selected.scan_indices),
        observed_tangent_altitude_m=(
            "los",
            selected.observed_tangent_altitude_m,
        ),
    )

    group_index = np.full(comparison.sizes["los"], -1, dtype=int)
    for diagnostic in engine.group_diagnostics:
        group_index[diagnostic["observation_indices"]] = diagnostic["group_index"]
    comparison = comparison.assign_coords(group_index=("los", group_index))
    comparison.attrs.update(
        l1g_file=str(l1g_path),
        ancillary_file=str(anc_path),
        ancillary_mapping=(
            "Nearest ordered OMPS 35-km tangent slice for every generated orbital row"
        ),
        modeled_measurement=(
            "Sun-normalized radiance calculated with unit incident solar irradiance"
        ),
        observed_normalized_dataset="GRIDDED_DATA/SunNormalizedRadiance",
        num_generated_orbital_positions=len(nearest_anc_scan),
    )
    return comparison


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l1g", type=Path, default=DEFAULT_L1G)
    parser.add_argument("--anc", type=Path, default=DEFAULT_ANC)
    parser.add_argument("--slit", type=int, default=1)
    parser.add_argument("--scan-start", type=int, default=0)
    parser.add_argument(
        "--num-scans",
        type=int,
        help="Number of scans to process; defaults to the rest of the orbit",
    )
    parser.add_argument("--scan-step", type=int, default=1)
    parser.add_argument("--tangent-altitude-min-km", type=float, default=10.0)
    parser.add_argument("--tangent-altitude-max-km", type=float, default=60.0)
    parser.add_argument(
        "--wavelengths-nm",
        type=float,
        nargs="+",
        default=[350.0, 500.0, 675.0],
    )
    parser.add_argument("--along-track-angle-deg", type=float, default=1.0)
    parser.add_argument("--path-padding-deg", type=float, default=10.0)
    parser.add_argument(
        "--group-padding-deg",
        type=float,
        default=10.0,
        help="Extra internal-grid margin before and after every time group",
    )
    parser.add_argument(
        "--time-group-duration-s",
        type=float,
        help=(
            "Time-bin width; defaults to 40 s for single scattering and "
            "240 s for successive orders"
        ),
    )
    parser.add_argument(
        "--multiple-scattering",
        choices=("none", "successive-orders"),
        default="none",
    )
    parser.add_argument(
        "--num-sza",
        type=int,
        default=5,
        help="Number of horizontal source columns for 2D successive orders",
    )
    parser.add_argument(
        "--successive-orders-altitude-points",
        type=int,
        default=25,
        help="Number of midpoint altitude source levels",
    )
    parser.add_argument(
        "--successive-orders-angular-points",
        type=int,
        default=110,
        help="Incoming and outgoing unit-sphere angular points",
    )
    parser.add_argument(
        "--successive-orders-iterations",
        type=int,
        default=50,
        help="Maximum successive-orders fixed-point iterations",
    )
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument(
        "--max-refraction-tangent-altitude-km",
        type=float,
        default=25.0,
        help="Trace higher-tangent-altitude LOS without refraction; use inf for no cutoff",
    )
    parser.add_argument(
        "--repeat-calculations",
        type=int,
        default=2,
        help="Additional optical-state calculations using the resident engine",
    )
    parser.add_argument(
        "--repeat-ozone-scale-step",
        type=float,
        default=1e-3,
        help="Fractional ozone VMR increment applied before each repeated calculation",
    )
    parser.add_argument(
        "--no-vjp",
        action="store_true",
        help="Skip derivative-atmosphere construction, linearization, and VJP",
    )
    parser.add_argument(
        "--derivative-execution",
        choices=("resident", "streaming"),
        default="resident",
        help=(
            "Keep derivative workspaces in every group for retrieval reuse, or "
            "construct one temporary derivative engine at a time"
        ),
    )
    parser.add_argument(
        "--group-details",
        action="store_true",
        help="Print one diagnostic line for every time group",
    )
    parser.add_argument(
        "--solar-source",
        choices=("astropy", "omps"),
        default="astropy",
        help=(
            "Use an Astropy solar ephemeris, or the OMPS 35-km solar-angle "
            "fields. The ephemeris avoids zero OMPS azimuths at dark orbit ends."
        ),
    )
    parser.add_argument(
        "--derivative-details",
        action="store_true",
        help="Print the storage used by every derivative mapping",
    )
    parser.add_argument("--no-refraction", action="store_true")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    profiler = ResourceProfiler()

    with profiler.phase("load OMPS inputs"):
        data = load_omps_inputs(args.l1g, args.anc, args.slit)
        if args.num_scans is None:
            scan_indices = np.arange(args.scan_start, len(data.times), args.scan_step)
        else:
            scan_indices = args.scan_start + args.scan_step * np.arange(args.num_scans)
        if np.any(scan_indices < 0) or np.any(scan_indices >= len(data.times)):
            raise ValueError("Requested scan indices lie outside the OMPS orbit")
        if args.tangent_altitude_min_km > args.tangent_altitude_max_km:
            raise ValueError("Minimum tangent altitude must not exceed the maximum")
        if args.repeat_calculations < 0:
            raise ValueError("repeat-calculations must be non-negative")
        if args.num_sza < 1:
            raise ValueError("num-sza must be positive")
        if args.successive_orders_altitude_points < 2:
            raise ValueError("successive-orders-altitude-points must be at least 2")
        if args.successive_orders_angular_points < 1:
            raise ValueError("successive-orders-angular-points must be positive")
        if args.successive_orders_iterations < 1:
            raise ValueError("successive-orders-iterations must be positive")
        if args.num_threads < 1:
            raise ValueError("num-threads must be positive")
        if (
            np.isnan(args.max_refraction_tangent_altitude_km)
            or args.max_refraction_tangent_altitude_km < 0
        ):
            raise ValueError(
                "max-refraction-tangent-altitude-km must be non-negative or inf"
            )
        if (
            not np.isfinite(args.repeat_ozone_scale_step)
            or args.repeat_ozone_scale_step < 0
        ):
            raise ValueError("repeat-ozone-scale-step must be finite and non-negative")

    with profiler.phase("construct viewing geometry"):
        requested_wavelengths_nm = np.asarray(args.wavelengths_nm)
        wavelength_indices = np.asarray(
            [
                int(np.argmin(np.abs(data.wavelengths_nm - wavelength)))
                for wavelength in requested_wavelengths_nm
            ]
        )
        if len(np.unique(wavelength_indices)) != len(wavelength_indices):
            raise ValueError("Requested wavelengths map to duplicate OMPS channels")
        wavelengths_nm = data.wavelengths_nm[wavelength_indices]

        geoid = sk.WGS84()
        selected = select_viewing_geometry(
            data,
            scan_indices,
            args.tangent_altitude_min_km * 1e3,
            args.tangent_altitude_max_km * 1e3,
            wavelength_indices,
            geoid,
        )

    with profiler.phase("construct atmosphere geometry"):
        geometry = selected.viewing.construct_atmosphere_geometry(
            data.atmosphere_altitude_m,
            np.deg2rad(args.along_track_angle_deg),
            path_padding_angle=np.deg2rad(args.path_padding_deg),
        )

    with profiler.phase("map ancillary state"):
        pressure_pa, temperature_k, ozone_vmr, nearest_anc_scan = map_ancillary_state(
            data, geometry, geoid
        )

    with profiler.phase("construct atmosphere"):
        config = sk.Config()
        config.num_threads = args.num_threads
        config.threading_model = sk.ThreadingModel.Wavelength
        config.single_scatter_source = sk.SingleScatterSource.Exact
        if args.multiple_scattering == "successive-orders":
            config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
            config.num_sza = args.num_sza
            altitude_edges_m = np.linspace(
                data.atmosphere_altitude_m[0],
                data.atmosphere_altitude_m[-1],
                args.successive_orders_altitude_points + 1,
            )
            config.successive_orders_altitude_grid_m = 0.5 * (
                altitude_edges_m[:-1] + altitude_edges_m[1:]
            )
            config.num_successive_orders_incoming = (
                args.successive_orders_angular_points
            )
            config.num_successive_orders_outgoing = (
                args.successive_orders_angular_points
            )
            config.num_successive_orders_iterations = args.successive_orders_iterations
        else:
            config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
        config.occultation_source = sk.OccultationSource.NoSource
        config.los_refraction = not args.no_refraction
        config.los_refraction_max_tangent_altitude_m = (
            args.max_refraction_tangent_altitude_km * 1.0e3
        )

        atmosphere, ozone = make_atmosphere(
            geometry,
            config,
            wavelengths_nm,
            pressure_pa,
            temperature_k,
            ozone_vmr,
            calculate_derivatives=False,
        )

    with profiler.phase("construct orbital engine"):
        time_group_duration_s = args.time_group_duration_s
        if time_group_duration_s is None:
            time_group_duration_s = (
                240.0 if args.multiple_scattering == "successive-orders" else 40.0
            )
        if args.solar_source == "astropy":
            solar_handler = sk.solar.SolarGeometryHandlerAstropy()
        else:
            solar_handler = OmpsSolarHandler(
                data.times,
                data.solar_zenith_35km_deg,
                data.solar_azimuth_35km_deg,
            )
        engine = sk.OrbitalPlaneEngine(
            config,
            geometry,
            selected.viewing,
            time_group_duration_s=time_group_duration_s,
            group_padding_angle=np.deg2rad(args.group_padding_deg),
            solar_handler=solar_handler,
            derivative_execution=args.derivative_execution,
        )

    setup_phase_count = len(profiler.measurements)
    with profiler.phase("first radiance calculation"):
        result = engine.calculate_radiance(atmosphere)
    repeated_changes = []
    for repeat_index in range(args.repeat_calculations):
        ozone.vmr = ozone_vmr * (
            1.0 + args.repeat_ozone_scale_step * (repeat_index + 1)
        )
        with profiler.phase(f"repeat calculation {repeat_index + 1}"):
            repeated_result = engine.calculate_radiance(atmosphere)
        repeated_changes.append(
            float(
                np.max(
                    np.abs(
                        repeated_result["radiance"].values - result["radiance"].values
                    )
                )
            )
        )

    linearization = None
    vjp = None
    if not args.no_vjp:
        with profiler.phase("construct derivative atmosphere"):
            derivative_atmosphere, _ = make_atmosphere(
                geometry,
                config,
                wavelengths_nm,
                pressure_pa,
                temperature_k,
                ozone_vmr,
                calculate_derivatives=True,
            )
            derivative_atmosphere.internal_object()
        if args.derivative_details:
            print("Derivative mapping storage:")
            for name in derivative_atmosphere.storage.derivative_mapping_names():
                mapping = derivative_atmosphere.storage.get_derivative_mapping(name)
                interpolator = np.asarray(mapping.interpolator)
                print(
                    f"  {name}: output={mapping.assign_name or name}, "
                    f"interpolator_shape={interpolator.shape}, "
                    f"interpolator={interpolator.nbytes / 1024**2:.1f} MiB"
                )
        with profiler.phase("create linearization"):
            linearization = engine.linearize(
                derivative_atmosphere, prepare_parameters=("ozone_vmr",)
            )
        cotangent = xr.ones_like(linearization.value)
        with profiler.phase("one VJP"):
            vjp = linearization.vjp(cotangent, parameters=("ozone_vmr",))

    comparison = make_comparison(
        result,
        selected,
        nearest_anc_scan,
        engine,
        args.l1g,
        args.anc,
    )

    print(f"L1G: {args.l1g}")
    print(f"ANC: {args.anc}")
    print(
        f"Slit: {args.slit}; {len(scan_indices)} scans "
        f"({scan_indices[0]} through {scan_indices[-1]}); "
        f"{len(selected.scan_indices)} LOS between "
        f"{args.tangent_altitude_min_km:g} and "
        f"{args.tangent_altitude_max_km:g} km"
    )
    print(f"OMPS wavelengths [nm]: {wavelengths_nm.tolist()}")
    print(f"Solar geometry source: {args.solar_source}")
    print(
        f"Multiple scattering: {args.multiple_scattering}; "
        f"time groups={time_group_duration_s:g} s"
    )
    if args.multiple_scattering == "successive-orders":
        print(
            "Successive orders: "
            f"num_sza={args.num_sza}; "
            f"source altitudes={args.successive_orders_altitude_points}; "
            f"angular points={args.successive_orders_angular_points}; "
            f"maximum iterations={args.successive_orders_iterations}; "
            f"threads={args.num_threads}; max refracted tangent altitude="
            f"{args.max_refraction_tangent_altitude_km:g} km"
        )
    print(
        "Geometry: "
        f"{geometry.shape[0]} orbital positions x {geometry.shape[1]} altitudes; "
        f"{np.rad2deg(geometry.cumulative_angles[-1]):.3f} deg track span; "
        f"{args.along_track_angle_deg:g} deg maximum spacing; "
        f"{args.path_padding_deg:g} deg master padding per end; "
        f"{args.group_padding_deg:g} deg group margin; "
        f"{len(engine.group_diagnostics)} time groups"
    )
    profiler.print_summary(setup_phase_count)
    if repeated_changes:
        print(
            "  Repeats changed ozone VMR by "
            f"{args.repeat_ozone_scale_step:.3g} per iteration; "
            "pressure, temperature, and refracted paths were unchanged."
        )
        print(
            "  Maximum sun-normalized radiance changes from the initial state: "
            f"{repeated_changes}"
        )
    if linearization is not None and vjp is not None:
        print(f"  VJP backend: {linearization.backends['vjp'].value}")
        print(f"  VJP parameter shapes: { {name: vjp[name].shape for name in vjp} }")
    if args.group_details:
        for diagnostic in engine.group_diagnostics:
            print(
                f"  group {diagnostic['group_index']}: "
                f"{len(diagnostic['observation_indices'])} LOS, "
                f"{len(diagnostic['grid_indices'])} atmosphere rows, "
                f"earth radius={diagnostic['earth_radius_m'] / 1e3:.3f} km, "
                "plane residual="
                f"{np.rad2deg(diagnostic['max_out_of_plane_angle']):.4f} deg, "
                "horizontal scale residual="
                f"{100 * diagnostic['maximum_relative_horizontal_scale_residual']:.4f}%, "
                f"edge clipping={diagnostic['edge_clipping']}"
            )
    for wavelength in wavelengths_nm:
        observed = comparison.observed_sun_normalized_radiance.sel(
            wavelength=wavelength
        )
        modeled = comparison.modeled_sun_normalized_radiance.sel(wavelength=wavelength)
        print(
            f"  {wavelength:.3f} nm: "
            f"mean observed={float(observed.mean(skipna=True)):.6e}, "
            f"mean modeled={float(modeled.mean(skipna=True)):.6e} "
            "sun-normalized radiance"
        )

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        comparison.to_netcdf(args.output)
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
