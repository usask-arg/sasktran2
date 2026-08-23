from __future__ import annotations

import operator
from collections.abc import Iterable, Sequence
from typing import Literal

import numpy as np
import pandas as pd
import xarray as xr

import sasktran2 as sk
from sasktran2._core_rust import (
    PyOrbitalPlaneEngine,
    PyOrbitalPlaneGeometry,
    PyOrbitalPlaneViewingGeometry,
    orbital_group_reference_data,
)
from sasktran2.engine import Engine
from sasktran2.geometry import Geometry2D
from sasktran2.linearization import _selected_parameters
from sasktran2.viewinggeo.base import ViewingGeometryContainer

_DEFAULT_MAX_EAGER_JACOBIAN_BYTES = 2 * 1024**3


def _readonly_view(value: np.ndarray) -> np.ndarray:
    """Return a view whose private read-only base prevents flag escalation."""
    result = value.view()
    result.setflags(write=False)
    return result


def _times_ns_array(times: Sequence) -> np.ndarray:
    """Return a validated, owned one-dimensional nanosecond time array."""
    try:
        result = np.asarray(times, dtype="datetime64[ns]")
    except (TypeError, ValueError) as error:
        msg = "times must be convertible to a one-dimensional datetime64 array"
        raise ValueError(msg) from error
    if result.ndim != 1 or result.size == 0:
        msg = (
            f"times must have shape (los,) with at least one sample; got {result.shape}"
        )
        raise ValueError(msg)
    if np.any(np.isnat(result)):
        msg = "times must not contain NaT"
        raise ValueError(msg)
    return np.array(result, dtype="datetime64[ns]", order="C", copy=True)


def _ecef_array(value: np.ndarray, name: str, num_los: int) -> np.ndarray:
    """Return validated, owned ECEF vectors with one row per LOS."""
    try:
        result = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as error:
        msg = f"{name} must be convertible to floating-point ECEF vectors"
        raise ValueError(msg) from error
    expected = (num_los, 3)
    if result.shape != expected:
        msg = f"{name} must have shape {expected}; got {result.shape}"
        raise ValueError(msg)
    if np.any(~np.isfinite(result)):
        msg = f"{name} must contain only finite values"
        raise ValueError(msg)
    norms = np.linalg.norm(result, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms == 0):
        msg = f"{name} must contain non-zero vectors"
        raise ValueError(msg)
    return np.array(result, dtype=np.float64, order="C", copy=True)


def _vertical_slice_array(
    value: Sequence, num_los: int
) -> tuple[np.ndarray, np.ndarray]:
    """Validate slice labels and return labels plus contiguous native codes."""
    labels = np.asarray(value)
    if labels.ndim != 1 or len(labels) != num_los:
        msg = f"vertical_slice must have shape ({num_los},); got {labels.shape}"
        raise ValueError(msg)
    native_labels = np.empty(num_los, dtype=np.int64)
    for index, label in enumerate(labels.tolist()):
        if isinstance(label, bool | np.bool_):
            msg = "vertical_slice must contain non-negative integer identifiers"
            raise TypeError(msg)
        try:
            native_labels[index] = operator.index(label)
        except (TypeError, ValueError, OverflowError) as error:
            msg = "vertical_slice must contain non-negative integer identifiers"
            raise TypeError(msg) from error
    if np.any(native_labels < 0):
        msg = "vertical_slice identifiers must be non-negative"
        raise ValueError(msg)

    # A vertical slice is an atomic image/scan in the flattened LOS order. A
    # label that reappears later would make complete-slice chunking ambiguous.
    closed: set[int] = set()
    previous = int(native_labels[0])
    for label in native_labels[1:]:
        current = int(label)
        if current != previous:
            closed.add(previous)
            if current in closed:
                msg = (
                    "Each vertical_slice identifier must occupy one contiguous "
                    "block of LOS"
                )
                raise ValueError(msg)
            previous = current

    code_by_label: dict[int, int] = {}
    codes = np.empty(num_los, dtype=np.int64)
    for index, label in enumerate(native_labels):
        codes[index] = code_by_label.setdefault(int(label), len(code_by_label))
    return native_labels, codes


def _structural_config_signature(config: sk.Config) -> tuple:
    """Configuration fields cached by native group-engine construction."""
    fields = (
        "num_threads",
        "wavelength_batch_size",
        "threading_model",
        "num_stokes",
        "num_streams",
        "single_scatter_source",
        "occultation_source",
        "emission_source",
        "multiple_scatter_source",
        "multiple_scatter_refraction",
        "stokes_basis",
        "delta_m_scaling",
        "los_refraction",
        "los_refraction_max_tangent_altitude_m",
        "solar_refraction",
        "output_los_optical_depth",
        "num_singlescatter_moments",
        "num_sza",
        "num_successive_orders_iterations",
        "successive_orders_relative_tolerance",
        "successive_orders_absolute_tolerance",
        "successive_orders_anderson_depth",
        "successive_orders_damping",
        "successive_orders_altitude_grid_m",
        "successive_orders_horizontal_angle_grid_radians",
        "init_successive_orders_with_discrete_ordinates",
        "num_successive_order_points",
        "num_successive_orders_incoming",
        "num_successive_orders_outgoing",
        "spectral_grid_mode",
    )
    values = []
    for field in fields:
        value = getattr(config, field)
        if isinstance(value, np.ndarray):
            value = tuple(value.tolist())
        values.append(value)
    return tuple(values)


class OrbitalPlaneGeometry(Geometry2D):
    """A full ordered ground-track grid used by rolling local 2D engines.

    Unlike :class:`Geometry2D`, the master grid may cover an entire orbit. Each
    calculation group selects a contiguous window whose projected span is less
    than pi.
    """

    _geometry: PyOrbitalPlaneGeometry

    def __init__(
        self,
        earth_radius_m: float,
        altitude_grid_m: np.ndarray,
        ground_track_ecef_m: np.ndarray,
        interpolation_method: sk.InterpolationMethod = sk.InterpolationMethod.LinearInterpolation,
        *,
        surface_radii_m: np.ndarray | None = None,
    ) -> None:
        """Construct an ordered master atmosphere grid.

        ``ground_track_ecef_m`` supplies directions and does not need to lie
        exactly on the requested spherical surface. It is normalized and then
        scaled by ``surface_radii_m`` when provided, otherwise by
        ``earth_radius_m``.
        """
        if isinstance(earth_radius_m, bool | np.bool_):
            msg = "earth_radius_m must be a finite positive number"
            raise TypeError(msg)
        try:
            earth_radius_m = float(earth_radius_m)
        except (TypeError, ValueError, OverflowError) as error:
            msg = "earth_radius_m must be a finite positive number"
            raise TypeError(msg) from error
        altitude_grid_m = np.asarray(altitude_grid_m, dtype=np.float64)
        if altitude_grid_m.ndim != 1:
            msg = f"altitude_grid_m must have shape (altitude,); got {altitude_grid_m.shape}"
            raise ValueError(msg)
        ground_track_ecef_m = np.asarray(ground_track_ecef_m, dtype=np.float64)
        if ground_track_ecef_m.ndim != 2 or ground_track_ecef_m.shape[1:] != (3,):
            msg = (
                "ground_track_ecef_m must have shape (orbital_position, 3); "
                f"got {ground_track_ecef_m.shape}"
            )
            raise ValueError(msg)
        if surface_radii_m is not None:
            surface_radii_m = np.asarray(surface_radii_m, dtype=np.float64)
            if surface_radii_m.ndim != 1:
                msg = (
                    "surface_radii_m must have shape (orbital_position,); "
                    f"got {surface_radii_m.shape}"
                )
                raise ValueError(msg)
            surface_radii_m = np.ascontiguousarray(surface_radii_m)
        self._geometry = PyOrbitalPlaneGeometry(
            earth_radius_m,
            np.ascontiguousarray(altitude_grid_m),
            np.ascontiguousarray(ground_track_ecef_m),
            interpolation_method,
            surface_radii_m,
        )
        self._anchor_vertical_slice_labels: np.ndarray | None = None

    def altitudes(self) -> np.ndarray:
        return self._geometry.altitudes()

    def horizontal_angles(self) -> np.ndarray:
        """Return the cumulative master along-track angle."""
        return self._geometry.cumulative_angles()

    @property
    def orbital_positions(self) -> np.ndarray:
        return np.arange(self.shape[0])

    @property
    def ground_track_ecef_m(self) -> np.ndarray:
        return self._geometry.ground_track_ecef()

    @property
    def cumulative_angles(self) -> np.ndarray:
        return self._geometry.cumulative_angles()

    @property
    def geodetic_latitude_deg(self) -> np.ndarray | None:
        """Grid geodetic latitude, or ``None`` for an explicit spherical grid."""
        coordinates = np.asarray(self._geometry.geodetic_coordinates_degrees())
        if coordinates.size == 0:
            return None
        return coordinates[:, 0].copy()

    @property
    def geodetic_longitude_deg(self) -> np.ndarray | None:
        """Grid geodetic longitude, or ``None`` for an explicit spherical grid."""
        coordinates = np.asarray(self._geometry.geodetic_coordinates_degrees())
        if coordinates.size == 0:
            return None
        return coordinates[:, 1].copy()

    @property
    def reference_times(self) -> np.ndarray | None:
        """Slice-time interpolation at every atmosphere-grid position.

        Values are clamped to the first and last vertical-slice mean timestamp
        through the path-padding portions of a viewing-constructed grid.
        Explicitly constructed spherical grids have no time provenance.
        """
        values = np.asarray(self._geometry.grid_reference_times_ns(), dtype=np.int64)
        if values.size == 0:
            return None
        return values.astype("datetime64[ns]")

    @property
    def vertical_slice_anchors(self) -> xr.Dataset:
        """Construction anchors linking vertical slices to the master track."""
        codes = np.asarray(
            self._geometry.anchor_vertical_slice_indices(), dtype=np.int64
        )
        if codes.size == 0:
            return xr.Dataset(coords={"vertical_slice_anchor": np.arange(0)})
        labels = (
            codes
            if self._anchor_vertical_slice_labels is None
            else self._anchor_vertical_slice_labels
        )
        return xr.Dataset(
            coords={
                "vertical_slice_anchor": np.arange(codes.size),
                "xyz": ["x", "y", "z"],
                "vertical_slice": ("vertical_slice_anchor", labels.copy()),
                "reference_time": (
                    "vertical_slice_anchor",
                    np.asarray(
                        self._geometry.anchor_reference_times_ns(), dtype=np.int64
                    ).astype("datetime64[ns]"),
                ),
                "along_track_angle_rad": (
                    "vertical_slice_anchor",
                    np.asarray(self._geometry.anchor_along_track_angles()).copy(),
                ),
            },
            data_vars={
                "ground_track_ecef_m": (
                    ("vertical_slice_anchor", "xyz"),
                    np.asarray(self._geometry.anchor_ground_track_ecef()).copy(),
                )
            },
        )

    def _xarray_orbital_coordinates(self) -> dict:
        coordinates: dict = {
            "orbital_position": self.orbital_positions,
            "along_track_angle_rad": (
                "orbital_position",
                self.cumulative_angles.copy(),
            ),
        }
        latitude = self.geodetic_latitude_deg
        longitude = self.geodetic_longitude_deg
        reference_time = self.reference_times
        if latitude is not None:
            coordinates["geodetic_latitude_deg"] = (
                "orbital_position",
                latitude,
            )
            coordinates["geodetic_longitude_deg"] = (
                "orbital_position",
                longitude,
            )
        if reference_time is not None:
            coordinates["reference_time"] = (
                "orbital_position",
                reference_time,
            )
        return coordinates

    @property
    def grid_dataset(self) -> xr.Dataset:
        """Labeled master atmosphere grid and its construction provenance."""
        coordinates = self._xarray_orbital_coordinates()
        coordinates.update(
            {
                "altitude": self.altitudes().copy(),
                "xyz": ["x", "y", "z"],
            }
        )
        return xr.Dataset(
            coords=coordinates,
            data_vars={
                "ground_track_ecef_m": (
                    ("orbital_position", "xyz"),
                    self.ground_track_ecef_m.copy(),
                ),
                "surface_radius_m": (
                    "orbital_position",
                    self.surface_radii_m.copy(),
                ),
            },
        )

    def project_ground_track(
        self,
        ecef_positions_m: np.ndarray,
        *,
        order: Literal["independent", "increasing", "decreasing"] = "independent",
    ) -> xr.Dataset:
        """Project ECEF positions onto the ordered master ground track.

        ``order="increasing"`` or ``"decreasing"`` constrains every search to
        the remaining portion of the track. This disambiguates repeated ECEF
        locations on full-orbit or self-crossing grids and is the appropriate
        mode for an ancillary product already ordered along the orbit.
        """
        try:
            positions = np.asarray(ecef_positions_m, dtype=np.float64)
        except (TypeError, ValueError) as error:
            msg = "ecef_positions_m must be floating-point ECEF vectors"
            raise ValueError(msg) from error
        if positions.ndim == 1:
            positions = positions[np.newaxis, :]
        if positions.ndim != 2 or positions.shape[1:] != (3,):
            msg = (
                "ecef_positions_m must have shape (position, 3) or (3,); "
                f"got {positions.shape}"
            )
            raise ValueError(msg)
        if np.any(~np.isfinite(positions)) or np.any(
            np.linalg.norm(positions, axis=1) == 0
        ):
            msg = "ecef_positions_m must contain finite, non-zero vectors"
            raise ValueError(msg)
        if order not in ("independent", "increasing", "decreasing"):
            msg = "order must be 'independent', 'increasing', or 'decreasing'"
            raise ValueError(msg)
        segment, fraction, angle, residual, at_edge = (
            self._geometry.project_ground_track(np.ascontiguousarray(positions), order)
        )
        return xr.Dataset(
            coords={"position": np.arange(len(positions))},
            data_vars={
                "segment_index": ("position", np.asarray(segment)),
                "segment_fraction": ("position", np.asarray(fraction)),
                "along_track_angle_rad": ("position", np.asarray(angle)),
                "cross_track_angle_rad": ("position", np.asarray(residual)),
                "at_track_edge": ("position", np.asarray(at_edge)),
            },
        )

    @property
    def earth_radius_m(self) -> float:
        """Mean geocentric surface radius of the master orbital grid."""
        return self._geometry.earth_radius_m

    @property
    def surface_radii_m(self) -> np.ndarray:
        """Geocentric geoid radius at each orbital position."""
        return self._geometry.surface_radii()

    @property
    def shape(self) -> tuple[int, int]:
        return self._geometry.location_shape()

    @property
    def horizontal_dimension(self) -> str:
        return "orbital_position"


class OrbitalPlaneViewingGeometry(ViewingGeometryContainer):
    """Timestamped ECEF limb measurements for an orbital-plane calculation.

    Parameters
    ----------
    times
        One timestamp per line of sight. Values are converted to nanosecond
        ``datetime64`` and must not contain ``NaT``.
    observer_positions_ecef_m
        Finite, non-zero observer positions with shape ``(los, 3)``.
    look_directions_ecef
        Finite, non-zero look vectors with shape ``(los, 3)``. They are copied
        and normalized on construction.
    vertical_slice
        Required non-negative integer image/scan identifier for every LOS.
        Identical identifiers must form one contiguous block. Track construction
        and time grouping treat each vertical slice atomically while retaining
        every LOS for ray placement and coverage.
    geoid
        Geoid used by :meth:`construct_atmosphere_geometry`; defaults to WGS84.
    """

    _viewing_geometry: PyOrbitalPlaneViewingGeometry

    def __init__(
        self,
        times: Sequence,
        observer_positions_ecef_m: np.ndarray,
        look_directions_ecef: np.ndarray,
        *,
        vertical_slice: Sequence,
        geoid: sk.Geodetic | None = None,
    ) -> None:
        self._initialize(
            times,
            observer_positions_ecef_m,
            look_directions_ecef,
            vertical_slice,
            geoid,
            None,
        )

    def _initialize(
        self,
        times: Sequence,
        observer_positions_ecef_m: np.ndarray,
        look_directions_ecef: np.ndarray,
        vertical_slice: Sequence,
        geoid: sk.Geodetic | None,
        time_bin_origin_ns: int | None,
    ) -> None:
        times = _times_ns_array(times)
        observers = _ecef_array(
            observer_positions_ecef_m, "observer_positions_ecef_m", len(times)
        )
        looks = _ecef_array(look_directions_ecef, "look_directions_ecef", len(times))
        looks /= np.linalg.norm(looks, axis=1)[:, np.newaxis]
        vertical_slice, vertical_slice_codes = _vertical_slice_array(
            vertical_slice, len(times)
        )
        earliest_time_ns = int(times.astype(np.int64).min())
        if time_bin_origin_ns is None:
            time_bin_origin_ns = earliest_time_ns
        else:
            try:
                time_bin_origin_ns = operator.index(time_bin_origin_ns)
            except TypeError as error:
                msg = "time_bin_origin_ns must be an integer nanosecond timestamp"
                raise TypeError(msg) from error
            if time_bin_origin_ns < np.iinfo(np.int64).min + 1:
                msg = "time_bin_origin_ns must be a valid integer nanosecond timestamp"
                raise ValueError(msg)
            if time_bin_origin_ns > earliest_time_ns:
                msg = (
                    "time_bin_origin_ns must not be later than the earliest observation"
                )
                raise ValueError(msg)
        self._geoid = sk.WGS84() if geoid is None else geoid
        if not isinstance(self._geoid, sk.Geodetic):
            msg = "geoid must be a sasktran2.Geodetic object"
            raise TypeError(msg)
        self._viewing_geometry = PyOrbitalPlaneViewingGeometry(
            np.ascontiguousarray(times.astype(np.int64)),
            observers,
            looks,
            vertical_slice_codes,
            time_bin_origin_ns,
        )
        geometry_ds = xr.Dataset(
            coords={
                "time": ("los", times.copy()),
                "vertical_slice": ("los", vertical_slice.copy()),
                "xyz": ["x", "y", "z"],
            },
            data_vars={
                "observer_position_ecef_m": (("los", "xyz"), observers.copy()),
                "look_direction_ecef": (("los", "xyz"), looks.copy()),
            },
        )
        times.setflags(write=False)
        observers.setflags(write=False)
        looks.setflags(write=False)
        vertical_slice.setflags(write=False)
        vertical_slice_codes.setflags(write=False)
        self._times = times
        self._observer_positions_ecef_m = observers
        self._look_directions_ecef = looks
        self._vertical_slice = vertical_slice
        self._vertical_slice_codes = vertical_slice_codes
        self._time_bin_origin_ns = time_bin_origin_ns
        super().__init__(geometry_ds)

    @classmethod
    def from_tangent_locations(
        cls,
        times: Sequence,
        observer_positions_ecef_m: np.ndarray,
        tangent_locations_ecef_m: np.ndarray,
        *,
        vertical_slice: Sequence,
        geoid: sk.Geodetic | None = None,
    ) -> OrbitalPlaneViewingGeometry:
        """Construct grouped unit look vectors from observer and tangent points."""
        times = _times_ns_array(times)
        observers = _ecef_array(
            observer_positions_ecef_m, "observer_positions_ecef_m", len(times)
        )
        tangents = _ecef_array(
            tangent_locations_ecef_m, "tangent_locations_ecef_m", len(times)
        )
        looks = tangents - observers
        norms = np.linalg.norm(looks, axis=1)
        if np.any(~np.isfinite(norms)) or np.any(norms == 0):
            msg = "observer and tangent locations must define finite, non-zero look directions"
            raise ValueError(msg)
        return cls(
            times,
            observers,
            looks / norms[:, np.newaxis],
            vertical_slice=vertical_slice,
            geoid=geoid,
        )

    def construct_atmosphere_geometry(
        self,
        altitude_grid_m: np.ndarray,
        along_track_angle_delta: float,
        interpolation_method: sk.InterpolationMethod = sk.InterpolationMethod.LinearInterpolation,
        *,
        path_padding_angle: float = np.deg2rad(10.0),
        geoid: sk.Geodetic | None = None,
        max_orbital_positions: int | None = 100_000,
    ) -> OrbitalPlaneGeometry:
        """Construct a geoid-following atmosphere grid from the viewing track.

        The spherical-mean tangent location of each time-ordered vertical slice
        is projected to the geoid and the resulting piecewise great-circle track
        is resampled. Individual LOS set the boundary-slice coverage but their
        order within a slice does not affect the grid. Angular arguments are in
        radians. End points are extended by ``path_padding_angle`` beyond the
        full observed boundary-slice extent.

        Parameters
        ----------
        altitude_grid_m
            Strictly increasing atmospheric altitudes above the geoid.
        along_track_angle_delta
            Maximum angular separation of adjacent orbital positions, in
            radians.
        interpolation_method
            Vertical interpolation used by each local ``Geometry2D`` engine.
        path_padding_angle
            Along-track angular padding placed before the first and after the
            last nominal tangent location, in radians. The default is 10
            degrees. A group that requires atmosphere beyond this margin uses
            the existing clamped edge behavior. ``group_diagnostics`` reports
            master-grid truncation as ``edge_clipping`` and actual traced LOS
            use of an extended edge cell as ``path_edge_clipping``.
        geoid
            Optional override for the geoid stored on this viewing geometry.
        max_orbital_positions
            Safety cap on the generated along-track grid size. The default is
            100,000 positions. Set to ``None`` to explicitly disable the cap.
        """
        selected_geoid = self._geoid if geoid is None else geoid
        if not isinstance(selected_geoid, sk.Geodetic):
            msg = "geoid must be a sasktran2.Geodetic object"
            raise TypeError(msg)
        if max_orbital_positions is not None:
            if isinstance(max_orbital_positions, bool | np.bool_):
                msg = "max_orbital_positions must be a positive integer or None"
                raise TypeError(msg)
            try:
                max_orbital_positions = operator.index(max_orbital_positions)
            except TypeError as error:
                msg = "max_orbital_positions must be a positive integer or None"
                raise TypeError(msg) from error
            if max_orbital_positions < 2:
                msg = "max_orbital_positions must be at least two when supplied"
                raise ValueError(msg)
        altitude_grid_m = np.asarray(altitude_grid_m, dtype=np.float64)
        if altitude_grid_m.ndim != 1:
            msg = f"altitude_grid_m must have shape (altitude,); got {altitude_grid_m.shape}"
            raise ValueError(msg)
        if isinstance(along_track_angle_delta, bool | np.bool_):
            msg = "along_track_angle_delta must be a finite angle in radians"
            raise TypeError(msg)
        try:
            along_track_angle_delta = float(along_track_angle_delta)
        except (TypeError, ValueError, OverflowError) as error:
            msg = "along_track_angle_delta must be a finite angle in radians"
            raise TypeError(msg) from error
        if isinstance(path_padding_angle, bool | np.bool_):
            msg = "path_padding_angle must be a finite angle in radians"
            raise TypeError(msg)
        try:
            path_padding_angle = float(path_padding_angle)
        except (TypeError, ValueError, OverflowError) as error:
            msg = "path_padding_angle must be a finite angle in radians"
            raise TypeError(msg) from error
        result = OrbitalPlaneGeometry.__new__(OrbitalPlaneGeometry)
        result._geometry = self._viewing_geometry.construct_atmosphere_geometry(
            np.ascontiguousarray(altitude_grid_m),
            along_track_angle_delta,
            path_padding_angle,
            interpolation_method,
            max_orbital_positions,
            selected_geoid._internal,
        )
        anchor_codes = np.asarray(
            result._geometry.anchor_vertical_slice_indices(), dtype=np.int64
        )
        label_by_code = {
            int(code): int(label)
            for code, label in zip(
                self._vertical_slice_codes, self._vertical_slice, strict=True
            )
        }
        result._anchor_vertical_slice_labels = np.asarray(
            [label_by_code[int(code)] for code in anchor_codes], dtype=np.int64
        )
        return result

    @property
    def times(self) -> np.ndarray:
        """Observation timestamps as a read-only ``datetime64[ns]`` array."""
        return _readonly_view(self._times)

    @property
    def observer_positions_ecef_m(self) -> np.ndarray:
        """Observer positions as a read-only ``(los, 3)`` array."""
        return _readonly_view(self._observer_positions_ecef_m)

    @property
    def look_directions_ecef(self) -> np.ndarray:
        """Normalized look directions as a read-only ``(los, 3)`` array."""
        return _readonly_view(self._look_directions_ecef)

    @property
    def vertical_slice(self) -> np.ndarray:
        """Read-only vertical image/scan identifier for every LOS."""
        return _readonly_view(self._vertical_slice)

    @property
    def num_vertical_slices(self) -> int:
        """Number of complete vertical images/scans."""
        return int(self._vertical_slice_codes[-1]) + 1

    @property
    def geoid(self) -> sk.Geodetic:
        """Geoid used when constructing the atmosphere geometry."""
        return self._geoid

    @property
    def flux_observers(self) -> tuple:
        return ()

    def __len__(self) -> int:
        return len(self._times)

    def isel(
        self, los: int | slice | Sequence[int] | np.ndarray
    ) -> OrbitalPlaneViewingGeometry:
        """Return a positional LOS subset as an independent viewing geometry."""
        if isinstance(los, bool | np.bool_):
            msg = "los must be an integer, slice, integer indices, or a boolean mask"
            raise TypeError(msg)
        try:
            selected = np.arange(len(self))[los]
        except (IndexError, TypeError, ValueError) as error:
            msg = "los must be valid positional indices for this viewing geometry"
            raise ValueError(msg) from error
        selected = np.atleast_1d(selected)
        if selected.ndim != 1 or selected.size == 0:
            msg = "An orbital viewing-geometry subset must contain at least one LOS"
            raise ValueError(msg)
        result = type(self).__new__(type(self))
        result._initialize(
            self._times[selected],
            self._observer_positions_ecef_m[selected],
            self._look_directions_ecef[selected],
            self._vertical_slice[selected],
            self._geoid,
            self._time_bin_origin_ns,
        )
        return result

    def split(
        self, num_chunks: int, *, time_group_duration_s: float
    ) -> tuple[OrbitalPlaneViewingGeometry, ...]:
        """Split complete native time groups into contiguous retrieval chunks.

        The duration must match the ``OrbitalPlaneEngine`` duration. This keeps
        every vertical slice and resulting native group in one chunk and retains
        the full viewing geometry's time-bin origin.
        """
        if isinstance(num_chunks, bool | np.bool_):
            msg = "num_chunks must be a positive integer"
            raise TypeError(msg)
        try:
            num_chunks = operator.index(num_chunks)
        except TypeError as error:
            msg = "num_chunks must be a positive integer"
            raise TypeError(msg) from error
        if isinstance(time_group_duration_s, bool | np.bool_):
            msg = "time_group_duration_s must be a finite positive number"
            raise TypeError(msg)
        try:
            duration_ns_float = float(time_group_duration_s) * 1e9
        except (TypeError, ValueError, OverflowError) as error:
            msg = "time_group_duration_s must be a finite positive number"
            raise TypeError(msg) from error
        if (
            not np.isfinite(duration_ns_float)
            or duration_ns_float < 0.5
            or duration_ns_float > np.iinfo(np.int64).max
        ):
            msg = "time_group_duration_s must be finite, positive, and at least 1 ns"
            raise ValueError(msg)
        duration_ns = int(np.floor(duration_ns_float + 0.5))
        slice_runs = np.flatnonzero(
            np.r_[True, np.diff(self._vertical_slice_codes) != 0, True]
        )
        times_ns = self._times.astype(np.int64)
        slice_bins = []
        for slice_index in range(self.num_vertical_slices):
            start = slice_runs[slice_index]
            stop = slice_runs[slice_index + 1]
            reference_time_ns = sum(int(value) for value in times_ns[start:stop]) // (
                stop - start
            )
            slice_bins.append(
                (reference_time_ns - self._time_bin_origin_ns) // duration_ns
            )
        if np.any(np.diff(slice_bins) < 0):
            msg = "Vertical slices must be ordered by non-decreasing reference time to split"
            raise ValueError(msg)
        group_slice_runs = np.flatnonzero(np.r_[True, np.diff(slice_bins) != 0, True])
        num_groups = len(group_slice_runs) - 1
        if num_chunks < 1 or num_chunks > num_groups:
            msg = f"num_chunks must be between 1 and {num_groups}; got {num_chunks}"
            raise ValueError(msg)
        chunks = []
        for group_indices in np.array_split(np.arange(num_groups), num_chunks):
            first_slice = group_slice_runs[group_indices[0]]
            final_slice = group_slice_runs[group_indices[-1] + 1]
            start = slice_runs[first_slice]
            stop = slice_runs[final_slice]
            chunks.append(self.isel(slice(start, stop)))
        return tuple(chunks)

    @property
    def time_bin_origin(self) -> np.datetime64:
        """Global half-open time-bin origin retained by subsets and chunks."""
        return np.datetime64(self._time_bin_origin_ns, "ns")


def _solar_vector_from_angles(
    up: np.ndarray,
    north: np.ndarray,
    east: np.ndarray,
    solar_zenith_degrees: float,
    solar_azimuth_degrees: float,
) -> np.ndarray:
    zenith = np.deg2rad(solar_zenith_degrees)
    azimuth = np.deg2rad(solar_azimuth_degrees)
    result = np.cos(zenith) * up + np.sin(zenith) * (
        np.cos(azimuth) * north + np.sin(azimuth) * east
    )
    return result / np.linalg.norm(result)


class OrbitalPlaneEngine(Engine):
    """Composite of persistent local Geometry2D engines along an orbit.

    Parameters
    ----------
    config
        Radiative-transfer configuration used by every local group engine.
    geometry
        Master orbital atmosphere grid.
    viewing_geometry
        Timestamped ECEF observer and look vectors.
    time_group_duration_s
        Width of the half-open time bins used to form local calculations. A
        complete vertical slice is assigned by its mean sample time and is never
        split between groups.
    group_padding_angle
        Extra angular margin, in radians, added before the first and after the
        last tangent location in every internal group. No separate horizon
        margin is added. The default is ten degrees on each side.
    max_horizontal_scale_residual
        Maximum permitted relative difference between a master-grid horizontal
        segment length and its locally spherical 2D representation. The default
        is 1%. Set to ``None`` to disable this validation.
    solar_handler
        Solar-position handler evaluated at each group's mean sample time.
        Supply this or ``sun_vectors_ecef`` when single scattering is enabled.
    sun_vectors_ecef
        Optional normalized or non-normalized ECEF Sun vector for every group.
    refraction_wavelength_nm
        Reference wavelength for automatic Ciddor refractivity.
    refraction_co2_ppm
        CO2 concentration for automatic Ciddor refractivity.
    derivative_execution
        ``"resident"`` (the default) reuses each group's engine and derivative
        atmosphere during repeated retrieval operations. ``"streaming"``
        constructs one temporary derivative engine at a time during a VJP,
        reducing retained derivative memory at the cost of repeated setup.
        Forward calculations always reuse the persistent group engines.
    max_eager_jacobian_bytes
        Safety limit for the estimated stitched full-Jacobian allocation used
        by ``calculate_radiance`` with a derivative-enabled atmosphere or by
        ``linearize(...).jacobian``. The default is 2 GiB. Set to ``None`` to
        explicitly permit a larger eager Jacobian; JVP and VJP are unaffected.
    """

    _engine: PyOrbitalPlaneEngine

    def __init__(
        self,
        config: sk.Config,
        geometry: OrbitalPlaneGeometry,
        viewing_geometry: OrbitalPlaneViewingGeometry,
        *,
        time_group_duration_s: float,
        group_padding_angle: float = np.deg2rad(10.0),
        max_horizontal_scale_residual: float | None = 0.01,
        solar_handler=None,
        sun_vectors_ecef: np.ndarray | None = None,
        refraction_wavelength_nm: float = 600.0,
        refraction_co2_ppm: float = 400.0,
        derivative_execution: Literal["resident", "streaming"] = "resident",
        max_eager_jacobian_bytes: int | None = _DEFAULT_MAX_EAGER_JACOBIAN_BYTES,
    ) -> None:
        if not isinstance(config, sk.Config):
            msg = "config must be a sasktran2.Config object"
            raise TypeError(msg)
        if not isinstance(geometry, OrbitalPlaneGeometry):
            msg = "OrbitalPlaneEngine requires OrbitalPlaneGeometry"
            raise TypeError(msg)
        if not isinstance(viewing_geometry, OrbitalPlaneViewingGeometry):
            msg = "OrbitalPlaneEngine requires OrbitalPlaneViewingGeometry"
            raise TypeError(msg)
        if isinstance(time_group_duration_s, bool | np.bool_):
            msg = "time_group_duration_s must be a finite positive number"
            raise TypeError(msg)
        try:
            duration_s = float(time_group_duration_s)
        except (TypeError, ValueError, OverflowError) as error:
            msg = "time_group_duration_s must be a finite positive number"
            raise TypeError(msg) from error
        if (
            not np.isfinite(duration_s)
            or duration_s <= 0
            or duration_s * 1e9 < 0.5
            or duration_s * 1e9 > np.iinfo(np.int64).max
        ):
            msg = (
                "time_group_duration_s must be finite, positive, and representable "
                "as at least one integer nanosecond"
            )
            raise ValueError(msg)
        if isinstance(group_padding_angle, bool | np.bool_):
            msg = "group_padding_angle must be a finite angle in radians"
            raise TypeError(msg)
        try:
            group_padding_angle = float(group_padding_angle)
        except (TypeError, ValueError, OverflowError) as error:
            msg = "group_padding_angle must be a finite angle in radians"
            raise TypeError(msg) from error
        if (
            not np.isfinite(group_padding_angle)
            or group_padding_angle < 0
            or group_padding_angle >= np.pi
        ):
            msg = "group_padding_angle must be finite and in [0, pi) radians"
            raise ValueError(msg)
        if max_horizontal_scale_residual is not None:
            if isinstance(max_horizontal_scale_residual, bool | np.bool_):
                msg = (
                    "max_horizontal_scale_residual must be a finite non-negative "
                    "number or None"
                )
                raise TypeError(msg)
            try:
                max_horizontal_scale_residual = float(max_horizontal_scale_residual)
            except (TypeError, ValueError, OverflowError) as error:
                msg = (
                    "max_horizontal_scale_residual must be a finite non-negative "
                    "number or None"
                )
                raise TypeError(msg) from error
            if (
                not np.isfinite(max_horizontal_scale_residual)
                or max_horizontal_scale_residual < 0
            ):
                msg = (
                    "max_horizontal_scale_residual must be finite and non-negative "
                    "when supplied"
                )
                raise ValueError(msg)
        if isinstance(refraction_wavelength_nm, bool | np.bool_):
            msg = "refraction_wavelength_nm must be a finite positive number"
            raise TypeError(msg)
        try:
            refraction_wavelength_nm = float(refraction_wavelength_nm)
        except (TypeError, ValueError, OverflowError) as error:
            msg = "refraction_wavelength_nm must be a finite positive number"
            raise TypeError(msg) from error
        if not np.isfinite(refraction_wavelength_nm) or refraction_wavelength_nm <= 0:
            msg = "refraction_wavelength_nm must be finite and positive"
            raise ValueError(msg)
        if isinstance(refraction_co2_ppm, bool | np.bool_):
            msg = "refraction_co2_ppm must be a finite non-negative number"
            raise TypeError(msg)
        try:
            refraction_co2_ppm = float(refraction_co2_ppm)
        except (TypeError, ValueError, OverflowError) as error:
            msg = "refraction_co2_ppm must be a finite non-negative number"
            raise TypeError(msg) from error
        if not np.isfinite(refraction_co2_ppm) or refraction_co2_ppm < 0:
            msg = "refraction_co2_ppm must be finite and non-negative"
            raise ValueError(msg)
        if derivative_execution not in ("resident", "streaming"):
            msg = (
                "derivative_execution must be 'resident' or 'streaming'; "
                f"got {derivative_execution!r}"
            )
            raise ValueError(msg)
        if max_eager_jacobian_bytes is not None:
            if isinstance(max_eager_jacobian_bytes, bool | np.bool_):
                msg = "max_eager_jacobian_bytes must be a positive integer or None"
                raise TypeError(msg)
            try:
                max_eager_jacobian_bytes = operator.index(max_eager_jacobian_bytes)
            except TypeError as error:
                msg = "max_eager_jacobian_bytes must be a positive integer or None"
                raise TypeError(msg) from error
            if max_eager_jacobian_bytes <= 0:
                msg = "max_eager_jacobian_bytes must be positive when supplied"
                raise ValueError(msg)
        if (
            config.single_scatter_source
            not in (sk.SingleScatterSource.NoSource, sk.SingleScatterSource.Exact)
            or config.multiple_scatter_source
            not in (
                sk.MultipleScatterSource.NoSource,
                sk.MultipleScatterSource.SuccessiveOrders,
            )
            or config.emission_source
            not in (
                sk.EmissionSource.NoSource,
                sk.EmissionSource.Standard,
                sk.EmissionSource.VolumeEmissionRate,
            )
        ):
            msg = (
                "OrbitalPlaneEngine supports exact single scattering, "
                "successive-orders multiple scattering, occultation, standard "
                "emission, and volume emission rate sources"
            )
            raise NotImplementedError(msg)
        if (
            config.multiple_scatter_source == sk.MultipleScatterSource.SuccessiveOrders
            and config.multiple_scatter_refraction
        ):
            msg = (
                "OrbitalPlaneEngine successive orders does not support "
                "diffuse-ray refraction"
            )
            raise NotImplementedError(msg)
        if config.solar_refraction:
            msg = (
                "OrbitalPlaneEngine supports line-of-sight refraction only; "
                "solar_refraction must be False"
            )
            raise NotImplementedError(msg)
        (
            group_times_ns,
            group_positions,
            group_geodetic_coordinates,
            group_local_up,
            group_local_north,
            group_local_east,
            group_horizontal_angle_bounds,
        ) = orbital_group_reference_data(
            geometry._geometry,
            viewing_geometry._viewing_geometry,
            duration_s,
            group_padding_angle,
            max_horizontal_scale_residual,
        )
        horizontal_source_angles = (
            config.successive_orders_horizontal_angle_grid_radians
        )
        if (
            config.multiple_scatter_source == sk.MultipleScatterSource.SuccessiveOrders
            and horizontal_source_angles is not None
        ):
            outside = (
                horizontal_source_angles[np.newaxis, :]
                < group_horizontal_angle_bounds[:, :1]
            ) | (
                horizontal_source_angles[np.newaxis, :]
                > group_horizontal_angle_bounds[:, 1:]
            )
            if np.any(outside):
                group_index = int(np.flatnonzero(np.any(outside, axis=1))[0])
                lower, upper = group_horizontal_angle_bounds[group_index]
                msg = (
                    "successive_orders_horizontal_angle_grid_radians must lie "
                    f"inside every local group grid; group {group_index} spans "
                    f"[{lower}, {upper}] radians"
                )
                raise ValueError(msg)
        scattering = (
            config.single_scatter_source != sk.SingleScatterSource.NoSource
            or config.multiple_scatter_source != sk.MultipleScatterSource.NoSource
        )
        if solar_handler is not None and sun_vectors_ecef is not None:
            msg = "Specify either solar_handler or sun_vectors_ecef, not both"
            raise ValueError(msg)
        if scattering and solar_handler is None and sun_vectors_ecef is None:
            msg = "Scattering orbital calculations require solar_handler or sun_vectors_ecef"
            raise ValueError(msg)
        if sun_vectors_ecef is None and solar_handler is not None:
            target_solar_angles = getattr(solar_handler, "target_solar_angles", None)
            if not callable(target_solar_angles):
                msg = "solar_handler must provide a callable target_solar_angles method"
                raise TypeError(msg)
            vectors = []
            for group_index, (time_ns, coordinates, up, north, east) in enumerate(
                zip(
                    group_times_ns,
                    group_geodetic_coordinates,
                    group_local_up,
                    group_local_north,
                    group_local_east,
                    strict=True,
                )
            ):
                latitude, longitude, altitude = coordinates
                angles = np.asarray(
                    target_solar_angles(
                        latitude,
                        longitude,
                        altitude,
                        pd.Timestamp(np.datetime64(int(time_ns), "ns")),
                    ),
                    dtype=np.float64,
                )
                if angles.shape != (2,) or np.any(~np.isfinite(angles)):
                    msg = (
                        "solar_handler.target_solar_angles must return two finite "
                        f"angles for group {group_index}; got shape {angles.shape}"
                    )
                    raise ValueError(msg)
                solar_zenith, solar_azimuth = angles
                vectors.append(
                    _solar_vector_from_angles(
                        up, north, east, solar_zenith, solar_azimuth
                    )
                )
            sun_vectors_ecef = np.asarray(vectors)
        elif sun_vectors_ecef is None:
            sun_vectors_ecef = np.asarray(group_positions)
        try:
            sun_vectors_ecef = np.array(
                sun_vectors_ecef, dtype=np.float64, order="C", copy=True
            )
        except (TypeError, ValueError) as error:
            msg = "sun_vectors_ecef must be convertible to floating-point vectors"
            raise ValueError(msg) from error
        expected_sun_shape = (len(group_times_ns), 3)
        if sun_vectors_ecef.shape != expected_sun_shape:
            msg = (
                f"sun_vectors_ecef must have shape {expected_sun_shape}; "
                f"got {sun_vectors_ecef.shape}"
            )
            raise ValueError(msg)
        norms = np.linalg.norm(sun_vectors_ecef, axis=1)
        if np.any(~np.isfinite(norms)) or np.any(norms == 0):
            msg = "sun_vectors_ecef must contain finite, non-zero vectors"
            raise ValueError(msg)
        sun_vectors_ecef /= norms[:, np.newaxis]

        self._engine = PyOrbitalPlaneEngine(
            config._config,
            geometry._geometry,
            viewing_geometry._viewing_geometry,
            duration_s,
            group_padding_angle,
            sun_vectors_ecef,
            refraction_wavelength_nm,
            refraction_co2_ppm,
            derivative_execution,
            max_horizontal_scale_residual,
        )
        self._config = config
        self._geometry = geometry
        self._viewing_geometry = viewing_geometry
        self._time_group_duration_s = duration_s
        self._group_padding_angle = group_padding_angle
        self._max_horizontal_scale_residual = max_horizontal_scale_residual
        self._refraction_wavelength_nm = refraction_wavelength_nm
        self._refraction_co2_ppm = refraction_co2_ppm
        self._derivative_execution = derivative_execution
        self._max_eager_jacobian_bytes = max_eager_jacobian_bytes
        self._num_groups = len(group_times_ns)
        self._sun_vectors_ecef = sun_vectors_ecef
        self._sun_vectors_ecef.setflags(write=False)
        self._config_signature = _structural_config_signature(config)

    def _validate_config_unchanged(self) -> None:
        if _structural_config_signature(self._config) != self._config_signature:
            msg = (
                "The Config used to construct OrbitalPlaneEngine was modified in a "
                "way that changes native engine structure; construct a new engine"
            )
            raise ValueError(msg)

    def _validate_orbital_atmosphere(self, atmosphere: sk.Atmosphere) -> None:
        if not isinstance(atmosphere, sk.Atmosphere):
            msg = "atmosphere must be a sasktran2.Atmosphere object"
            raise TypeError(msg)
        self._validate_config_unchanged()
        self._validate_atmosphere_geometry(atmosphere)
        if atmosphere.nstokes != self._config.num_stokes:
            msg = (
                "Atmosphere and OrbitalPlaneEngine Stokes dimensions do not match: "
                f"{atmosphere.nstokes} != {self._config.num_stokes}"
            )
            raise ValueError(msg)

    def _estimated_eager_jacobian_bytes(self, atmosphere: sk.Atmosphere) -> int:
        if not atmosphere.calculate_derivatives:
            return 0
        output_values = (
            atmosphere.num_wavel
            * len(self._viewing_geometry.times)
            * atmosphere.nstokes
        )
        parameter_values = 0
        for name in atmosphere.storage.derivative_mapping_names():
            mapping = atmosphere.storage.get_derivative_mapping(name)
            interpolator = np.asarray(mapping.interpolator)
            parameter_values += (
                int(interpolator.shape[1])
                if interpolator.size
                else atmosphere.num_locations
            )
        for name in atmosphere.surface.derivative_mapping_names():
            layout = atmosphere.surface_derivative_output_layout(name)
            if layout is not None:
                parameter_values += int(np.prod(layout[1]))
                continue
            mapping = atmosphere.surface.get_derivative_mapping(name)
            interpolator = np.asarray(mapping.interpolator)
            parameter_values += int(interpolator.shape[1]) if interpolator.size else 1
        return int(np.dtype(np.float64).itemsize * output_values * parameter_values)

    def estimate_eager_jacobian_bytes(self, atmosphere: sk.Atmosphere) -> int:
        """Estimate the stitched eager-Jacobian payload for ``atmosphere``.

        This materializes constituent derivative mappings but performs no
        radiative-transfer calculation. Temporary native group output and
        xarray overhead are not included, so the result is a lower bound.
        """
        self._validate_orbital_atmosphere(atmosphere)
        atmosphere.internal_object()
        return self._estimated_eager_jacobian_bytes(atmosphere)

    def _enforce_eager_jacobian_limit(self, atmosphere: sk.Atmosphere) -> None:
        estimate = self._estimated_eager_jacobian_bytes(atmosphere)
        limit = self._max_eager_jacobian_bytes
        if limit is not None and estimate > limit:
            msg = (
                "The requested eager orbital Jacobian is estimated to require at "
                f"least {estimate / 1024**3:.2f} GiB, exceeding the configured "
                f"{limit / 1024**3:.2f} GiB limit. Use "
                "calculate_radiance(..., derivatives=False) for radiance only, use linearize() "
                "for JVP/VJP retrieval operations, or explicitly set "
                "max_eager_jacobian_bytes=None when constructing the engine."
            )
            raise MemoryError(msg)

    def _prepare_refraction(self, atmosphere: sk.Atmosphere) -> None:
        if not self._config.los_refraction:
            return

        def native_state(name: str):
            value = atmosphere._native_state(name)
            if value is None:
                return None
            return np.ascontiguousarray(value.reshape(self._geometry.shape))

        refractive_index = atmosphere.refractive_index
        if refractive_index is not None:
            refractive_index = atmosphere._spatial_layout.native_state(
                refractive_index, "refractive_index"
            ).reshape(self._geometry.shape)
            refractive_index = np.ascontiguousarray(refractive_index)
        self._engine._prepare_refraction(
            native_state("pressure_pa"),
            native_state("temperature_k"),
            native_state("specific_humidity"),
            refractive_index,
        )

    def _prepare_surface(self, atmosphere: sk.Atmosphere) -> None:
        field = atmosphere._orbital_lambertian_surface
        if field is not None:
            field = np.ascontiguousarray(field, dtype=np.float64)
        spectral_interpolator = atmosphere._orbital_lambertian_spectral_interpolator
        if spectral_interpolator is not None:
            spectral_interpolator = np.ascontiguousarray(
                spectral_interpolator, dtype=np.float64
            )
        self._engine._set_lambertian_surface(
            field,
            atmosphere._orbital_lambertian_derivative_name,
            spectral_interpolator,
            atmosphere._orbital_lambertian_spatial_parameters,
        )

    def _validate_linearization_session(self, state_generation: int) -> None:
        self._validate_config_unchanged()
        if self._engine.state_generation() != state_generation:
            msg = (
                "The orbital engine's refractive geometry or spatial surface state "
                "has changed since this linearization was created. Call "
                "engine.linearize(atmosphere) again."
            )
            raise sk.StaleLinearizationError(msg)

    def _calculate_radiance(
        self,
        atmosphere,
        internal_atmosphere=None,
        *,
        include_derivatives: bool = True,
    ):
        self._validate_orbital_atmosphere(atmosphere)
        if internal_atmosphere is None:
            internal_atmosphere = atmosphere.internal_object()
            if include_derivatives:
                self._enforce_eager_jacobian_limit(atmosphere)
            self._prepare_refraction(atmosphere)
            self._prepare_surface(atmosphere)
        elif include_derivatives:
            self._enforce_eager_jacobian_limit(atmosphere)
        return super()._calculate_radiance(
            atmosphere,
            internal_atmosphere=internal_atmosphere,
            include_derivatives=include_derivatives,
        )

    def linearize(
        self,
        atmosphere: sk.Atmosphere,
        *,
        prepare_parameters: Iterable[str] | None = None,
    ):
        """Construct a retrieval linearization for the orbital atmosphere.

        ``prepare_parameters`` optionally prepares the named resident group
        workspaces during the primal calculation. This avoids repeating the
        successive-orders forward solve in the first matching JVP or VJP;
        other registered parameters remain available lazily. Streaming mode
        validates but does not retain the preparation hint.
        """
        self._validate_orbital_atmosphere(atmosphere)
        internal_atmosphere = atmosphere.internal_object()
        self._prepare_refraction(atmosphere)
        self._prepare_surface(atmosphere)
        state_generation = self._engine.state_generation()
        prepared_parameters = prepare_parameters
        if self._derivative_execution == "streaming" and prepare_parameters is not None:
            _selected_parameters(
                prepare_parameters,
                self._linearization_registry(atmosphere).specs,
                operation="linearization preparation",
            )
            prepared_parameters = None

        def validate_session() -> None:
            self._validate_linearization_session(state_generation)

        return super()._linearize(
            atmosphere,
            internal_atmosphere=internal_atmosphere,
            validate_session=validate_session,
            prepare_parameters=prepared_parameters,
        )

    @property
    def time_group_duration_s(self) -> float:
        return self._time_group_duration_s

    @property
    def group_padding_angle(self) -> float:
        return self._group_padding_angle

    @property
    def max_horizontal_scale_residual(self) -> float | None:
        return self._max_horizontal_scale_residual

    @property
    def refraction_wavelength_nm(self) -> float:
        return self._refraction_wavelength_nm

    @property
    def refraction_co2_ppm(self) -> float:
        return self._refraction_co2_ppm

    @property
    def derivative_execution(self) -> Literal["resident", "streaming"]:
        return self._derivative_execution

    @property
    def max_eager_jacobian_bytes(self) -> int | None:
        return self._max_eager_jacobian_bytes

    @property
    def sun_vectors_ecef(self) -> np.ndarray:
        """Normalized group Sun vectors as a read-only array."""
        return _readonly_view(self._sun_vectors_ecef)

    @property
    def num_groups(self) -> int:
        return self._num_groups

    @property
    def group_diagnostics(self) -> list[dict]:
        """Return detached diagnostic snapshots for every non-empty group.

        The local ``earth_radius_m`` is the least-squares constant radius at
        the group's actual tangent locations. Each ECEF LOS is transformed
        through a geometry-relative ray policy that preserves its geoid tangent
        altitude and surface location. ``surface_radius_residuals_m`` and
        ``maximum_relative_observation_radius_residual`` describe the radius fit.
        ``maximum_relative_horizontal_scale_residual`` includes both that radius
        approximation and plane-projection distortion across the full selected
        grid window. ``ray_horizontal_angles`` and
        ``ray_out_of_plane_angles`` show the transformed tangent locations;
        ``horizontal_angles`` and ``horizontal_distances_m`` describe the
        projected internal grid. ``edge_clipping`` reports a requested group
        window truncated by the master orbital grid. ``path_edge_clipping``
        and ``path_edge_clipping_per_los`` report traced atmospheric path
        segments that actually use a local grid's extended horizontal edge
        cells, including after a refractive geometry refresh.
        """
        diagnostics = self._engine.group_diagnostics()
        for item in diagnostics:
            item["reference_time"] = np.datetime64(item["reference_time_ns"], "ns")
            item["reference_position_ecef_m"] = item["reference_position_ecef"]
            item["max_out_of_plane_angle_rad"] = item["max_out_of_plane_angle"]
        return diagnostics
