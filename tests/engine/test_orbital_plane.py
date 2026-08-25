from __future__ import annotations

import gc

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2._core_rust import PyKokhanovsky, orbital_ciddor_index
from sasktran2.optical.refraction import ciddor_index_of_refraction

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.array([0.0, 10_000.0, 20_000.0, 35_000.0, 60_000.0])


def orbital_geometry(*, warped: bool = False) -> sk.OrbitalPlaneGeometry:
    angles = np.linspace(-0.8, 0.8, 65)
    across = 0.025 * np.sin(3 * angles) if warped else np.zeros_like(angles)
    directions = np.column_stack((np.sin(angles), across, np.cos(angles)))
    directions /= np.linalg.norm(directions, axis=1)[:, np.newaxis]
    return sk.OrbitalPlaneGeometry(
        EARTH_RADIUS_M,
        ALTITUDES_M,
        EARTH_RADIUS_M * directions,
    )


def limb_viewing(
    geometry: sk.OrbitalPlaneGeometry,
    angles: np.ndarray,
    times: np.ndarray | None = None,
    vertical_slice: np.ndarray | None = None,
    tangent_altitude_m: float | np.ndarray = 20_000.0,
) -> sk.OrbitalPlaneViewingGeometry:
    if times is None:
        times = np.datetime64("2026-01-01T00:00:00", "ns") + np.arange(
            len(angles)
        ) * np.timedelta64(30, "s")
    if vertical_slice is None:
        vertical_slice = np.arange(len(angles))
    observers = []
    tangents = []
    track = geometry.ground_track_ecef_m
    track /= np.linalg.norm(track, axis=1)[:, np.newaxis]
    cumulative = geometry.cumulative_angles
    requested = angles + 0.8
    tangent_altitudes = np.broadcast_to(
        np.asarray(tangent_altitude_m, dtype=np.float64), requested.shape
    )
    for coordinate, tangent_altitude in zip(requested, tangent_altitudes, strict=True):
        upper = np.searchsorted(cumulative, coordinate, side="right")
        upper = np.clip(upper, 1, len(cumulative) - 1)
        fraction = (coordinate - cumulative[upper - 1]) / (
            cumulative[upper] - cumulative[upper - 1]
        )
        up = (1 - fraction) * track[upper - 1] + fraction * track[upper]
        up /= np.linalg.norm(up)
        trial = np.array([0.0, 1.0, 0.0])
        forward = np.cross(trial, up)
        if np.linalg.norm(forward) < 1e-8:
            forward = np.cross(np.array([1.0, 0.0, 0.0]), up)
        forward /= np.linalg.norm(forward)
        tangent_radius = EARTH_RADIUS_M + tangent_altitude
        observer_radius = EARTH_RADIUS_M + 700_000.0
        distance = np.sqrt(observer_radius**2 - tangent_radius**2)
        tangent = tangent_radius * up
        observers.append(tangent - distance * forward)
        tangents.append(tangent)
    return sk.OrbitalPlaneViewingGeometry.from_tangent_locations(
        times,
        np.asarray(observers),
        np.asarray(tangents),
        vertical_slice=vertical_slice,
    )


def transmission_config(*, refraction: bool = False) -> sk.Config:
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.occultation_source = sk.OccultationSource.Standard
    config.los_refraction = refraction
    return config


def raw_atmosphere(
    geometry: sk.OrbitalPlaneGeometry,
    config: sk.Config,
    *,
    calculate_derivatives: bool = True,
) -> sk.Atmosphere:
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        calculate_derivatives=calculate_derivatives,
        legendre_derivative=False,
    )
    atmosphere.storage.total_extinction[:] = 1e-5
    atmosphere.storage.ssa[:] = 0.0
    return atmosphere


def test_orbital_geometry_normalizes_track_and_uses_altitude_fastest_layout():
    geometry = orbital_geometry(warped=True)
    np.testing.assert_allclose(
        np.linalg.norm(geometry.ground_track_ecef_m, axis=1), EARTH_RADIUS_M
    )
    assert geometry.shape == (65, len(ALTITUDES_M))
    assert np.all(np.diff(geometry.cumulative_angles) > 0)
    atmosphere = sk.Atmosphere(geometry, transmission_config(), numwavel=1)
    assert atmosphere.storage.total_extinction.shape[0] == np.prod(geometry.shape)


@pytest.mark.parametrize(
    ("earth_radius_m", "altitudes", "track", "surface_radii", "match"),
    [
        (0.0, ALTITUDES_M, np.ones((2, 3)), None, "earth_radius_m"),
        (
            EARTH_RADIUS_M,
            np.array([0.0, 10_000.0, 10_000.0]),
            np.ones((2, 3)),
            None,
            "strictly increasing",
        ),
        (EARTH_RADIUS_M, ALTITUDES_M, np.ones(3), None, "ground_track_ecef_m"),
        (
            EARTH_RADIUS_M,
            ALTITUDES_M,
            np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            None,
            "non-zero",
        ),
        (
            EARTH_RADIUS_M,
            ALTITUDES_M,
            np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
            np.array([EARTH_RADIUS_M]),
            "surface_radii_m",
        ),
    ],
)
def test_orbital_geometry_rejects_invalid_public_inputs(
    earth_radius_m, altitudes, track, surface_radii, match
):
    with pytest.raises(ValueError, match=match):
        sk.OrbitalPlaneGeometry(
            earth_radius_m,
            altitudes,
            track,
            surface_radii_m=surface_radii,
        )


def test_orbital_geometry_rejects_boolean_earth_radius():
    with pytest.raises(TypeError, match="earth_radius_m"):
        sk.OrbitalPlaneGeometry(True, ALTITUDES_M, np.ones((2, 3)))


def test_viewing_geometry_constructs_geoid_following_atmosphere_grid():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([-0.3, 0.0, 0.3]))
    delta = 0.04
    geometry = viewing.construct_atmosphere_geometry(ALTITUDES_M, delta)

    assert geometry.shape[0] > 3
    assert np.all(np.diff(geometry.cumulative_angles) > 0)
    assert np.max(np.diff(geometry.cumulative_angles)) <= delta + 1e-12
    np.testing.assert_allclose(
        np.linalg.norm(geometry.ground_track_ecef_m, axis=1),
        geometry.surface_radii_m,
    )
    assert np.ptp(geometry.surface_radii_m) > 100.0
    geoid = sk.WGS84()
    for location in geometry.ground_track_ecef_m:
        geoid.from_xyz(location)
        assert geoid.altitude == pytest.approx(0.0, abs=1e-5)


def test_constructed_grid_exposes_geoid_and_vertical_slice_provenance():
    source_geometry = orbital_geometry()
    times = np.array(
        [
            "2026-01-01T00:00:00",
            "2026-01-01T00:00:02",
            "2026-01-01T00:00:30",
            "2026-01-01T00:01:00",
        ],
        dtype="datetime64[ns]",
    )
    vertical_slice = np.array([10, 10, 20, 30])
    viewing = limb_viewing(
        source_geometry,
        np.array([-0.31, -0.29, 0.0, 0.3]),
        times=times,
        vertical_slice=vertical_slice,
    )
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M, 0.02, path_padding_angle=0.1
    )

    grid = geometry.grid_dataset
    anchors = geometry.vertical_slice_anchors
    assert grid.sizes == {
        "orbital_position": geometry.shape[0],
        "xyz": 3,
        "altitude": len(ALTITUDES_M),
    }
    assert {
        "along_track_angle_rad",
        "geodetic_latitude_deg",
        "geodetic_longitude_deg",
        "reference_time",
    } <= set(grid.coords)
    np.testing.assert_array_equal(grid.altitude, ALTITUDES_M)
    np.testing.assert_allclose(grid.ground_track_ecef_m, geometry.ground_track_ecef_m)
    np.testing.assert_array_equal(anchors.vertical_slice, [10, 20, 30])
    np.testing.assert_array_equal(
        anchors.reference_time,
        np.array(
            [
                "2026-01-01T00:00:01",
                "2026-01-01T00:00:30",
                "2026-01-01T00:01:00",
            ],
            dtype="datetime64[ns]",
        ),
    )
    assert np.all(np.diff(anchors.along_track_angle_rad) >= 0)
    assert grid.reference_time[0] == anchors.reference_time[0]
    assert grid.reference_time[-1] == anchors.reference_time[-1]

    geoid = sk.WGS84()
    for index, location in enumerate(geometry.ground_track_ecef_m):
        geoid.from_xyz(location)
        assert grid.geodetic_latitude_deg[index] == pytest.approx(
            geoid.latitude, abs=2e-11
        )
        longitude_difference = (
            grid.geodetic_longitude_deg[index] - geoid.longitude + 180
        ) % 360 - 180
        assert longitude_difference == pytest.approx(0.0, abs=2e-11)

    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=120,
        max_horizontal_scale_residual=None,
    )
    template = engine.linearize(raw_atmosphere(geometry, config)).tangent_template
    for coordinate in (
        "along_track_angle_rad",
        "geodetic_latitude_deg",
        "geodetic_longitude_deg",
        "reference_time",
    ):
        assert coordinate in template.extinction.coords
        xr.testing.assert_equal(
            template.extinction.coords[coordinate], grid.coords[coordinate]
        )


def test_project_ground_track_returns_interpolation_coordinates_and_order_policy():
    angles = np.linspace(0.0, 2 * np.pi, 73)
    directions = np.column_stack(
        (np.sin(angles), np.zeros_like(angles), np.cos(angles))
    )
    geometry = sk.OrbitalPlaneGeometry(
        EARTH_RADIUS_M,
        ALTITUDES_M,
        EARTH_RADIUS_M * directions,
    )

    increasing = geometry.project_ground_track(
        geometry.ground_track_ecef_m[[0, 1]], order="increasing"
    )
    decreasing = geometry.project_ground_track(
        geometry.ground_track_ecef_m[[-1, -2]], order="decreasing"
    )

    np.testing.assert_allclose(
        increasing.along_track_angle_rad, geometry.cumulative_angles[:2], atol=2e-14
    )
    np.testing.assert_allclose(
        decreasing.along_track_angle_rad,
        geometry.cumulative_angles[[-1, -2]],
        atol=2e-14,
    )
    np.testing.assert_allclose(increasing.cross_track_angle_rad, 0.0, atol=2e-8)
    np.testing.assert_allclose(decreasing.cross_track_angle_rad, 0.0, atol=2e-8)
    assert bool(increasing.at_track_edge[0])
    assert bool(decreasing.at_track_edge[0])

    with pytest.raises(ValueError, match="order"):
        geometry.project_ground_track(directions[:1], order="sideways")
    with pytest.raises(ValueError, match="shape"):
        geometry.project_ground_track(np.ones((2, 2)))


def test_geoid_resampling_preserves_the_selected_radial_directions():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([-0.35, 0.27]))
    expected_directions = []
    geoid = sk.WGS84()
    for observer, look in zip(
        viewing.observer_positions_ecef_m,
        viewing.look_directions_ecef,
        strict=True,
    ):
        geoid.from_tangent_point(observer, look)
        latitude = geoid.latitude
        longitude = geoid.longitude
        geoid.from_lat_lon_alt(latitude, longitude, 0.0)
        expected_directions.append(geoid.location / np.linalg.norm(geoid.location))

    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M, 0.04, path_padding_angle=0.0
    )
    actual = geometry.ground_track_ecef_m[[0, -1]]
    actual /= np.linalg.norm(actual, axis=1, keepdims=True)
    np.testing.assert_allclose(actual, expected_directions, atol=2e-14)


def test_track_construction_is_invariant_to_los_order_within_vertical_slices():
    source_geometry = orbital_geometry()
    angles = np.array([-0.32, -0.30, -0.28, 0.08, 0.10, 0.12])
    vertical_slice = np.repeat(np.arange(2), 3)
    times = np.datetime64("2026-01-01", "ns") + vertical_slice * np.timedelta64(30, "s")
    viewing = limb_viewing(
        source_geometry,
        angles,
        times=times,
        vertical_slice=vertical_slice,
    )
    reordered = viewing.isel([2, 1, 0, 5, 4, 3])

    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M, 0.02, path_padding_angle=0.1
    )
    reordered_geometry = reordered.construct_atmosphere_geometry(
        ALTITUDES_M, 0.02, path_padding_angle=0.1
    )

    assert geometry.shape == reordered_geometry.shape
    np.testing.assert_allclose(
        geometry.ground_track_ecef_m,
        reordered_geometry.ground_track_ecef_m,
        atol=1e-8,
    )
    np.testing.assert_allclose(
        geometry.cumulative_angles, reordered_geometry.cumulative_angles, atol=2e-14
    )
    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=60
    )
    reordered_engine = sk.OrbitalPlaneEngine(
        transmission_config(),
        reordered_geometry,
        reordered,
        time_group_duration_s=60,
    )
    np.testing.assert_array_equal(
        engine.group_diagnostics[0]["grid_indices"],
        reordered_engine.group_diagnostics[0]["grid_indices"],
    )


def test_equal_time_image_uses_tangent_altitude_gradient_for_track_direction():
    geoid = sk.WGS84()
    geoid.from_lat_lon_alt(45.0, 10.0, 700_000.0)
    observer = geoid.location.copy()
    tangent_altitudes_m = np.array([10_000.0, 20_000.0, 40_000.0, 60_000.0])
    looks = []
    tangent_surface_directions = []
    for altitude_m in tangent_altitudes_m:
        look = geoid.from_tangent_altitude(
            altitude_m, observer, np.array([0.0, 0.0, 1.0])
        )
        looks.append(look)
        geoid.from_tangent_point(observer, look)
        latitude = geoid.latitude
        longitude = geoid.longitude
        geoid.from_lat_lon_alt(latitude, longitude, 0.0)
        tangent_surface_directions.append(
            geoid.location / np.linalg.norm(geoid.location)
        )

    times = np.full(
        len(tangent_altitudes_m),
        np.datetime64("2026-01-01T00:00:00", "ns"),
    )
    viewing = sk.OrbitalPlaneViewingGeometry(
        times,
        np.broadcast_to(observer, (len(times), 3)),
        np.asarray(looks),
        vertical_slice=np.zeros(len(times), dtype=int),
    )
    reversed_viewing = viewing.isel(np.arange(len(viewing) - 1, -1, -1))

    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M, 0.01, path_padding_angle=0.1
    )
    reversed_geometry = reversed_viewing.construct_atmosphere_geometry(
        ALTITUDES_M, 0.01, path_padding_angle=0.1
    )

    np.testing.assert_allclose(
        geometry.ground_track_ecef_m,
        reversed_geometry.ground_track_ecef_m,
        atol=2e-8,
    )
    track = geometry.ground_track_ecef_m
    track = track / np.linalg.norm(track, axis=1, keepdims=True)
    tangent_gradient = tangent_surface_directions[-1] - tangent_surface_directions[0]
    assert np.dot(track[-1] - track[0], tangent_gradient) > 0


def test_one_time_varying_vertical_scan_defines_one_track_slice():
    source_geometry = orbital_geometry()
    angles = np.array([-0.12, -0.04, 0.05, 0.13])
    times = np.datetime64("2026-01-01", "ns") + np.arange(4) * np.timedelta64(8, "s")
    viewing = limb_viewing(
        source_geometry,
        angles,
        times=times,
        vertical_slice=np.zeros(len(times), dtype=int),
    )
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M, 0.02, path_padding_angle=0.1
    )

    assert viewing.num_vertical_slices == 1
    assert geometry.cumulative_angles[-1] >= angles[-1] - angles[0] + 0.2 - 1e-10
    with pytest.raises(ValueError, match="between 1 and 1"):
        viewing.split(2, time_group_duration_s=10)
    engine = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=10,
        max_horizontal_scale_residual=None,
    )
    assert sum(
        len(group["observation_indices"]) for group in engine.group_diagnostics
    ) == len(viewing)
    assert [group["observation_indices"] for group in engine.group_diagnostics] == [
        [0, 1, 2, 3]
    ]


def test_viewing_geometry_normalizes_reported_looks_without_mutating_inputs():
    times = np.array(["2026-01-01", "2026-01-01T00:00:01"], dtype="datetime64[ns]")
    observers = np.array(
        [[EARTH_RADIUS_M + 700_000.0, 0.0, 0.0], [0.0, 0.0, EARTH_RADIUS_M + 700_000.0]]
    )
    looks = np.array([[-2.0, 0.0, 0.0], [0.0, 0.0, -3.0]])
    original_looks = looks.copy()

    viewing = sk.OrbitalPlaneViewingGeometry(
        times, observers, looks, vertical_slice=np.arange(len(times))
    )

    np.testing.assert_array_equal(looks, original_looks)
    np.testing.assert_allclose(
        np.linalg.norm(viewing.look_directions_ecef, axis=1), 1.0
    )
    np.testing.assert_array_equal(
        viewing.geometry_ds.look_direction_ecef, viewing.look_directions_ecef
    )
    assert not viewing.times.flags.writeable
    assert not viewing.observer_positions_ecef_m.flags.writeable
    assert not viewing.look_directions_ecef.flags.writeable
    with pytest.raises(ValueError, match="read-only"):
        viewing.look_directions_ecef[0, 0] = 0.0
    with pytest.raises(ValueError, match="WRITEABLE"):
        viewing.times.setflags(write=True)


def test_viewing_geometry_requires_valid_vertical_slice_groups():
    times = np.datetime64("2026-01-01", "ns") + np.arange(3) * np.timedelta64(1, "s")
    observers = np.tile(np.array([0.0, 0.0, EARTH_RADIUS_M + 700_000.0]), (3, 1))
    looks = np.tile(np.array([0.0, 0.0, -1.0]), (3, 1))

    with pytest.raises(TypeError, match="vertical_slice"):
        sk.OrbitalPlaneViewingGeometry(times, observers, looks)
    with pytest.raises(ValueError, match="shape"):
        sk.OrbitalPlaneViewingGeometry(
            times, observers, looks, vertical_slice=np.array([0, 1])
        )
    with pytest.raises(TypeError, match="non-negative integer"):
        sk.OrbitalPlaneViewingGeometry(
            times, observers, looks, vertical_slice=np.array([False, False, False])
        )
    with pytest.raises(ValueError, match="non-negative"):
        sk.OrbitalPlaneViewingGeometry(
            times, observers, looks, vertical_slice=np.array([0, -1, -1])
        )
    with pytest.raises(ValueError, match="contiguous block"):
        sk.OrbitalPlaneViewingGeometry(
            times, observers, looks, vertical_slice=np.array([0, 1, 0])
        )


def test_vertical_slice_metadata_is_read_only_and_split_keeps_slices_complete():
    source_geometry = orbital_geometry()
    vertical_slice = np.array([10, 10, 20, 30, 30, 30])
    viewing = limb_viewing(
        source_geometry,
        np.linspace(-0.3, 0.3, len(vertical_slice)),
        vertical_slice=vertical_slice,
    )

    np.testing.assert_array_equal(viewing.vertical_slice, vertical_slice)
    np.testing.assert_array_equal(viewing.geometry_ds.vertical_slice, vertical_slice)
    assert not viewing.vertical_slice.flags.writeable
    assert viewing.num_vertical_slices == 3
    chunks = viewing.split(2, time_group_duration_s=60)
    assert [len(chunk) for chunk in chunks] == [3, 3]
    np.testing.assert_array_equal(chunks[0].vertical_slice, [10, 10, 20])
    np.testing.assert_array_equal(chunks[1].vertical_slice, [30, 30, 30])


@pytest.mark.parametrize(
    ("times", "observers", "looks", "match"),
    [
        (np.datetime64("2026-01-01"), np.ones((1, 3)), np.ones((1, 3)), "times"),
        (
            np.array(["NaT"], dtype="datetime64[ns]"),
            np.ones((1, 3)),
            np.ones((1, 3)),
            "NaT",
        ),
        (
            np.array(["2026-01-01"], dtype="datetime64[ns]"),
            np.ones(3),
            np.ones((1, 3)),
            "observer_positions_ecef_m",
        ),
        (
            np.array(["2026-01-01"], dtype="datetime64[ns]"),
            np.ones((1, 3)),
            np.zeros((1, 3)),
            "look_directions_ecef",
        ),
        (
            np.array(["2026-01-01"], dtype="datetime64[ns]"),
            np.array([[np.nan, 0.0, 1.0]]),
            np.ones((1, 3)),
            "finite",
        ),
    ],
)
def test_viewing_geometry_rejects_invalid_public_inputs(times, observers, looks, match):
    with pytest.raises(ValueError, match=match):
        sk.OrbitalPlaneViewingGeometry(
            times, observers, looks, vertical_slice=np.array([0])
        )


def test_from_tangent_locations_reports_shape_errors_at_the_public_boundary():
    times = np.array(["2026-01-01"], dtype="datetime64[ns]")
    with pytest.raises(ValueError, match="observer_positions_ecef_m"):
        sk.OrbitalPlaneViewingGeometry.from_tangent_locations(
            times,
            np.ones(3),
            np.ones((1, 3)),
            vertical_slice=np.array([0]),
        )
    with pytest.raises(ValueError, match="tangent_locations_ecef_m"):
        sk.OrbitalPlaneViewingGeometry.from_tangent_locations(
            times,
            np.ones((1, 3)),
            np.ones((2, 3)),
            vertical_slice=np.array([0]),
        )


def test_viewing_geometry_subsets_and_splits_preserve_observation_order():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.linspace(-0.3, 0.3, 7))

    subset = viewing.isel([5, 1, 3])
    np.testing.assert_array_equal(subset.times, viewing.times[[5, 1, 3]])
    np.testing.assert_array_equal(
        subset.observer_positions_ecef_m,
        viewing.observer_positions_ecef_m[[5, 1, 3]],
    )

    chunks = viewing.split(4, time_group_duration_s=20)
    assert [len(chunk) for chunk in chunks] == [2, 2, 2, 1]
    np.testing.assert_array_equal(
        np.concatenate([chunk.times for chunk in chunks]), viewing.times
    )
    with pytest.raises(ValueError, match="at least one LOS"):
        viewing.isel([])
    with pytest.raises(ValueError, match="between 1 and 7"):
        viewing.split(8, time_group_duration_s=20)
    with pytest.raises(TypeError, match="positive integer"):
        viewing.split(1.5, time_group_duration_s=20)


def test_viewing_geometry_uses_ten_degree_default_and_user_path_padding():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([-0.3, 0.3]))
    delta = 0.01
    padding = np.deg2rad(10.0)
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M,
        delta,
        path_padding_angle=padding,
        geoid=sk.SphericalGeoid(EARTH_RADIUS_M),
    )
    default_geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M,
        delta,
        geoid=sk.SphericalGeoid(EARTH_RADIUS_M),
    )

    expected_span = 0.6 + 2 * padding
    assert geometry.cumulative_angles[-1] == pytest.approx(expected_span)
    assert geometry.shape[0] == int(np.ceil(expected_span / delta)) + 1
    np.testing.assert_allclose(
        default_geometry.cumulative_angles, geometry.cumulative_angles
    )


def test_viewing_geometry_supports_zero_path_padding():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([-0.3, 0.3]))
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M,
        0.01,
        path_padding_angle=0.0,
        geoid=sk.SphericalGeoid(EARTH_RADIUS_M),
    )

    assert geometry.cumulative_angles[-1] == pytest.approx(0.6)
    assert np.all(np.diff(geometry.cumulative_angles) > 0)


def test_viewing_geometry_caps_accidentally_huge_atmosphere_grids():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([-0.3, 0.3]))
    with pytest.raises(ValueError, match="max_orbital_positions=10"):
        viewing.construct_atmosphere_geometry(
            ALTITUDES_M,
            0.01,
            path_padding_angle=0.0,
            geoid=sk.SphericalGeoid(EARTH_RADIUS_M),
            max_orbital_positions=10,
        )
    with pytest.raises(TypeError, match="max_orbital_positions"):
        viewing.construct_atmosphere_geometry(
            ALTITUDES_M,
            0.01,
            max_orbital_positions=1.5,
        )
    with pytest.raises(TypeError, match="along_track_angle_delta"):
        viewing.construct_atmosphere_geometry(ALTITUDES_M, True)
    with pytest.raises(TypeError, match="path_padding_angle"):
        viewing.construct_atmosphere_geometry(
            ALTITUDES_M, 0.01, path_padding_angle=np.bool_(False)
        )


def test_each_group_radius_and_ray_policy_are_anchored_to_observations():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([-0.25, 0.25]))
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M,
        0.025,
        path_padding_angle=np.deg2rad(30.0),
        geoid=sk.WGS84(),
    )
    unpadded = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=20,
        group_padding_angle=0.0,
    )
    padded = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=20,
        group_padding_angle=np.deg2rad(5.0),
    )

    for unpadded_diagnostic, padded_diagnostic in zip(
        unpadded.group_diagnostics, padded.group_diagnostics, strict=True
    ):
        observation_index = unpadded_diagnostic["observation_indices"][0]
        geoid = sk.WGS84()
        geoid.from_tangent_point(
            viewing.observer_positions_ecef_m[observation_index],
            viewing.look_directions_ecef[observation_index],
        )
        expected_tangent_altitude_m = geoid.altitude
        latitude = geoid.latitude
        longitude = geoid.longitude
        geoid.from_lat_lon_alt(latitude, longitude, 0.0)
        expected = np.linalg.norm(geoid.location)
        assert unpadded_diagnostic["earth_radius_m"] == pytest.approx(expected)
        assert padded_diagnostic["earth_radius_m"] == pytest.approx(expected)
        assert unpadded_diagnostic["surface_radius_residuals_m"] == [
            pytest.approx(0.0, abs=1e-9)
        ]
        assert unpadded_diagnostic["tangent_altitudes_m"] == [
            pytest.approx(expected_tangent_altitude_m, abs=1e-6)
        ]
        assert unpadded_diagnostic[
            "maximum_tangent_direction_error_radians"
        ] == pytest.approx(0.0, abs=2e-8)
        theta = unpadded_diagnostic["ray_horizontal_angles"][0]
        phi = unpadded_diagnostic["ray_out_of_plane_angles"][0]
        reference_z = np.asarray(unpadded_diagnostic["reference_z"])
        reference_x = np.asarray(unpadded_diagnostic["reference_x"])
        reference_y = np.asarray(unpadded_diagnostic["reference_y"])
        reconstructed_direction = (
            np.cos(phi) * (np.cos(theta) * reference_z + np.sin(theta) * reference_x)
            - np.sin(phi) * reference_y
        )
        np.testing.assert_allclose(
            reconstructed_direction,
            geoid.location / np.linalg.norm(geoid.location),
            atol=2e-14,
        )
        assert unpadded_diagnostic[
            "maximum_relative_observation_radius_residual"
        ] == pytest.approx(0.0, abs=1e-15)
        assert (
            padded_diagnostic["grid_indices"][0]
            < unpadded_diagnostic["grid_indices"][0]
        )
        assert (
            padded_diagnostic["grid_indices"][-1]
            > unpadded_diagnostic["grid_indices"][-1]
        )
        assert padded_diagnostic["edge_clipping"] == (False, False)


def test_ray_policy_uses_the_geoid_selected_for_atmosphere_geometry():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([-0.2, 0.2]))
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M,
        0.025,
        path_padding_angle=np.deg2rad(5.0),
        geoid=sk.SphericalGeoid(EARTH_RADIUS_M),
    )
    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=60
    )

    diagnostic = engine.group_diagnostics[0]
    np.testing.assert_allclose(diagnostic["tangent_altitudes_m"], 20_000.0, atol=1e-7)
    np.testing.assert_allclose(
        diagnostic["observation_surface_radii_m"], EARTH_RADIUS_M, atol=1e-7
    )


def test_ray_policy_conserves_tangent_altitude_and_horizontal_location():
    angles = np.linspace(-0.8, 0.8, 65)
    directions = np.column_stack(
        (np.sin(angles), np.zeros_like(angles), np.cos(angles))
    )
    surface_radii_m = EARTH_RADIUS_M + 12_000.0 * np.sin(2.0 * angles + 0.4)
    geometry = sk.OrbitalPlaneGeometry(
        EARTH_RADIUS_M,
        ALTITUDES_M,
        surface_radii_m[:, np.newaxis] * directions,
        surface_radii_m=surface_radii_m,
    )

    grid_indices = np.array([28, 29])
    tangent_altitude_m = 19_500.0
    tangent_radii_m = surface_radii_m[grid_indices] + tangent_altitude_m
    tangent_directions = directions[grid_indices]
    forward = np.column_stack(
        (
            np.cos(angles[grid_indices]),
            np.zeros(len(grid_indices)),
            -np.sin(angles[grid_indices]),
        )
    )
    observer_radii_m = tangent_radii_m + 700_000.0
    observer_distances_m = np.sqrt(observer_radii_m**2 - tangent_radii_m**2)
    tangents = tangent_radii_m[:, np.newaxis] * tangent_directions
    observers = tangents - observer_distances_m[:, np.newaxis] * forward
    times = np.datetime64("2026-01-01", "ns") + np.arange(2) * np.timedelta64(1, "s")
    viewing = sk.OrbitalPlaneViewingGeometry.from_tangent_locations(
        times, observers, tangents, vertical_slice=np.arange(len(times))
    )
    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=60
    )
    diagnostic = engine.group_diagnostics[0]

    expected_radius_m = np.mean(surface_radii_m[grid_indices])
    expected_residuals_m = surface_radii_m[grid_indices] - expected_radius_m
    assert diagnostic["earth_radius_m"] == pytest.approx(expected_radius_m)
    np.testing.assert_allclose(
        diagnostic["observation_surface_radii_m"],
        surface_radii_m[grid_indices],
        atol=1e-8,
    )
    np.testing.assert_allclose(
        diagnostic["surface_radius_residuals_m"],
        expected_residuals_m,
        atol=1e-8,
    )
    assert diagnostic["maximum_absolute_surface_radius_residual_m"] == pytest.approx(
        np.max(np.abs(expected_residuals_m))
    )
    assert diagnostic["rms_surface_radius_residual_m"] == pytest.approx(
        np.sqrt(np.mean(expected_residuals_m**2))
    )
    np.testing.assert_allclose(
        diagnostic["tangent_altitudes_m"], tangent_altitude_m, atol=1e-7
    )
    assert diagnostic["maximum_tangent_direction_error_radians"] == pytest.approx(
        0.0, abs=2e-8
    )

    theta = np.asarray(diagnostic["ray_horizontal_angles"])
    phi = np.asarray(diagnostic["ray_out_of_plane_angles"])
    reference_z = np.asarray(diagnostic["reference_z"])
    reference_x = np.asarray(diagnostic["reference_x"])
    reference_y = np.asarray(diagnostic["reference_y"])
    reconstructed_directions = (
        np.cos(phi)[:, np.newaxis]
        * (
            np.cos(theta)[:, np.newaxis] * reference_z
            + np.sin(theta)[:, np.newaxis] * reference_x
        )
        - np.sin(phi)[:, np.newaxis] * reference_y
    )
    np.testing.assert_allclose(reconstructed_directions, tangent_directions, atol=2e-14)
    np.testing.assert_allclose(
        diagnostic["horizontal_location_residuals_radians"], 0.0, atol=2e-14
    )

    local_lookup = {
        global_index: local_index
        for local_index, global_index in enumerate(diagnostic["grid_indices"])
    }
    local_indices = [local_lookup[index] for index in grid_indices]
    horizontal_angles = np.asarray(diagnostic["horizontal_angles"])
    horizontal_distances_m = np.asarray(diagnostic["horizontal_distances_m"])
    internal_distance_m = np.diff(horizontal_distances_m[local_indices])[0]
    geoid_distance_m = expected_radius_m * np.diff(angles[grid_indices])[0]
    assert internal_distance_m == pytest.approx(geoid_distance_m, rel=1e-12)
    np.testing.assert_allclose(
        np.diff(horizontal_angles),
        np.diff(geometry.cumulative_angles[diagnostic["grid_indices"]]),
        atol=2e-14,
    )


def test_ray_policy_preserves_locations_outside_the_fitted_group_plane():
    geometry = orbital_geometry(warped=True)
    viewing = limb_viewing(geometry, np.array([-0.3, 0.0, 0.3]))
    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=120
    )
    diagnostic = engine.group_diagnostics[0]

    theta = np.asarray(diagnostic["ray_horizontal_angles"])
    phi = np.asarray(diagnostic["ray_out_of_plane_angles"])
    reference_z = np.asarray(diagnostic["reference_z"])
    reference_x = np.asarray(diagnostic["reference_x"])
    reference_y = np.asarray(diagnostic["reference_y"])
    reconstructed = (
        np.cos(phi)[:, np.newaxis]
        * (
            np.cos(theta)[:, np.newaxis] * reference_z
            + np.sin(theta)[:, np.newaxis] * reference_x
        )
        - np.sin(phi)[:, np.newaxis] * reference_y
    )
    observers = viewing.observer_positions_ecef_m
    looks = viewing.look_directions_ecef
    tangents = observers - np.sum(observers * looks, axis=1)[:, np.newaxis] * looks
    tangent_directions = tangents / np.linalg.norm(tangents, axis=1)[:, np.newaxis]

    assert np.max(np.abs(phi)) > 1e-5
    np.testing.assert_allclose(reconstructed, tangent_directions, atol=2e-14)
    np.testing.assert_allclose(diagnostic["tangent_altitudes_m"], 20_000.0, atol=1e-7)


def test_group_reference_position_is_the_mean_sample_track_location():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.33, 0.10]))
    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=120
    )

    reference = np.asarray(engine.group_diagnostics[0]["reference_position_ecef_m"])
    reference /= np.linalg.norm(reference)
    observers = viewing.observer_positions_ecef_m
    looks = viewing.look_directions_ecef
    distances = -np.einsum("ij,ij->i", observers, looks)
    tangents = observers + distances[:, np.newaxis] * looks
    tangent_directions = tangents / np.linalg.norm(tangents, axis=1)[:, np.newaxis]
    mean_angle = np.arctan2(tangent_directions[:, 0], tangent_directions[:, 2]).mean()
    expected = np.array([np.sin(mean_angle), 0.0, np.cos(mean_angle)])

    np.testing.assert_allclose(reference, expected, atol=2e-12)


def test_single_los_constructs_a_padded_spherical_geoid_track():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([0.0]))
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M,
        0.03,
        path_padding_angle=np.deg2rad(20.0),
        geoid=sk.SphericalGeoid(EARTH_RADIUS_M),
    )

    assert geometry.shape[0] >= 3
    np.testing.assert_allclose(geometry.surface_radii_m, EARTH_RADIUS_M, atol=1e-8)
    engine = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=60,
        group_padding_angle=0.0,
    )
    assert engine.group_diagnostics[0]["edge_clipping"] == (False, False)


def test_full_orbit_uses_separate_contiguous_local_windows():
    master_angles = np.linspace(0.0, 2 * np.pi, 361)
    directions = np.column_stack(
        (np.sin(master_angles), np.zeros_like(master_angles), np.cos(master_angles))
    )
    geometry = sk.OrbitalPlaneGeometry(
        EARTH_RADIUS_M, ALTITUDES_M, EARTH_RADIUS_M * directions
    )
    # The first and last observations have the same ECEF direction. Their
    # vertical-slice time order must select opposite ends of the non-periodic
    # master coordinate rather than mapping both to the first segment.
    sample_angles = np.array([0.0, 3.2, 2 * np.pi])
    tangent_radius = EARTH_RADIUS_M + 20_000.0
    observer_radius = EARTH_RADIUS_M + 700_000.0
    distance = np.sqrt(observer_radius**2 - tangent_radius**2)
    tangent_directions = np.column_stack(
        (
            np.sin(sample_angles),
            np.zeros_like(sample_angles),
            np.cos(sample_angles),
        )
    )
    forward = np.column_stack(
        (
            np.cos(sample_angles),
            np.zeros_like(sample_angles),
            -np.sin(sample_angles),
        )
    )
    tangents = tangent_radius * tangent_directions
    observers = tangents - distance * forward
    times = np.datetime64("2026-01-01", "ns") + np.arange(3) * np.timedelta64(60, "s")
    viewing = sk.OrbitalPlaneViewingGeometry.from_tangent_locations(
        times, observers, tangents, vertical_slice=np.arange(len(times))
    )

    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=20
    )
    diagnostics = engine.group_diagnostics
    centers = [int(np.median(item["grid_indices"])) for item in diagnostics]

    assert centers[0] < 60
    assert 150 < centers[1] < 220
    assert centers[2] > 300
    for item in diagnostics:
        indices = item["grid_indices"]
        assert indices == list(range(indices[0], indices[-1] + 1))
        assert (
            geometry.cumulative_angles[indices[-1]]
            - geometry.cumulative_angles[indices[0]]
            < np.pi
        )


def test_warped_track_scale_guard_and_local_refraction_coordinates():
    angles = np.linspace(-1.0, 1.0, 101)
    across = 0.08 * np.sin(4.0 * angles)
    directions = np.column_stack((np.sin(angles), across, np.cos(angles)))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    geometry = sk.OrbitalPlaneGeometry(
        EARTH_RADIUS_M, ALTITUDES_M, EARTH_RADIUS_M * directions
    )
    lower_indices = np.array([30, 40, 50, 60, 70])
    fraction = 0.37
    separation = np.arccos(
        np.sum(directions[lower_indices] * directions[lower_indices + 1], axis=1)
    )
    tangent_directions = (
        np.sin((1.0 - fraction) * separation)[:, np.newaxis]
        / np.sin(separation)[:, np.newaxis]
        * directions[lower_indices]
        + np.sin(fraction * separation)[:, np.newaxis]
        / np.sin(separation)[:, np.newaxis]
        * directions[lower_indices + 1]
    )
    along = np.gradient(directions, axis=0)[lower_indices]
    along -= (
        np.sum(along * tangent_directions, axis=1)[:, np.newaxis] * tangent_directions
    )
    along /= np.linalg.norm(along, axis=1, keepdims=True)
    tangent_radius = EARTH_RADIUS_M + 20_000.0
    observer_radius = EARTH_RADIUS_M + 700_000.0
    distance = np.sqrt(observer_radius**2 - tangent_radius**2)
    tangents = tangent_radius * tangent_directions
    observers = tangents - distance * along
    times = np.datetime64("2026-01-01", "ns") + np.arange(5) * np.timedelta64(5, "s")
    viewing = sk.OrbitalPlaneViewingGeometry.from_tangent_locations(
        times,
        observers,
        tangents,
        vertical_slice=np.arange(len(times)),
        geoid=sk.SphericalGeoid(EARTH_RADIUS_M),
    )

    with pytest.raises(ValueError, match="max_horizontal_scale_residual"):
        sk.OrbitalPlaneEngine(
            transmission_config(), geometry, viewing, time_group_duration_s=60
        )
    engine = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=60,
        max_horizontal_scale_residual=None,
    )
    diagnostic = engine.group_diagnostics[0]
    assert diagnostic["maximum_relative_horizontal_scale_residual"] > 0.05
    residuals = np.asarray(diagnostic["horizontal_distance_relative_residuals"])
    assert len(residuals) == len(diagnostic["grid_indices"]) - 1
    assert np.max(np.abs(residuals)) == pytest.approx(
        diagnostic["maximum_relative_horizontal_scale_residual"]
    )

    horizontal_angles = np.asarray(diagnostic["horizontal_angles"])
    grid_indices = np.asarray(diagnostic["grid_indices"])
    for theta, segment, local_fraction in zip(
        diagnostic["ray_horizontal_angles"],
        diagnostic["refractive_profile_segments"],
        diagnostic["refractive_profile_fractions"],
        strict=True,
    ):
        local_cell = int(np.flatnonzero(grid_indices == segment)[0])
        expected_fraction = (theta - horizontal_angles[local_cell]) / (
            horizontal_angles[local_cell + 1] - horizontal_angles[local_cell]
        )
        assert local_fraction == pytest.approx(expected_fraction, abs=2e-14)


def test_engine_rejects_away_pointing_rays():
    geometry = orbital_geometry()
    times = np.array(["2026-01-01"], dtype="datetime64[ns]")
    away = sk.OrbitalPlaneViewingGeometry(
        times,
        np.array([[0.0, 0.0, EARTH_RADIUS_M + 700_000.0]]),
        np.array([[0.0, 0.0, 1.0]]),
        vertical_slice=np.array([0]),
    )
    with pytest.raises(ValueError, match="looking toward Earth"):
        sk.OrbitalPlaneEngine(
            transmission_config(), geometry, away, time_group_duration_s=60
        )


def test_group_padding_extends_selected_internal_window():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    unpadded = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=60,
        group_padding_angle=0.0,
    )
    padding = np.deg2rad(10.0)
    padded = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=60,
        group_padding_angle=padding,
    )
    default_padded = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=60,
    )

    unpadded_indices = unpadded.group_diagnostics[0]["grid_indices"]
    padded_diagnostic = padded.group_diagnostics[0]
    padded_indices = padded_diagnostic["grid_indices"]
    assert len(unpadded_indices) == 2
    assert padded_indices[0] < unpadded_indices[0]
    assert padded_indices[-1] > unpadded_indices[-1]
    assert padded_diagnostic["padding_angle"] == pytest.approx(padding)
    assert padded.group_padding_angle == pytest.approx(padding)
    assert default_padded.group_padding_angle == pytest.approx(padding)
    assert default_padded.group_diagnostics[0]["grid_indices"] == padded_indices
    assert padded_diagnostic["path_edge_clipping"] == (False, False)
    assert padded_diagnostic["path_edge_clipping_per_los"] == [(False, False)]
    assert padded_diagnostic["uses_extended_horizontal_edge"] is False


def test_group_diagnostics_report_actual_extended_horizontal_edge_usage():
    source_geometry = orbital_geometry()
    viewing = limb_viewing(source_geometry, np.array([0.0]))
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M,
        0.025,
        path_padding_angle=0.0,
        geoid=sk.SphericalGeoid(EARTH_RADIUS_M),
    )
    engine = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=60,
        group_padding_angle=0.0,
    )

    diagnostic = engine.group_diagnostics[0]
    assert diagnostic["edge_clipping"] == (False, False)
    assert diagnostic["path_edge_clipping"] == (True, True)
    assert diagnostic["path_edge_clipping_per_los"] == [(True, True)]
    assert diagnostic["uses_extended_horizontal_edge"] is True


@pytest.mark.parametrize("value", [True, -0.1, np.pi, np.inf])
def test_group_padding_rejects_invalid_values(value):
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    error = TypeError if isinstance(value, bool) else ValueError
    with pytest.raises(error, match="group_padding_angle"):
        sk.OrbitalPlaneEngine(
            transmission_config(),
            geometry,
            viewing,
            time_group_duration_s=60,
            group_padding_angle=value,
        )


def test_time_groups_are_half_open_and_stitch_original_observation_order():
    geometry = orbital_geometry()
    offsets_s = np.array([120, 0, 60, 59])
    times = np.datetime64("2026-01-01", "ns") + offsets_s * np.timedelta64(1, "s")
    viewing = limb_viewing(geometry, np.array([-0.2, -0.1, 0.0, 0.1]), times)
    engine = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=60,
    )
    diagnostics = engine.group_diagnostics
    assert [item["observation_indices"] for item in diagnostics] == [
        [1, 3],
        [2],
        [0],
    ]
    selected_ns = times[[1, 3]].astype(np.int64)
    expected_mean = np.datetime64(int(selected_ns.sum() // selected_ns.size), "ns")
    assert diagnostics[0]["reference_time"] == expected_mean
    result = engine.calculate_radiance(
        raw_atmosphere(geometry, engine._config, calculate_derivatives=False)
    )
    np.testing.assert_array_equal(result.time, times)
    assert result.radiance.shape == (1, 4, 1)


def test_time_groups_keep_complete_vertical_slices():
    geometry = orbital_geometry()
    offsets_s = np.array([0, 15, 20, 25])
    times = np.datetime64("2026-01-01", "ns") + offsets_s * np.timedelta64(1, "s")
    viewing = limb_viewing(
        geometry,
        np.array([-0.2, -0.18, 0.1, 0.12]),
        times,
        vertical_slice=np.array([0, 0, 1, 1]),
    )

    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=10
    )

    assert [item["observation_indices"] for item in engine.group_diagnostics] == [
        [0, 1],
        [2, 3],
    ]


def test_split_preserves_global_slice_atomic_time_groups():
    geometry = orbital_geometry()
    offsets_s = np.arange(10) * 13
    times = np.datetime64("2026-01-01", "ns") + offsets_s * np.timedelta64(1, "s")
    viewing = limb_viewing(geometry, np.linspace(-0.3, 0.3, 10), times)
    duration_s = 40

    full_engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=duration_s
    )
    chunks = viewing.split(2, time_group_duration_s=duration_s)
    chunk_engines = [
        sk.OrbitalPlaneEngine(
            transmission_config(), geometry, chunk, time_group_duration_s=duration_s
        )
        for chunk in chunks
    ]

    full_group_times = [
        item["reference_time"] for item in full_engine.group_diagnostics
    ]
    chunk_group_times = [
        item["reference_time"]
        for engine in chunk_engines
        for item in engine.group_diagnostics
    ]
    assert chunk_group_times == full_group_times
    assert all(chunk.time_bin_origin == viewing.time_bin_origin for chunk in chunks)


def test_group_reference_time_uses_floor_for_sub_nanosecond_means():
    geometry = orbital_geometry()
    times = np.array([np.datetime64(-1, "ns"), np.datetime64(0, "ns")])
    viewing = limb_viewing(geometry, np.array([-0.1, 0.1]), times)
    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=1
    )

    assert engine.group_diagnostics[0]["reference_time"] == np.datetime64(-1, "ns")


def test_los_optical_depth_is_stitched_and_matches_transmission():
    geometry = orbital_geometry()
    offsets_s = np.array([60, 0])
    times = np.datetime64("2026-01-01", "ns") + offsets_s * np.timedelta64(1, "s")
    viewing = limb_viewing(geometry, np.array([-0.2, 0.2]), times)
    config = transmission_config()
    config.output_los_optical_depth = True
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=20)
    result = engine.calculate_radiance(
        raw_atmosphere(geometry, config, calculate_derivatives=False)
    )

    assert result.los_optical_depth.shape == (1, 2)
    assert np.all(result.los_optical_depth > 0)
    np.testing.assert_allclose(
        result.los_optical_depth.values,
        -np.log(result.radiance.sel(stokes="I").values),
        rtol=2e-12,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    "emission_source",
    [sk.EmissionSource.Standard, sk.EmissionSource.VolumeEmissionRate],
)
def test_supported_emission_sources_run_across_multiple_groups(emission_source):
    geometry = orbital_geometry()
    offsets_s = np.array([60, 0])
    times = np.datetime64("2026-01-01", "ns") + offsets_s * np.timedelta64(1, "s")
    viewing = limb_viewing(geometry, np.array([-0.2, 0.2]), times)
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = emission_source
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=20)
    atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)
    atmosphere.storage.emission_source[:] = 2e-6

    result = engine.calculate_radiance(atmosphere)

    np.testing.assert_array_equal(result.time, times)
    assert np.all(np.isfinite(result.radiance))
    assert np.all(result.radiance.sel(stokes="I") > 0)


def test_exact_single_scattering_runs_across_multiple_groups():
    class OverheadSolarHandler:
        def target_solar_angles(self, *_args):
            return 0.0, 0.0

    geometry = orbital_geometry()
    offsets_s = np.array([60, 0])
    times = np.datetime64("2026-01-01", "ns") + offsets_s * np.timedelta64(1, "s")
    viewing = limb_viewing(geometry, np.array([-0.2, 0.2]), times)
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=20,
        solar_handler=OverheadSolarHandler(),
    )
    atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)
    atmosphere.storage.ssa[:] = 0.5
    atmosphere.storage.solar_irradiance[:] = 1.0
    atmosphere.leg_coeff.a1[0] = 1.0

    result = engine.calculate_radiance(atmosphere)

    np.testing.assert_array_equal(result.time, times)
    assert np.all(np.isfinite(result.radiance))
    assert np.all(result.radiance.sel(stokes="I") > 0)


def test_successive_orders_runs_across_multiple_groups_without_single_scatter():
    class OverheadSolarHandler:
        def target_solar_angles(self, *_args):
            return 0.0, 0.0

    geometry = orbital_geometry()
    offsets_s = np.array([60, 0])
    times = np.datetime64("2026-01-01", "ns") + offsets_s * np.timedelta64(1, "s")
    viewing = limb_viewing(geometry, np.array([-0.2, 0.2]), times)
    config = sk.Config()
    config.num_threads = 1
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
    config.num_sza = 2
    config.successive_orders_horizontal_angle_grid_radians = np.array([-0.1, 0.1])
    config.successive_orders_altitude_grid_m = np.array([5_000.0, 25_000.0, 55_000.0])
    config.num_successive_orders_incoming = 6
    config.num_successive_orders_outgoing = 6
    config.num_successive_orders_iterations = 2
    config.successive_orders_relative_tolerance = 0.0
    config.successive_orders_absolute_tolerance = 0.0
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=20,
        solar_handler=OverheadSolarHandler(),
    )
    atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)
    atmosphere.storage.ssa[:] = 0.5
    atmosphere.storage.solar_irradiance[:] = 1.0
    atmosphere.leg_coeff.a1[0] = 1.0

    result = engine.calculate_radiance(atmosphere)

    np.testing.assert_array_equal(result.time, times)
    assert np.all(np.isfinite(result.radiance))
    assert np.all(result.radiance.sel(stokes="I") > 0)

    config.successive_orders_horizontal_angle_grid_radians = np.array([-0.3, 0.1])
    with pytest.raises(ValueError, match="inside every local group grid"):
        sk.OrbitalPlaneEngine(
            config,
            geometry,
            viewing,
            time_group_duration_s=20,
            solar_handler=OverheadSolarHandler(),
        )


def test_engine_exposes_normalized_solar_and_execution_settings_without_mutation():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.1, 0.1]))
    vectors = np.array([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
    original = vectors.copy()
    engine = sk.OrbitalPlaneEngine(
        transmission_config(),
        geometry,
        viewing,
        time_group_duration_s=20,
        sun_vectors_ecef=vectors,
        refraction_wavelength_nm=550.0,
        refraction_co2_ppm=425.0,
        derivative_execution="streaming",
        max_horizontal_scale_residual=0.02,
        max_eager_jacobian_bytes=123_456,
    )

    np.testing.assert_array_equal(vectors, original)
    np.testing.assert_allclose(np.linalg.norm(engine.sun_vectors_ecef, axis=1), 1.0)
    assert not engine.sun_vectors_ecef.flags.writeable
    with pytest.raises(ValueError, match="WRITEABLE"):
        engine.sun_vectors_ecef.setflags(write=True)
    assert engine.time_group_duration_s == 20
    assert engine.refraction_wavelength_nm == 550
    assert engine.refraction_co2_ppm == 425
    assert engine.derivative_execution == "streaming"
    assert engine.max_horizontal_scale_residual == pytest.approx(0.02)
    assert engine.max_eager_jacobian_bytes == 123_456
    assert engine.num_groups == 2


@pytest.mark.parametrize(
    ("kwargs", "error", "match"),
    [
        ({"time_group_duration_s": True}, TypeError, "time_group_duration_s"),
        ({"time_group_duration_s": 0.0}, ValueError, "time_group_duration_s"),
        ({"time_group_duration_s": 1e-12}, ValueError, "nanosecond"),
        (
            {"time_group_duration_s": 60, "refraction_wavelength_nm": 0.0},
            ValueError,
            "refraction_wavelength_nm",
        ),
        (
            {"time_group_duration_s": 60, "refraction_wavelength_nm": np.bool_(True)},
            TypeError,
            "refraction_wavelength_nm",
        ),
        (
            {"time_group_duration_s": 60, "refraction_co2_ppm": -1.0},
            ValueError,
            "refraction_co2_ppm",
        ),
        (
            {"time_group_duration_s": 60, "max_horizontal_scale_residual": True},
            TypeError,
            "max_horizontal_scale_residual",
        ),
        (
            {"time_group_duration_s": 60, "max_horizontal_scale_residual": -1.0},
            ValueError,
            "max_horizontal_scale_residual",
        ),
        (
            {"time_group_duration_s": 60, "max_eager_jacobian_bytes": 0},
            ValueError,
            "max_eager_jacobian_bytes",
        ),
        (
            {"time_group_duration_s": 60, "max_eager_jacobian_bytes": True},
            TypeError,
            "max_eager_jacobian_bytes",
        ),
        (
            {"time_group_duration_s": 60, "max_eager_jacobian_bytes": 1.5},
            TypeError,
            "max_eager_jacobian_bytes",
        ),
    ],
)
def test_engine_rejects_invalid_scalar_settings(kwargs, error, match):
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    with pytest.raises(error, match=match):
        sk.OrbitalPlaneEngine(transmission_config(), geometry, viewing, **kwargs)


def test_engine_validates_solar_vector_shape_and_handler_results():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.1, 0.1]))
    with pytest.raises(ValueError, match=r"shape \(2, 3\)"):
        sk.OrbitalPlaneEngine(
            transmission_config(),
            geometry,
            viewing,
            time_group_duration_s=20,
            sun_vectors_ecef=np.ones(3),
        )

    class InvalidSolarHandler:
        def target_solar_angles(self, *_args):
            return np.array([45.0, np.nan])

    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    with pytest.raises(ValueError, match="two finite angles"):
        sk.OrbitalPlaneEngine(
            config,
            geometry,
            viewing,
            time_group_duration_s=20,
            solar_handler=InvalidSolarHandler(),
        )


def test_engine_rejects_unsupported_config_before_calling_solar_handler():
    class FailingSolarHandler:
        def target_solar_angles(self, *_args):
            msg = "handler must not be called"
            raise AssertionError(msg)

    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    with pytest.raises(NotImplementedError, match="successive-orders"):
        sk.OrbitalPlaneEngine(
            config,
            geometry,
            viewing,
            time_group_duration_s=60,
            solar_handler=FailingSolarHandler(),
        )


def test_warped_ground_track_fits_oriented_local_planes():
    geometry = orbital_geometry(warped=True)
    viewing = limb_viewing(geometry, np.array([-0.3, 0.0, 0.3]))
    engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=120
    )
    diagnostics = engine.group_diagnostics
    assert any(item["max_out_of_plane_angle"] > 0 for item in diagnostics)
    for item in diagnostics:
        assert len(item["grid_indices"]) >= 2
        assert item["window_expanded"] is False
        assert np.isfinite(item["reference_z"]).all()
        assert np.isfinite(item["reference_x"]).all()
        assert item["reference_position_ecef_m"] == item["reference_position_ecef"]
        assert item["max_out_of_plane_angle_rad"] == item["max_out_of_plane_angle"]

    diagnostics[0]["observation_indices"].append(999)
    assert 999 not in engine.group_diagnostics[0]["observation_indices"]


def test_ciddor_rust_parity_for_humidity_co2_and_orbital_variation():
    temperature = np.array([[288.15, 250.0], [280.0, 240.0]])
    pressure = np.array([[101_325.0, 20_000.0], [80_000.0, 10_000.0]])
    humidity = np.array([[0.0, 0.001], [0.01, 0.02]])
    for supplied_humidity in (humidity, None):
        expected_humidity = (
            np.zeros_like(pressure) if supplied_humidity is None else supplied_humidity
        )
        expected = ciddor_index_of_refraction(
            temperature, pressure, expected_humidity, 425.0, 600.0
        )
        actual = orbital_ciddor_index(
            temperature, pressure, supplied_humidity, 425.0, 600.0
        )
        np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-16)


def test_refraction_cache_reuses_engines_and_retraces_changed_profiles():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.1, 0.1]))
    config = transmission_config(refraction=True)
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=120)
    atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)
    profile = 1 + 2.7e-4 * np.exp(-ALTITUDES_M / 8_000)
    atmosphere.refractive_index = profile

    first = engine.calculate_radiance(atmosphere).radiance.values.copy()
    first_diagnostics = engine.group_diagnostics
    identities = [item["engine_identity"] for item in first_diagnostics]
    workspaces = [item["atmosphere_workspace_identity"] for item in first_diagnostics]
    assert all(item["geometry_refresh_count"] == 1 for item in first_diagnostics)
    assert all(value != 0 for value in workspaces)

    atmosphere.storage.total_extinction[:] *= 1.01
    engine.calculate_radiance(atmosphere)
    second_diagnostics = engine.group_diagnostics
    assert [item["engine_identity"] for item in second_diagnostics] == identities
    assert [
        item["atmosphere_workspace_identity"] for item in second_diagnostics
    ] == workspaces
    assert all(item["geometry_refresh_count"] == 1 for item in second_diagnostics)

    atmosphere.refractive_index = profile + np.linspace(0, 2e-8, profile.size)
    changed = engine.calculate_radiance(atmosphere).radiance.values
    assert all(item["geometry_refresh_count"] == 2 for item in engine.group_diagnostics)
    assert not np.array_equal(first, changed)


def test_refraction_tangent_altitude_limit_traces_high_los_straight():
    geometry = orbital_geometry()
    viewing = limb_viewing(
        geometry,
        np.array([-0.05, 0.05]),
        tangent_altitude_m=np.array([20_000.0, 35_000.0]),
    )
    profile = 1.0 + 2.7e-4 * np.exp(-ALTITUDES_M / 8_000.0)

    def calculate(*, refraction: bool, limit_m: float) -> tuple[np.ndarray, list[dict]]:
        config = transmission_config(refraction=refraction)
        config.los_refraction_max_tangent_altitude_m = limit_m
        atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)
        atmosphere.refractive_index = profile
        engine = sk.OrbitalPlaneEngine(
            config,
            geometry,
            viewing,
            time_group_duration_s=120,
        )
        return (
            engine.calculate_radiance(atmosphere).radiance.values,
            engine.group_diagnostics,
        )

    limited, diagnostics = calculate(refraction=True, limit_m=25_000.0)
    fully_refracted, _ = calculate(refraction=True, limit_m=np.inf)
    straight, _ = calculate(refraction=False, limit_m=25_000.0)

    np.testing.assert_allclose(limited[:, 0, :], fully_refracted[:, 0, :])
    np.testing.assert_allclose(limited[:, 1, :], straight[:, 1, :])
    assert not np.array_equal(limited[:, 0, :], straight[:, 0, :])
    assert sum(item["refracted_observation_count"] for item in diagnostics) == 1
    assert all(
        item["max_refraction_tangent_altitude_m"] == 25_000.0 for item in diagnostics
    )


def test_parallel_group_refraction_refresh_matches_serial():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.45, -0.15, 0.15, 0.45]))
    horizontal_scale = np.linspace(0.9, 1.1, geometry.shape[0])[:, np.newaxis]
    base_profile = 1.0 + horizontal_scale * (2.7e-4 * np.exp(-ALTITUDES_M / 8_000.0))
    changed_profile = base_profile + horizontal_scale * np.linspace(
        0.0, 2.0e-8, geometry.shape[1]
    )

    def calculate(num_threads: int) -> tuple[np.ndarray, np.ndarray, list[dict]]:
        config = transmission_config(refraction=True)
        config.num_threads = num_threads
        config.threading_lib = sk.ThreadingLib.Rayon
        config.threading_model = sk.ThreadingModel.Wavelength
        atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)
        atmosphere.refractive_index = base_profile
        engine = sk.OrbitalPlaneEngine(
            config,
            geometry,
            viewing,
            time_group_duration_s=20,
        )
        initial = engine.calculate_radiance(atmosphere).radiance.values.copy()
        atmosphere.refractive_index = changed_profile
        changed = engine.calculate_radiance(atmosphere).radiance.values.copy()
        return initial, changed, engine.group_diagnostics

    serial_initial, serial_changed, serial_diagnostics = calculate(1)
    parallel_initial, parallel_changed, parallel_diagnostics = calculate(2)

    np.testing.assert_array_equal(parallel_initial, serial_initial)
    np.testing.assert_array_equal(parallel_changed, serial_changed)
    assert not np.array_equal(serial_initial, serial_changed)
    assert len(serial_diagnostics) == len(parallel_diagnostics) == 4
    assert all(item["geometry_refresh_count"] == 2 for item in serial_diagnostics)
    assert all(item["geometry_refresh_count"] == 2 for item in parallel_diagnostics)


def test_native_refractive_index_change_retraces_only_affected_group():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.3, 0.3]))
    config = transmission_config(refraction=True)
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=20)
    atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)
    refractive_index = np.ones(geometry.shape)
    atmosphere.refractive_index = refractive_index

    engine.calculate_radiance(atmosphere)
    assert [item["geometry_refresh_count"] for item in engine.group_diagnostics] == [
        1,
        1,
    ]

    # The first nominal tangent lies exactly at master row 20. Change both
    # columns adjoining that coordinate to make the expectation independent
    # of which segment owns an exact grid boundary.
    refractive_index[19:21] += 1e-8
    engine.calculate_radiance(atmosphere)

    assert [item["geometry_refresh_count"] for item in engine.group_diagnostics] == [
        2,
        1,
    ]


def test_explicit_refractive_index_takes_precedence_over_missing_state():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config(refraction=True)
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)

    with pytest.raises(ValueError, match="pressure_pa and temperature_k"):
        engine.calculate_radiance(atmosphere)

    atmosphere.refractive_index = 1 + 2.7e-4 * np.exp(-ALTITUDES_M / 8_000)
    result = engine.calculate_radiance(atmosphere)
    assert np.all(np.isfinite(result.radiance))


def test_automatic_refraction_and_state_derivative_dimensions():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config(refraction=True)
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.array([500.0]), legendre_derivative=False
    )
    pressure_profile = np.array([101_325.0, 26_000.0, 5_500.0, 600.0, 10.0])
    temperature_profile = np.array([288.0, 250.0, 225.0, 230.0, 250.0])
    atmosphere.pressure_pa = pressure_profile
    atmosphere.temperature_k = temperature_profile
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    linearization = engine.linearize(atmosphere)
    assert linearization.parameter_dims["pressure_pa"] == ("altitude",)
    assert linearization.parameter_dims["temperature_k"] == ("altitude",)
    refresh_count = engine.group_diagnostics[0]["geometry_refresh_count"]
    linearization.jvp(linearization.tangent_template)
    linearization.vjp(xr.ones_like(linearization.value))
    assert engine.group_diagnostics[0]["geometry_refresh_count"] == refresh_count

    native = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.array([500.0]), legendre_derivative=False
    )
    native.pressure_pa = np.broadcast_to(pressure_profile, geometry.shape).copy()
    native.temperature_k = np.broadcast_to(temperature_profile, geometry.shape).copy()
    native["rayleigh"] = sk.constituent.Rayleigh()
    native_linearization = engine.linearize(native)
    assert native_linearization.parameter_dims["pressure_pa"] == (
        "orbital_position",
        "altitude",
    )


def test_stitched_jvp_vjp_are_adjoint_and_lazy_jacobian_agrees():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.2, 0.15]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=20)
    atmosphere = raw_atmosphere(geometry, config)
    linearization = engine.linearize(atmosphere)
    tangent = linearization.tangent_template[["extinction"]]
    rng = np.random.default_rng(1234)
    tangent["extinction"].values[:] = rng.normal(size=geometry.shape)
    cotangent = xr.DataArray(
        rng.normal(size=linearization.value.shape),
        dims=linearization.value.dims,
        coords=linearization.value.coords,
    )
    jvp = linearization.jvp(tangent)
    vjp = linearization.vjp(cotangent, parameters=["extinction"])
    np.testing.assert_allclose(
        (jvp * cotangent).sum(),
        (tangent.extinction * vjp.extinction).sum(),
        rtol=2e-11,
        atol=1e-12,
    )
    expected = np.einsum(
        "ij...,ij->...", linearization.jacobian.extinction, tangent.extinction
    )
    np.testing.assert_allclose(jvp, expected, rtol=2e-11, atol=1e-12)


def test_eager_jacobian_memory_guard_leaves_native_vjp_available():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.2, 0.15]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=20,
        max_eager_jacobian_bytes=1_024,
    )
    atmosphere = raw_atmosphere(geometry, config)

    estimate = engine.estimate_eager_jacobian_bytes(atmosphere)
    assert estimate > engine.max_eager_jacobian_bytes
    with pytest.raises(MemoryError, match="derivatives=False"):
        engine.calculate_radiance(atmosphere)

    radiance_only = engine.calculate_radiance(atmosphere, derivatives=False)
    assert "radiance" in radiance_only
    assert not any(name.startswith("wf_") for name in radiance_only.data_vars)
    for diagnostics in engine.group_diagnostics:
        assert diagnostics["resident_volume_derivative_mappings"] == []
        assert diagnostics["resident_surface_derivative_mappings"] == []

    linearization = engine.linearize(atmosphere)
    result = linearization.vjp(xr.ones_like(linearization.value))
    assert "extinction" in result
    with pytest.raises(MemoryError, match="max_eager_jacobian_bytes=None"):
        _ = linearization.jacobian


def test_eager_jacobian_limit_can_be_explicitly_disabled():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=60,
        max_eager_jacobian_bytes=None,
    )
    result = engine.calculate_radiance(raw_atmosphere(geometry, config))
    assert "wf_extinction" in result


def test_engine_rejects_structural_config_mutation_and_stokes_mismatch():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)

    config.output_los_optical_depth = not config.output_los_optical_depth
    with pytest.raises(ValueError, match=r"Config.*modified"):
        engine.calculate_radiance(atmosphere)

    config3 = transmission_config()
    config3.num_stokes = 3
    engine3_atmosphere = raw_atmosphere(geometry, config3, calculate_derivatives=False)
    fresh_engine = sk.OrbitalPlaneEngine(
        transmission_config(), geometry, viewing, time_group_duration_s=60
    )
    with pytest.raises(ValueError, match="Stokes dimensions"):
        fresh_engine.calculate_radiance(engine3_atmosphere)


def test_linearization_rejects_structural_config_mutation():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    linearization = engine.linearize(raw_atmosphere(geometry, config))
    config.output_los_optical_depth = not config.output_los_optical_depth

    with pytest.raises(ValueError, match=r"Config.*modified"):
        linearization.vjp(xr.ones_like(linearization.value))


def test_linearization_rejects_later_refractive_geometry_state():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config(refraction=True)
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    atmosphere = raw_atmosphere(geometry, config)
    atmosphere.refractive_index = np.ones(len(ALTITUDES_M))
    linearization = engine.linearize(atmosphere)
    tangent = linearization.tangent_template[["extinction"]]
    tangent.extinction.values[:] = 1.0
    linearization.jvp(tangent)

    other = raw_atmosphere(geometry, config, calculate_derivatives=False)
    other.refractive_index = 1.0 + 2.7e-4 * np.exp(-ALTITUDES_M / 8_000.0)
    engine.calculate_radiance(other)

    with pytest.raises(sk.StaleLinearizationError, match="orbital engine"):
        linearization.jvp(tangent)
    with pytest.raises(sk.StaleLinearizationError, match="orbital engine"):
        linearization.vjp(xr.ones_like(linearization.value))


def test_linearization_rejects_later_spatial_surface_state():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)

    def atmosphere_with_albedo(albedo):
        atmosphere = sk.Atmosphere(
            geometry,
            config,
            wavelengths_nm=np.array([500.0]),
            legendre_derivative=False,
        )
        atmosphere["optics"] = sk.constituent.Manual(
            np.full((*geometry.shape, 1), 1e-5),
            np.zeros((*geometry.shape, 1)),
        )
        atmosphere["surface"] = sk.constituent.LambertianSurface2D(albedo)
        return atmosphere

    atmosphere = atmosphere_with_albedo(np.linspace(0.1, 0.3, geometry.shape[0]))
    linearization = engine.linearize(atmosphere)
    other = atmosphere_with_albedo(np.linspace(0.2, 0.4, geometry.shape[0]))
    engine.calculate_radiance(other)

    with pytest.raises(sk.StaleLinearizationError, match="spatial surface"):
        linearization.vjp(xr.ones_like(linearization.value))


def test_engine_rejects_equal_but_distinct_atmosphere_geometry_early():
    geometry = orbital_geometry()
    other_geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config(refraction=True)
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    atmosphere = raw_atmosphere(other_geometry, config, calculate_derivatives=False)
    atmosphere.refractive_index = np.ones(len(ALTITUDES_M))

    with pytest.raises(ValueError, match="same Geometry2D object"):
        engine.calculate_radiance(atmosphere)


def test_engine_rejects_one_argument_non_lambertian_surface():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    atmosphere = raw_atmosphere(geometry, config, calculate_derivatives=False)
    atmosphere.surface.brdf = PyKokhanovsky(atmosphere.nstokes)

    with pytest.raises(ValueError, match="Lambertian surfaces only"):
        engine.calculate_radiance(atmosphere)


def test_orbital_linearize_materializes_constituents_once():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        legendre_derivative=False,
    )
    atmosphere.pressure_pa = np.array([101_325.0, 26_000.0, 5_500.0, 600.0, 10.0])
    atmosphere.temperature_k = np.array([288.0, 250.0, 225.0, 230.0, 250.0])
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    revision = atmosphere.revision

    engine.linearize(atmosphere)

    assert atmosphere.revision == revision + 1


def test_resident_and_streaming_vjp_agree_and_keep_group_engines():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.2, 0.15]))
    config = transmission_config()
    cotangent = None
    results = {}

    for execution in ("resident", "streaming"):
        engine = sk.OrbitalPlaneEngine(
            config,
            geometry,
            viewing,
            time_group_duration_s=20,
            derivative_execution=execution,
        )
        identities = [item["engine_identity"] for item in engine.group_diagnostics]
        linearization = engine.linearize(
            raw_atmosphere(geometry, config), prepare_parameters=("extinction",)
        )
        if execution == "streaming":
            assert all(
                not item["resident_volume_derivative_mappings"]
                for item in engine.group_diagnostics
            )
        if cotangent is None:
            cotangent = xr.ones_like(linearization.value)
        results[execution] = linearization.vjp(
            cotangent, parameters=["extinction"]
        ).extinction
        assert [
            item["engine_identity"] for item in engine.group_diagnostics
        ] == identities

    np.testing.assert_allclose(
        results["resident"], results["streaming"], rtol=2e-11, atol=1e-12
    )


def test_resident_group_workspaces_contain_only_selected_derivative_mappings():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.2, 0.15]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=20,
        derivative_execution="resident",
    )
    atmosphere = raw_atmosphere(geometry, config)
    linearization = engine.linearize(atmosphere, prepare_parameters=["extinction"])
    cotangent = xr.ones_like(linearization.value)
    first_identities = [
        item["atmosphere_workspace_identity"] for item in engine.group_diagnostics
    ]
    first_update_counts = [
        item["atmosphere_update_count"] for item in engine.group_diagnostics
    ]
    for diagnostics in engine.group_diagnostics:
        assert diagnostics["resident_volume_derivative_mappings"] == ["wf_extinction"]
        assert diagnostics["resident_surface_derivative_mappings"] == []

    linearization.vjp(cotangent, parameters=["extinction"])
    linearization.vjp(cotangent, parameters=["extinction"])
    assert [
        item["atmosphere_workspace_identity"] for item in engine.group_diagnostics
    ] == first_identities
    assert [
        item["atmosphere_update_count"] for item in engine.group_diagnostics
    ] == first_update_counts

    linearization.vjp(cotangent, parameters=["ssa"])
    for diagnostics in engine.group_diagnostics:
        assert diagnostics["resident_volume_derivative_mappings"] == ["wf_ssa"]

    tangent = linearization.tangent_template[["extinction"]]
    tangent.extinction.values[:] = 1.0
    linearization.jvp(tangent)
    for diagnostics in engine.group_diagnostics:
        assert diagnostics["resident_volume_derivative_mappings"] == ["wf_extinction"]

    changed_counts = [
        item["atmosphere_update_count"] for item in engine.group_diagnostics
    ]
    atmosphere.mark_changed()
    engine.linearize(atmosphere, prepare_parameters=["extinction"])
    assert [item["atmosphere_update_count"] for item in engine.group_diagnostics] == [
        count + 1 for count in changed_counts
    ]


def test_surface_only_update_preserves_group_volume_state():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.1, 0.1]))
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
    config.num_sza = 2
    config.successive_orders_altitude_grid_m = np.array([5_000.0, 25_000.0, 55_000.0])
    config.num_successive_orders_incoming = 6
    config.num_successive_orders_outgoing = 6
    config.num_successive_orders_iterations = 40
    config.successive_orders_relative_tolerance = 1.0e-11
    config.successive_orders_absolute_tolerance = 1.0e-13

    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([600.0]),
        legendre_derivative=False,
    )
    extinction = np.full((*geometry.shape, 1), 1.0e-5)
    ssa = np.full_like(extinction, 0.9)
    legendre = np.zeros((atmosphere.storage.leg_coeff.shape[0], *extinction.shape))
    legendre[0] = 1.0
    optics = sk.constituent.Manual(extinction, ssa, legendre)
    atmosphere["optics"] = optics
    surface = sk.constituent.LambertianSurface2D(
        np.linspace(0.1, 0.3, geometry.shape[0])
    )
    atmosphere["surface"] = surface
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=60,
        sun_vectors_ecef=np.array([[0.0, 0.0, 1.0]]),
    )

    first = engine.linearize(atmosphere, prepare_parameters=("surface_albedo",))
    first_diagnostics = engine.group_diagnostics
    assert all(item["volume_update_count"] == 1 for item in first_diagnostics)
    assert all(item["surface_only_update_count"] == 0 for item in first_diagnostics)

    surface.albedo = np.linspace(0.2, 0.4, geometry.shape[0])
    second = engine.linearize(atmosphere, prepare_parameters=("surface_albedo",))
    second_diagnostics = engine.group_diagnostics
    assert [item["volume_update_count"] for item in second_diagnostics] == [
        item["volume_update_count"] for item in first_diagnostics
    ]
    assert [item["surface_only_update_count"] for item in second_diagnostics] == [
        item["surface_only_update_count"] + 1 for item in first_diagnostics
    ]
    assert not np.array_equal(first.value.values, second.value.values)

    optics.extinction[...] *= 1.2
    third = engine.linearize(atmosphere, prepare_parameters=("surface_albedo",))
    third_diagnostics = engine.group_diagnostics
    assert [item["volume_update_count"] for item in third_diagnostics] == [
        item["volume_update_count"] + 1 for item in second_diagnostics
    ]
    assert [item["surface_only_update_count"] for item in third_diagnostics] == [
        item["surface_only_update_count"] for item in second_diagnostics
    ]
    assert not np.array_equal(second.value.values, third.value.values)

    fresh_engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=60,
        sun_vectors_ecef=np.array([[0.0, 0.0, 1.0]]),
    )
    fresh_third = fresh_engine.linearize(
        atmosphere, prepare_parameters=("surface_albedo",)
    )
    xr.testing.assert_allclose(third.value, fresh_third.value, rtol=2.0e-11)


def test_parallel_group_construction_and_scheduler_match_serial_products_repeatedly():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.2, 0.0, 0.2]))

    def calculate(num_threads: int):
        config = transmission_config()
        config.num_threads = num_threads
        config.threading_lib = sk.ThreadingLib.Rayon
        config.threading_model = sk.ThreadingModel.Wavelength
        config.wavelength_batch_size = 2
        atmosphere = sk.Atmosphere(
            geometry,
            config,
            wavelengths_nm=np.array([450.0, 500.0, 550.0]),
            legendre_derivative=False,
        )
        atmosphere.storage.total_extinction[:] = np.linspace(
            0.8e-5,
            1.2e-5,
            atmosphere.num_locations,
        )[:, np.newaxis]
        atmosphere.storage.ssa[:] = 0.0
        engine = sk.OrbitalPlaneEngine(
            config,
            geometry,
            viewing,
            time_group_duration_s=20,
        )
        linearization = engine.linearize(atmosphere, prepare_parameters=("extinction",))
        tangent = linearization.tangent_template[["extinction"]]
        tangent.extinction.values[:] = np.linspace(
            -0.4, 0.7, tangent.extinction.size
        ).reshape(tangent.extinction.shape)
        cotangent = xr.ones_like(linearization.value)
        cotangent.values[:] = np.linspace(-0.5, 0.8, cotangent.size).reshape(
            cotangent.shape
        )
        return (
            linearization.value.copy(),
            linearization.jvp(tangent),
            linearization.vjp(cotangent, parameters=("extinction",)).extinction,
            engine.group_diagnostics,
        )

    serial_value, serial_jvp, serial_vjp, serial_diagnostics = calculate(1)
    assert all(
        not diagnostics["composite_wavelength_scheduler"]
        for diagnostics in serial_diagnostics
    )

    # Repeatedly transfer the independently constructed native group engines
    # from Rayon workers, exercise every native product, and destroy the whole
    # object graph before constructing the next engine. Four threads with only
    # three groups also exercises a construction pool smaller than the shared
    # calculation pool.
    for num_threads in (2, 4, 2, 4, 2, 4):
        threaded_value, threaded_jvp, threaded_vjp, threaded_diagnostics = calculate(
            num_threads
        )
        xr.testing.assert_allclose(threaded_value, serial_value, rtol=1e-13, atol=1e-14)
        xr.testing.assert_allclose(threaded_jvp, serial_jvp, rtol=1e-13, atol=1e-14)
        xr.testing.assert_allclose(threaded_vjp, serial_vjp, rtol=1e-13, atol=1e-14)
        assert all(
            diagnostics["composite_wavelength_scheduler"]
            and diagnostics["composite_wavelength_threads"] == num_threads
            for diagnostics in threaded_diagnostics
        )
        del threaded_value, threaded_jvp, threaded_vjp, threaded_diagnostics
        gc.collect()


def test_orbital_aerosol_altitude_auxiliary_parameter_is_adjoint_and_selected():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([-0.15, 0.12]))
    config = transmission_config()
    optical_property = sk.optical.HenyeyGreenstein(
        db=xr.Dataset(
            {
                "xs_total": (
                    ("radius", "wavelength_nm"),
                    np.array([[2.0e-12, 3.0e-12], [6.0e-12, 8.0e-12]]),
                ),
                "ssa": (
                    ("radius", "wavelength_nm"),
                    np.array([[1.0e-12, 1.5e-12], [3.0e-12, 4.0e-12]]),
                ),
                "asymmetry_parameter": (
                    ("radius", "wavelength_nm"),
                    np.array([[0.4, 0.45], [0.6, 0.65]]),
                ),
            },
            coords={"radius": [1.0, 3.0], "wavelength_nm": [500.0, 600.0]},
        ),
        max_num_moments=16,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([550.0]),
        legendre_derivative=False,
    )
    atmosphere["aerosol"] = sk.constituent.NumberDensityScatterer2D(
        optical_property,
        np.full(geometry.shape, 2.0e6),
        radius=np.linspace(1.2, 2.0, geometry.shape[1]),
    )
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=20)
    linearization = engine.linearize(atmosphere)

    assert linearization.parameter_dims["aerosol_radius"] == ("altitude",)
    tangent = linearization.tangent_template[["aerosol_radius"]]
    tangent.aerosol_radius.values[:] = np.linspace(-0.2, 0.3, geometry.shape[1])
    cotangent = xr.ones_like(linearization.value)
    jvp = linearization.jvp(tangent)
    gradient = linearization.vjp(cotangent, parameters=["aerosol_radius"])
    np.testing.assert_allclose(
        (jvp * cotangent).sum(),
        (tangent.aerosol_radius * gradient.aerosol_radius).sum(),
        rtol=2e-11,
        atol=1e-14,
    )
    for diagnostics in engine.group_diagnostics:
        assert diagnostics["resident_volume_derivative_mappings"] == [
            "wf_aerosol_radius"
        ]
        assert diagnostics["resident_surface_derivative_mappings"] == []


@pytest.mark.parametrize("altitude_m", [20_000.0, 35_000.0])
def test_orbital_aerosol_extinction_jvp_matches_finite_difference_with_phase_mixing(
    altitude_m,
):
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.los_refraction = False
    config.num_singlescatter_moments = 16

    optical_wavelengths = np.array([499.0, 501.0])
    aerosol_optical = sk.optical.HenyeyGreenstein.from_parameters(
        wavelength_nm=optical_wavelengths,
        xs_total=np.full(2, 2.0e-12),
        ssa=np.full(2, 1.8e-12),
        g=np.full(2, 0.7),
        max_num_moments=16,
    )
    background_optical = sk.optical.HenyeyGreenstein.from_parameters(
        wavelength_nm=optical_wavelengths,
        xs_total=np.full(2, 1.0e-12),
        ssa=np.full(2, 0.7e-12),
        g=np.full(2, -0.2),
        max_num_moments=16,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        legendre_derivative=False,
    )
    aerosol = sk.constituent.ExtinctionScatterer2D(
        aerosol_optical,
        np.full(geometry.shape, 1.0e-7),
        500.0,
    )
    atmosphere["aerosol"] = aerosol
    atmosphere["background"] = sk.constituent.NumberDensityScatterer2D(
        background_optical,
        np.full(geometry.shape, 2.0e6),
    )
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=60,
        sun_vectors_ecef=np.array([[0.0, 0.0, 1.0]]),
    )

    horizontal_index = geometry.shape[0] // 2
    altitude_index = int(np.flatnonzero(altitude_m == ALTITUDES_M)[0])
    linearization = engine.linearize(atmosphere)
    tangent = linearization.tangent_template[["aerosol_extinction"]] * 0.0
    tangent.aerosol_extinction.values[horizontal_index, altitude_index] = 1.0
    analytic = linearization.jvp(tangent)
    cotangent = xr.ones_like(linearization.value)
    gradient = linearization.vjp(cotangent, parameters=("aerosol_extinction",))

    delta = 1.0e-10
    original = aerosol.extinction_per_m[horizontal_index, altitude_index]
    aerosol.extinction_per_m[horizontal_index, altitude_index] = original + delta
    above = engine.calculate_radiance(atmosphere).radiance
    aerosol.extinction_per_m[horizontal_index, altitude_index] = original - delta
    below = engine.calculate_radiance(atmosphere).radiance
    aerosol.extinction_per_m[horizontal_index, altitude_index] = original
    finite_difference = (above - below) / (2.0 * delta)

    assert np.any(np.abs(finite_difference) > 0.0)
    np.testing.assert_allclose(analytic, finite_difference, rtol=2.0e-6, atol=1.0e-10)
    np.testing.assert_allclose(
        gradient.aerosol_extinction.values[horizontal_index, altitude_index],
        (finite_difference * cotangent).sum().item(),
        rtol=2.0e-6,
        atol=1.0e-10,
    )


def test_derivative_execution_must_be_resident_or_streaming():
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    with pytest.raises(ValueError, match="derivative_execution"):
        sk.OrbitalPlaneEngine(
            transmission_config(),
            geometry,
            viewing,
            time_group_duration_s=60,
            derivative_execution="invalid",  # type: ignore[arg-type]
        )


@pytest.mark.parametrize(
    ("albedo", "expected_dims"),
    [
        (0.2, ()),
        (np.linspace(0.1, 0.4, 65), ("orbital_position",)),
        (
            np.linspace(0.1, 0.4, 130).reshape(65, 2),
            ("orbital_position", "surface_wavelength"),
        ),
    ],
)
def test_lambertian_surface_2d_parameter_layout(albedo, expected_dims):
    geometry = orbital_geometry()
    viewing = limb_viewing(geometry, np.array([0.0]))
    config = transmission_config()
    engine = sk.OrbitalPlaneEngine(config, geometry, viewing, time_group_duration_s=60)
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0, 600.0]),
        legendre_derivative=False,
    )
    atmosphere["optics"] = sk.constituent.Manual(
        np.full((*geometry.shape, 2), 1e-5), np.zeros((*geometry.shape, 2))
    )
    atmosphere["surface"] = sk.constituent.LambertianSurface2D(albedo)
    linearization = engine.linearize(atmosphere)
    assert linearization.parameter_dims["surface_albedo"] == expected_dims
    if expected_dims:
        assert linearization.tangent_template.surface_albedo.shape == np.shape(albedo)


def test_lambertian_surface_2d_translation_preserves_spatial_spectral_interpolation():
    geometry = orbital_geometry()
    lower_indices = np.array([24, 39])
    upper_indices = lower_indices + 1
    fractions = np.array([0.35, 0.65])
    lower_directions = geometry.ground_track_ecef_m[lower_indices]
    lower_directions /= np.linalg.norm(lower_directions, axis=1)[:, np.newaxis]
    upper_directions = geometry.ground_track_ecef_m[upper_indices]
    upper_directions /= np.linalg.norm(upper_directions, axis=1)[:, np.newaxis]
    separation = np.arccos(np.sum(lower_directions * upper_directions, axis=1))
    target_directions = (
        np.sin((1.0 - fractions) * separation)[:, np.newaxis]
        / np.sin(separation)[:, np.newaxis]
        * lower_directions
        + np.sin(fractions * separation)[:, np.newaxis]
        / np.sin(separation)[:, np.newaxis]
        * upper_directions
    )
    target_positions = EARTH_RADIUS_M * target_directions
    observer_angles = np.array([-0.23, 0.23])
    observer_directions = np.column_stack(
        (
            np.sin(observer_angles),
            np.zeros_like(observer_angles),
            np.cos(observer_angles),
        )
    )
    observers = (EARTH_RADIUS_M + 100_000.0) * observer_directions
    looks = target_positions - observers
    looks /= np.linalg.norm(looks, axis=1)[:, np.newaxis]
    viewing = sk.OrbitalPlaneViewingGeometry(
        np.array(
            ["2026-01-01T00:00:00", "2026-01-01T00:00:10"],
            dtype="datetime64[ns]",
        ),
        observers,
        looks,
        vertical_slice=np.arange(len(observers)),
    )

    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    atmosphere_wavelengths = np.array([450.0, 500.0, 550.0])
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=atmosphere_wavelengths,
        legendre_derivative=False,
    )
    atmosphere["optics"] = sk.constituent.Manual(
        np.full((*geometry.shape, len(atmosphere_wavelengths)), 1.0e-8),
        np.zeros((*geometry.shape, len(atmosphere_wavelengths))),
    )
    source_albedo = np.column_stack(
        (
            np.linspace(0.08, 0.38, geometry.shape[0]),
            np.linspace(0.58, 0.88, geometry.shape[0]),
        )
    )
    surface = sk.constituent.LambertianSurface2D(
        source_albedo,
        wavelengths_nm=np.array([400.0, 600.0]),
    )
    atmosphere["surface"] = surface
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=5,
        group_padding_angle=np.deg2rad(25.0),
        sun_vectors_ecef=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
    )

    assert len(engine.group_diagnostics) == 2
    assert set(engine.group_diagnostics[0]["grid_indices"]) & set(
        engine.group_diagnostics[1]["grid_indices"]
    )
    spatial_radiance = engine.calculate_radiance(atmosphere).radiance.values[:, :, 0]
    surface.albedo = np.ones_like(source_albedo)
    unit_albedo_radiance = engine.calculate_radiance(atmosphere).radiance.values[
        :, :, 0
    ]

    spectral_weights = np.array(
        [
            [0.75, 0.25],
            [0.50, 0.50],
            [0.25, 0.75],
        ]
    )
    native_albedo = source_albedo @ spectral_weights.T
    expected = (1.0 - fractions[:, np.newaxis]) * native_albedo[
        lower_indices
    ] + fractions[:, np.newaxis] * native_albedo[upper_indices]
    np.testing.assert_allclose(
        (spatial_radiance / unit_albedo_radiance).T,
        expected,
        rtol=2.0e-13,
        atol=2.0e-14,
    )

    surface.albedo = source_albedo
    linearization = engine.linearize(atmosphere)
    direction = np.zeros_like(source_albedo)
    direction[lower_indices] = np.array([[0.3, -0.1], [-0.2, 0.4]])
    direction[upper_indices] = np.array([[-0.15, 0.2], [0.1, -0.3]])
    tangent = linearization.tangent_template[["surface_albedo"]] * 0.0
    tangent.surface_albedo.values[:] = direction
    jvp = linearization.jvp(tangent)
    cotangent = xr.ones_like(linearization.value)
    gradient = linearization.vjp(cotangent, parameters=("surface_albedo",))
    np.testing.assert_allclose(
        (jvp * cotangent).sum(),
        (tangent.surface_albedo * gradient.surface_albedo).sum(),
        rtol=2.0e-12,
        atol=2.0e-14,
    )

    step = 1.0e-5
    surface.albedo = source_albedo + step * direction
    above = engine.calculate_radiance(atmosphere).radiance
    surface.albedo = source_albedo - step * direction
    below = engine.calculate_radiance(atmosphere).radiance
    surface.albedo = source_albedo
    finite_difference = (above - below) / (2.0 * step)
    np.testing.assert_allclose(jvp, finite_difference, rtol=2.0e-8, atol=2.0e-12)


@pytest.mark.parametrize(
    ("refraction", "multiple_scatter"),
    [
        (False, sk.MultipleScatterSource.NoSource),
        (True, sk.MultipleScatterSource.NoSource),
        (False, sk.MultipleScatterSource.SuccessiveOrders),
    ],
)
def test_lambertian_surface_2d_varies_at_ground_intersections_and_linearizes(
    refraction,
    multiple_scatter,
):
    geometry = orbital_geometry()
    lower_indices = np.array([24, 39])
    upper_indices = lower_indices + 1
    fractions = np.array([0.35, 0.65])
    lower_directions = geometry.ground_track_ecef_m[lower_indices]
    lower_directions /= np.linalg.norm(lower_directions, axis=1)[:, np.newaxis]
    upper_directions = geometry.ground_track_ecef_m[upper_indices]
    upper_directions /= np.linalg.norm(upper_directions, axis=1)[:, np.newaxis]
    separation = np.arccos(np.sum(lower_directions * upper_directions, axis=1))
    directions = (
        np.sin((1.0 - fractions) * separation)[:, np.newaxis]
        / np.sin(separation)[:, np.newaxis]
        * lower_directions
        + np.sin(fractions * separation)[:, np.newaxis]
        / np.sin(separation)[:, np.newaxis]
        * upper_directions
    )
    target_positions = EARTH_RADIUS_M * directions
    observer_angles = np.array([-0.23, 0.23])
    observer_directions = np.column_stack(
        (
            np.sin(observer_angles),
            np.zeros_like(observer_angles),
            np.cos(observer_angles),
        )
    )
    observers = (EARTH_RADIUS_M + 100_000.0) * observer_directions
    looks = target_positions - observers
    looks /= np.linalg.norm(looks, axis=1)[:, np.newaxis]
    viewing = sk.OrbitalPlaneViewingGeometry(
        np.array(
            ["2026-01-01T00:00:00", "2026-01-01T00:00:10"],
            dtype="datetime64[ns]",
        ),
        observers,
        looks,
        vertical_slice=np.arange(len(observers)),
    )
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = multiple_scatter
    config.los_refraction = refraction
    if multiple_scatter == sk.MultipleScatterSource.SuccessiveOrders:
        config.num_threads = 1
        config.num_sza = 2
        config.successive_orders_altitude_grid_m = np.array(
            [5_000.0, 25_000.0, 55_000.0]
        )
        config.num_successive_orders_incoming = 6
        config.num_successive_orders_outgoing = 6
        config.num_successive_orders_iterations = 60
        config.successive_orders_relative_tolerance = 1.0e-11
        config.successive_orders_absolute_tolerance = 1.0e-13
        config.successive_orders_anderson_depth = 3
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0, 600.0]),
        legendre_derivative=False,
    )
    extinction = np.full(
        (*geometry.shape, 2),
        (
            1.0e-5
            if multiple_scatter == sk.MultipleScatterSource.SuccessiveOrders
            else 1.0e-8
        ),
    )
    if multiple_scatter == sk.MultipleScatterSource.SuccessiveOrders:
        legendre_moments = np.zeros(
            (atmosphere.storage.leg_coeff.shape[0], *extinction.shape)
        )
        legendre_moments[0] = 1.0
        ssa = np.full_like(extinction, 0.9)
    else:
        legendre_moments = None
        ssa = np.zeros_like(extinction)
    atmosphere["optics"] = sk.constituent.Manual(extinction, ssa, legendre_moments)
    spatial_albedo = np.linspace(0.1, 0.7, geometry.shape[0])
    albedo = np.column_stack((spatial_albedo, 0.8 - 0.5 * spatial_albedo))
    surface = sk.constituent.LambertianSurface2D(albedo)
    atmosphere["surface"] = surface
    if refraction:
        refractive_profile = 1.0 + 3.0e-4 * np.exp(-ALTITUDES_M / 7_000.0)
        atmosphere.refractive_index = np.broadcast_to(
            refractive_profile, geometry.shape
        ).copy()
    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=60,
        sun_vectors_ecef=np.array([[0.0, 0.0, 1.0]]),
    )

    result = engine.calculate_radiance(atmosphere)
    measured = result.radiance.values[:, :, 0]
    expected_albedo = (1.0 - fractions[:, np.newaxis]) * albedo[
        lower_indices, :
    ] + fractions[:, np.newaxis] * albedo[upper_indices, :]
    expected_ratio = expected_albedo[1, :] / expected_albedo[0, :]
    if not refraction and multiple_scatter == sk.MultipleScatterSource.NoSource:
        np.testing.assert_allclose(
            measured[:, 1] / measured[:, 0], expected_ratio, rtol=2.0e-5
        )
    assert result.wf_surface_albedo.dims == (
        "orbital_position",
        "surface_wavelength",
        "wavelength",
        "los",
        "stokes",
    )

    linearization = engine.linearize(atmosphere, prepare_parameters=("surface_albedo",))
    prepared_update_counts = [
        item["atmosphere_update_count"] for item in engine.group_diagnostics
    ]
    tangent = xr.zeros_like(linearization.tangent_template[["surface_albedo"]])
    tangent["surface_albedo"].values[lower_indices, :] = [
        [0.3, -0.1],
        [-0.2, 0.4],
    ]
    tangent["surface_albedo"].values[upper_indices, :] = [
        [-0.15, 0.2],
        [0.1, -0.3],
    ]
    jvp = linearization.jvp(tangent)
    eager_jvp = (result.wf_surface_albedo * tangent.surface_albedo).sum(
        ("orbital_position", "surface_wavelength")
    )
    np.testing.assert_allclose(jvp, eager_jvp, rtol=2.0e-11, atol=1.0e-13)
    if multiple_scatter == sk.MultipleScatterSource.SuccessiveOrders:
        direct_nodes = {*lower_indices, *upper_indices}
        remote_sensitivity = np.abs(result.wf_surface_albedo.values[:, 0, 0, 0, 0])
        remote_sensitivity[list(direct_nodes)] = 0.0
        remote_index = int(np.argmax(remote_sensitivity))
        assert remote_sensitivity[remote_index] > 1.0e-8
        remote_tangent = xr.zeros_like(
            linearization.tangent_template[["surface_albedo"]]
        )
        remote_tangent["surface_albedo"].values[remote_index, 0] = 1.0
        remote_jvp = linearization.jvp(remote_tangent)
        remote_eager_jvp = (
            result.wf_surface_albedo * remote_tangent.surface_albedo
        ).sum(("orbital_position", "surface_wavelength"))
        np.testing.assert_allclose(
            remote_jvp, remote_eager_jvp, rtol=1.0e-9, atol=1.0e-12
        )
    cotangent = xr.zeros_like(linearization.value)
    cotangent.values[:] = np.array([[[0.7], [-0.4]], [[-0.2], [0.9]]])
    gradient = linearization.vjp(cotangent, parameters=("surface_albedo",))
    assert [
        item["atmosphere_update_count"] for item in engine.group_diagnostics
    ] == prepared_update_counts
    lhs = float((jvp * cotangent).sum())
    rhs = float((tangent.surface_albedo * gradient.surface_albedo).sum())
    assert lhs == pytest.approx(rhs, rel=2.0e-11, abs=1.0e-13)
    if multiple_scatter == sk.MultipleScatterSource.SuccessiveOrders:
        step = 1.0e-4
        direction = tangent.surface_albedo.values
        surface.albedo = albedo + step * direction
        above = engine.calculate_radiance(atmosphere).radiance
        surface.albedo = albedo - step * direction
        below = engine.calculate_radiance(atmosphere).radiance
        surface.albedo = albedo
        finite_difference = (above - below) / (2.0 * step)
        # Reusing the converged successive-orders primal changes the warm-start
        # sequence relative to the two perturbed forward solves. Their
        # absolute solver tolerance is 1e-13, which is amplified by the 1e-4
        # finite-difference step; retain a tight derivative check while
        # allowing that numerical floor.
        np.testing.assert_allclose(jvp, finite_difference, rtol=3.0e-6, atol=1.0e-9)


def test_solar_handler_is_evaluated_at_group_mean_sample_times():
    class RecordingSolarHandler:
        def __init__(self):
            self.times = []

        def target_solar_angles(self, _latitude, _longitude, _altitude, time):
            self.times.append(time.to_datetime64())
            return 45.0, 120.0

    geometry = orbital_geometry()
    times = np.array(
        [
            "2026-01-01T00:00:00",
            "2026-01-01T00:00:20",
            "2026-01-01T00:02:00",
        ],
        dtype="datetime64[ns]",
    )
    viewing = limb_viewing(geometry, np.array([-0.2, 0.0, 0.2]), times)
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    handler = RecordingSolarHandler()
    sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=60,
        solar_handler=handler,
    )
    np.testing.assert_array_equal(
        np.asarray(handler.times),
        np.array(
            [
                np.datetime64(
                    int(times[:2].astype(np.int64).sum() // times[:2].size), "ns"
                ),
                times[2],
            ],
            dtype="datetime64[ns]",
        ),
    )


def test_solar_handler_uses_geodetic_coordinates_and_local_basis():
    class RecordingOverheadSolarHandler:
        def __init__(self):
            self.coordinates = []

        def target_solar_angles(self, latitude, longitude, altitude, _time):
            self.coordinates.append((latitude, longitude, altitude))
            return 0.0, 0.0

    geoid = sk.WGS84()
    geoid.from_lat_lon_alt(45.0, 10.0, 700_000.0)
    observer = geoid.location.copy()
    look = geoid.from_tangent_altitude(20_000.0, observer, np.array([0.0, 0.0, 1.0]))
    viewing = sk.OrbitalPlaneViewingGeometry(
        np.array(["2026-01-01T00:00:00"], dtype="datetime64[ns]"),
        observer[np.newaxis, :],
        np.asarray(look)[np.newaxis, :],
        vertical_slice=np.array([0]),
        geoid=geoid,
    )
    geometry = viewing.construct_atmosphere_geometry(
        ALTITUDES_M, 0.02, path_padding_angle=0.1
    )
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    handler = RecordingOverheadSolarHandler()

    engine = sk.OrbitalPlaneEngine(
        config,
        geometry,
        viewing,
        time_group_duration_s=60,
        solar_handler=handler,
    )

    latitude, longitude, altitude = handler.coordinates[0]
    reference_geoid = sk.WGS84()
    reference_geoid.from_xyz(engine.group_diagnostics[0]["reference_position_ecef_m"])
    assert latitude == pytest.approx(reference_geoid.latitude, abs=1e-8)
    assert longitude == pytest.approx(reference_geoid.longitude, abs=1e-8)
    assert altitude == pytest.approx(0.0, abs=1e-12)
    reference_geoid.from_lat_lon_alt(latitude, longitude, 0.0)
    np.testing.assert_allclose(
        engine.sun_vectors_ecef[0], reference_geoid.local_up, atol=2e-14
    )
