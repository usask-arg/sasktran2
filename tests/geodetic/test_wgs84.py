from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _wgs84() -> sk.Geodetic:
    return sk.WGS84()


def test_geodetic_construction():
    _ = _wgs84()


def test_geodetic_from_lat_lon_alt():
    geo = _wgs84()

    # Test case from sasktran1 geodetic
    geo.from_lat_lon_alt(latitude=60, longitude=30, altitude=0)

    np.testing.assert_almost_equal(
        geo.location, np.array([2768773.79083189, 1598552.29346197, 5500477.13393864])
    )
    np.testing.assert_almost_equal(geo.local_south, np.array([0.75, 0.4330127, -0.5]))
    np.testing.assert_almost_equal(geo.local_up, np.array([0.4330127, 0.25, 0.8660254]))
    np.testing.assert_almost_equal(
        geo.local_west, np.array([5.00000000e-01, -8.66025404e-01, -2.77555756e-17])
    )

    # Check internal consistency
    lats = [-80, 0, 80]
    lons = [0, 180, 330]
    alts = [0, 10000, 80000]

    for lat, lon, alt in zip(lats, lons, alts, strict=True):
        geo.from_lat_lon_alt(latitude=lat, longitude=lon, altitude=alt)

        np.testing.assert_almost_equal(geo.latitude, lat)
        np.testing.assert_almost_equal(geo.longitude, lon)
        np.testing.assert_almost_equal(geo.altitude, alt)


def test_geodetic_from_xyz():
    geo = _wgs84()

    geo.from_xyz([5797230.47518212, -2110019.3341472, -1642001.16317228])

    np.testing.assert_almost_equal(geo.latitude, -14.999999973747736)
    np.testing.assert_almost_equal(geo.longitude, 340.00000000000006)
    np.testing.assert_almost_equal(geo.altitude, 7344.999610390202)


def test_geodetic_from_tangent_altitude():
    geo = _wgs84()

    look = geo.from_tangent_altitude(
        15322,
        [3.676013154788849600e005, 1.009976313640051500e006, -6.871601202127538600e006],
        [0, 0, 1],
    )

    np.testing.assert_almost_equal(look, np.array([0.28880556, 0.79348676, 0.53569591]))
    np.testing.assert_almost_equal(
        geo.location, np.array([1174200.8589215, 3226090.34579105, -5375466.33023017])
    )


def test_geodetic_from_tangent_point():
    geo = _wgs84()

    geo.from_tangent_point(
        [3.676013154788849600e005, 1.009976313640051500e006, -6.871601202127538600e006],
        [
            2.884568631765662100e-001,
            7.925287180643269000e-001,
            5.372996083468238900e-001,
        ],
    )

    np.testing.assert_almost_equal(
        geo.location, np.array([1176731.78951674, 3233044.02045169, -5364459.08730841])
    )
