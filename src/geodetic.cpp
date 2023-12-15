#include "sasktran2/math/geodetic.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

void init_geodetic(py::module_& m) {
    py::class_<sasktran2::math::geodetic::Geodetic>(m, "Geodetic")
        .def(py::init<double, double>(), R"(
            A geodetic object based on a given equatorial (semi-major) radius and flattening factor.

            Standard usage is to create a geodetic object, and then initialize it through one of the
            `from_*` methods.

            Parameters
            ----------
            equatorial_radius: float
                Radius at the equator (semi-major) of the ellipsoid.
            flattening_factor: float
                Flattening factor of the ellipsoid. This is defined as (a-b)/a, where a is the semi-major axis
                and b is the semi-minor radius

        )",
             "equatorial_radius"_a, "flattening_factor"_a)
        .def_property_readonly(
            "altitude", &sasktran2::math::geodetic::Geodetic::altitude, R"(
                Altitude in [m] above the surface of the ellipsoid.
            )")
        .def_property_readonly(
            "latitude", &sasktran2::math::geodetic::Geodetic::latitude, R"(
                Geodetic latitude in degrees
            )")
        .def_property_readonly(
            "longitude", &sasktran2::math::geodetic::Geodetic::longitude, R"(
                Geodetic longitude in degrees
            )")
        .def_property_readonly(
            "location", &sasktran2::math::geodetic::Geodetic::location, R"(
                Geocentric location in cartesian coordinates
            )")
        .def_property_readonly(
            "local_south", &sasktran2::math::geodetic::Geodetic::local_south,
            R"(
                A unit vector pointing in the local south direction
            )")
        .def_property_readonly(
            "local_up", &sasktran2::math::geodetic::Geodetic::local_up, R"(
                A unit vector pointing up (perpindicular to the ellipsoidal surface)
            )")
        .def_property_readonly("local_west",
                               &sasktran2::math::geodetic::Geodetic::local_west,
                               R"(
                                A unit vector pointing in the local west direction
                               )")
        .def("altitude_intercepts",
             &sasktran2::math::geodetic::Geodetic::altitude_intercepts, R"(
                Calculate the two intersections of a line of sight and an altitude.

                Parameters
                ----------
                altitude : float
                    Altitude in meters.
                observer : np.ndarray
                    Three element array containing the obsever position in geocentric coordinates.
                look_vector : np.ndarray
                    Three element array containing a normalized look vector.

                Returns
                -------
                np.ndarray
                    Three element array containing the first (entering) intercept in geocentric coordinates.
                np.ndarray
                    Three element array containing the second (exiting) intercept in geocentric coordinates.

                Examples
                --------
                >>> import sasktran2 as sk
                >>> import numpy as np
                >>> geodetic = sk.WGS84()
                >>> look = geodetic.from_tangent_altitude(15322, [3.676013154788849600e+005, 1.009976313640051500e+006, \
                                                            -6.871601202127538600e+006], [0, 0, 1])
                >>> obs = geodetic.location
                >>> intercept1, intercept2 = geodetic.altitude_intercepts(16000, obs, look)
                >>> print(np.array_str(intercept1, precision=3))
                [ 1147302.059  3152186.5   -5425360.027]
                >>> print(np.array_str(intercept2, precision=3))
                [ 1201098.489  3299990.978 -5325574.803]
             )",
             "altitude"_a, "observer"_a, "look_vector"_a)
        .def("from_lat_lon_alt",
             &sasktran2::math::geodetic::Geodetic::from_lat_lon_alt, R"(
                Initializes the Geodetic based on a specifiec latitude, longitude, and altitude.

                Parameters
                ----------
                latitude : float
                    Latitude in degrees (-90 to 90)
                longitude : float
                    Longitude in degrees (0 to 360 or -180 to 180)
                altitude : float
                    Altitude above the geoid in metres

                Examples
                --------
                >>> import sasktran2 as sk
                >>> geodetic = sk.WGS84()
                >>> geodetic.from_lat_lon_alt(latitude=-15, longitude=-20, altitude=7342)
                >>> print(geodetic)
                WGS84 Location:
                Latitude: -15.0, Longitude: 340.0, Altitude: 7342.0
             )",
             "latitude"_a, "longitude"_a, "altitude"_a)
        .def("from_tangent_altitude",
             &sasktran2::math::geodetic::Geodetic::from_tangent_altitude,
             R"(
                Initialized the Geodetic from a specified tangent altitude, obsever location, and bore sight plane.

                Parameters
                ----------
                altitude : float
                    Tangent altitude in meters
                observer : np.ndarray
                    Three element array containing the obsever position in geocentric coordinates
                boresight : np.ndarray
                    Three element array containing a normalized look vector that is within the bore sight plane.

                Returns
                -------
                np.ndarray
                    Three element array containing the normalized look vector to the tangent point.

                Examples
                --------
                >>> import sasktran2 as sk
                >>> geodetic = sk.WGS84()
                >>> look = geodetic.from_tangent_altitude(15322, [ 3.676013154788849600e+005, 1.009976313640051500e+006,\
                                                                  -6.871601202127538600e+006], [0, 0, 1])
                >>> print(look)
                [0.28880556 0.79348676 0.53569591]
                >>> print(geodetic)
                WGS84 Location:
                Latitude: -57.60888188776806, Longitude: 70.00000000000001, Altitude: 15321.971935882739
             )",
             "altitude"_a, "observer"_a, "boresight"_a)
        .def("from_tangent_point",
             &sasktran2::math::geodetic::Geodetic::from_tangent_point,
             R"(
                Initializes  the Geodetic by calculating the tangent point from an observer position and look vector

                Parameters
                ----------
                observer : np.ndarray
                    Three element array containing the observer position in geocentric coordinates
                look_vector : np.ndarray
                    Three element array containing a normalized look vector

                Examples
                --------
                >>> import sasktran2 as sk
                >>> geodetic = sk.WGS84()
                >>> geodetic.from_tangent_point([ 3.676013154788849600e+005, 1.009976313640051500e+006,\
                                                 -6.871601202127538600e+006], [ 2.884568631765662100e-001,\
                                                  7.925287180643269000e-001,  5.372996083468238900e-001])
                >>> print(geodetic)
                WGS84 Location:
                Latitude: -57.500000192733594, Longitude: 70.0, Altitude: 10002.99586173162
             )",
             "observer"_a, "look_vector"_a)
        .def("from_xyz", &sasktran2::math::geodetic::Geodetic::from_xyz,
             R"(
                Initializes the Geodetic from a geocentric location

                Parameters
                ----------
                location : np.ndarray
                    Three element vector containing a location in geocentric coordinates

                Examples
                --------
                >>> import sasktran2 as sk
                >>> geodetic = sk.WGS84()
                >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
                >>> print(geodetic)
                WGS84 Location:
                Latitude: -14.999999973747736, Longitude: 340.00000000000006, Altitude: 7344.999610390202
             )",
             "location"_a)
        .def_property_readonly(
            "valid", &sasktran2::math::geodetic::Geodetic::is_valid, R"(
                True if the geodetic object has been initialized, False otherwise.
            )")

        ;
}
