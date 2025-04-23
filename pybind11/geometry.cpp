#include "sasktran2/geometry.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

void init_geometry(py::module_& m) {
    py::class_<sasktran2::Geometry1D>(m, "Geometry1D")
        .def(py::init<double, double, double, Eigen::VectorXd&&,
                      sasktran2::grids::interpolation,
                      sasktran2::geometrytype>(),
             R"(
                   Initializes a geometry where the atmosphere varies only in 1 dimension (altitude).  The reference point
                   is defined by solar angles at the reference point.

                   Parameters
                   ----------
                   cos_sza: float
                       Cosine of solar zenith angle at the reference point
                   solar_azimuth: float
                       Solar azimuth angle at the reference point.
                   earth_radius_m: float
                       Radius of the earth.  Only has an effect if geometry_type is not set to PlaneParallel
                   altitude_grid_m: np.array
                       One dimensional altitude grid
                   interpolation_method: sasktran2.InterpolationMethod
                       The interpolation method to use in-between geometry grid points

                       `sasktran2.InterpolationMethod.LinearInterpolation`
                           In-between grid points, linear interpolation is assumed.  This means that Atmospheric quantities
                           such as extinction, single scatter albedo, should be thought of as sampled on the geometry grid points.

                       `sasktran2.InterpolationMethod.ShellInterpolation`
                           Atmospheric quantities such as extinction, single scatter albedo, are assumed to be constant in-between
                           geometry grid points.

                   geometry_type: sasktran2.GeometryType
                       The global geometry type used inside the radiative transfer calculation.

                       `sasktran2.GeometryType.Spherical`
                           All aspects of the calculation are done using spherical geometry.

                       `sasktran2.GeometryType.PlaneParallel`
                           All aspects of the calculation are done using plane-parallel geometry.

                       `sasktran2.GeometryType.PseudoSpherical`
                           Line of sight integration and the multiple scatter calculation is done using
                           plane parallel geometry, however the initial solar source function is calculated
                           using a spherical geometry.
                 )",
             "cos_sza"_a, "solar_azimuth"_a, "earth_radius_m"_a,
             "altitude_grid_m"_a, "interpolation_method"_a, "geometry_type"_a)
        .def("altitudes",
             [](const sasktran2::Geometry1D& geo) {
                 return geo.altitude_grid().grid();
             })
        .def_property(
            "refractive_index",
            [](sasktran2::Geometry1D& storage) -> Eigen::VectorXd& {
                return storage.refractive_index();
            },
            [](sasktran2::Geometry1D& storage, const Eigen::VectorXd& ssa) {
                storage.refractive_index() = ssa;
            },
            R"(
                The refractive index of the atmosphere.  This is used to calculate refraction in the radiative transfer calculation.
                Defaults to 1.0 which indicates no refractive effects.  Only has an effect if the refraction configuration options are
                enabled in the `sasktran2.Config` object.
            )");
}
