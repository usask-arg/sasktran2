#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;


void init_geometry(py::module_ &  m) {
    py::class_<sasktran2::Geometry1D>(m, "Geometry1D")
            .def(py::init<double, double, double, Eigen::VectorXd&&, sasktran2::grids::interpolation, sasktran2::geometrytype>(),
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
                       One of sasktran2.InterpolationMethod.LinearInterpolation or sasktran2.InterpolationMethod.ShellInterpolation
                   geometry_type: sasktran2.GeometryType
                       One of sasktran2.GeometryType.Spherical or sasktran2.GeometryType.PlaneParallel
                 )",
                 "cos_sza"_a,
                 "solar_azimuth"_a,
                 "earth_radius_m"_a,
                 "altitude_grid_m"_a,
                 "interpolation_method"_a,
                 "geometry_type"_a
                 )
            .def("altitudes",
                 [](const sasktran2::Geometry1D& geo) {
                      return geo.altitude_grid().grid();
            })
            ;
}

