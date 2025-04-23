#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;

void init_grids(py::module_& m) {
    py::enum_<sasktran2::grids::interpolation>(m, "InterpolationMethod")
        .value("ShellInterpolation", sasktran2::grids::interpolation::shell)
        .value("LinearInterpolation", sasktran2::grids::interpolation::linear)
        .value("LowerInterpolation", sasktran2::grids::interpolation::lower)
        .export_values();

    py::enum_<sasktran2::grids::gridspacing>(m, "GridSpacing")
        .value("ConstantSpacing", sasktran2::grids::gridspacing::constant)
        .value("LinearSpacing", sasktran2::grids::gridspacing::variable)
        .export_values();

    py::enum_<sasktran2::grids::outofbounds>(m, "OutOfBoundsPolicy")
        .value("OutOfBoundsExtend", sasktran2::grids::outofbounds::extend)
        .value("OutOfBoundsSetZero", sasktran2::grids::outofbounds::setzero)
        .export_values();

    py::class_<sasktran2::grids::AltitudeGrid>(m, "AltitudeGrid")
        .def(py::init<Eigen::VectorXd&&, sasktran2::grids::gridspacing,
                      sasktran2::grids::outofbounds,
                      sasktran2::grids::interpolation>());
}
