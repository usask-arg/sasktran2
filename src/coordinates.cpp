#include "sasktran2/geometry.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;

void init_coordinates(py::module_& m) {
    py::enum_<sasktran2::geometrytype>(m, "GeometryType")
        .value("PlaneParallel", sasktran2::geometrytype::planeparallel)
        .value("Spherical", sasktran2::geometrytype::spherical)
        .value("PseudoSpherical", sasktran2::geometrytype::pseudospherical)
        .value("Ellipsoidal", sasktran2::geometrytype::ellipsoidal)
        .export_values();

    py::class_<sasktran2::Coordinates>(m, "Coordinates")
        .def(py::init<double, double, double, sasktran2::geometrytype, bool>());
}
