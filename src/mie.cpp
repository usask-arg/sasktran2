#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

void init_mie(py::module_& m) {
    // MieData
    py::class_<sasktran2::mie::MieData>(m, "MieData")
        .def_readonly("Qext", &sasktran2::mie::MieData::Qext)
        .def_readonly("Qsca", &sasktran2::mie::MieData::Qsca)
        .def_readonly("S1", &sasktran2::mie::MieData::S1)
        .def_readonly("S2", &sasktran2::mie::MieData::S2);

    // MieOutput
    py::class_<sasktran2::mie::MieOutput>(m, "MieOutput")
        .def_readonly("size_parameter", &sasktran2::mie::MieOutput::size_param)
        .def_readonly("cos_angles", &sasktran2::mie::MieOutput::cos_angles)
        .def_readonly("refractive_index",
                      &sasktran2::mie::MieOutput::refractive_index)
        .def_readonly("values", &sasktran2::mie::MieOutput::values);

    py::class_<sasktran2::mie::LinearizedMie>(m, "LinearizedMie")
        .def(py::init<>())
        .def("calculate", &sasktran2::mie::LinearizedMie::calculate,
             "size_param"_a, "refractive_index"_a, "cos_angles"_a,
             "calculate_derivative"_a);
}
