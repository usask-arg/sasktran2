#include <pybind11/pybind11.h>
#include <sasktran2.h>

namespace py = pybind11;



void init_config(py::module_ &  m) {
    py::class_<sasktran2::Config>(m, "Config")
            .def(py::init<>())
            .def_property("num_threads", &sasktran2::Config::num_threads, &sasktran2::Config::set_num_threads);
}

