#include <pybind11/pybind11.h>
#include <sasktran2.h>

namespace py = pybind11;



void init_config(py::module_ &  m) {
    py::class_<sasktran2::Config>(m, "Config")
            .def(py::init<>(),
            R"(
                Object which stores all of the configuration settings for the radiative transfer calculation.
            )"
            )
            .def_property("num_threads", &sasktran2::Config::num_threads, &sasktran2::Config::set_num_threads,
            R"(
                Controls the number of threads used in the calculation.  For maximum performance it is
                recommended to set this to the number of physical cores on your machine.  Defaults to
                1
            )"
            )
            .def_property("num_stokes", &sasktran2::Config::num_stokes, &sasktran2::Config::set_num_stokes,
            R"(
                Sets the number of Stokes parameters used in the calculation. 1 is equivalent to the scalar approximation.
                Currently the only supported values are 1, and 3.
            )")
            ;
}

