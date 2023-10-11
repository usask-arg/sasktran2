#include <pybind11/pybind11.h>
#include <sasktran2.h>

namespace py = pybind11;



void init_config(py::module_ &  m) {
    py::enum_<sasktran2::Config::MultipleScatterSource>(m, "MultipleScatterSource")
            .value("DiscreteOrdinates", sasktran2::Config::MultipleScatterSource::discrete_ordinates)
            .value("SuccessiveOrders", sasktran2::Config::MultipleScatterSource::hr)
            .value("NoSource", sasktran2::Config::MultipleScatterSource::none)
            .export_values();


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
            .def_property("multiple_scatter_source", &sasktran2::Config::multiple_scatter_source, &sasktran2::Config::set_multiple_scatter_source,
            R"(

            )")
            .def_property("num_streams", &sasktran2::Config::num_do_streams, &sasktran2::Config::set_num_do_streams)
            ;
}

