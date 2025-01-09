#include "sasktran2/output.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <sasktran2.h>

namespace py = pybind11;

template <int NSTOKES>
void declareOutputBase(py::module_& m, const std::string& suffix) {
    using Output = sasktran2::Output<NSTOKES>;

    py::class_<Output>(m, ("Output" + suffix).c_str());
}

template <int NSTOKES>
void declareOutputIdeal(py::module_& m, const std::string& suffix) {
    using Output = sasktran2::OutputIdealDense<NSTOKES>;

    py::class_<Output, sasktran2::Output<NSTOKES>>(
        m, ("OutputIdeal" + suffix).c_str())
        .def(py::init<>())
        .def_property(
            "radiance",
            [](Output& output) -> Eigen::VectorXd& {
                return output.radiance().value;
            },
            nullptr)
        .def_property(
            "d_radiance",
            [](Output& output) -> Eigen::MatrixXd& {
                return output.radiance().deriv;
            },
            nullptr);
}

template <int NSTOKES>
void declareOutputDerivMapped(py::module_& m, const std::string& suffix) {
    using Output = sasktran2::OutputDerivMapped<NSTOKES>;

    py::class_<Output, sasktran2::Output<NSTOKES>>(
        m, ("OutputDerivMapped" + suffix).c_str())
        .def(py::init<>())
        .def_property(
            "radiance",
            [](Output& output) -> Eigen::VectorXd& {
                return output.radiance().value;
            },
            nullptr)
        .def_property(
            "deriv_map",
            [](Output& output)
                -> const std::map<std::string, Eigen::MatrixXd>& {
                return output.derivatives();
            },
            nullptr)
        .def_property(
            "surface_deriv_map",
            [](Output& output)
                -> const std::map<std::string, Eigen::MatrixXd>& {
                return output.surface_derivatives();
            },
            nullptr);
}

void init_output(py::module_& m) {
    declareOutputBase<1>(m, "Stokes_1");
    declareOutputBase<3>(m, "Stokes_3");

    declareOutputIdeal<1>(m, "Stokes_1");
    declareOutputIdeal<3>(m, "Stokes_3");

    declareOutputDerivMapped<1>(m, "Stokes_1");
    declareOutputDerivMapped<3>(m, "Stokes_3");
}
