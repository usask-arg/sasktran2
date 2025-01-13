#include "sasktran2/derivative_mapping.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/eigen/tensor.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

void init_derivative_mappings(py::module_& m) {
    py::class_<sasktran2::DerivativeMapping>(m, "DerivativeMapping")
        .def(py::init<int, int, int>())
        .def_property(
            "d_extinction",
            [](sasktran2::DerivativeMapping& mapping) -> Eigen::MatrixXd& {
                if (!mapping.native_mapping().d_extinction.has_value()) {
                    mapping.allocate_extinction_derivatives();
                }
                return mapping.native_mapping().d_extinction.value();
            },
            nullptr)
        .def_property(
            "scat_factor",
            [](sasktran2::DerivativeMapping& mapping) -> Eigen::MatrixXd& {
                if (!mapping.native_mapping().scat_factor.has_value()) {
                    mapping.allocate_legendre_derivatives();
                }
                return mapping.native_mapping().scat_factor.value();
            },
            nullptr)
        .def_property("scat_deriv_index",
                      &sasktran2::DerivativeMapping::get_scattering_index,
                      &sasktran2::DerivativeMapping::set_scattering_index)
        .def_property("interp_dim",
                      &sasktran2::DerivativeMapping::get_interp_dim,
                      &sasktran2::DerivativeMapping::set_interp_dim)
        .def_property("assign_name",
                      &sasktran2::DerivativeMapping::get_assign_name,
                      &sasktran2::DerivativeMapping::set_assign_name)
        .def_property("log_radiance_space",
                      &sasktran2::DerivativeMapping::log_radiance_space,
                      &sasktran2::DerivativeMapping::set_log_radiance_space)
        .def_property(
            "d_ssa",
            [](sasktran2::DerivativeMapping& mapping) -> Eigen::MatrixXd& {
                if (!mapping.native_mapping().d_ssa.has_value()) {
                    mapping.allocate_ssa_derivatives();
                }
                return mapping.native_mapping().d_ssa.value();
            },
            nullptr)
        .def_property(
            "d_emission",
            [](sasktran2::DerivativeMapping& mapping) -> Eigen::MatrixXd& {
                if (!mapping.native_mapping().d_emission.has_value()) {
                    mapping.allocate_emission_derivatives();
                }
                return mapping.native_mapping().d_emission.value();
            },
            nullptr)
        .def_property_readonly(
            "is_scattering_derivative",
            &sasktran2::DerivativeMapping::is_scattering_derivative)
        .def_property_readonly("num_output",
                               &sasktran2::DerivativeMapping::num_output)
        .def_property(
            "d_leg_coeff",
            [](sasktran2::DerivativeMapping& mapping)
                -> Eigen::Tensor<double, 3>& {
                if (!mapping.native_mapping().d_legendre.has_value()) {
                    mapping.allocate_legendre_derivatives();
                }
                return mapping.native_mapping().d_legendre.value();
            },
            nullptr)
        .def("set_zero", &sasktran2::DerivativeMapping::set_zero)
        .def_property("interpolator",
                      &sasktran2::DerivativeMapping::get_interpolator,
                      &sasktran2::DerivativeMapping::set_interpolator);

    py::class_<sasktran2::SurfaceDerivativeMapping>(m,
                                                    "SurfaceDerivativeMapping")
        .def(py::init<int, int>())
        .def_property(
            "d_brdf",
            [](sasktran2::SurfaceDerivativeMapping& mapping)
                -> Eigen::MatrixXd& {
                if (!mapping.native_surface_mapping().d_brdf.has_value()) {
                    mapping.allocate_brdf_derivatives();
                }
                return mapping.native_surface_mapping().d_brdf.value();
            },
            nullptr)
        .def_property(
            "d_emission",
            [](sasktran2::SurfaceDerivativeMapping& mapping)
                -> Eigen::MatrixXd& {
                if (!mapping.native_surface_mapping().d_emission.has_value()) {
                    mapping.allocate_emission_derivatives();
                }
                return mapping.native_surface_mapping().d_emission.value();
            },
            nullptr)
        .def_property("interpolator",
                      &sasktran2::SurfaceDerivativeMapping::get_interpolator,
                      &sasktran2::SurfaceDerivativeMapping::set_interpolator)
        .def_property("interp_dim",
                      &sasktran2::SurfaceDerivativeMapping::get_interp_dim,
                      &sasktran2::SurfaceDerivativeMapping::set_interp_dim)
        .def("set_zero", &sasktran2::SurfaceDerivativeMapping::set_zero);
}
