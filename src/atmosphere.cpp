#include <pybind11/pybind11.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/eigen/tensor.h>
#include <pybind11/stl.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

template <int NSTOKES>
void declareAtmosphere(py::module_& m, const std::string& suffix) {
    using Atmosphere = sasktran2::atmosphere::Atmosphere<NSTOKES>;

    py::class_<Atmosphere>(m, ("Atmosphere" + suffix).c_str())
        .def(py::init<int, const sasktran2::Geometry1D&,
                      const sasktran2::Config&, bool>())
        .def("apply_delta_m_scaling", &Atmosphere::apply_delta_m_scaling)
        .def_property(
            "storage",
            [](Atmosphere& atmo)
                -> sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>& {
                return atmo.storage();
            },
            [](Atmosphere& atmo,
               const sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>&
                   storage) { return atmo.storage() = storage; })
        .def_property(
            "surface",
            [](Atmosphere& atmo) -> sasktran2::atmosphere::Surface<NSTOKES>& {
                return atmo.surface();
            },
            nullptr);
}

template <int NSTOKES>
void declareAtmosphereStorage(py::module_& m, const std::string& suffix) {
    using AtmosphereGridStorage =
        sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>;

    py::class_<AtmosphereGridStorage>(m, ("AtmosphereStorage" + suffix).c_str())
        .def(py::init<int, int, int>())
        .def("resize_derivatives", &AtmosphereGridStorage::resize_derivatives,
             "num_deriv"_a)
        .def("get_derivative_mapping",
             &AtmosphereGridStorage::get_derivative_mapping, "name"_a,
             py::return_value_policy::reference)
        .def_property_readonly("derivative_mappings",
                               &AtmosphereGridStorage::derivative_mappings)
        .def_property(
            "ssa",
            [](AtmosphereGridStorage& storage) -> Eigen::MatrixXd& {
                return storage.ssa;
            },
            [](AtmosphereGridStorage& storage, const Eigen::MatrixXd& ssa) {
                storage.ssa = ssa;
            })
        .def_property(
            "total_extinction",
            [](AtmosphereGridStorage& storage) -> Eigen::MatrixXd& {
                return storage.total_extinction;
            },
            [](AtmosphereGridStorage& storage,
               const Eigen::MatrixXd& total_extinction) {
                storage.total_extinction = total_extinction;
            })
        .def_property(
            "emission_source",
            [](AtmosphereGridStorage& storage) -> Eigen::MatrixXd& {
                return storage.emission_source;
            },
            [](AtmosphereGridStorage& storage,
               const Eigen::MatrixXd& emission_source) {
                storage.emission_source = emission_source;
            })
        .def_property(
            "f",
            [](AtmosphereGridStorage& storage) -> const Eigen::MatrixXd& {
                return storage.f;
            },
            nullptr)
        .def_property(
            "d_f",
            [](AtmosphereGridStorage& storage)
                -> const Eigen::Tensor<double, 3>& { return storage.d_f; },
            nullptr)
        .def_property(
            "leg_coeff",
            [](AtmosphereGridStorage& storage) -> Eigen::Tensor<double, 3>& {
                return storage.leg_coeff;
            },
            [](AtmosphereGridStorage& storage,
               const Eigen::Tensor<double, 3>& leg_coeff) {
                storage.leg_coeff = leg_coeff;
            })
        .def_property(
            "d_leg_coeff",
            [](AtmosphereGridStorage& storage) -> Eigen::Tensor<double, 4>& {
                return storage.d_leg_coeff;
            },
            [](AtmosphereGridStorage& storage,
               const Eigen::Tensor<double, 4>& d_leg_coeff) {
                storage.d_leg_coeff = d_leg_coeff;
            })
        .def_property(
            "solar_irradiance",
            [](AtmosphereGridStorage& storage) -> Eigen::VectorXd& {
                return storage.solar_irradiance;
            },
            [](AtmosphereGridStorage& storage,
               const Eigen::VectorXd& irradiance) {
                storage.solar_irradiance = irradiance;
            });
}

template <int NSTOKES>
void declareSurface(py::module_& m, const std::string& suffix) {
    using Surface = sasktran2::atmosphere::Surface<NSTOKES>;

    py::class_<Surface>(m, ("Surface" + suffix).c_str())
        .def_property(
            "max_azimuthal_order",
            [](Surface& surface) -> int {
                return surface.max_azimuthal_order();
            },
            nullptr)
        .def("get_derivative_mapping", &Surface::get_derivative_mapping,
             "name"_a, py::return_value_policy::reference)
        .def_property_readonly("derivative_mappings",
                               &Surface::derivative_mappings)
        .def_property("brdf", &Surface::brdf_object, &Surface::set_brdf_object)
        .def_property(
            "brdf_args",
            [](Surface& surface) -> Eigen::MatrixXd& {
                return surface.brdf_args();
            },
            nullptr)
        .def_property(
            "d_brdf_args",
            [](Surface& surface) -> std::vector<Eigen::MatrixXd>& {
                return surface.d_brdf_args();
            },
            nullptr)
        .def_property(
            "albedo",
            [](Surface& surface) -> Eigen::MatrixXd& {
                return surface.brdf_args();
            },
            nullptr)
        .def_property(
            "emission",
            [](Surface& surface) -> Eigen::VectorXd& {
                return surface.emission();
            },
            nullptr);
}

void init_atmosphere(py::module_& m) {
    declareAtmosphere<1>(m, "Stokes_1");
    declareAtmosphere<3>(m, "Stokes_3");

    declareAtmosphereStorage<1>(m, "Stokes_1");
    declareAtmosphereStorage<3>(m, "Stokes_3");

    declareSurface<1>(m, "Stokes_1");
    declareSurface<3>(m, "Stokes_3");
}
