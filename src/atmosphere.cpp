#include <pybind11/pybind11.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <pybind11/eigen.h>
#include <pybind11/eigen/tensor.h>
#include <pybind11/stl.h>
#include <sasktran2.h>


namespace py = pybind11;

template<int NSTOKES>
void declareAtmosphere(py::module_ & m, const std::string & suffix) {
    using Atmosphere = sasktran2::atmosphere::Atmosphere<NSTOKES>;

    py::class_<Atmosphere>(m, ("Atmosphere" + suffix).c_str())
            .def(py::init<int, const sasktran2::Geometry1D&, const sasktran2::Config&, bool>())
            .def_property("storage",
                          [](Atmosphere& atmo) -> sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>& { return atmo.storage(); },
                          [](Atmosphere& atmo, const sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>& storage) { return atmo.storage() = storage; }
            )
            .def_property("surface",
                          [](Atmosphere& atmo) -> sasktran2::atmosphere::Surface& { return atmo.surface(); },
                          nullptr
            )
            ;

}

template<int NSTOKES>
void declareAtmosphereStorage(py::module_ & m, const std::string & suffix) {
    using AtmosphereGridStorage = sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>;

    py::class_<AtmosphereGridStorage>(m, ("AtmosphereStorage" + suffix).c_str())
            .def(py::init<int, int, int>())
            .def_property("ssa",
                          [](AtmosphereGridStorage& storage) -> Eigen::MatrixXd& { return storage.ssa; },
                          [](AtmosphereGridStorage& storage, const Eigen::MatrixXd& ssa) { storage.ssa = ssa; }
            )
            .def_property("total_extinction",
                          [](AtmosphereGridStorage& storage) -> Eigen::MatrixXd& { return storage.total_extinction; },
                          [](AtmosphereGridStorage& storage, const Eigen::MatrixXd& total_extinction) { storage.total_extinction = total_extinction; }
            )
            .def_property("f",
                          [](AtmosphereGridStorage& storage) -> const Eigen::MatrixXd& { return storage.f; },
                          nullptr
            )
            .def_property("leg_coeff",
                          [](AtmosphereGridStorage& storage) -> Eigen::Tensor<double, 3>& { return storage.leg_coeff; },
                          [](AtmosphereGridStorage& storage, const Eigen::Tensor<double, 3>& leg_coeff) { storage.leg_coeff = leg_coeff; }
            )
            ;

}



void init_atmosphere(py::module_ &  m) {
    declareAtmosphere<1>(m, "Stokes_1");
    declareAtmosphere<3>(m, "Stokes_3");

    declareAtmosphereStorage<1>(m, "Stokes_1");
    declareAtmosphereStorage<3>(m, "Stokes_3");

    py::class_<sasktran2::atmosphere::Surface>(m, "Surface")
        .def(py::init<>())
        .def_property("albedo",
                      [](sasktran2::atmosphere::Surface& surface) -> Eigen::VectorXd& { return surface.albedo(); },
                      [](sasktran2::atmosphere::Surface& surface, const Eigen::VectorXd& albedo) { surface.albedo() = albedo; }
        );
}

