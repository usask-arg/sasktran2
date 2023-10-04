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
            .def_readwrite("phase", &AtmosphereGridStorage::phase);

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


    // For Phase Storage we declare manually because leg_coeff has a different number of dimensions
    py::class_<sasktran2::atmosphere::PhaseStorage<1>>(m, "PhaseStorageStokes_1")
        .def(py::init<>())
        .def_property("leg_coeff",
                      [](sasktran2::atmosphere::PhaseStorage<1>& storage) -> Eigen::MatrixXd& { return storage.storage(); },
                      [](sasktran2::atmosphere::PhaseStorage<1>& storage, const Eigen::MatrixXd& lc) { storage.storage() = lc; }
        )
        ;

    py::class_<sasktran2::atmosphere::PhaseStorage<3>>(m, "PhaseStorageStokes_3")
            .def(py::init<>())
            .def_property("leg_coeff",
                          [](sasktran2::atmosphere::PhaseStorage<3>& storage) -> Eigen::TensorMap<Eigen::Tensor<double, 3>> { return Eigen::TensorMap<Eigen::Tensor<double, 3>>(storage.storage().data(), storage.storage().rows() / 4, 4, storage.storage().cols()); },
                          nullptr
            )
            ;
}

