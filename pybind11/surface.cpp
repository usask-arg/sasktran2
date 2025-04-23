#include "sasktran2/atmosphere/surface.h"
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;

template <int NSTOKES>
void declareBRDF(py::module_& m, const std::string& suffix) {
    using BRDF = sasktran2::atmosphere::brdf::BRDF<NSTOKES>;

    py::class_<BRDF, std::shared_ptr<BRDF>>(m, ("BRDF" + suffix).c_str());
}

template <int NSTOKES>
void declareLambertian(py::module_& m, const std::string& suffix) {
    using Lambertian = sasktran2::atmosphere::brdf::Lambertian<NSTOKES>;
    using BRDF = sasktran2::atmosphere::brdf::BRDF<NSTOKES>;

    py::class_<Lambertian, std::shared_ptr<Lambertian>, BRDF>(
        m, ("Lambertian" + suffix).c_str())
        .def(py::init<>());
}

template <typename Derived, int NSTOKES>
void declareDerived(py::module_& m, const std::string& name) {
    using BRDF = sasktran2::atmosphere::brdf::BRDF<NSTOKES>;

    py::class_<Derived, std::shared_ptr<Derived>, BRDF>(m, (name).c_str())
        .def(py::init<>())
        .def_property_readonly("num_deriv", &BRDF::num_deriv, R"(
                Number of derivatives this BRDF will calculate.
            )");
    ;
}

void init_surface(py::module_& m) {
    declareBRDF<1>(m, "Stokes_1");
    declareBRDF<3>(m, "Stokes_3");

    // declareLambertian<1>(m, "Stokes_1");
    // declareLambertian<3>(m, "Stokes_3");

    declareDerived<sasktran2::atmosphere::brdf::Lambertian<1>, 1>(
        m, "LambertianStokes_1");
    declareDerived<sasktran2::atmosphere::brdf::Lambertian<3>, 3>(
        m, "LambertianStokes_3");

    declareDerived<sasktran2::atmosphere::brdf::SnowKokhanovsky<1>, 1>(
        m, "SnowKokhanovskyStokes_1");
    declareDerived<sasktran2::atmosphere::brdf::SnowKokhanovsky<3>, 3>(
        m, "SnowKokhanovskyStokes_3");

    declareDerived<sasktran2::atmosphere::brdf::MODIS<1>, 1>(m,
                                                             "MODISStokes_1");
    declareDerived<sasktran2::atmosphere::brdf::MODIS<3>, 3>(m,
                                                             "MODISStokes_3");
}
