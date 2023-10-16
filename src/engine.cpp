#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

template <int NSTOKES>
void declareEngine(py::module_& m, const std::string& suffix) {
    using Engine = Sasktran2<NSTOKES>;

    py::class_<Engine>(m, ("Engine" + suffix).c_str())
        .def(py::init<
                 const sasktran2::Config&, const sasktran2::Geometry1D*,
                 const sasktran2::viewinggeometry::ViewingGeometryContainer&>(),
             R"(
                 Internal SASKTRAN2 object which handles the radiative transfer calculation.

                 Parameters
                 ----------
                 config: sasktran2.Config
                     Configuration settings
                 model_geometry: sasktran2.Geometry1D
                     The model geometry
                 viewing_geometry: sasktan2.ViewingGeometry
                     The viewing geometry
                 )",
             "config"_a, "model_geometry"_a, "viewing_geometry"_a)
        .def("calculate_radiance", &Engine::calculate_radiance,
             R"(
                    Performs the radiative transfer calculation for the given atmosphere, placing the result in output

                    Parameters
                    ----------
                    atmosphere: sasktran2.Atmosphere
                        The atmosphere object

                    output: sasktran2.Output
                        The result to place the output inside
                 )",
             "atmosphere"_a, "output"_a);
}

void init_engine(py::module_& m) {
    declareEngine<1>(m, "Stokes_1");
    declareEngine<3>(m, "Stokes_3");
}
