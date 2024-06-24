#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

void init_viewing_geometry(py::module_& m) {
    py::class_<sasktran2::viewinggeometry::ViewingGeometryBase>(
        m, "ViewingGeometryBase")
        .def("__repr__",
             &sasktran2::viewinggeometry::ViewingGeometryBase::to_string);

    py::class_<sasktran2::viewinggeometry::TangentAltitudeSolar,
               sasktran2::viewinggeometry::ViewingGeometryBase>(
        m, "TangentAltitudeSolar")
        .def(py::init<double, double, double, double>(),
             R"(
                    Defines a viewing ray from the observer altitude, and tangent point parameters. Note that all of
                    these parameters assume straight line paths (i.e. no atmospheric refraction)

                    Parameters
                    ----------
                    tangent_altitude_m: float
                        Tangent altitude in [m]
                    relative_azimuth: float
                        Relative azimuth angle to the sun. An angle of 0 degrees corresponds to the forward scattering plane. [rad]
                    observer_altitude_m: float
                        Observer altitude relative to the earth [m]
                    cos_sza: float
                        Cosine of the solar zenith angle at the tangent point [unitless]

                 )",
             "tangent_altitude_m"_a, "relative_azimuth"_a,
             "observer_altitude_m"_a, "cos_sza"_a);

    py::class_<sasktran2::viewinggeometry::GroundViewingSolar,
               sasktran2::viewinggeometry::ViewingGeometryBase>(
        m, "GroundViewingSolar")
        .def(py::init<double, double, double, double>(),
             R"(
                Defines a viewing ray that is looking at the ground from angles defined at the ground location. Note that
                all of these parameters assumes straight line paths (i.e. no atmospheric refraction)

                Parameters
                ----------
                cos_sza: float
                    Cosine of solar zenith angle at the ground point [unitless]
                relative_azimuth: float
                    Relative azimuth angle to the sun [rad] at the ground point. An angle of 0 degrees corresponds to the forward scattering plane.
                observer_altitude_m: float
                    Observer altitude relative to the earth [m]
                cos_viewing_zenith: float
                    Cosine of the viewing zenith angle at the ground point [unitless]
            )",
             "cos_sza"_a, "relative_azimuth"_a, "cos_viewing_zenith"_a,
             "observer_altitude_m"_a);

    py::class_<sasktran2::viewinggeometry::ViewingGeometryContainer>(
        m, "ViewingGeometry")
        .def(py::init<>())
        .def_property(
            "observer_rays",
            [](sasktran2::viewinggeometry::ViewingGeometryContainer&
                   container) {
                auto pylist = py::list();
                for (auto& ptr : container.observer_rays()) {
                    auto pyobj =
                        py::cast(*ptr, py::return_value_policy::reference);
                    pylist.append(pyobj);
                }
                return pylist;
            },
            nullptr)
        .def("add_ray",
             &sasktran2::viewinggeometry::ViewingGeometryContainer::add_ray);
}
