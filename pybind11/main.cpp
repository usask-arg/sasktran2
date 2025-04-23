#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_coordinates(py::module_&);
void init_geometry(py::module_&);
void init_geodetic(py::module_&);
void init_grids(py::module_&);
void init_config(py::module_&);
void init_surface(py::module_&);
void init_atmosphere(py::module_&);
void init_viewing_geometry(py::module_&);
void init_output(py::module_&);
void init_engine(py::module_&);
void init_mie(py::module_&);
void init_math(py::module_&);
void init_derivative_mappings(py::module_&);

PYBIND11_MODULE(_core, m) {
    init_config(m);
    init_grids(m);
    init_derivative_mappings(m);
    init_coordinates(m);
    init_geometry(m);
    init_geodetic(m);
    init_surface(m);
    init_atmosphere(m);
    init_viewing_geometry(m);
    init_output(m);
    init_engine(m);
    init_math(m);
    init_mie(m);
}
