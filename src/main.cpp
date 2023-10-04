#include <pybind11/pybind11.h>

namespace py = pybind11;


void init_coordinates(py::module_ &);
void init_geometry(py::module_ &);
void init_grids(py::module_ &);
void init_config(py::module_ &);
void init_atmosphere(py::module_ &);
void init_viewing_geometry(py::module_ &);
void init_output(py::module_ &);
void init_engine(py::module_ &);


PYBIND11_MODULE(_core, m) {
    init_config(m);
    init_grids(m);
    init_coordinates(m);
    init_geometry(m);
    init_atmosphere(m);
    init_viewing_geometry(m);
    init_output(m);
    init_engine(m);
}

