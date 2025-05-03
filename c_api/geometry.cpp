#include "geometry.h"
#include <sasktran2.h>
#include "internal_types.h"

Geometry1D::Geometry1D(double cos_sza, double saa, double earth_radius,
                       double* grid_values, int ngrid_values, int interp_method,
                       int geotype) {
    Eigen::VectorXd grid_values_vec =
        Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(grid_values, ngrid_values));

    try {
        impl = std::make_unique<sasktran2::Geometry1D>(
            cos_sza, saa, earth_radius, std::move(grid_values_vec),
            static_cast<sasktran2::grids::interpolation>(interp_method),
            static_cast<sasktran2::geometrytype>(geotype));
    } catch (const std::exception& e) {
        impl = nullptr;
    }
}

extern "C" {
Geometry1D* sk_geometry1d_create(double cos_sza, double saa,
                                 double earth_radius, double* grid_values,
                                 int ngrid_values, int interp_method,
                                 int geotype) {
    Geometry1D* geometry = nullptr;
    geometry = new Geometry1D(cos_sza, saa, earth_radius, grid_values,
                              ngrid_values, interp_method, geotype);
    if (geometry->impl == nullptr) {
        delete geometry;
        return nullptr; // Error: failed to create Geometry1D
    }
    return geometry;
}

void sk_geometry1d_destroy(Geometry1D* geometry) { delete geometry; }

int sk_geometry1d_get_num_altitudes(const Geometry1D* geometry) {
    return static_cast<int>(geometry->impl->altitude_grid().grid().size());
}

int sk_geometry1d_get_altitudes(const Geometry1D* geometry, double* altitudes) {
    if (geometry == nullptr || altitudes == nullptr) {
        return -1; // Error: null pointer
    }
    const auto& grid = geometry->impl->altitude_grid().grid();
    for (size_t i = 0; i < grid.size(); ++i) {
        altitudes[i] = grid[i];
    }
    return 0; // Success
}

int sk_geometry1d_get_refractive_index_ptr(const Geometry1D* geometry,
                                           double** refractive_index) {
    if (geometry == nullptr || refractive_index == nullptr) {
        return -1; // Error: null pointer
    }
    auto& n = geometry->impl->refractive_index();
    *refractive_index = n.data();
    return 0; // Success
}
}
