#include "geometry.h"
#include <sasktran2.h>
#include "internal_types.h"
#include <algorithm>

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

Geometry2D::Geometry2D(double cos_sza, double saa, double earth_radius,
                       const double* altitude_grid_values, int num_altitudes,
                       const double* horizontal_angle_grid_values,
                       int num_horizontal_locations,
                       int altitude_interp_method) {
    try {
        const Eigen::VectorXd altitude_grid = Eigen::Map<const Eigen::VectorXd>(
            altitude_grid_values, num_altitudes);
        const Eigen::VectorXd horizontal_angle_grid =
            Eigen::Map<const Eigen::VectorXd>(horizontal_angle_grid_values,
                                              num_horizontal_locations);
        impl = std::make_unique<sasktran2::Geometry2D>(
            cos_sza, saa, earth_radius, Eigen::VectorXd(altitude_grid),
            Eigen::VectorXd(horizontal_angle_grid),
            static_cast<sasktran2::grids::interpolation>(
                altitude_interp_method));
    } catch (const std::exception&) {
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

Geometry2D*
sk_geometry2d_create(double cos_sza, double saa, double earth_radius,
                     const double* altitude_grid_values, int num_altitudes,
                     const double* horizontal_angle_grid_values,
                     int num_horizontal_locations, int altitude_interp_method) {
    if (altitude_grid_values == nullptr ||
        horizontal_angle_grid_values == nullptr || num_altitudes < 0 ||
        num_horizontal_locations < 0) {
        return nullptr;
    }
    auto* geometry =
        new Geometry2D(cos_sza, saa, earth_radius, altitude_grid_values,
                       num_altitudes, horizontal_angle_grid_values,
                       num_horizontal_locations, altitude_interp_method);
    if (geometry->impl == nullptr) {
        delete geometry;
        return nullptr;
    }
    return geometry;
}

void sk_geometry2d_destroy(Geometry2D* geometry) { delete geometry; }

int sk_geometry2d_get_location_shape(const Geometry2D* geometry,
                                     int* num_horizontal_locations,
                                     int* num_altitudes) {
    if (geometry == nullptr || num_horizontal_locations == nullptr ||
        num_altitudes == nullptr) {
        return -1;
    }
    const auto shape = geometry->impl->location_shape();
    *num_horizontal_locations = shape.first;
    *num_altitudes = shape.second;
    return 0;
}

int sk_geometry2d_get_altitudes(const Geometry2D* geometry, double* altitudes) {
    if (geometry == nullptr || altitudes == nullptr) {
        return -1;
    }
    const auto& grid = geometry->impl->altitude_grid().grid();
    std::copy(grid.begin(), grid.end(), altitudes);
    return 0;
}

int sk_geometry2d_get_horizontal_angles(const Geometry2D* geometry,
                                        double* horizontal_angles) {
    if (geometry == nullptr || horizontal_angles == nullptr) {
        return -1;
    }
    const auto& grid = geometry->impl->horizontal_angle_grid();
    std::copy(grid.begin(), grid.end(), horizontal_angles);
    return 0;
}

int sk_geometry2d_get_location_index(const Geometry2D* geometry,
                                     int altitude_index, int horizontal_index,
                                     int* location_index) {
    if (geometry == nullptr || location_index == nullptr) {
        return -1;
    }
    try {
        *location_index =
            geometry->impl->location_index(altitude_index, horizontal_index);
    } catch (const std::exception&) {
        return -1;
    }
    return 0;
}
}
