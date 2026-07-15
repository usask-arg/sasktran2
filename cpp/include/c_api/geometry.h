#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Geometry1D Geometry1D;
typedef struct Geometry2D Geometry2D;

Geometry1D* sk_geometry1d_create(double cos_sza, double saa,
                                 double earth_radius, double* grid_values,
                                 int ngrid_values, int interp_method,
                                 int geotype);

void sk_geometry1d_destroy(Geometry1D* geometry);

int sk_geometry1d_get_num_altitudes(const Geometry1D* geometry);

int sk_geometry1d_get_altitudes(const Geometry1D* geometry, double* altitudes);

int sk_geometry1d_get_refractive_index_ptr(const Geometry1D* geometry,
                                           double** refractive_index);

Geometry2D*
sk_geometry2d_create(double cos_sza, double saa, double earth_radius,
                     const double* altitude_grid_values, int num_altitudes,
                     const double* horizontal_angle_grid_values,
                     int num_horizontal_locations, int altitude_interp_method);

void sk_geometry2d_destroy(Geometry2D* geometry);

int sk_geometry2d_get_location_shape(const Geometry2D* geometry,
                                     int* num_horizontal_locations,
                                     int* num_altitudes);

int sk_geometry2d_get_altitudes(const Geometry2D* geometry, double* altitudes);

int sk_geometry2d_get_horizontal_angles(const Geometry2D* geometry,
                                        double* horizontal_angles);

int sk_geometry2d_get_location_index(const Geometry2D* geometry,
                                     int altitude_index, int horizontal_index,
                                     int* location_index);

#ifdef __cplusplus
}
#endif
