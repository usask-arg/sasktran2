#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Geometry1D Geometry1D;

Geometry1D* sk_geometry1d_create(double cos_sza, double saa,
                                 double earth_radius, double* grid_values,
                                 int ngrid_values, int interp_method,
                                 int geotype);

void sk_geometry1d_destroy(Geometry1D* geometry);

int sk_geometry1d_get_num_altitudes(const Geometry1D* geometry);

int sk_geometry1d_get_altitudes(const Geometry1D* geometry, double* altitudes);

int sk_geometry1d_get_refractive_index_ptr(const Geometry1D* geometry,
                                           double** refractive_index);

#ifdef __cplusplus
}
#endif
