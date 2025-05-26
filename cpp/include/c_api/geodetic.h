#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Geodetic Geodetic;

Geodetic* sk_geodetic_create(double equatorial_radius,
                             double flattening_factor);

void sk_geodetic_destroy(Geodetic* geodetic);

int sk_geodetic_get_altitude(const Geodetic* geodetic, double* altitude);
int sk_geodetic_get_latitude(const Geodetic* geodetic, double* latitude);
int sk_geodetic_get_longitude(const Geodetic* geodetic, double* longitude);
int sk_geodetic_get_location(const Geodetic* geodetic, double* x, double* y,
                             double* z);

int sk_geodetic_get_local_south(const Geodetic* geodetic, double* x, double* y,
                                double* z);

int sk_geodetic_get_local_up(const Geodetic* geodetic, double* x, double* y,
                             double* z);

int sk_geodetic_get_local_west(const Geodetic* geodetic, double* x, double* y,
                               double* z);

int sk_geodetic_get_altitude_intercepts(
    const Geodetic* geodetic, double altitude, double observer_x,
    double observer_y, double observer_z, double look_vector_x,
    double look_vector_y, double look_vector_z, double* x1, double* y1,
    double* z1, double* x2, double* y2, double* z2);
int sk_geodetic_from_lat_lon_altitude(const Geodetic* geodetic, double latitude,
                                      double longitude, double altitude);

int sk_geodetic_from_tangent_altitude(const Geodetic* geodetic, double altitude,
                                      double observer_x, double observer_y,
                                      double observer_z, double boresight_x,
                                      double boresight_y, double boresight_z,
                                      double* look_vector_x,
                                      double* look_vector_y,
                                      double* look_vector_z);

int sk_geodetic_from_tangent_point(const Geodetic* geodetic, double observer_x,
                                   double observer_y, double observer_z,
                                   double look_vector_x, double look_vector_y,
                                   double look_vector_z);

int sk_geodetic_from_xyz(const Geodetic* geodetic, double x, double y,
                         double z);

int sk_geodetic_is_valid(const Geodetic* geodetic, int* is_valid);

int sk_geodetic_get_osculating_spheroid(const Geodetic* geodetic,
                                        double* radius, double* offset_x,
                                        double* offset_y, double* offset_z);

#ifdef __cplusplus
}
#endif
