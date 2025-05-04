#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ViewingGeometry ViewingGeometry;

ViewingGeometry* sk_viewing_geometry_create();

void sk_viewing_geometry_destroy(ViewingGeometry* geometry);

void sk_viewing_geometry_add_ground_viewing_solar(ViewingGeometry* geometry,
                                                  double cos_sza,
                                                  double relative_azimuth_angle,
                                                  double observeraltitude,
                                                  double cos_viewing_zenith);

int sk_viewing_geometry_add_tangent_altitude_solar(
    ViewingGeometry* geometry, double tangent_altitude_m,
    double relative_azimuth_angle, double observeraltitude, double cos_sza);

int sk_viewing_geometry_add_solar_angles_observer_location(
    ViewingGeometry* geometry, double cos_sza, double relative_azimuth_angle,
    double cos_viewing_zenith, double observeraltitude);

int sk_viewing_geometry_num_rays(ViewingGeometry* geometry, int* num_rays);

#ifdef __cplusplus
}
#endif
