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

int sk_viewing_geometry_num_rays(ViewingGeometry* geometry, int* num_rays);

#ifdef __cplusplus
}
#endif
