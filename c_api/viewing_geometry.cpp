#include "viewing_geometry.h"
#include "sasktran2/viewinggeometry.h"
#include <sasktran2.h>
#include "internal_types.h"

extern "C" {
ViewingGeometry* sk_viewing_geometry_create() { return new ViewingGeometry(); }

void sk_viewing_geometry_destroy(ViewingGeometry* storage) { delete storage; }

void sk_viewing_geometry_add_ground_viewing_solar(ViewingGeometry* geometry,
                                                  double cos_sza,
                                                  double relative_azimuth_angle,
                                                  double observeraltitude,
                                                  double cos_viewing_zenith) {
    geometry->impl.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            cos_sza, relative_azimuth_angle, cos_viewing_zenith,
            observeraltitude));
}

int sk_viewing_geometry_add_tangent_altitude_solar(
    ViewingGeometry* geometry, double tangent_altitude_m,
    double relative_azimuth_angle, double observeraltitude, double cos_sza) {
    if (geometry == nullptr) {
        return -1; // Error: null pointer
    }
    geometry->impl.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(
            tangent_altitude_m, relative_azimuth_angle, observeraltitude,
            cos_sza));
    return 0; // Success
}

int sk_viewing_geometry_add_solar_angles_observer_location(
    ViewingGeometry* geometry, double cos_sza, double relative_azimuth_angle,
    double cos_viewing_zenith, double observeraltitude) {
    if (geometry == nullptr) {
        return -1; // Error: null pointer
    }
    geometry->impl.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::ViewingUpSolar>(
            cos_sza, relative_azimuth_angle, cos_viewing_zenith,
            observeraltitude));
    return 0; // Success
}

int sk_viewing_geometry_num_rays(ViewingGeometry* geometry, int* num_rays) {
    if (geometry == nullptr || num_rays == nullptr) {
        return -1; // Error: null pointer
    }
    *num_rays = static_cast<int>(geometry->impl.observer_rays().size());
    return 0; // Success
}
}
