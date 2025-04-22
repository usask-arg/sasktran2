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

int sk_viewing_geometry_num_rays(ViewingGeometry* geometry, int* num_rays) {
    if (geometry == nullptr || num_rays == nullptr) {
        return -1; // Error: null pointer
    }
    *num_rays = static_cast<int>(geometry->impl.observer_rays().size());
    return 0; // Success
}
}
