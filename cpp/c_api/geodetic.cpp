#include "geodetic.h"
#include "internal_types.h"

Geodetic* sk_geodetic_create(double equatorial_radius,
                             double flattening_factor) {
    return new Geodetic(equatorial_radius, flattening_factor);
}

void sk_geodetic_destroy(Geodetic* geodetic) { delete geodetic; }

int sk_geodetic_get_altitude(const Geodetic* geodetic, double* altitude) {
    if (geodetic == nullptr || altitude == nullptr) {
        return -1; // Error: null pointer
    }
    *altitude = geodetic->impl->altitude();
    return 0; // Success
}

int sk_geodetic_get_latitude(const Geodetic* geodetic, double* latitude) {
    if (geodetic == nullptr || latitude == nullptr) {
        return -1; // Error: null pointer
    }
    *latitude = geodetic->impl->latitude();
    return 0; // Success
}

int sk_geodetic_get_longitude(const Geodetic* geodetic, double* longitude) {
    if (geodetic == nullptr || longitude == nullptr) {
        return -1; // Error: null pointer
    }
    *longitude = geodetic->impl->longitude();
    return 0; // Success
}

int sk_geodetic_get_location(const Geodetic* geodetic, double* x, double* y,
                             double* z) {
    if (geodetic == nullptr || x == nullptr || y == nullptr || z == nullptr) {
        return -1; // Error: null pointer
    }
    auto location = geodetic->impl->location();
    *x = location.x();
    *y = location.y();
    *z = location.z();
    return 0; // Success
}

int sk_geodetic_get_local_south(const Geodetic* geodetic, double* x, double* y,
                                double* z) {
    if (geodetic == nullptr || x == nullptr || y == nullptr || z == nullptr) {
        return -1; // Error: null pointer
    }
    auto local_south = geodetic->impl->local_south();
    *x = local_south.x();
    *y = local_south.y();
    *z = local_south.z();
    return 0; // Success
}

int sk_geodetic_get_local_up(const Geodetic* geodetic, double* x, double* y,
                             double* z) {
    if (geodetic == nullptr || x == nullptr || y == nullptr || z == nullptr) {
        return -1; // Error: null pointer
    }
    auto local_up = geodetic->impl->local_up();
    *x = local_up.x();
    *y = local_up.y();
    *z = local_up.z();
    return 0; // Success
}

int sk_geodetic_get_local_west(const Geodetic* geodetic, double* x, double* y,
                               double* z) {
    if (geodetic == nullptr || x == nullptr || y == nullptr || z == nullptr) {
        return -1; // Error: null pointer
    }
    auto local_west = geodetic->impl->local_west();
    *x = local_west.x();
    *y = local_west.y();
    *z = local_west.z();
    return 0; // Success
}

int sk_geodetic_get_altitude_intercepts(
    const Geodetic* geodetic, double altitude, double observer_x,
    double observer_y, double observer_z, double look_vector_x,
    double look_vector_y, double look_vector_z, double* x1, double* y1,
    double* z1, double* x2, double* y2, double* z2) {
    if (geodetic == nullptr || x1 == nullptr || y1 == nullptr ||
        z1 == nullptr || x2 == nullptr || y2 == nullptr || z2 == nullptr) {
        return -1; // Error: null pointer
    }
    Eigen::Vector3d observer(observer_x, observer_y, observer_z);
    Eigen::Vector3d look_vector(look_vector_x, look_vector_y, look_vector_z);
    auto intercepts =
        geodetic->impl->altitude_intercepts(altitude, observer, look_vector);
    *x1 = intercepts.first.x();
    *y1 = intercepts.first.y();
    *z1 = intercepts.first.z();
    *x2 = intercepts.second.x();
    *y2 = intercepts.second.y();
    *z2 = intercepts.second.z();
    return 0; // Success
}

int sk_geodetic_from_lat_lon_altitude(const Geodetic* geodetic, double latitude,
                                      double longitude, double altitude) {
    if (geodetic == nullptr) {
        return -1; // Error: null pointer
    }
    geodetic->impl->from_lat_lon_alt(latitude, longitude, altitude);
    return 0; // Success
}

int sk_geodetic_from_tangent_altitude(const Geodetic* geodetic, double altitude,
                                      double observer_x, double observer_y,
                                      double observer_z, double boresight_x,
                                      double boresight_y, double boresight_z,
                                      double* look_vector_x,
                                      double* look_vector_y,
                                      double* look_vector_z) {
    if (geodetic == nullptr || look_vector_x == nullptr ||
        look_vector_y == nullptr || look_vector_z == nullptr) {
        return -1; // Error: null pointer
    }
    Eigen::Vector3d observer(observer_x, observer_y, observer_z);
    Eigen::Vector3d boresight(boresight_x, boresight_y, boresight_z);
    auto look_vector =
        geodetic->impl->from_tangent_altitude(altitude, observer, boresight);
    *look_vector_x = look_vector.x();
    *look_vector_y = look_vector.y();
    *look_vector_z = look_vector.z();
    return 0; // Success
}

int sk_geodetic_from_tangent_point(const Geodetic* geodetic, double observer_x,
                                   double observer_y, double observer_z,
                                   double look_vector_x, double look_vector_y,
                                   double look_vector_z) {
    if (geodetic == nullptr) {
        return -1; // Error: null pointer
    }
    Eigen::Vector3d observer(observer_x, observer_y, observer_z);
    Eigen::Vector3d look_vector(look_vector_x, look_vector_y, look_vector_z);
    geodetic->impl->from_tangent_point(observer, look_vector);
    return 0; // Success
}

int sk_geodetic_from_xyz(const Geodetic* geodetic, double x, double y,
                         double z) {
    if (geodetic == nullptr) {
        return -1; // Error: null pointer
    }
    Eigen::Vector3d location(x, y, z);
    geodetic->impl->from_xyz(location);
    return 0; // Success
}

int sk_geodetic_is_valid(const Geodetic* geodetic, int* is_valid) {
    if (geodetic == nullptr || is_valid == nullptr) {
        return -1; // Error: null pointer
    }
    *is_valid = geodetic->impl->is_valid() ? 1 : 0;
    return 0; // Success
}

int sk_geodetic_get_osculating_spheroid(const Geodetic* geodetic,
                                        double* radius, double* offset_x,
                                        double* offset_y, double* offset_z) {
    if (geodetic == nullptr || radius == nullptr || offset_x == nullptr ||
        offset_y == nullptr || offset_z == nullptr) {
        return -1; // Error: null pointer
    }
    Eigen::Vector3d offset;
    geodetic->impl->get_osculating_spheroid(radius, &offset);
    *offset_x = offset.x();
    *offset_y = offset.y();
    *offset_z = offset.z();
    return 0; // Success
}
