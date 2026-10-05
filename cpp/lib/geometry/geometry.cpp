#include <sasktran2/geometry.h>
#include <sasktran2/math/trig.h>
#include <sasktran2/validation/validation.h>

namespace sasktran2 {
    Coordinates::Coordinates(double cos_sza, double saa, double earth_radius,
                             geometrytype geotype, bool force_sun_z)
        : m_earth_radius(earth_radius), m_geotype(geotype),
          m_force_sun_z(force_sun_z) {

        if (!force_sun_z) {
            // Create the default x,y,z local coordinates
            m_x_unit << 1, 0, 0;
            m_y_unit << 0, 1, 0;
            m_z_unit << 0, 0, 1;

            Eigen::Vector3d sun_horiz =
                cos(saa) * m_x_unit + sin(saa) * m_y_unit;

            m_sun_unit =
                cos_sza * m_z_unit + sqrt(1 - cos_sza * cos_sza) * sun_horiz;
        } else {
            // We are forcing the sun unit vector to be z
            m_sun_unit << 0, 0, 1;

            // And the reference point unit vector (m_z_unit) then must be
            m_z_unit << sqrt(1 - cos_sza * cos_sza), 0, cos_sza;
            // y will be unchanged
            m_y_unit << 0, 1, 0;

            m_x_unit = m_y_unit.cross(m_z_unit);
        }
        validate();
    }

    Coordinates::Coordinates(Eigen::Vector3d ref_point_unit,
                             Eigen::Vector3d ref_plane_unit,
                             Eigen::Vector3d sun_unit, double earth_radius,
                             geometrytype geotype)
        : m_x_unit(ref_plane_unit),
          m_y_unit(ref_point_unit.cross(ref_plane_unit).normalized()),
          m_z_unit(ref_point_unit), m_sun_unit(sun_unit),
          m_earth_radius(earth_radius), m_geotype(geotype),
          m_force_sun_z(false) {
        validate();
    }

    void Coordinates::validate() const {
        switch (m_geotype) {
        case geometrytype::planeparallel:
        case geometrytype::pseudospherical:
        case geometrytype::spherical:
        case geometrytype::ellipsoidal:
            break;
        default:
            spdlog::critical("Invalid geometry type: {}",
                             static_cast<int>(m_geotype));
            sasktran2::validation::throw_configuration_error();
        }

        if (!std::isfinite(m_earth_radius) || m_earth_radius < 0) {
            spdlog::critical("Invalid earth radius: {}", m_earth_radius);
            sasktran2::validation::throw_configuration_error();
        }

        if (!m_x_unit.allFinite() || !m_y_unit.allFinite() ||
            !m_z_unit.allFinite() || !m_sun_unit.allFinite()) {
            spdlog::critical("Coordinate unit vectors must be finite");
            sasktran2::validation::throw_configuration_error();
        }

        constexpr double unit_tolerance = 1e-12;
        if (std::abs(m_x_unit.squaredNorm() - 1.0) > unit_tolerance ||
            std::abs(m_y_unit.squaredNorm() - 1.0) > unit_tolerance ||
            std::abs(m_z_unit.squaredNorm() - 1.0) > unit_tolerance ||
            std::abs(m_sun_unit.squaredNorm() - 1.0) > unit_tolerance ||
            std::abs(m_x_unit.dot(m_y_unit)) > unit_tolerance ||
            std::abs(m_x_unit.dot(m_z_unit)) > unit_tolerance ||
            std::abs(m_y_unit.dot(m_z_unit)) > unit_tolerance) {
            spdlog::critical(
                "Coordinate basis and sun vectors must be unit vectors, and "
                "the coordinate basis must be orthogonal");
            sasktran2::validation::throw_configuration_error();
        }

        if (m_geotype == geometrytype::planeparallel) {
            // Have to make sure that the SZA is not greater than 90, or that
            // cos_sza is positive
            double cos_sza = m_sun_unit.dot(m_z_unit);

            if (cos_sza < 0) {
                spdlog::critical(
                    "Invalid solar zenith angle for planeparallel geometry, "
                    "cos_sza: {}, and it should be greater than 0",
                    cos_sza);
                sasktran2::validation::throw_configuration_error();
            }
        }
    }

    Eigen::Vector3d Coordinates::unit_vector_from_angles(double theta,
                                                         double phi) const {
        // First calculate the unit vector in the primary tangent plane
        Eigen::AngleAxis<double> primary_transform(theta, m_y_unit);

        Eigen::Vector3d in_plane = primary_transform.matrix() * m_z_unit;

        // Then rotate the unitvector in the secondary plane, the rotation
        // vector is the x_unit vector rotated to the location of
        Eigen::AngleAxis<double> transform(phi, primary_transform.matrix() *
                                                    m_x_unit);

        return transform.matrix() * in_plane;
    }

    std::pair<Eigen::Vector3d, Eigen::Vector3d>
    Coordinates::local_x_y_from_angles(double theta, double phi) const {
        // First calculate the unit vector in the primary tangent plane
        Eigen::AngleAxis<double> primary_transform(theta, m_y_unit);

        // Then rotate the unitvector in the secondary plane, the rotation
        // vector is the x_unit vector rotated to the location of
        Eigen::AngleAxis<double> transform(phi, primary_transform.matrix() *
                                                    m_x_unit);

        return {transform.matrix() * primary_transform.matrix() * m_x_unit,
                transform.matrix() * primary_transform.matrix() * m_y_unit};
    }

    std::pair<double, double>
    Coordinates::angles_from_unit_vector(const Eigen::Vector3d& uv) const {
        // Start by calculating theta, which we can do by projecting the
        // unitvector into the primary plane
        Eigen::Vector3d projected =
            (uv - uv.dot(m_y_unit) * m_y_unit).normalized();

        double theta = atan2(projected.dot(m_x_unit), projected.dot(m_z_unit));

        // cos_phi is simply the angle between the projected unit vector and the
        // original unitvector
        double cos_phi = projected.dot(uv);
        double phi = acos(cos_phi);

        // We are still missing the sign of phi, we have to check what side of
        // y_unit we are on
        double sign_phi = uv.dot(m_y_unit);
        if (sign_phi > 0) {
            phi *= -1;
        }
        return std::pair<double, double>(theta, phi);
    }

    Eigen::Vector3d
    Coordinates::solar_coordinate_vector(double cos_sza, double saa,
                                         double altitude) const {
        if (m_geotype == sasktran2::geometrytype::planeparallel ||
            m_geotype == sasktran2::geometrytype::pseudospherical) {
            return m_z_unit * (altitude + m_earth_radius);
        }

        // First find the plane to rotate the sun vector around

        Eigen::Vector3d normal = m_sun_unit.cross(m_z_unit);

        constexpr double singular_squared_tolerance = 1.0e-24;
        if (normal.squaredNorm() <= singular_squared_tolerance) {
            // Special case where sun is parallel to z-axis
            normal = m_y_unit;
        } else {
            normal = normal.normalized();
        }

        Eigen::AngleAxis<double> zenith_transform(acos(cos_sza), normal);

        Eigen::AngleAxis<double> azimuth_transform(saa, m_sun_unit);

        return azimuth_transform.matrix() * zenith_transform.matrix() *
               m_sun_unit * (altitude + m_earth_radius);
    }

    std::pair<double, double> Coordinates::solar_angles_at_location(
        const Eigen::Vector3d& location) const {
        std::pair<double, double> result;

        result.first = location.normalized().dot(m_sun_unit);

        // Solar azimuth at location
        Eigen::Vector3d sun_horiz =
            m_sun_unit -
            location.normalized() * (location.normalized().dot(m_sun_unit));

        result.second = std::numeric_limits<double>::quiet_NaN();

        return result;
    }

    Eigen::Vector3d Coordinates::look_vector_from_azimuth(
        const Eigen::Vector3d& location, double saa, double cos_viewing) const {
        // Calculate the normalized sun vector at the location

        Eigen::Vector3d local_up;
        if (m_geotype == sasktran2::geometrytype::spherical) {
            local_up = location.normalized();
        } else {
            local_up = m_z_unit;
        }

        Eigen::Vector3d sun_horiz =
            m_sun_unit - local_up * (local_up.dot(m_sun_unit));

        if (sun_horiz.norm() == 0) {
            // sun azimuth is ambiguous, and so we define? the x vector as
            // horizontal sun
            // TODO: Check this
            sun_horiz = m_y_unit;
        }
        sun_horiz = sun_horiz.normalized();

        // Rotate the horizontal sun vector around SAA
        Eigen::AngleAxis<double> azimuth_transform(-saa, local_up);

        Eigen::Vector3d horiz_look = azimuth_transform.matrix() * sun_horiz;

        double viewing_angle = sasktran2::math::PiOver2 - acos(-cos_viewing);

        Eigen::AngleAxis<double> vertical_transform(viewing_angle,
                                                    local_up.cross(horiz_look));

        return vertical_transform.matrix() * horiz_look;
    }

    std::pair<double, double>
    Coordinates::stokes_rotation(const Eigen::Vector3d& look_vector,
                                 const Eigen::Vector3d& from_reference,
                                 const Eigen::Vector3d& to_reference) {
        // Project both references perpendicular to the look vector
        const Eigen::Vector3d perp_from =
            from_reference - from_reference.dot(look_vector) * look_vector;
        const Eigen::Vector3d perp_to =
            to_reference - to_reference.dot(look_vector) * look_vector;

        constexpr double parallel_tolerance = 1e-8;
        if (perp_from.norm() <= parallel_tolerance * from_reference.norm() ||
            perp_to.norm() <= parallel_tolerance * to_reference.norm()) {
            // A reference along the look vector does not define a basis
            return std::make_pair(1.0, 0.0);
        }

        // Signed angle from perp_from to perp_to about the propagation
        // direction, -look_vector.  The sign is required, the reference can
        // turn either way
        const double angle = atan2(-look_vector.dot(perp_from.cross(perp_to)),
                                   perp_from.dot(perp_to));

        return std::make_pair(cos(2 * angle), -sin(2 * angle));
    }

} // namespace sasktran2
