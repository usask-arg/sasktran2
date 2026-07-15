#include <sasktran2.h>
#include "sasktran2/viewinggeometry.h"

namespace sasktran2::viewinggeometry {
    TangentAltitude::TangentAltitude(double tangentaltitude,
                                     double viewing_azimuth,
                                     double observeraltitude,
                                     double horizontal_angle,
                                     double tangent_out_of_plane_angle)
        : m_observeraltitude(observeraltitude),
          m_tangentaltitude(tangentaltitude),
          m_viewing_azimuth(viewing_azimuth),
          m_horizontal_angle(horizontal_angle),
          m_tangent_out_of_plane_angle(tangent_out_of_plane_angle) {
        if (!std::isfinite(m_tangentaltitude)) {
            throw std::invalid_argument("Tangent altitude must be finite");
        }
        if (!std::isfinite(m_observeraltitude) ||
            m_observeraltitude < m_tangentaltitude) {
            throw std::invalid_argument(
                "Observer altitude must be finite and greater than or equal "
                "to the tangent altitude");
        }
        if (!std::isfinite(m_horizontal_angle) ||
            !std::isfinite(m_viewing_azimuth) ||
            !std::isfinite(m_tangent_out_of_plane_angle)) {
            throw std::invalid_argument("Viewing angles must be finite");
        }
    }

    ViewingRay
    TangentAltitude::construct_ray(const sasktran2::Coordinates& geometry) {
        if (geometry.geometry_type() != sasktran2::geometrytype::spherical) {
            auto msg =
                "Error constructing ray in TangentAltitude::construct_ray, "
                "TangentAltitude ray construction can only be used in "
                "spherical geometry mode.";
            spdlog::error(msg);
            throw std::invalid_argument(msg);
        }
        ViewingRay ray;

        // Get the unit vector pointing to the tangent altitude
        Eigen::Vector3d uv = geometry.unit_vector_from_angles(
            m_horizontal_angle, m_tangent_out_of_plane_angle);
        Eigen::Vector3d tangent_point =
            uv * (geometry.earth_radius() + m_tangentaltitude);

        // And the local coordinate system at the tangent point
        std::pair<Eigen::Vector3d, Eigen::Vector3d> x_y =
            geometry.local_x_y_from_angles(m_horizontal_angle,
                                           m_tangent_out_of_plane_angle);

        // Calculate the local look vector
        ray.look_away = cos(m_viewing_azimuth) * x_y.first +
                        sin(m_viewing_azimuth) * x_y.second;

        // Now we need to back calculate the observer position based upon
        // altitude
        double s =
            sqrt(math::sqr(geometry.earth_radius() + m_observeraltitude) -
                 math::sqr(geometry.earth_radius() + m_tangentaltitude));

        ray.observer.position = tangent_point - s * ray.look_away;

        ray.relative_azimuth = m_viewing_azimuth;

        return ray;
    }

    std::string TangentAltitude::to_string() const {
        return fmt::format(
            "Tangent Viewing Ray: tangent_altitude_m: {}, viewing_azimuth: "
            "{}, observer_altitude_m: {}, horizontal_angle: {}, "
            "tangent_out_of_plane_angle: {}",
            m_tangentaltitude, m_viewing_azimuth, m_observeraltitude,
            m_horizontal_angle, m_tangent_out_of_plane_angle);
    }
} // namespace sasktran2::viewinggeometry
