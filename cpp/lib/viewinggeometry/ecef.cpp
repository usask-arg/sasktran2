#include <sasktran2/viewinggeometry.h>

namespace sasktran2::viewinggeometry {
    ECEFViewingRay::ECEFViewingRay(Eigen::Vector3d observer_position,
                                   Eigen::Vector3d look_direction)
        : m_observer_position(std::move(observer_position)),
          m_look_direction(std::move(look_direction)) {
        if (!m_observer_position.allFinite() || !m_look_direction.allFinite() ||
            m_observer_position.norm() == 0.0 ||
            m_look_direction.norm() == 0.0) {
            throw std::invalid_argument(
                "ECEF viewing rays require finite, non-zero vectors");
        }
        m_look_direction.normalize();
    }

    ViewingRay
    ECEFViewingRay::construct_ray(const sasktran2::Coordinates& geometry) {
        ViewingRay result;
        result.observer.position = m_observer_position;
        result.look_away = m_look_direction;
        const double distance = -m_observer_position.dot(m_look_direction);
        const Eigen::Vector3d tangent =
            m_observer_position + distance * m_look_direction;
        if (distance > 0.0 && tangent.norm() > 0.0) {
            const auto angles =
                geometry.angles_from_unit_vector(tangent.normalized());
            const auto local_xy =
                geometry.local_x_y_from_angles(angles.first, angles.second);
            result.relative_azimuth =
                std::atan2(m_look_direction.dot(local_xy.second),
                           m_look_direction.dot(local_xy.first));
        } else {
            result.relative_azimuth = 0.0;
        }
        return result;
    }

    std::string ECEFViewingRay::to_string() const { return "ECEFViewingRay"; }
} // namespace sasktran2::viewinggeometry
