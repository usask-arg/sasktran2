#include <sasktran2/math/unitsphere.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace sasktran2::math {
    UnitSphereGround::UnitSphereGround(
        std::unique_ptr<const UnitSphere>&& sphere,
        const Eigen::Vector3d location)
        : m_full_sphere(std::move(sphere)), m_location(location) {
        Eigen::Vector3d quad;

        m_contributing_map.reserve(m_full_sphere->num_points() / 2);

        m_quadrature_normalization = 0.0;

        // We need to find the points that are looking up
        for (int i = 0; i < m_full_sphere->num_points(); ++i) {
            quad = m_full_sphere->get_quad_position(i);

            if (quad.dot(location) > 0) {
                // Looking up
                m_contributing_map.push_back(i);

                m_quadrature_normalization +=
                    m_full_sphere->quadrature_weight(i);
            }
        }
    }

    int UnitSphereGround::num_points() const {
        return m_contributing_map.size();
    }

    Eigen::Vector3d UnitSphereGround::get_quad_position(int index) const {
        return m_full_sphere->get_quad_position(m_contributing_map[index]);
    }

    double UnitSphereGround::quadrature_weight(int i) const {
        return m_full_sphere->quadrature_weight(m_contributing_map[i]) /
               m_quadrature_normalization * 0.5;
    }

    void UnitSphereGround::interpolate(
        const Eigen::Vector3d& direction,
        std::vector<std::pair<int, double>>& index_weights,
        int& num_interp) const {
        num_interp = std::min(3, num_points());
        index_weights.assign(static_cast<std::size_t>(num_interp),
                             {-1, std::numeric_limits<double>::infinity()});

        // Select neighbors directly from the upward hemisphere. Filtering a
        // full-sphere stencil after selection can discard every neighbor for
        // near-horizon directions.
        for (int local_index = 0; local_index < num_points(); ++local_index) {
            double squared_distance =
                (get_quad_position(local_index) - direction).squaredNorm();
            for (int slot = 0; slot < num_interp; ++slot) {
                if (squared_distance < index_weights[slot].second) {
                    for (int shifted = num_interp - 1; shifted > slot;
                         --shifted) {
                        index_weights[shifted] = index_weights[shifted - 1];
                    }
                    index_weights[slot] = {local_index, squared_distance};
                    break;
                }
            }
        }

        double total_inverse_distance = 0.0;
        for (int slot = 0; slot < num_interp; ++slot) {
            if (index_weights[slot].second < 1.0e-8) {
                for (auto& weight : index_weights) {
                    weight.second = 0.0;
                }
                index_weights[slot].second = 1.0;
                return;
            }
            total_inverse_distance +=
                1.0 / std::sqrt(index_weights[slot].second);
        }
        for (auto& weight : index_weights) {
            weight.second =
                (1.0 / std::sqrt(weight.second)) / total_inverse_distance;
        }
    }

} // namespace sasktran2::math
