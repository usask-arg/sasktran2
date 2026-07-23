#include <sasktran2/math/unitsphere.h>
#include <sasktran2/math/unitsphere/lebedev.h>

namespace sasktran2::math {
    LebedevSphere::LebedevSphere(int npoints) {
        // 1) Load points/weights + faces
        Eigen::MatrixXd xyzw;  // 4×N
        Eigen::MatrixXi faces; // 3×F
        unitsphere::lebedev::get_lebedev_data(npoints, xyzw, faces);

        // 2) Optional rotation (pure rotation preserves triangle orientation)
        const Eigen::AngleAxisd rot(0.01, Eigen::Vector3d::UnitY());
        const Eigen::AngleAxisd rot2(0.01, Eigen::Vector3d::UnitZ());

        const_cast<Eigen::MatrixXd&>(m_xyz) = xyzw.topRows<3>(); // 3×N
        const_cast<Eigen::VectorXd&>(m_weights) =
            xyzw.row(3).transpose();                   // N×1
        const_cast<Eigen::MatrixXi&>(m_faces) = faces; // 3×F
    }

    template <int N>
    void
    LebedevSphere::integrate_on_grid(const Eigen::Matrix<double, N, -1>& values,
                                     Eigen::Vector<double, N>& result) {
        result.setZero();
        for (int i = 0; i < m_weights.size(); i++) {
            result += m_weights(i) * values.col(i);
        }
        result *= 4 * EIGEN_PI;
    }

    Eigen::Vector3d LebedevSphere::get_quad_position(int index) const {
        return Eigen::Vector3d(m_xyz(0, index), m_xyz(1, index),
                               m_xyz(2, index));
    }

    void LebedevSphere::interpolate(
        const Eigen::Vector3d& direction,
        std::vector<std::pair<int, double>>& index_weights,
        int& num_interp) const {
        num_interp = 3;
        index_weights.resize(num_interp);
        if (!direction.allFinite()) {
            // Preserve the established zenith-degenerate fallback. Relative
            // azimuth is undefined for exactly vertical directions, and the
            // historical search consequently returned index zero with equal
            // weights.
            for (auto& index_weight : index_weights) {
                index_weight = {0, 1.0 / num_interp};
            }
            return;
        }

        std::array<int, 3> nearest_indices = {-1, -1, -1};
        std::array<double, 3> nearest_distances = {
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()};

        // Keep the three exact Euclidean nearest neighbours. Using scalar
        // coordinates avoids constructing an Eigen temporary for every
        // quadrature point; this is a hot geometry-setup path.
        for (int i = 0; i < m_weights.size(); ++i) {
            const double dx = m_xyz(0, i) - direction.x();
            const double dy = m_xyz(1, i) - direction.y();
            const double dz = m_xyz(2, i) - direction.z();
            double sqdist = dx * dx + dy * dy + dz * dz;
            for (int j = 0; j < num_interp; ++j) {
                if (sqdist < nearest_distances[j]) {
                    for (int k = num_interp - 1; k > j; --k) {
                        nearest_distances[k] = nearest_distances[k - 1];
                        nearest_indices[k] = nearest_indices[k - 1];
                    }
                    nearest_distances[j] = sqdist;
                    nearest_indices[j] = i;
                    break;
                }
            }
        }
        for (int i = 0; i < num_interp; ++i) {
            index_weights[i] = {nearest_indices[i], nearest_distances[i]};
        }

        double total_distance = 0;
        for (int i = 0; i < num_interp; ++i) {
            // first check for a special case, if any of the distances are ~0 we
            // just want to return that one back
            if (index_weights[i].second < 1e-8) {
                for (int j = 0; j < num_interp; ++j) {
                    index_weights[j].second = 0;
                }
                index_weights[i].second = 1;
                return;
            }

            total_distance += 1 / sqrt(index_weights[i].second);
        }

        // Else we weight each one by inverse distance
        for (int i = 0; i < num_interp; ++i) {
            index_weights[i].second =
                (1 / sqrt(index_weights[i].second)) / total_distance;
        }
        /*
        for(int i = 0; i < 3; ++i) {
            BOOST_LOG_TRIVIAL(warning) << "dir:" <<
        m_xyz(Eigen::placeholders::all, index_weights[i].first).transpose() << "
        sqdist: " << index_weights[i].second;

        }
         */
    }

    template void LebedevSphere::integrate_on_grid<1>(
        const Eigen::Matrix<double, 1, -1>& values,
        Eigen::Vector<double, 1>& result);
    template void LebedevSphere::integrate_on_grid<3>(
        const Eigen::Matrix<double, 3, -1>& values,
        Eigen::Vector<double, 3>& result);

} // namespace sasktran2::math
