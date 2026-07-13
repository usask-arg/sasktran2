#include <sasktran2/solartransmission.h>

namespace {
    bool solar_ray_hits_ground(const sasktran2::Location& location,
                               const sasktran2::Geometry2D& geometry) {
        const Eigen::Vector3d& direction = geometry.coordinates().sun_unit();
        const double ground_radius = geometry.coordinates().earth_radius() +
                                     geometry.altitude_grid().grid()[0];
        const double projected_distance = location.position.dot(direction);
        double discriminant =
            projected_distance * projected_distance -
            (location.position.squaredNorm() - ground_radius * ground_radius);
        const double discriminant_scale =
            std::max({1.0, projected_distance * projected_distance,
                      std::abs(location.position.squaredNorm() -
                               ground_radius * ground_radius)});
        const double discriminant_tolerance =
            128.0 * std::numeric_limits<double>::epsilon() * discriminant_scale;
        if (discriminant < -discriminant_tolerance) {
            return false;
        }
        discriminant = std::max(0.0, discriminant);

        const double far_intersection =
            -projected_distance + std::sqrt(discriminant);
        const double tolerance =
            64.0 * std::numeric_limits<double>::epsilon() * ground_radius;
        return far_intersection > tolerance;
    }
} // namespace

namespace sasktran2::solartransmission {
    void SolarTransmissionExact::generate_geometry_matrix(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        Eigen::MatrixXd& od_matrix, std::vector<bool>& ground_hit_flag) const {
        // First calculate the number of points we need to create the matrix for
        // We calculate solar transmission at the boundaries of layers, so it is
        // nlayer+1 for each ray
        int numpoints = 0;
        for (const auto& ray : rays) {
            numpoints += (int)ray.layers.size() + 1;
        }

        // od matrix is such that matrix @ extinction = od
        od_matrix.resize(numpoints, m_geometry.size());
        od_matrix.setZero();

        // Have to handle rays that hit the ground separately since they have no
        // solar transmission
        ground_hit_flag.resize(numpoints, false);

        sasktran2::viewinggeometry::ViewingRay ray_to_sun;

        ray_to_sun.look_away = m_geometry.coordinates().sun_unit();

        // common memory
        std::vector<std::pair<int, double>> index_weights;

        raytracing::TracedRay traced_ray;

        int row = 0;
        for (int i = 0; i < rays.size(); ++i) {
            const auto& ray = rays[i];
            for (int j = 0; j < ray.layers.size(); ++j) {
                const auto& layer = ray.layers[j];

                if (j == 0) {
                    // End layer at TOA, need to use layer exit
                    ray_to_sun.observer = layer.exit;

                    // Always don't use refraction for this
                    m_raytracer->trace_ray(ray_to_sun, traced_ray, false);

                    if (!traced_ray.ground_is_hit) {
                        assign_dense_matrix_column(row, traced_ray, m_geometry,
                                                   od_matrix, index_weights);
                    } else {
                        ground_hit_flag[row] = true;
                    }
                    ++row;
                }

                ray_to_sun.observer = layer.entrance;
                m_raytracer->trace_ray(ray_to_sun, traced_ray, false);

                if (!traced_ray.ground_is_hit) {
                    assign_dense_matrix_column(row, traced_ray, m_geometry,
                                               od_matrix, index_weights);
                } else {
                    ground_hit_flag[row] = true;
                }

                ++row;
            }
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    void SolarTransmissionExact::generate_geometry_matrix(
        const std::vector<sasktran2::raytracing::TracedRay2D>& rays,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& od_matrix,
        std::vector<bool>& ground_hit_flag) const {
        int numpoints = 0;
        for (const auto& ray : rays) {
            numpoints += static_cast<int>(ray.layers.size()) + 1;
        }

        od_matrix.resize(numpoints, m_geometry.size());
        ground_hit_flag.assign(numpoints, false);

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(static_cast<std::size_t>(numpoints) * 16);

        sasktran2::viewinggeometry::ViewingRay ray_to_sun;
        ray_to_sun.look_away = m_geometry.coordinates().sun_unit();
        raytracing::TracedRay2D traced_ray;

        const auto append_ray = [&](int row) {
            for (const auto& layer : traced_ray.layers) {
                const std::array<int, 4> indices = {
                    m_geometry_2d->location_index(layer.altitude_cell,
                                                  layer.horizontal_cell),
                    m_geometry_2d->location_index(layer.altitude_cell + 1,
                                                  layer.horizontal_cell),
                    m_geometry_2d->location_index(layer.altitude_cell,
                                                  layer.horizontal_cell + 1),
                    m_geometry_2d->location_index(layer.altitude_cell + 1,
                                                  layer.horizontal_cell + 1)};
                for (int local_index = 0; local_index < 4; ++local_index) {
                    const double weight =
                        layer.integrated_od.weights[local_index];
                    if (weight != 0.0) {
                        triplets.emplace_back(row, indices[local_index],
                                              weight);
                    }
                }
            }
        };

        int row = 0;
        for (const auto& ray : rays) {
            for (int layer_index = 0; layer_index < ray.layers.size();
                 ++layer_index) {
                const auto& layer = ray.layers[layer_index];
                if (layer_index == 0) {
                    ray_to_sun.observer = layer.exit;
                    if (solar_ray_hits_ground(layer.exit, *m_geometry_2d)) {
                        ground_hit_flag[row] = true;
                    } else {
                        m_raytracer_2d->trace_ray(ray_to_sun, traced_ray);
                        append_ray(row);
                    }
                    ++row;
                }

                ray_to_sun.observer = layer.entrance;
                if (solar_ray_hits_ground(layer.entrance, *m_geometry_2d)) {
                    ground_hit_flag[row] = true;
                } else {
                    m_raytracer_2d->trace_ray(ray_to_sun, traced_ray);
                    append_ray(row);
                }
                ++row;
            }
        }
        od_matrix.setFromTriplets(triplets.begin(), triplets.end());
    }
#endif
} // namespace sasktran2::solartransmission
