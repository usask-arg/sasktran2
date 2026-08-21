#include <sasktran2/solartransmission.h>

#include <spdlog/spdlog.h>

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
                        assign_dense_matrix_column(row, traced_ray, od_matrix);
                    } else {
                        ground_hit_flag[row] = true;
                    }
                    ++row;
                }

                ray_to_sun.observer = layer.entrance;
                m_raytracer->trace_ray(ray_to_sun, traced_ray, false);

                if (!traced_ray.ground_is_hit) {
                    assign_dense_matrix_column(row, traced_ray, od_matrix);
                } else {
                    ground_hit_flag[row] = true;
                }

                ++row;
            }
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    void SolarTransmissionExact::generate_geometry_matrix(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        SolarGeometryMatrix& od_matrix,
        std::vector<bool>& ground_hit_flag) const {
        int numpoints = 0;
        for (const auto& ray : rays) {
            numpoints += static_cast<int>(ray.layers.size()) + 1;
        }

        ground_hit_flag.assign(numpoints, false);
        // A straight path through a structured grid usually accumulates about
        // one sparse coefficient per point along its dominant grid dimension.
        // Reserving at that scale avoids Eigen's doubling growth for large 2D
        // matrices while the bounded-growth fallback below handles longer
        // paths without assuming a hard upper limit.
        const Eigen::Index dominant_grid_size =
            std::max(m_geometry_2d->altitude_grid().grid().size(),
                     m_geometry_2d->horizontal_angle_grid().size());
        const Eigen::Index initial_entries_per_row = std::max<Eigen::Index>(
            16, dominant_grid_size - dominant_grid_size / 16);
        od_matrix.initialize_exact(numpoints, m_geometry.size(),
                                   static_cast<Eigen::Index>(numpoints) *
                                       initial_entries_per_row);

        sasktran2::viewinggeometry::ViewingRay ray_to_sun;
        ray_to_sun.look_away = m_geometry.coordinates().sun_unit();
        raytracing::TracedRay traced_ray;
        std::vector<std::pair<int, double>> row_weights;
        row_weights.reserve(static_cast<std::size_t>(
                                m_geometry_2d->altitude_grid().grid().size() +
                                m_geometry_2d->horizontal_angle_grid().size()) *
                            4);

        const auto ensure_entry_capacity = [&](Eigen::Index additional) {
            // Eigen doubles compressed storage when insertBack exhausts its
            // reservation. Large 2D calculations can therefore retain almost
            // twice their final matrix size because shrinking reallocations are
            // commonly kept in place by the system allocator. Grow explicitly
            // by 50% instead, bounding unused retained capacity without adding
            // a counting/tracing pass.
            od_matrix.ensure_capacity(additional);
        };

        const auto append_ray = [&](int row) {
            row_weights.clear();
            for (std::size_t layer_index = 0;
                 layer_index < traced_ray.layers.size(); ++layer_index) {
                const auto weights =
                    traced_ray.optical_depth_weights(layer_index);
                for (std::size_t index = 0; index < weights.size(); ++index) {
                    const auto weight = weights[index];
                    if (weight.second != 0.0) {
                        row_weights.push_back(weight);
                    }
                }
            }

            std::sort(row_weights.begin(), row_weights.end(),
                      [](const auto& left, const auto& right) {
                          return left.first < right.first;
                      });
            std::size_t write = 0;
            for (std::size_t begin = 0; begin < row_weights.size();) {
                const int column = row_weights[begin].first;
                double value = 0.0;
                std::size_t end = begin;
                while (end < row_weights.size() &&
                       row_weights[end].first == column) {
                    value += row_weights[end].second;
                    ++end;
                }
                if (value != 0.0) {
                    row_weights[write++] = {column, value};
                }
                begin = end;
            }
            row_weights.resize(write);

            ensure_entry_capacity(
                static_cast<Eigen::Index>(row_weights.size()));
            od_matrix.start_row(row);
            for (const auto& [column, value] : row_weights) {
                od_matrix.insert_back(row, column, value);
            }
        };

        const auto append_empty_row = [&](int row) {
            od_matrix.start_row(row);
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
                        append_empty_row(row);
                    } else {
                        m_raytracer_2d->trace_ray(ray_to_sun, traced_ray);
                        append_ray(row);
                    }
                    ++row;
                }

                ray_to_sun.observer = layer.entrance;
                if (solar_ray_hits_ground(layer.entrance, *m_geometry_2d)) {
                    ground_hit_flag[row] = true;
                    append_empty_row(row);
                } else {
                    m_raytracer_2d->trace_ray(ray_to_sun, traced_ray);
                    append_ray(row);
                }
                ++row;
            }
        }
        const Eigen::Index construction_capacity = od_matrix.allocated_size();
        od_matrix.finalize();
        spdlog::debug(
            "Geometry2D exact solar matrix: {} rows, {} nonzeros, {} peak "
            "construction entries, {} ground-blocked rows, compact storage {}",
            od_matrix.rows(), od_matrix.non_zeros(), construction_capacity,
            std::count(ground_hit_flag.begin(), ground_hit_flag.end(), true),
            od_matrix.is_compact());
    }
#endif
} // namespace sasktran2::solartransmission
