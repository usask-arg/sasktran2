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
    SolarTransmissionExact::Query
    SolarTransmissionExact::query(const sasktran2::Location& location) const {
        Query result;
        sasktran2::viewinggeometry::ViewingRay ray_to_sun;
        ray_to_sun.observer = location;
        ray_to_sun.look_away = m_geometry.coordinates().sun_unit();

        if (m_geometry_1d != nullptr) {
            m_raytracer->trace_ray(ray_to_sun, result.ray, false);
            result.shadowed = result.ray.ground_is_hit;
            return result;
        }

#ifdef SKTRAN_RUST_SUPPORT
        if (m_geometry_2d != nullptr) {
            if (solar_ray_hits_ground(location, *m_geometry_2d)) {
                result.shadowed = true;
            } else {
                m_raytracer_2d->trace_ray(ray_to_sun, result.ray);
            }
            return result;
        }
#endif

        throw std::logic_error(
            "Exact solar transmission has no compatible ray tracer");
    }

    void SolarTransmissionExact::generate_geometry_matrix(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        Eigen::MatrixXd& od_matrix, std::vector<bool>& ground_hit_flag,
        const std::vector<bool>* required_rows) const {
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
        ground_hit_flag.assign(numpoints, false);
        if (required_rows != nullptr &&
            required_rows->size() != static_cast<std::size_t>(numpoints)) {
            throw std::invalid_argument(
                "Solar geometry row mask has an invalid size");
        }

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
                    if (required_rows == nullptr || required_rows->at(row)) {
                        ray_to_sun.observer = layer.exit;

                        // Always don't use refraction for this
                        m_raytracer->trace_ray(ray_to_sun, traced_ray, false);

                        if (!traced_ray.ground_is_hit) {
                            assign_dense_matrix_column(row, traced_ray,
                                                       od_matrix);
                        } else {
                            ground_hit_flag[row] = true;
                        }
                    }
                    ++row;
                }

                if (required_rows == nullptr || required_rows->at(row)) {
                    ray_to_sun.observer = layer.entrance;
                    m_raytracer->trace_ray(ray_to_sun, traced_ray, false);

                    if (!traced_ray.ground_is_hit) {
                        assign_dense_matrix_column(row, traced_ray, od_matrix);
                    } else {
                        ground_hit_flag[row] = true;
                    }
                }

                ++row;
            }
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    void SolarTransmissionExact::generate_geometry_matrix(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
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
        raytracing::TracedRay traced_ray;

        const auto append_ray = [&](int row) {
            for (std::size_t layer_index = 0;
                 layer_index < traced_ray.layers.size(); ++layer_index) {
                const auto weights =
                    traced_ray.optical_depth_weights(layer_index);
                for (std::size_t index = 0; index < weights.size(); ++index) {
                    const auto weight = weights[index];
                    if (weight.second != 0.0) {
                        triplets.emplace_back(row, weight.first, weight.second);
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
