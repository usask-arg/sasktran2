#include <sasktran2/solartransmission.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace sasktran2::solartransmission {
    namespace {
        constexpr double IMPACT_PARAMETER_MERGE_TOLERANCE_M = 1e-6;
        constexpr double LOW_ALTITUDE_IMPACT_SPACING_M = 500.0;
        constexpr double HIGH_ALTITUDE_IMPACT_SPACING_M = 1'000.0;
        constexpr double IMPACT_SPACING_TRANSITION_ALTITUDE_M = 30'000.0;
        constexpr int BELOW_GROUND_IMPACT_COUNT = 128;

        void add_unique_impact_parameter(std::vector<double>& values,
                                         double value) {
            if (values.empty() || std::abs(values.back() - value) >
                                      IMPACT_PARAMETER_MERGE_TOLERANCE_M) {
                values.push_back(value);
            }
        }
    } // namespace

    void SolarRayTable::initialize() {
        m_sun_unit = m_geometry.coordinates().sun_unit().normalized();
        m_perpendicular_unit = m_sun_unit.unitOrthogonal().normalized();
        m_ground_radius = m_geometry.coordinates().earth_radius() +
                          m_geometry.altitude_grid().grid()(0);
        m_top_radius =
            m_geometry.coordinates().earth_radius() +
            m_geometry.altitude_grid().grid()(Eigen::placeholders::last);

        std::vector<double> candidates;
        const int below_ground_count = BELOW_GROUND_IMPACT_COUNT;
        candidates.reserve(below_ground_count + 128);
        for (int index = 0; index < below_ground_count; ++index) {
            const double fraction =
                static_cast<double>(index) / (below_ground_count - 1.0);
            const double theta =
                0.5 * EIGEN_PI * (1.0 - std::pow(1.0 - fraction, 2.5));
            candidates.push_back(m_ground_radius * std::sin(theta));
        }
        const double transition_radius =
            std::min(m_top_radius,
                     m_ground_radius + IMPACT_SPACING_TRANSITION_ALTITUDE_M);
        for (double impact = m_ground_radius; impact < transition_radius;
             impact += LOW_ALTITUDE_IMPACT_SPACING_M) {
            candidates.push_back(impact);
        }
        for (double impact = transition_radius; impact < m_top_radius;
             impact += HIGH_ALTITUDE_IMPACT_SPACING_M) {
            candidates.push_back(impact);
        }
        candidates.push_back(m_top_radius);
        std::sort(candidates.begin(), candidates.end());

        m_impact_parameters.clear();
        m_impact_parameters.reserve(candidates.size());
        for (const double candidate : candidates) {
            add_unique_impact_parameter(m_impact_parameters, candidate);
        }

        m_rays.resize(m_impact_parameters.size());
        m_cumulative_row_offsets.assign(m_impact_parameters.size(), -1);
        m_num_cumulative_rows = 0;
    }

    void SolarRayTable::ensure_ray(int ray_index) {
        if (m_cumulative_row_offsets.at(ray_index) >= 0) {
            return;
        }

        sasktran2::viewinggeometry::ViewingRay solar_ray;
        solar_ray.look_away = -m_sun_unit;
        solar_ray.relative_azimuth = 0.0;
        const double impact_parameter = m_impact_parameters[ray_index];
        const double top_coordinate =
            std::sqrt(std::max(0.0, m_top_radius * m_top_radius -
                                        impact_parameter * impact_parameter));
        solar_ray.observer.position = impact_parameter * m_perpendicular_unit +
                                      top_coordinate * m_sun_unit;
        solar_ray.observer.on_exact_altitude = true;
        solar_ray.observer.lower_alt_index =
            static_cast<int>(m_geometry.altitude_grid().grid().size()) - 1;
        m_raytracer.trace_ray(solar_ray, m_rays[ray_index], false);

        m_cumulative_row_offsets[ray_index] = m_num_cumulative_rows;
        m_num_cumulative_rows +=
            static_cast<int>(m_rays[ray_index].layers.size()) + 1;
    }

    SolarRayTable::SingleRayQuery
    SolarRayTable::query_ray(int ray_index, double ray_coordinate) const {
        SingleRayQuery result;
        result.sample.ray_index = ray_index;
        const auto& ray = m_rays.at(ray_index);
        const double impact_parameter = m_impact_parameters.at(ray_index);
        const double top_coordinate =
            std::sqrt(std::max(0.0, m_top_radius * m_top_radius -
                                        impact_parameter * impact_parameter));
        const double coordinate_tolerance =
            128.0 * std::numeric_limits<double>::epsilon() * m_top_radius;

        if (ray_coordinate >= top_coordinate - coordinate_tolerance) {
            result.sample.boundary_index = static_cast<int>(ray.layers.size());
            result.valid = true;
            return result;
        }

        for (int layer_index = 0; layer_index < ray.layers.size();
             ++layer_index) {
            const auto& layer = ray.layers[layer_index];
            const double entrance_coordinate =
                layer.entrance.position.dot(m_sun_unit);
            const double exit_coordinate = layer.exit.position.dot(m_sun_unit);
            const double high = std::max(entrance_coordinate, exit_coordinate) +
                                coordinate_tolerance;
            const double low = std::min(entrance_coordinate, exit_coordinate) -
                               coordinate_tolerance;
            if (ray_coordinate < low || ray_coordinate > high) {
                continue;
            }

            // Ray storage is ordered from the far end back toward the
            // observer.  The boundary immediately sunward of this layer is
            // therefore layer_index + 1.
            result.sample.boundary_index = layer_index + 1;
            const double sunward_coordinate =
                std::max(entrance_coordinate, exit_coordinate);
            const double distance = sunward_coordinate - ray_coordinate;
            if (distance <= coordinate_tolerance) {
                result.valid = true;
                return result;
            }

            sasktran2::raytracing::TracedRay partial_ray;
            partial_ray.layers.resize(1);
            auto& partial = partial_ray.layers.front();
            partial.entrance.position =
                impact_parameter * m_perpendicular_unit +
                sunward_coordinate * m_sun_unit;
            partial.exit.position = impact_parameter * m_perpendicular_unit +
                                    ray_coordinate * m_sun_unit;
            partial.r_entrance = partial.entrance.radius();
            partial.r_exit = partial.exit.radius();
            partial.layer_distance = distance;
            partial.curvature_factor = 1.0;
            partial.type = sasktran2::raytracing::LayerType::partial;
            sasktran2::raytracing::add_od_quadrature(
                partial, sasktran2::geometrytype::spherical,
                m_geometry.altitude_grid().interpolation_method());
            std::vector<std::pair<int, double>> interpolation_workspace;
            sasktran2::raytracing::add_interpolation_weights(
                partial_ray, 0, m_geometry, interpolation_workspace);
            const auto weights = partial_ray.optical_depth_weights(0);
            if (weights.size() > MAX_PARTIAL_WEIGHTS) {
                throw std::runtime_error(
                    "Solar ray partial shell exceeds fixed stencil capacity");
            }
            result.sample.partial.count = weights.size();
            for (std::size_t index = 0; index < weights.size(); ++index) {
                const auto weight = weights[index];
                result.sample.partial.index[index] = weight.first;
                result.sample.partial.weight[index] = weight.second;
            }
            result.valid = true;
            return result;
        }

        if (ray.ground_is_hit) {
            result.shadowed = true;
            return result;
        }

        // A point beyond the far-side TOA sees the full limb chord between it
        // and the sun-facing boundary.
        result.sample.boundary_index = 0;
        result.valid = true;
        return result;
    }

    SolarRayTable::Query
    SolarRayTable::query(const sasktran2::Location& location) {
        if (m_rays.empty()) {
            throw std::logic_error("Solar ray table is not initialized");
        }

        Query result;
        const double ray_coordinate = location.position.dot(m_sun_unit);
        const double radius = location.radius();
        const double impact_parameter =
            location.position.cross(m_sun_unit).norm();
        if (impact_parameter < m_ground_radius && ray_coordinate < 0.0) {
            result.shadowed = true;
            return result;
        }

        auto upper =
            std::lower_bound(m_impact_parameters.begin(),
                             m_impact_parameters.end(), impact_parameter);
        int upper_index =
            static_cast<int>(std::distance(m_impact_parameters.begin(), upper));
        if (upper_index == m_impact_parameters.size()) {
            upper_index = static_cast<int>(m_impact_parameters.size()) - 1;
        }
        int lower_index = std::max(0, upper_index - 1);
        if (upper_index == 0 ||
            std::abs(m_impact_parameters[upper_index] - impact_parameter) <=
                IMPACT_PARAMETER_MERGE_TOLERANCE_M) {
            lower_index = upper_index;
        }

        const double coordinate_sign = ray_coordinate < 0.0 ? -1.0 : 1.0;
        const auto query_at_same_radius = [&](int ray_index) {
            const double table_impact = m_impact_parameters[ray_index];
            if (table_impact > radius) {
                return SingleRayQuery{};
            }
            ensure_ray(ray_index);
            return query_ray(
                ray_index,
                coordinate_sign *
                    std::sqrt(std::max(0.0, radius * radius -
                                                table_impact * table_impact)));
        };
        const auto lower_sample = query_at_same_radius(lower_index);
        const auto upper_sample = upper_index == lower_index
                                      ? lower_sample
                                      : query_at_same_radius(upper_index);
        if (!lower_sample.valid && !upper_sample.valid) {
            result.shadowed = lower_sample.shadowed || upper_sample.shadowed;
            return result;
        }
        if (!lower_sample.valid) {
            result.samples[0] = upper_sample.sample;
            result.samples[0].interpolation_weight = 1.0;
            result.count = 1;
            return result;
        }
        if (!upper_sample.valid || upper_index == lower_index) {
            result.samples[0] = lower_sample.sample;
            result.samples[0].interpolation_weight = 1.0;
            result.count = 1;
            return result;
        }

        const double interval =
            m_impact_parameters[upper_index] - m_impact_parameters[lower_index];
        const double upper_weight =
            interval > 0.0
                ? (impact_parameter - m_impact_parameters[lower_index]) /
                      interval
                : 0.0;
        result.samples[0] = lower_sample.sample;
        result.samples[0].interpolation_weight = 1.0 - upper_weight;
        result.samples[1] = upper_sample.sample;
        result.samples[1].interpolation_weight = upper_weight;
        result.count = 2;
        return result;
    }
} // namespace sasktran2::solartransmission
