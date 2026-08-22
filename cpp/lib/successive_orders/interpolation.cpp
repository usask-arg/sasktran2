#include "interpolation.h"

#include "geometry.h"

#include <sasktran2/grids.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace sasktran2::successive_orders {
    namespace {
        std::uint32_t checked_u32(std::size_t value, const char* description) {
            if (value > std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(std::string(description) +
                                        " exceeds the compact offset range");
            }
            return static_cast<std::uint32_t>(value);
        }

        void sort_and_combine(const std::vector<std::pair<int, double>>& input,
                              std::vector<InterpolationWeight>& output) {
            output.clear();
            output.reserve(input.size());
            for (const auto& [index, weight] : input) {
                if (index < 0 || !std::isfinite(weight)) {
                    throw std::runtime_error(
                        "Invalid successive-orders interpolation weight");
                }
                if (weight != 0.0) {
                    output.push_back({index, weight});
                }
            }
            std::stable_sort(output.begin(), output.end(),
                             [](const auto& left, const auto& right) {
                                 return left.index < right.index;
                             });

            std::size_t write = 0;
            for (const auto& weight : output) {
                if (write != 0 && output[write - 1].index == weight.index) {
                    output[write - 1].weight += weight.weight;
                } else {
                    output[write++] = weight;
                }
            }
            output.resize(write);
            output.erase(std::remove_if(output.begin(), output.end(),
                                        [](const auto& weight) {
                                            return weight.weight == 0.0;
                                        }),
                         output.end());
        }

        void sort_and_combine(std::vector<SourceInterpolationWeight>& values) {
            std::stable_sort(values.begin(), values.end(),
                             [](const auto& left, const auto& right) {
                                 return left.source_index < right.source_index;
                             });

            std::size_t write = 0;
            for (const auto& value : values) {
                if (value.source_index < 0 || !std::isfinite(value.weight)) {
                    throw std::runtime_error(
                        "Invalid successive-orders source interpolation");
                }
                if (value.weight == 0.0) {
                    continue;
                }
                if (write != 0 &&
                    values[write - 1].source_index == value.source_index) {
                    values[write - 1].weight += value.weight;
                } else {
                    values[write++] = value;
                }
            }
            values.resize(write);
            values.erase(std::remove_if(values.begin(), values.end(),
                                        [](const auto& value) {
                                            return value.weight == 0.0;
                                        }),
                         values.end());
        }

        Eigen::Vector3d
        deterministic_horizontal(const Eigen::Vector3d& local_up,
                                 const sasktran2::Geometry& geometry) {
            constexpr double singular_squared_tolerance = 1.0e-24;
            Eigen::Vector3d horizontal =
                geometry.coordinates().reference_x() -
                local_up * local_up.dot(geometry.coordinates().reference_x());
            if (horizontal.squaredNorm() <= singular_squared_tolerance) {
                horizontal =
                    geometry.coordinates().reference_y() -
                    local_up *
                        local_up.dot(geometry.coordinates().reference_y());
            }
            if (horizontal.squaredNorm() <= singular_squared_tolerance) {
                horizontal = local_up.unitOrthogonal();
            }
            return horizontal.normalized();
        }

        Eigen::Vector3d
        horizontal_solar_reference(const Eigen::Vector3d& local_up,
                                   const sasktran2::Geometry& geometry) {
            constexpr double singular_squared_tolerance = 1.0e-24;
            Eigen::Vector3d horizontal =
                geometry.coordinates().sun_unit() -
                local_up * local_up.dot(geometry.coordinates().sun_unit());
            if (horizontal.squaredNorm() <= singular_squared_tolerance) {
                return deterministic_horizontal(local_up, geometry);
            }
            return horizontal.normalized();
        }

        Eigen::Vector3d
        rotate_unit_vector(const Eigen::Vector3d& vector,
                           const Eigen::Vector3d& initial_position,
                           const Eigen::Vector3d& new_position,
                           const sasktran2::Geometry& geometry) {
            const double vector_norm = vector.norm();
            const double initial_norm = initial_position.norm();
            const double new_norm = new_position.norm();
            if (!std::isfinite(vector_norm) || vector_norm == 0.0 ||
                !std::isfinite(initial_norm) || initial_norm == 0.0 ||
                !std::isfinite(new_norm) || new_norm == 0.0) {
                throw std::runtime_error(
                    "Cannot rotate a non-finite successive-orders direction "
                    "or position");
            }

            const Eigen::Vector3d initial_up = initial_position / initial_norm;
            const Eigen::Vector3d new_up = new_position / new_norm;
            const Eigen::Vector3d unit_vector = vector / vector_norm;
            const double mu =
                std::clamp(initial_up.dot(unit_vector), -1.0, 1.0);

            // Azimuth is undefined at exact zenith and nadir. Preserve the
            // radial direction directly instead of allowing a normalized zero
            // horizontal projection to introduce NaNs.
            constexpr double radial_tolerance = 1.0e-12;
            if (std::abs(mu) >= 1.0 - radial_tolerance) {
                return std::copysign(1.0, mu) * new_up;
            }

            const Eigen::Vector3d initial_horizontal =
                (unit_vector - mu * initial_up).normalized();
            const Eigen::Vector3d initial_solar =
                horizontal_solar_reference(initial_up, geometry);
            const double solar_azimuth = std::atan2(
                initial_up.cross(initial_solar).dot(initial_horizontal),
                initial_solar.dot(initial_horizontal));

            const Eigen::Vector3d new_solar =
                horizontal_solar_reference(new_up, geometry);
            const Eigen::Vector3d new_horizontal =
                Eigen::AngleAxis<double>(solar_azimuth, new_up) * new_solar;
            const double sin_zenith = std::sqrt(std::max(0.0, 1.0 - mu * mu));
            return (mu * new_up + sin_zenith * new_horizontal).normalized();
        }

        void compile_source_weights(
            const Eigen::Vector3d& direction,
            const sasktran2::Location& interpolation_location, bool ground,
            const sasktran2::Geometry& geometry,
            sasktran2::grids::SourceLocationInterpolator& location_interpolator,
            const std::vector<SourcePoint>& source_points,
            std::vector<SourceInterpolationWeight>& result,
            InterpolationScratch& scratch) {
            int num_locations = 0;
            if (ground) {
                location_interpolator.ground_interpolation_weights(
                    geometry.coordinates(), interpolation_location,
                    scratch.location, num_locations);
            } else {
                location_interpolator.interior_interpolation_weights(
                    geometry.coordinates(), interpolation_location,
                    scratch.location, num_locations);
            }

            result.clear();
            for (int location_index = 0; location_index < num_locations;
                 ++location_index) {
                const auto& location_weight = scratch.location[location_index];
                if (location_weight.first < 0 ||
                    location_weight.first >=
                        static_cast<int>(source_points.size())) {
                    throw std::runtime_error(
                        "Source-location interpolation returned an invalid "
                        "point index");
                }

                const auto& point = source_points[location_weight.first];
                Eigen::Vector3d rotated_direction;
                if (geometry.coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    rotated_direction = rotate_unit_vector(
                        direction, interpolation_location.position,
                        point.location().position, geometry);
                } else {
                    rotated_direction = direction;
                }

                int num_directions = 0;
                point.outgoing_sphere().interpolate(
                    rotated_direction, scratch.direction, num_directions);
                for (int direction_index = 0; direction_index < num_directions;
                     ++direction_index) {
                    const auto& direction_weight =
                        scratch.direction[direction_index];
                    if (direction_weight.first < 0 ||
                        direction_weight.first >= point.num_outgoing()) {
                        throw std::runtime_error(
                            "Unit-sphere interpolation returned an invalid "
                            "direction index");
                    }
                    result.push_back(
                        {point.outgoing_offset() + direction_weight.first,
                         location_weight.second * direction_weight.second, 0});
                }
            }
            sort_and_combine(result);
        }

        void compile_column_map(RayInterpolation& result,
                                std::vector<int>& columns) {
            const std::size_t num_weights =
                result.source_weights.size() + result.ground_weights.size();
            columns.clear();
            columns.reserve(num_weights);
            for (const auto& weight : result.source_weights) {
                columns.push_back(weight.source_index);
            }
            for (const auto& weight : result.ground_weights) {
                columns.push_back(weight.source_index);
            }
            std::sort(columns.begin(), columns.end());
            columns.erase(std::unique(columns.begin(), columns.end()),
                          columns.end());

            if (columns.size() > std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(
                    "Successive-orders transport row exceeds its compact "
                    "index range");
            }
            const auto assign_inner_index = [&](auto& weights) {
                for (auto& weight : weights) {
                    const auto iterator = std::lower_bound(
                        columns.begin(), columns.end(), weight.source_index);
                    if (iterator == columns.end() ||
                        *iterator != weight.source_index) {
                        throw std::logic_error(
                            "Successive-orders transport column map is "
                            "inconsistent");
                    }
                    weight.row_inner_index =
                        static_cast<std::uint32_t>(iterator - columns.begin());
                }
            };
            assign_inner_index(result.source_weights);
            assign_inner_index(result.ground_weights);
            result.transport_row_nnz =
                checked_u32(columns.size(), "Transport row size");
        }

        const Eigen::Vector3d&
        usable_layer_direction(const Eigen::Vector3d& average_direction,
                               const sasktran2::raytracing::TracedRay& ray) {
            if (average_direction.allFinite() &&
                average_direction.squaredNorm() > 0.0) {
                return average_direction;
            }
            // A ray tracer can emit a zero-length tangent segment whose
            // endpoint-derived average direction is undefined. Its transport
            // contribution is zero, but compiling deterministic topology is
            // still useful; the unrefracted viewing direction is the natural
            // limiting direction for that segment.
            return ray.observer_and_look.look_away;
        }
    } // namespace

    void compile_ray_interpolation(
        const sasktran2::raytracing::TracedRay& ray,
        const sasktran2::Geometry& geometry,
        sasktran2::grids::SourceLocationInterpolator& location_interpolator,
        const std::vector<SourcePoint>& source_points, RayInterpolation& result,
        InterpolationScratch& scratch) {
        result = {};
        result.traced_ray = &ray;
        result.ground_hit = ray.ground_is_hit;
        result.layers.resize(ray.layers.size());
        result.atmosphere_weights.reserve(ray.layers.size() * 2);
        // The common one-SZA grid uses at most two altitude locations and
        // three Lebedev interpolation directions. Multi-SZA grids can grow
        // this vector naturally when their four-location stencil is needed.
        result.source_weights.reserve(ray.layers.size() * 6);

        for (std::size_t layer_index = 0; layer_index < ray.layers.size();
             ++layer_index) {
            const auto& traced_layer = ray.layers[layer_index];
            auto& interpolation = result.layers[layer_index];

            interpolation.optical_depth_offset =
                traced_layer.grid_weight_offset;
            interpolation.optical_depth_count = traced_layer.grid_weight_count;

            sasktran2::Location midpoint;
            midpoint.position = 0.5 * (traced_layer.entrance.position +
                                       traced_layer.exit.position);
            geometry.assign_interpolation_weights(midpoint, scratch.atmosphere);
            sort_and_combine(scratch.atmosphere, scratch.compiled_atmosphere);
            interpolation.atmosphere_offset = checked_u32(
                result.atmosphere_weights.size(), "Atmosphere weight offset");
            interpolation.atmosphere_count = checked_u32(
                scratch.compiled_atmosphere.size(), "Atmosphere weight count");
            result.atmosphere_weights.insert(
                result.atmosphere_weights.end(),
                scratch.compiled_atmosphere.begin(),
                scratch.compiled_atmosphere.end());

            compile_source_weights(
                usable_layer_direction(traced_layer.average_look_away, ray),
                midpoint, false, geometry, location_interpolator, source_points,
                scratch.compiled_source, scratch);
            interpolation.source_offset = checked_u32(
                result.source_weights.size(), "Source weight offset");
            interpolation.source_count = checked_u32(
                scratch.compiled_source.size(), "Source weight count");
            result.source_weights.insert(result.source_weights.end(),
                                         scratch.compiled_source.begin(),
                                         scratch.compiled_source.end());
        }

        if (ray.ground_is_hit && !ray.layers.empty()) {
            const auto& ground_layer = ray.layers.front();
            sasktran2::Location ground_location;
            ground_location.position = ground_layer.exit.position;
            if (const auto* geometry_2d =
                    dynamic_cast<const sasktran2::Geometry2D*>(&geometry);
                geometry_2d != nullptr) {
                geometry_2d->assign_horizontal_interpolation_weights(
                    ground_location, result.ground_horizontal_weights);
            }
            compile_source_weights(
                -usable_layer_direction(ground_layer.average_look_away, ray),
                ground_location, true, geometry, location_interpolator,
                source_points, result.ground_weights, scratch);
        }
    }

    void adopt_optical_depth_storage(sasktran2::raytracing::TracedRay& ray,
                                     RayInterpolation& interpolation) {
        if (interpolation.traced_ray != &ray ||
            !interpolation.optical_depth_indices.empty() ||
            !interpolation.optical_depth_weights.empty()) {
            throw std::invalid_argument(
                "Invalid successive-orders OD stencil ownership transfer");
        }
        ray.release_optical_depth_storage(interpolation.optical_depth_indices,
                                          interpolation.optical_depth_weights);
        if (interpolation.optical_depth_indices.size() !=
            interpolation.optical_depth_weights.size()) {
            throw std::logic_error(
                "Successive-orders OD stencil arrays have different sizes");
        }
        interpolation.traced_ray = nullptr;
    }

    void compile_transport_row(RayInterpolation& interpolation,
                               std::vector<int>& columns) {
        compile_column_map(interpolation, columns);
    }

} // namespace sasktran2::successive_orders
