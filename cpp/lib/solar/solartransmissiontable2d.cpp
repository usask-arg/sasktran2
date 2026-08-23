#include <sasktran2/solartransmission.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#ifdef SKTRAN_RUST_SUPPORT
namespace sasktran2::solartransmission {
    namespace {
        constexpr double pi = 3.141592653589793238462643383279502884;
        constexpr double two_pi = 2.0 * pi;
        constexpr int maximum_segment_weights = 4;

        int number_of_endpoint_rows(
            const std::vector<sasktran2::raytracing::TracedRay>& rays) {
            std::size_t rows = 0;
            for (const auto& ray : rays) {
                rows += ray.layers.size() + 1;
            }
            if (rows >
                static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::length_error(
                    "Solar-table endpoint count exceeds the public index "
                    "range");
            }
            return static_cast<int>(rows);
        }

        void consolidate_weights(std::vector<std::pair<int, double>>& weights) {
            std::sort(weights.begin(), weights.end(),
                      [](const auto& left, const auto& right) {
                          return left.first < right.first;
                      });
            std::size_t output = 0;
            for (std::size_t begin = 0; begin < weights.size();) {
                const int index = weights[begin].first;
                double value = 0.0;
                std::size_t end = begin;
                while (end < weights.size() && weights[end].first == index) {
                    value += weights[end].second;
                    ++end;
                }
                if (value != 0.0) {
                    weights[output++] = {index, value};
                }
                begin = end;
            }
            weights.resize(output);
        }
    } // namespace

    class SolarTransmissionTable2D::Impl {
      private:
        struct SliceEntry {
            double cos_sza = 0.0;
            std::uint32_t impact_index = 0;
            std::uint8_t side = 0;
        };

        const Geometry2D& m_geometry;
        const sasktran2::raytracing::RustRayTracer2D& m_raytracer;
        const sasktran2::Config* m_config = nullptr;

        Eigen::Vector3d m_sun;
        Eigen::Vector3d m_azimuth_zero;
        Eigen::Vector3d m_azimuth_quarter;
        std::vector<double> m_radii;
        std::vector<double> m_impact_parameters;
        // Each azimuth plane owns its sampled characteristic slice. A
        // structured trace can occasionally coalesce a shell crossing with an
        // angular-boundary event in one plane but not its neighbors, so the
        // valid interpolation topology is not guaranteed to be identical in
        // every plane.
        std::vector<std::vector<SliceEntry>> m_azimuth_altitude_slices;
        // Union topology used as the common impact/SZA interpolation axis.
        // Missing plane-local nodes are completed with direct rays before the
        // temporary per-plane slices are released.
        std::vector<std::vector<SliceEntry>> m_altitude_slices;

        int m_num_azimuths = 0;
        int m_num_impacts = 0;
        int m_num_nodes = 0;
        bool m_compact_indices = false;
        bool m_initialized = false;

        // Characteristic rays are stored as a compact collection of
        // incremental atmosphere-cell stencils.
        std::vector<std::uint32_t> m_ray_segment_offsets;
        std::vector<std::uint8_t> m_segment_weight_counts;
        std::vector<std::uint16_t> m_segment_indices16;
        std::vector<std::uint32_t> m_segment_indices32;
        std::vector<double> m_segment_weights;
        std::vector<std::int32_t> m_segment_nodes;

        // (azimuth, impact, altitude, side), where side 0 is the first
        // sun-inward crossing and side 1 is the crossing after the tangent.
        std::vector<std::int32_t> m_node_lookup;

        std::size_t ray_index(int azimuth, int impact) const {
            return static_cast<std::size_t>(azimuth) * m_num_impacts + impact;
        }

        std::size_t lookup_index(int azimuth, int impact, int altitude,
                                 int side) const {
            return (
                ((static_cast<std::size_t>(azimuth) * m_num_impacts + impact) *
                     m_radii.size() +
                 altitude) *
                    2 +
                side);
        }

        std::size_t slice_index(int azimuth, int altitude) const {
            return static_cast<std::size_t>(azimuth) * m_radii.size() +
                   altitude;
        }

        std::vector<SliceEntry>& slice(int azimuth, int altitude) {
            return m_azimuth_altitude_slices[slice_index(azimuth, altitude)];
        }

        const std::vector<SliceEntry>& slice(int azimuth, int altitude) const {
            return m_azimuth_altitude_slices[slice_index(azimuth, altitude)];
        }

        bool is_shell_tangent(const SliceEntry& entry, int altitude) const {
            const double shell_parameter =
                m_radii[altitude] *
                (m_config->solar_refraction()
                     ? m_geometry.refractive_index()[altitude]
                     : 1.0);
            const double impact = m_impact_parameters[entry.impact_index];
            const double parameter_tolerance =
                64.0 * std::numeric_limits<double>::epsilon() *
                std::max({1.0, std::abs(shell_parameter), std::abs(impact)});
            constexpr double tangent_cos_tolerance = 1.0e-10;
            return std::abs(impact - shell_parameter) <= parameter_tolerance &&
                   std::abs(entry.cos_sza) <= tangent_cos_tolerance;
        }

        double surface_parameter() const {
            return m_radii.front() * (m_config->solar_refraction()
                                          ? m_geometry.refractive_index()[0]
                                          : 1.0);
        }

        bool is_surface_blocked_far_side(const SliceEntry& entry) const {
            if (entry.side == 0) {
                return false;
            }
            const double impact = m_impact_parameters[entry.impact_index];
            const double surface = surface_parameter();
            const double tolerance =
                64.0 * std::numeric_limits<double>::epsilon() *
                std::max({1.0, std::abs(surface), std::abs(impact)});
            return impact <= surface + tolerance;
        }

        int node(int azimuth, const SliceEntry& entry, int altitude) const {
            int result = m_node_lookup[lookup_index(
                azimuth, static_cast<int>(entry.impact_index), altitude,
                static_cast<int>(entry.side))];
            if (result >= 0) {
                return result;
            }

            if (is_shell_tangent(entry, altitude)) {
                // Exact shell tangencies have one physical node. Cartesian
                // roundoff can label it side 0 in one azimuth plane and side 1
                // in another, so use whichever label that plane recorded.
                result = m_node_lookup[lookup_index(
                    azimuth, static_cast<int>(entry.impact_index), altitude,
                    1 - static_cast<int>(entry.side))];
            }
            return result;
        }

        std::uint32_t segment_atmosphere_index(std::size_t segment,
                                               int slot) const {
            const std::size_t position =
                segment * maximum_segment_weights + slot;
            return m_compact_indices ? m_segment_indices16[position]
                                     : m_segment_indices32[position];
        }

        std::uint32_t current_segment_offset() const {
            if (m_segment_weight_counts.size() >
                std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(
                    "Solar-table segment count exceeds the compact offset "
                    "range");
            }
            return static_cast<std::uint32_t>(m_segment_weight_counts.size());
        }

        int altitude_index(double radius) const {
            const auto upper =
                std::lower_bound(m_radii.begin(), m_radii.end(), radius);
            int candidate;
            if (upper == m_radii.end()) {
                candidate = static_cast<int>(m_radii.size()) - 1;
            } else if (upper == m_radii.begin()) {
                candidate = 0;
            } else {
                const int upper_index =
                    static_cast<int>(upper - m_radii.begin());
                candidate = std::abs(radius - m_radii[upper_index - 1]) <
                                    std::abs(radius - m_radii[upper_index])
                                ? upper_index - 1
                                : upper_index;
            }
            const double tolerance = std::max(
                1.0e-3, 128.0 * std::numeric_limits<double>::epsilon() *
                            m_radii[candidate]);
            return std::abs(radius - m_radii[candidate]) <= tolerance
                       ? candidate
                       : -1;
        }

        Eigen::Vector3d impact_direction(double azimuth) const {
            return std::cos(azimuth) * m_azimuth_zero +
                   std::sin(azimuth) * m_azimuth_quarter;
        }

        double position_azimuth(const Eigen::Vector3d& position,
                                double cos_sza) const {
            Eigen::Vector3d transverse =
                position.normalized() - cos_sza * m_sun;
            if (transverse.squaredNorm() < 1.0e-28) {
                return 0.0;
            }
            transverse.normalize();
            double result = std::atan2(transverse.dot(m_azimuth_quarter),
                                       transverse.dot(m_azimuth_zero));
            if (result < 0.0) {
                result += two_pi;
            }
            return result;
        }

        void construct_axes() {
            m_sun = m_geometry.coordinates().sun_unit().normalized();
            const Eigen::Vector3d reference =
                m_geometry.coordinates().reference_z();
            m_azimuth_zero = reference - reference.dot(m_sun) * m_sun;
            if (m_azimuth_zero.squaredNorm() < 1.0e-24) {
                const Eigen::Vector3d fallback =
                    m_geometry.coordinates().reference_x();
                m_azimuth_zero = fallback - fallback.dot(m_sun) * m_sun;
            }
            m_azimuth_zero.normalize();
            m_azimuth_quarter = m_sun.cross(m_azimuth_zero).normalized();
        }

        void construct_impact_grid() {
            const bool refracted = m_config->solar_refraction();
            const auto& refractive_index = m_geometry.refractive_index();
            const double top_radius = m_radii.back();
            const double top_parameter =
                top_radius *
                (refracted ? refractive_index[m_radii.size() - 1] : 1.0);
            const double surface = surface_parameter();
            // Keep the horizon-limit characteristic unambiguously above the
            // surface. At exact tangency, roundoff can make one azimuth plane
            // stop at the ground while another continues to the far side.
            // The exact surface characteristic is still added below from the
            // altitude grid and supplies the illuminated tangent node.
            const double grazing_clearance = std::max(
                1.0e-3, 256.0 * std::numeric_limits<double>::epsilon() *
                            std::abs(surface));
            const double grazing_parameter =
                std::min(std::nextafter(top_parameter, 0.0),
                         surface + grazing_clearance);

            m_impact_parameters.clear();
            const int interior_count =
                std::max(16, 2 * m_geometry.num_horizontal_locations());
            for (int index = 0; index < interior_count; ++index) {
                const double angle =
                    0.5 * pi * index / static_cast<double>(interior_count - 1);
                m_impact_parameters.push_back(grazing_parameter *
                                              std::sin(angle));
            }
            for (int altitude = 0; altitude < m_geometry.num_altitudes();
                 ++altitude) {
                const double impact =
                    m_radii[altitude] *
                    (refracted ? refractive_index[altitude] : 1.0);
                m_impact_parameters.push_back(
                    std::min(impact, std::nextafter(top_parameter, 0.0)));
            }
            std::sort(m_impact_parameters.begin(), m_impact_parameters.end());
            const auto new_end = std::unique(
                m_impact_parameters.begin(), m_impact_parameters.end(),
                [](double left, double right) {
                    return std::abs(left - right) <=
                           64.0 * std::numeric_limits<double>::epsilon() *
                               std::max({1.0, std::abs(left), std::abs(right)});
                });
            m_impact_parameters.erase(new_end, m_impact_parameters.end());
            m_num_impacts = static_cast<int>(m_impact_parameters.size());
        }

        void reserve_characteristics() {
            const std::size_t rays =
                static_cast<std::size_t>(m_num_azimuths) * m_num_impacts;
            const std::size_t estimated_segments =
                rays * static_cast<std::size_t>(
                           3 *
                           (m_geometry.num_altitudes() +
                            m_geometry.num_horizontal_locations()) /
                           2);
            m_ray_segment_offsets.assign(rays + 1, 0);
            m_segment_weight_counts.clear();
            m_segment_nodes.clear();
            m_segment_weights.clear();
            m_segment_indices16.clear();
            m_segment_indices32.clear();
            m_segment_weight_counts.reserve(estimated_segments);
            m_segment_nodes.reserve(estimated_segments);
            m_segment_weights.reserve(estimated_segments *
                                      maximum_segment_weights);
            if (m_compact_indices) {
                m_segment_indices16.reserve(estimated_segments *
                                            maximum_segment_weights);
            } else {
                m_segment_indices32.reserve(estimated_segments *
                                            maximum_segment_weights);
            }
        }

        int allocate_node(int azimuth, int impact, int altitude, int side) {
            const auto lookup = lookup_index(azimuth, impact, altitude, side);
            if (m_node_lookup[lookup] >= 0) {
                return m_node_lookup[lookup];
            }
            if (m_num_nodes == std::numeric_limits<int>::max()) {
                throw std::length_error(
                    "Solar-table node count exceeds the public index range");
            }
            const int result = m_num_nodes++;
            m_node_lookup[lookup] = result;
            return result;
        }

        void append_segment(const sasktran2::raytracing::TracedRay& ray,
                            std::size_t layer_index, std::int32_t node_after) {
            const auto weights = ray.optical_depth_weights(layer_index);
            if (weights.size() > maximum_segment_weights) {
                throw std::logic_error(
                    "Structured 2D solar sweep encountered a non-cell "
                    "optical-depth stencil");
            }
            m_segment_weight_counts.push_back(
                static_cast<std::uint8_t>(weights.size()));
            m_segment_nodes.push_back(node_after);
            for (int slot = 0; slot < maximum_segment_weights; ++slot) {
                std::uint32_t atmosphere_index = 0;
                double weight = 0.0;
                if (slot < weights.size()) {
                    const auto value = weights[slot];
                    atmosphere_index = static_cast<std::uint32_t>(value.first);
                    weight = value.second;
                }
                if (m_compact_indices) {
                    m_segment_indices16.push_back(
                        static_cast<std::uint16_t>(atmosphere_index));
                } else {
                    m_segment_indices32.push_back(atmosphere_index);
                }
                m_segment_weights.push_back(weight);
            }
        }

        void trace_characteristic(int azimuth_index, int impact_index) {
            const double azimuth =
                two_pi * azimuth_index / static_cast<double>(m_num_azimuths);
            const double impact = m_impact_parameters[impact_index];
            const double top_radius = m_radii.back();
            const double top_refractive_index =
                m_config->solar_refraction()
                    ? m_geometry.refractive_index()[m_radii.size() - 1]
                    : 1.0;
            const double geometric_impact = impact / top_refractive_index;
            const double toa_along_sun = std::sqrt(
                std::max(0.0, top_radius * top_radius -
                                  geometric_impact * geometric_impact));
            const double observer_radius = top_radius + 1.0;
            const double observer_along_sun = std::sqrt(
                std::max(0.0, observer_radius * observer_radius -
                                  geometric_impact * geometric_impact));
            const Eigen::Vector3d transverse =
                geometric_impact * impact_direction(azimuth);
            const Eigen::Vector3d toa_position =
                transverse + toa_along_sun * m_sun;

            sasktran2::viewinggeometry::ViewingRay incoming;
            incoming.observer.position =
                transverse + observer_along_sun * m_sun;
            incoming.look_away = -m_sun;

            sasktran2::raytracing::TracedRay traced;
            if (m_config->solar_refraction()) {
                m_raytracer.trace_ray(incoming, m_geometry.refractive_index(),
                                      traced);
            } else {
                m_raytracer.trace_ray(incoming, traced);
            }

            const std::size_t ray = ray_index(azimuth_index, impact_index);
            m_ray_segment_offsets[ray] = current_segment_offset();
            std::vector<std::uint8_t> crossings(m_radii.size(), 0);

            const int top_altitude = static_cast<int>(m_radii.size()) - 1;
            allocate_node(azimuth_index, impact_index, top_altitude, 0);
            crossings[top_altitude] = 1;
            slice(azimuth_index, top_altitude)
                .push_back({toa_position.normalized().dot(m_sun),
                            static_cast<std::uint32_t>(impact_index), 0});

            for (std::size_t reverse_index = traced.layers.size();
                 reverse_index-- > 0;) {
                const auto& layer = traced.layers[reverse_index];
                const int altitude = altitude_index(layer.exit.radius());
                std::int32_t node_after = -1;
                if (altitude >= 0 && crossings[altitude] < 2) {
                    const int side = crossings[altitude]++;
                    node_after = allocate_node(azimuth_index, impact_index,
                                               altitude, side);
                    slice(azimuth_index, altitude)
                        .push_back({layer.exit.position.normalized().dot(m_sun),
                                    static_cast<std::uint32_t>(impact_index),
                                    static_cast<std::uint8_t>(side)});
                }
                append_segment(traced, reverse_index, node_after);
            }
            m_ray_segment_offsets[ray + 1] = current_segment_offset();
        }

        void finalize_slices() {
            m_altitude_slices.assign(m_radii.size(), {});
            std::vector<std::uint8_t> seen(
                static_cast<std::size_t>(m_num_impacts) * 2, 0);
            for (auto& current_slice : m_azimuth_altitude_slices) {
                std::sort(current_slice.begin(), current_slice.end(),
                          [](const SliceEntry& left, const SliceEntry& right) {
                              if (left.cos_sza != right.cos_sza) {
                                  return left.cos_sza < right.cos_sza;
                              }
                              if (left.impact_index != right.impact_index) {
                                  return left.impact_index < right.impact_index;
                              }
                              return left.side < right.side;
                          });
                if (current_slice.empty()) {
                    throw std::runtime_error(
                        "Solar characteristic table has an empty altitude "
                        "slice");
                }
            }

            for (int altitude = 0; altitude < static_cast<int>(m_radii.size());
                 ++altitude) {
                std::fill(seen.begin(), seen.end(), 0);
                auto& union_slice = m_altitude_slices[altitude];
                for (int azimuth = 0; azimuth < m_num_azimuths; ++azimuth) {
                    for (auto entry : slice(azimuth, altitude)) {
                        if (is_surface_blocked_far_side(entry)) {
                            continue;
                        }
                        if (!m_config->solar_refraction()) {
                            const double ratio = std::clamp(
                                m_impact_parameters[entry.impact_index] /
                                    m_radii[altitude],
                                0.0, 1.0);
                            const double magnitude =
                                std::sqrt(std::max(0.0, 1.0 - ratio * ratio));
                            entry.cos_sza =
                                entry.side == 0 ? magnitude : -magnitude;
                        }
                        if (is_shell_tangent(entry, altitude)) {
                            // A shell tangent is one physical node even when
                            // Cartesian roundoff gives different side labels
                            // in different azimuth planes.
                            entry.side = 0;
                        }
                        const std::size_t key =
                            static_cast<std::size_t>(entry.impact_index) * 2 +
                            entry.side;
                        if (seen[key] == 0) {
                            seen[key] = 1;
                            union_slice.push_back(entry);
                        }
                    }
                }
                std::sort(union_slice.begin(), union_slice.end(),
                          [](const SliceEntry& left, const SliceEntry& right) {
                              if (left.cos_sza != right.cos_sza) {
                                  return left.cos_sza < right.cos_sza;
                              }
                              if (left.impact_index != right.impact_index) {
                                  return left.impact_index < right.impact_index;
                              }
                              return left.side < right.side;
                          });
            }
        }

        bool straight_ray_hits_ground(const Location& location) const {
            const double ground_radius = m_radii.front();
            const double projected_distance = location.position.dot(m_sun);
            double discriminant = projected_distance * projected_distance -
                                  (location.position.squaredNorm() -
                                   ground_radius * ground_radius);
            const double scale =
                std::max({1.0, projected_distance * projected_distance,
                          std::abs(location.position.squaredNorm() -
                                   ground_radius * ground_radius)});
            const double tolerance =
                128.0 * std::numeric_limits<double>::epsilon() * scale;
            if (discriminant < -tolerance) {
                return false;
            }
            discriminant = std::max(0.0, discriminant);
            const double far_intersection =
                -projected_distance + std::sqrt(discriminant);
            return far_intersection >
                   64.0 * std::numeric_limits<double>::epsilon() *
                       ground_radius;
        }

        std::array<std::pair<int, double>, 2>
        altitude_weights(double radius) const {
            if (radius <= m_radii.front()) {
                return {{{0, 1.0}, {0, 0.0}}};
            }
            const int last = static_cast<int>(m_radii.size()) - 1;
            if (radius >= m_radii.back()) {
                return {{{last, 1.0}, {last, 0.0}}};
            }
            const auto upper =
                std::upper_bound(m_radii.begin(), m_radii.end(), radius);
            const int upper_index = static_cast<int>(upper - m_radii.begin());
            const int lower_index = upper_index - 1;
            const double upper_weight =
                (radius - m_radii[lower_index]) /
                (m_radii[upper_index] - m_radii[lower_index]);
            return {{{lower_index, 1.0 - upper_weight},
                     {upper_index, upper_weight}}};
        }

        std::array<std::pair<int, double>, 2>
        azimuth_weights(double azimuth) const {
            const double scaled =
                azimuth * m_num_azimuths / static_cast<double>(two_pi);
            int lower = static_cast<int>(std::floor(scaled));
            double upper_weight = scaled - lower;
            lower %= m_num_azimuths;
            if (lower < 0) {
                lower += m_num_azimuths;
            }
            const int upper = (lower + 1) % m_num_azimuths;
            return {{{lower, 1.0 - upper_weight}, {upper, upper_weight}}};
        }

        std::array<std::pair<const SliceEntry*, double>, 2>
        slice_weights(const std::vector<SliceEntry>& slice,
                      double cos_sza) const {
            if (cos_sza <= slice.front().cos_sza) {
                return {{{&slice.front(), 1.0}, {&slice.front(), 0.0}}};
            }
            if (cos_sza >= slice.back().cos_sza) {
                return {{{&slice.back(), 1.0}, {&slice.back(), 0.0}}};
            }
            const auto upper =
                std::upper_bound(slice.begin(), slice.end(), cos_sza,
                                 [](double value, const SliceEntry& entry) {
                                     return value < entry.cos_sza;
                                 });
            const auto lower = upper - 1;
            const double width = upper->cos_sza - lower->cos_sza;
            if (std::abs(width) <=
                32.0 * std::numeric_limits<double>::epsilon()) {
                return {{{&*lower, 1.0}, {&*lower, 0.0}}};
            }
            const double upper_weight = (cos_sza - lower->cos_sza) / width;
            return {{{&*lower, 1.0 - upper_weight}, {&*upper, upper_weight}}};
        }

        double visibility_limit(
            const std::array<std::pair<int, double>, 2>& altitude) const {
            double result = 0.0;
            for (const auto& [index, weight] : altitude) {
                result += weight * m_altitude_slices[index].front().cos_sza;
            }
            return result;
        }

        Eigen::Vector3d node_propagation_direction(int azimuth_index,
                                                   const SliceEntry& entry,
                                                   int altitude_index) const {
            if (!m_config->solar_refraction()) {
                return -m_sun;
            }
            const double cos_sza = std::clamp(entry.cos_sza, -1.0, 1.0);
            const double sin_sza =
                std::sqrt(std::max(0.0, 1.0 - cos_sza * cos_sza));
            const double azimuth =
                two_pi * azimuth_index / static_cast<double>(m_num_azimuths);
            const Eigen::Vector3d transverse = impact_direction(azimuth);
            const Eigen::Vector3d radial =
                cos_sza * m_sun + sin_sza * transverse;
            const Eigen::Vector3d angular =
                -sin_sza * m_sun + cos_sza * transverse;
            const double refractive_index =
                m_geometry.refractive_index()[altitude_index];
            const double sin_alpha =
                std::clamp(m_impact_parameters[entry.impact_index] /
                               (refractive_index * m_radii[altitude_index]),
                           0.0, 1.0);
            const double cos_alpha =
                std::sqrt(std::max(0.0, 1.0 - sin_alpha * sin_alpha));
            const double radial_sign = entry.side == 0 ? -1.0 : 1.0;
            return (radial_sign * cos_alpha * radial + sin_alpha * angular)
                .normalized();
        }

        bool append_completed_node(int azimuth_index, const SliceEntry& entry,
                                   int altitude_index) {
            const auto exact_lookup = lookup_index(
                azimuth_index, static_cast<int>(entry.impact_index),
                altitude_index, static_cast<int>(entry.side));
            const int existing = node(azimuth_index, entry, altitude_index);
            if (existing >= 0) {
                // Materialize tangent-side aliases so subsequent lookup does
                // not depend on which side label survived in this plane.
                m_node_lookup[exact_lookup] = existing;
                return false;
            }

            const double azimuth =
                two_pi * azimuth_index / static_cast<double>(m_num_azimuths);
            Location location;
            location.position =
                m_geometry.coordinates().solar_coordinate_vector(
                    entry.cos_sza, azimuth,
                    m_radii[altitude_index] -
                        m_geometry.coordinates().earth_radius());

            sasktran2::viewinggeometry::ViewingRay ray_to_sun;
            ray_to_sun.observer = location;
            ray_to_sun.look_away =
                m_config->solar_refraction()
                    ? -node_propagation_direction(azimuth_index, entry,
                                                  altitude_index)
                    : m_sun;
            sasktran2::raytracing::TracedRay traced;
            if (m_config->solar_refraction()) {
                m_raytracer.trace_ray_optical_depth(
                    ray_to_sun, m_geometry.refractive_index(), traced);
            } else {
                m_raytracer.trace_ray_optical_depth(ray_to_sun, traced);
            }
            // Keep the traced layers even if the tracer labels this direct
            // completion as a ground hit. Entries on the union interpolation
            // topology are known to exist in another azimuth plane; the
            // remaining discrepancy is precisely the plane-local tangent or
            // boundary classification that this completion resolves.

            const int completed_node = allocate_node(
                azimuth_index, static_cast<int>(entry.impact_index),
                altitude_index, static_cast<int>(entry.side));
            for (std::size_t layer = 0; layer < traced.layers.size(); ++layer) {
                append_segment(
                    traced, layer,
                    layer + 1 == traced.layers.size() ? completed_node : -1);
            }
            if (!traced.layers.empty()) {
                m_ray_segment_offsets.push_back(current_segment_offset());
            }
            return true;
        }

        int complete_azimuth_topology() {
            int completed = 0;
            for (int altitude = 0; altitude < static_cast<int>(m_radii.size());
                 ++altitude) {
                for (const auto& entry : m_altitude_slices[altitude]) {
                    for (int azimuth = 0; azimuth < m_num_azimuths; ++azimuth) {
                        completed +=
                            append_completed_node(azimuth, entry, altitude) ? 1
                                                                            : 0;
                    }
                }
            }
            for (int altitude = 0; altitude < static_cast<int>(m_radii.size());
                 ++altitude) {
                for (const auto& entry : m_altitude_slices[altitude]) {
                    for (int azimuth = 0; azimuth < m_num_azimuths; ++azimuth) {
                        if (node(azimuth, entry, altitude) < 0) {
                            throw std::logic_error(
                                "Solar characteristic table could not "
                                "complete its azimuth topology");
                        }
                    }
                }
            }
            m_azimuth_altitude_slices.clear();
            m_azimuth_altitude_slices.shrink_to_fit();
            return completed;
        }

      public:
        Impl(const Geometry2D& geometry,
             const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : m_geometry(geometry), m_raytracer(raytracer) {}

        void initialize_config(const sasktran2::Config& config) {
            if (m_initialized && m_config != nullptr &&
                m_config->solar_refraction() != config.solar_refraction()) {
                throw std::logic_error(
                    "A shared Geometry2D solar table cannot be reconfigured "
                    "after geometry initialization");
            }
            m_config = &config;
        }

        void initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>&) {
            if (m_config == nullptr) {
                throw std::logic_error(
                    "SolarTransmissionTable2D requires configuration before "
                    "geometry initialization");
            }
            if (m_initialized) {
                return;
            }
            construct_axes();
            const auto& altitudes = m_geometry.altitude_grid().grid();
            m_radii.resize(static_cast<std::size_t>(altitudes.size()));
            for (Eigen::Index index = 0; index < altitudes.size(); ++index) {
                m_radii[index] =
                    m_geometry.coordinates().earth_radius() + altitudes[index];
            }
            // The off-plane axis still scales with horizontal atmospheric
            // resolution, but half-resolution is a substantially better
            // setup/memory tradeoff for the smoothly interpolated OD field.
            m_num_azimuths =
                std::max(8, (m_geometry.num_horizontal_locations() + 1) / 2);
            construct_impact_grid();
            m_compact_indices =
                m_geometry.size() <=
                static_cast<int>(std::numeric_limits<std::uint16_t>::max()) + 1;
            m_num_nodes = 0;
            m_azimuth_altitude_slices.assign(
                static_cast<std::size_t>(m_num_azimuths) * m_radii.size(), {});
            m_node_lookup.assign(static_cast<std::size_t>(m_num_azimuths) *
                                     m_num_impacts * m_radii.size() * 2,
                                 -1);
            reserve_characteristics();

            for (int azimuth = 0; azimuth < m_num_azimuths; ++azimuth) {
                for (int impact = 0; impact < m_num_impacts; ++impact) {
                    trace_characteristic(azimuth, impact);
                }
            }
            finalize_slices();
            const int completed_nodes = complete_azimuth_topology();
            m_initialized = true;

            spdlog::debug(
                "Geometry2D solar characteristic table: {} azimuths, {} "
                "impact rays, {} nodes ({} completed), {} segments, compact "
                "indices {}",
                m_num_azimuths, m_num_impacts, m_num_nodes, completed_nodes,
                m_segment_weight_counts.size(), m_compact_indices);
        }

        void interpolation_weights(
            const Location& location,
            std::vector<std::pair<int, double>>& result, bool& ground_hit,
            Eigen::Vector3d* solar_propagation_direction = nullptr) const {
            result.clear();
            const double radius = location.radius();
            const Eigen::Vector3d radial = location.position / radius;
            const double cos_sza = std::clamp(radial.dot(m_sun), -1.0, 1.0);
            const auto altitude = altitude_weights(radius);

            ground_hit = m_config->solar_refraction()
                             ? cos_sza < visibility_limit(altitude) - 1.0e-12
                             : straight_ray_hits_ground(location);
            if (ground_hit) {
                if (solar_propagation_direction != nullptr) {
                    *solar_propagation_direction = -m_sun;
                }
                return;
            }

            const double azimuth = position_azimuth(location.position, cos_sza);
            const auto azimuth_interpolation = azimuth_weights(azimuth);
            result.reserve(8);
            Eigen::Vector3d direction = Eigen::Vector3d::Zero();
            for (const auto& [altitude_index, altitude_weight] : altitude) {
                if (altitude_weight == 0.0) {
                    continue;
                }
                const auto characteristic =
                    slice_weights(m_altitude_slices[altitude_index], cos_sza);
                for (const auto& [azimuth_index, azimuth_weight] :
                     azimuth_interpolation) {
                    if (azimuth_weight == 0.0) {
                        continue;
                    }
                    for (const auto& [entry, characteristic_weight] :
                         characteristic) {
                        const double weight = altitude_weight * azimuth_weight *
                                              characteristic_weight;
                        if (weight == 0.0) {
                            continue;
                        }
                        const int table_node =
                            node(azimuth_index, *entry, altitude_index);
                        if (table_node < 0) {
                            spdlog::error(
                                "Missing solar table node: azimuth {}, impact "
                                "{}, altitude {}, side {}, cos_sza {}",
                                azimuth_index, entry->impact_index,
                                altitude_index, entry->side, entry->cos_sza);
                            throw std::logic_error(
                                "Solar characteristic table topology differs "
                                "between azimuth planes");
                        }
                        result.emplace_back(table_node, weight);
                        if (solar_propagation_direction != nullptr) {
                            direction += weight * node_propagation_direction(
                                                      azimuth_index, *entry,
                                                      altitude_index);
                        }
                    }
                }
            }
            consolidate_weights(result);
            if (solar_propagation_direction != nullptr) {
                *solar_propagation_direction = direction.squaredNorm() > 1.0e-28
                                                   ? direction.normalized()
                                                   : -m_sun;
            }
        }

        void generate_interpolation(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            SolarTableInterpolation& interpolator,
            std::vector<bool>& ground_hit_flag,
            std::vector<Eigen::Vector3d>* solar_propagation_directions) const {
            const int rows = number_of_endpoint_rows(rays);
            interpolator.initialize(rows, m_num_nodes,
                                    static_cast<Eigen::Index>(rows) * 8);
            ground_hit_flag.assign(static_cast<std::size_t>(rows), false);
            if (solar_propagation_directions != nullptr) {
                solar_propagation_directions->assign(
                    static_cast<std::size_t>(rows), -m_sun);
            }
            std::vector<std::pair<int, double>> weights;
            int row = 0;
            for (const auto& ray : rays) {
                if (ray.layers.empty()) {
                    weights.clear();
                    interpolator.append_row(weights);
                    ++row;
                    continue;
                }
                for (int layer = 0; layer < ray.layers.size(); ++layer) {
                    if (layer == 0) {
                        bool ground_hit = false;
                        interpolation_weights(
                            ray.layers[layer].exit, weights, ground_hit,
                            solar_propagation_directions == nullptr
                                ? nullptr
                                : &(*solar_propagation_directions)[row]);
                        ground_hit_flag[row] = ground_hit;
                        interpolator.append_row(weights);
                        ++row;
                    }
                    bool ground_hit = false;
                    interpolation_weights(
                        ray.layers[layer].entrance, weights, ground_hit,
                        solar_propagation_directions == nullptr
                            ? nullptr
                            : &(*solar_propagation_directions)[row]);
                    ground_hit_flag[row] = ground_hit;
                    interpolator.append_row(weights);
                    ++row;
                }
            }
            interpolator.finalize();
            spdlog::debug(
                "Geometry2D solar table interpolation: {} rows, {} entries, "
                "{} ground-blocked rows",
                interpolator.rows(), interpolator.non_zeros(),
                std::count(ground_hit_flag.begin(), ground_hit_flag.end(),
                           true));
        }

        void generate_solar_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            std::vector<bool>& ground_hit_flag,
            std::vector<Eigen::Vector3d>& solar_propagation_directions) const {
            const int rows = number_of_endpoint_rows(rays);
            ground_hit_flag.assign(static_cast<std::size_t>(rows), false);
            solar_propagation_directions.assign(static_cast<std::size_t>(rows),
                                                -m_sun);
            std::vector<std::pair<int, double>> unused_weights;
            int row = 0;
            for (const auto& ray : rays) {
                if (ray.layers.empty()) {
                    ++row;
                    continue;
                }
                for (int layer = 0; layer < ray.layers.size(); ++layer) {
                    if (layer == 0) {
                        bool ground_hit = false;
                        interpolation_weights(
                            ray.layers[layer].exit, unused_weights, ground_hit,
                            &solar_propagation_directions[row]);
                        ground_hit_flag[row] = ground_hit;
                        ++row;
                    }
                    bool ground_hit = false;
                    interpolation_weights(ray.layers[layer].entrance,
                                          unused_weights, ground_hit,
                                          &solar_propagation_directions[row]);
                    ground_hit_flag[row] = ground_hit;
                    ++row;
                }
            }
        }

        int table_size() const { return m_num_nodes; }
        int atmosphere_size() const { return m_geometry.size(); }

        void apply(Eigen::Ref<const Eigen::VectorXd> extinction,
                   Eigen::Ref<Eigen::VectorXd> table_od) const {
            if (extinction.size() != atmosphere_size() ||
                table_od.size() != table_size()) {
                throw std::invalid_argument(
                    "Invalid Geometry2D solar characteristic product "
                    "dimensions");
            }
            table_od.setZero();
            const std::size_t num_rays = m_ray_segment_offsets.size() - 1;
            for (std::size_t ray = 0; ray < num_rays; ++ray) {
                double cumulative = 0.0;
                const std::size_t begin = m_ray_segment_offsets[ray];
                const std::size_t end = m_ray_segment_offsets[ray + 1];
                for (std::size_t segment = begin; segment < end; ++segment) {
                    const int count = m_segment_weight_counts[segment];
                    const std::size_t weight_offset =
                        segment * maximum_segment_weights;
                    for (int slot = 0; slot < count; ++slot) {
                        cumulative +=
                            m_segment_weights[weight_offset + slot] *
                            extinction(segment_atmosphere_index(segment, slot));
                    }
                    const int node_after = m_segment_nodes[segment];
                    if (node_after >= 0) {
                        table_od(node_after) = cumulative;
                    }
                }
            }
        }

        void
        accumulate_transpose(Eigen::Ref<const Eigen::VectorXd> table_cotangent,
                             Eigen::Ref<Eigen::VectorXd> extinction_cotangent,
                             double scale) const {
            if (table_cotangent.size() != table_size() ||
                extinction_cotangent.size() != atmosphere_size()) {
                throw std::invalid_argument(
                    "Invalid Geometry2D solar characteristic transpose "
                    "dimensions");
            }
            const std::size_t num_rays = m_ray_segment_offsets.size() - 1;
            for (std::size_t ray = 0; ray < num_rays; ++ray) {
                double cumulative = 0.0;
                const std::size_t begin = m_ray_segment_offsets[ray];
                const std::size_t end = m_ray_segment_offsets[ray + 1];
                for (std::size_t segment = end; segment-- > begin;) {
                    const int node_after = m_segment_nodes[segment];
                    if (node_after >= 0) {
                        cumulative += table_cotangent(node_after);
                    }
                    if (cumulative == 0.0) {
                        continue;
                    }
                    const int count = m_segment_weight_counts[segment];
                    const std::size_t weight_offset =
                        segment * maximum_segment_weights;
                    for (int slot = 0; slot < count; ++slot) {
                        extinction_cotangent(
                            segment_atmosphere_index(segment, slot)) +=
                            scale * cumulative *
                            m_segment_weights[weight_offset + slot];
                    }
                }
            }
        }

        std::size_t storage_bytes() const {
            std::size_t result =
                m_radii.capacity() * sizeof(double) +
                m_impact_parameters.capacity() * sizeof(double) +
                m_ray_segment_offsets.capacity() * sizeof(std::uint32_t) +
                m_segment_weight_counts.capacity() * sizeof(std::uint8_t) +
                m_segment_indices16.capacity() * sizeof(std::uint16_t) +
                m_segment_indices32.capacity() * sizeof(std::uint32_t) +
                m_segment_weights.capacity() * sizeof(double) +
                m_segment_nodes.capacity() * sizeof(std::int32_t) +
                m_node_lookup.capacity() * sizeof(std::int32_t);
            for (const auto& current_slice : m_azimuth_altitude_slices) {
                result += current_slice.capacity() * sizeof(SliceEntry);
            }
            for (const auto& current_slice : m_altitude_slices) {
                result += current_slice.capacity() * sizeof(SliceEntry);
            }
            return result;
        }
    };

    SolarTransmissionTable2D::SolarTransmissionTable2D(
        const Geometry2D& geometry,
        const sasktran2::raytracing::RustRayTracer2D& raytracer)
        : m_impl(std::make_unique<Impl>(geometry, raytracer)) {}

    SolarTransmissionTable2D::~SolarTransmissionTable2D() = default;

    void SolarTransmissionTable2D::initialize_config(
        const sasktran2::Config& config) {
        m_impl->initialize_config(config);
    }

    void SolarTransmissionTable2D::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& integration_rays) {
        m_impl->initialize_geometry(integration_rays);
    }

    void SolarTransmissionTable2D::generate_interpolation(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        SolarTableInterpolation& interpolator,
        std::vector<bool>& ground_hit_flag,
        std::vector<Eigen::Vector3d>* solar_propagation_directions) const {
        m_impl->generate_interpolation(rays, interpolator, ground_hit_flag,
                                       solar_propagation_directions);
    }

    void SolarTransmissionTable2D::generate_solar_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        std::vector<bool>& ground_hit_flag,
        std::vector<Eigen::Vector3d>& solar_propagation_directions) const {
        m_impl->generate_solar_geometry(rays, ground_hit_flag,
                                        solar_propagation_directions);
    }

    Eigen::Index SolarTransmissionTable2D::table_size() const {
        return m_impl->table_size();
    }

    Eigen::Index SolarTransmissionTable2D::atmosphere_size() const {
        return m_impl->atmosphere_size();
    }

    void SolarTransmissionTable2D::apply(
        Eigen::Ref<const Eigen::VectorXd> extinction,
        Eigen::Ref<Eigen::VectorXd> table_od) const {
        m_impl->apply(extinction, table_od);
    }

    void SolarTransmissionTable2D::accumulate_transpose(
        Eigen::Ref<const Eigen::VectorXd> table_cotangent,
        Eigen::Ref<Eigen::VectorXd> extinction_cotangent, double scale) const {
        m_impl->accumulate_transpose(table_cotangent, extinction_cotangent,
                                     scale);
    }

    std::size_t SolarTransmissionTable2D::storage_bytes() const {
        return m_impl->storage_bytes();
    }
} // namespace sasktran2::solartransmission
#endif
