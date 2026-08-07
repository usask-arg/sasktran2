#include "transport_geometry.h"

#include <sasktran2/math/scattering.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace sasktran2::successive_orders {
    template <int NSTOKES>
    void TransportGeometry<NSTOKES>::initialize(
        const std::vector<sasktran2::raytracing::TracedRay>& traced_rays) {
        m_traced_rays = &traced_rays;
        m_use_compact_geometry = false;
        m_released = false;
        m_compact_rays.clear();
    }

    template <int NSTOKES>
    void TransportGeometry<NSTOKES>::begin_compact_2d(
        std::vector<sasktran2::raytracing::TracedRay>& traced_rays,
        const sasktran2::Geometry& geometry) {
        if (geometry.num_atmosphere_dimensions() != 2) {
            throw std::logic_error(
                "Compact successive-orders geometry requires Geometry2D");
        }
        m_traced_rays = &traced_rays;
        m_use_compact_geometry = false;
        m_released = false;
        m_compact_rays.clear();
        m_compact_rays.resize(traced_rays.size());
    }

    template <int NSTOKES>
    void TransportGeometry<NSTOKES>::compact_2d_ray(
        int rayidx, sasktran2::raytracing::TracedRay& ray) {
        if (m_traced_rays == nullptr || m_released) {
            throw std::logic_error(
                "Compact successive-orders geometry was not initialized");
        }
        if (!ray.is_straight) {
            throw std::invalid_argument(
                "Geometry2D successive orders requires straight diffuse "
                "rays");
        }

        auto& compact = m_compact_rays.at(rayidx);
        compact.layers.clear();
        compact.layers.reserve(ray.layers.size());

        std::vector<int> first_indices;
        std::vector<double> first_entrance;
        std::vector<double> first_exit;
        std::vector<double> first_optical_depth;
        if (!ray.layers.empty()) {
            const auto entrance = ray.entrance_weights(0);
            const auto exit = ray.exit_weights(0);
            const auto optical_depth = ray.optical_depth_weights(0);
            first_indices.assign(entrance.indices_data(),
                                 entrance.indices_data() + entrance.size());
            first_entrance.assign(entrance.weights_data(),
                                  entrance.weights_data() + entrance.size());
            first_exit.assign(exit.weights_data(),
                              exit.weights_data() + exit.size());
            first_optical_depth.assign(optical_depth.weights_data(),
                                       optical_depth.weights_data() +
                                           optical_depth.size());
        }

        for (const auto& layer : ray.layers) {
            if (layer.grid_weight_count >
                std::numeric_limits<std::uint8_t>::max()) {
                throw std::length_error(
                    "Compact 2D layer stencil exceeds 8-bit count");
            }
            compact.layers.push_back(
                {layer.layer_distance, layer.od_quad_start, layer.od_quad_end,
                 layer.od_quad_start_fraction, layer.od_quad_end_fraction,
                 layer.grid_weight_offset, layer.grid_weight_count});
        }

        ray.move_grid_weights_to(compact.weights.indices,
                                 compact.weights.entrance, compact.weights.exit,
                                 compact.weights.optical_depth);

        sasktran2::raytracing::TracedRay retained;
        retained.observer_and_look = ray.observer_and_look;
        retained.is_straight = ray.is_straight;
        retained.ground_is_hit = ray.ground_is_hit;
        retained.tangent_radius = ray.tangent_radius;
        if (!ray.layers.empty()) {
            retained.layers.push_back(ray.layers.front());
            retained.set_layer_weights(0, first_indices.data(),
                                       first_entrance.data(), first_exit.data(),
                                       first_optical_depth.data(),
                                       first_indices.size());
        }
        ray = std::move(retained);
    }

    template <int NSTOKES>
    void TransportGeometry<NSTOKES>::finalize_compact_2d() {
        if (m_traced_rays == nullptr || m_released) {
            throw std::logic_error(
                "Compact successive-orders geometry was not initialized");
        }
        m_use_compact_geometry = true;
    }

    template <int NSTOKES>
    int TransportGeometry<NSTOKES>::num_layers(int rayidx) const {
        if (m_released || m_traced_rays == nullptr) {
            throw std::logic_error(
                "Successive-orders transport geometry has been released");
        }
        return m_use_compact_geometry
                   ? static_cast<int>(m_compact_rays.at(rayidx).layers.size())
                   : static_cast<int>(m_traced_rays->at(rayidx).layers.size());
    }

    template <int NSTOKES>
    const sasktran2::raytracing::TracedLayer& TransportGeometry<NSTOKES>::layer(
        int rayidx, int layeridx,
        sasktran2::raytracing::TracedLayer& scratch) const {
        if (!m_use_compact_geometry) {
            return m_traced_rays->at(rayidx).layers.at(layeridx);
        }
        const auto& compact = m_compact_rays.at(rayidx).layers.at(layeridx);
        scratch = {};
        scratch.layer_distance = compact.layer_distance;
        scratch.curvature_factor = 1.0;
        scratch.od_quad_start = compact.od_quad_start;
        scratch.od_quad_end = compact.od_quad_end;
        scratch.od_quad_start_fraction = compact.od_quad_start_fraction;
        scratch.od_quad_end_fraction = compact.od_quad_end_fraction;
        return scratch;
    }

    template <int NSTOKES>
    sasktran2::raytracing::GridWeightStencilView
    TransportGeometry<NSTOKES>::entrance_weights(int rayidx,
                                                 int layeridx) const {
        if (!m_use_compact_geometry) {
            return m_traced_rays->at(rayidx).entrance_weights(layeridx);
        }
        const auto& layer = m_compact_rays.at(rayidx).layers.at(layeridx);
        const auto& weights = m_compact_rays.at(rayidx).weights;
        return {weights.indices.data() + layer.grid_weight_offset,
                weights.entrance.data() + layer.grid_weight_offset,
                layer.grid_weight_count};
    }

    template <int NSTOKES>
    sasktran2::raytracing::GridWeightStencilView
    TransportGeometry<NSTOKES>::exit_weights(int rayidx, int layeridx) const {
        if (!m_use_compact_geometry) {
            return m_traced_rays->at(rayidx).exit_weights(layeridx);
        }
        const auto& layer = m_compact_rays.at(rayidx).layers.at(layeridx);
        const auto& weights = m_compact_rays.at(rayidx).weights;
        return {weights.indices.data() + layer.grid_weight_offset,
                weights.exit.data() + layer.grid_weight_offset,
                layer.grid_weight_count};
    }

    template <int NSTOKES>
    sasktran2::raytracing::GridWeightStencilView
    TransportGeometry<NSTOKES>::optical_depth_weights(int rayidx,
                                                      int layeridx) const {
        if (!m_use_compact_geometry) {
            return m_traced_rays->at(rayidx).optical_depth_weights(layeridx);
        }
        const auto& layer = m_compact_rays.at(rayidx).layers.at(layeridx);
        const auto& weights = m_compact_rays.at(rayidx).weights;
        return {weights.indices.data() + layer.grid_weight_offset,
                weights.optical_depth.data() + layer.grid_weight_offset,
                layer.grid_weight_count};
    }

    template <int NSTOKES>
    RayTransportGeometry TransportGeometry<NSTOKES>::pack(
        const SInterpolator& source_interpolator,
        const sasktran2::Geometry& model_geometry) const {
        if (m_traced_rays == nullptr || m_released ||
            source_interpolator.size() != m_traced_rays->size()) {
            throw std::invalid_argument(
                "Rust transport geometry ray count mismatch");
        }

        const auto as_u32 = [](std::size_t value, const char* description) {
            if (value > std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(std::string(description) +
                                        " exceeds 32-bit packed storage");
            }
            return static_cast<std::uint32_t>(value);
        };

        RayTransportGeometry geometry;
        geometry.ray_layer_offsets.reserve(source_interpolator.size() + 1);
        geometry.ray_ground_offsets.reserve(source_interpolator.size() + 1);
        geometry.ray_transport_value_offsets.reserve(
            source_interpolator.size());
        geometry.ray_transport_row_nnz.reserve(source_interpolator.size());
        geometry.ray_layer_offsets.push_back(0);
        geometry.layer_atmosphere_offsets.push_back(0);
        geometry.layer_source_offsets.push_back(0);
        geometry.ray_ground_offsets.push_back(0);

        std::size_t total_layers = 0;
        std::size_t total_atmosphere_weights = 0;
        std::size_t total_source_weights = 0;
        std::size_t total_ground_weights = 0;
        for (std::size_t rayidx = 0; rayidx < source_interpolator.size();
             ++rayidx) {
            total_layers += num_layers(static_cast<int>(rayidx));
            total_atmosphere_weights +=
                source_interpolator[rayidx].atmosphere_weights.size();
            total_source_weights +=
                source_interpolator[rayidx].source_weights.size();
            total_ground_weights +=
                source_interpolator[rayidx].ground_source_weights.size();
        }
        geometry.layer_atmosphere_offsets.reserve(total_layers + 1);
        geometry.layer_source_offsets.reserve(total_layers + 1);
        geometry.atmosphere_indices.reserve(total_atmosphere_weights);
        geometry.optical_depth_weights.reserve(total_atmosphere_weights);
        geometry.albedo_weights.reserve(total_atmosphere_weights);
        geometry.entrance_weights.reserve(total_atmosphere_weights);
        geometry.exit_weights.reserve(total_atmosphere_weights);
        geometry.layer_distance.reserve(total_layers);
        geometry.layer_start_fraction.reserve(total_layers);
        geometry.layer_end_fraction.reserve(total_layers);
        geometry.ray_scattering_cosine.reserve(source_interpolator.size());
        geometry.ray_phase_q_factor.reserve(source_interpolator.size());
        geometry.ray_phase_u_factor.reserve(source_interpolator.size());
        geometry.source_value_inner_indices.reserve(total_source_weights);
        geometry.source_weights.reserve(total_source_weights);
        geometry.ground_value_inner_indices.reserve(total_ground_weights);
        geometry.ground_weights.reserve(total_ground_weights);

        for (std::size_t rayidx = 0; rayidx < source_interpolator.size();
             ++rayidx) {
            const auto& interpolator = source_interpolator[rayidx];
            geometry.ray_transport_value_offsets.push_back(
                interpolator.accumulation_row_offset);
            geometry.ray_transport_row_nnz.push_back(
                interpolator.accumulation_row_nnz);
            const auto& ray = m_traced_rays->at(rayidx);
            if (!ray.is_straight) {
                throw std::invalid_argument(
                    "Rust first-order forcing requires straight rays");
            }
            if (ray.layers.empty()) {
                geometry.ray_scattering_cosine.push_back(1.0);
                if constexpr (NSTOKES == 3) {
                    geometry.ray_phase_q_factor.push_back(0.0);
                    geometry.ray_phase_u_factor.push_back(0.0);
                }
            } else {
                const auto& look_away = ray.layers.front().average_look_away;
                double theta;
                double C1;
                double C2;
                double S1;
                double S2;
                int negation;
                sasktran2::math::stokes_scattering_factors(
                    -model_geometry.coordinates().sun_unit(), -look_away, theta,
                    C1, C2, S1, S2, negation);
                geometry.ray_scattering_cosine.push_back(std::cos(theta));
                if constexpr (NSTOKES == 3) {
                    const auto observer_rotation =
                        model_geometry.coordinates()
                            .stokes_standard_to_observer_z(
                                look_away,
                                ray.observer_and_look.observer.position);
                    geometry.ray_phase_q_factor.push_back(
                        C2 * observer_rotation.first -
                        S2 * observer_rotation.second);
                    geometry.ray_phase_u_factor.push_back(
                        C2 * observer_rotation.second +
                        S2 * observer_rotation.first);
                } else {
                    static_assert(NSTOKES == 1);
                }
            }

            const int ray_num_layers = num_layers(static_cast<int>(rayidx));
            if (interpolator.interior_weights.size() !=
                static_cast<std::size_t>(ray_num_layers)) {
                throw std::invalid_argument(
                    "Rust transport layer interpolation mismatch");
            }

            for (int layeridx = 0; layeridx < ray_num_layers; ++layeridx) {
                const auto entrance =
                    entrance_weights(static_cast<int>(rayidx), layeridx);
                const auto exit =
                    exit_weights(static_cast<int>(rayidx), layeridx);
                const auto optical_depth =
                    optical_depth_weights(static_cast<int>(rayidx), layeridx);
                const auto& layer_interpolator =
                    interpolator.interior_weights[layeridx];
                if (entrance.size() != exit.size() ||
                    entrance.size() != optical_depth.size() ||
                    layer_interpolator.atmosphere_count != entrance.size()) {
                    throw std::invalid_argument(
                        "Rust transport atmosphere stencil mismatch");
                }

                for (std::size_t index = 0; index < entrance.size(); ++index) {
                    const auto entrance_weight = entrance[index];
                    const auto optical_depth_weight = optical_depth[index];
                    if (entrance_weight.first != optical_depth_weight.first ||
                        entrance_weight.first < 0) {
                        throw std::invalid_argument(
                            "Rust transport atmospheric indices do not "
                            "match");
                    }
                    geometry.atmosphere_indices.push_back(
                        as_u32(static_cast<std::size_t>(entrance_weight.first),
                               "Atmosphere index"));
                    geometry.optical_depth_weights.push_back(
                        optical_depth_weight.second);
                    geometry.albedo_weights.push_back(
                        interpolator.atmosphere_weights
                            [layer_interpolator.atmosphere_offset + index]);
                    geometry.entrance_weights.push_back(entrance_weight.second);
                    geometry.exit_weights.push_back(exit[index].second);
                }
                geometry.layer_atmosphere_offsets.push_back(
                    as_u32(geometry.atmosphere_indices.size(),
                           "Atmosphere stencil offset"));

                sasktran2::raytracing::TracedLayer layer_scratch;
                const auto& traced_layer =
                    layer(static_cast<int>(rayidx), layeridx, layer_scratch);
                geometry.layer_distance.push_back(traced_layer.layer_distance);
                geometry.layer_start_fraction.push_back(
                    traced_layer.od_quad_start_fraction);
                geometry.layer_end_fraction.push_back(
                    traced_layer.od_quad_end_fraction);

                const std::size_t source_end =
                    static_cast<std::size_t>(layer_interpolator.source_offset) +
                    layer_interpolator.source_count;
                for (std::size_t index = layer_interpolator.source_offset;
                     index < source_end; ++index) {
                    geometry.source_value_inner_indices.push_back(
                        interpolator.source_accumulation_inner_indices[index]);
                    geometry.source_weights.push_back(
                        interpolator.source_weights[index]);
                }
                geometry.layer_source_offsets.push_back(
                    as_u32(geometry.source_value_inner_indices.size(),
                           "Source stencil offset"));
            }
            geometry.ray_layer_offsets.push_back(as_u32(
                geometry.layer_source_offsets.size() - 1, "Layer offset"));

            for (std::size_t index = 0;
                 index < interpolator.ground_source_weights.size(); ++index) {
                geometry.ground_value_inner_indices.push_back(
                    interpolator.ground_accumulation_inner_indices[index]);
                geometry.ground_weights.push_back(
                    interpolator.ground_source_weights[index]);
            }
            geometry.ray_ground_offsets.push_back(
                as_u32(geometry.ground_value_inner_indices.size(),
                       "Ground stencil offset"));
        }
        return geometry;
    }

    template <int NSTOKES> void TransportGeometry<NSTOKES>::release() {
        m_compact_rays.clear();
        m_compact_rays.shrink_to_fit();
        m_traced_rays = nullptr;
        m_use_compact_geometry = false;
        m_released = true;
    }

    template class TransportGeometry<1>;
    template class TransportGeometry<3>;
} // namespace sasktran2::successive_orders
