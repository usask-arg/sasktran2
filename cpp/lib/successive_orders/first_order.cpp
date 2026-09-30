#include "first_order.h"
#include "first_order_scratch.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#ifdef SKTRAN_OPENMP_SUPPORT
#include <omp.h>
#endif

namespace sasktran2::successive_orders {
    namespace {
        constexpr double minimum_layer_distance_m = 1.0e-4;
        constexpr double inverse_four_pi = 1.0 / (4.0 * EIGEN_PI);

        void profile_product_scratch(const void* provider,
                                     const FirstOrderScratchLease& scratch,
                                     const char* operation, int wavelength,
                                     int wavelength_thread) {
            if (std::getenv("SASKTRAN2_PROFILE_MEMORY") == nullptr) {
                return;
            }
            std::fprintf(
                stderr,
                "SASKTRAN2_MEMORY {\"kind\":\"first_order_ephemeral_arena\","
                "\"provider\":\"%p\",\"operation\":\"%s\","
                "\"wavelength\":%d,\"worker\":%d,\"allocation_id\":%llu,"
                "\"shared\":%s,\"solar_product_bytes\":%zu,"
                "\"solar_table_product_bytes\":%zu,\"endpoint_jvp_bytes\":%zu,"
                "\"payload_bytes\":%zu}\n",
                provider, operation, wavelength, wavelength_thread,
                static_cast<unsigned long long>(scratch.allocation_id()),
                scratch.shared() ? "true" : "false", scratch.solar_bytes(),
                scratch.table_bytes(), scratch.endpoint_bytes(),
                scratch.storage_bytes());
        }

        using EndpointKey = std::vector<std::pair<int, std::uint64_t>>;

        std::pair<int, double>
        interpolation_pair(const InterpolationWeight& value) {
            return {value.index, value.weight()};
        }

        std::pair<int, double>
        interpolation_pair(const std::pair<int, double>& value) {
            return value;
        }

        struct EndpointKeyHash {
            std::size_t operator()(const EndpointKey& key) const noexcept {
                std::size_t result = key.size();
                for (const auto& [index, weight] : key) {
                    const std::size_t index_hash = std::hash<int>{}(index);
                    const std::size_t weight_hash =
                        std::hash<std::uint64_t>{}(weight);
                    result ^= index_hash + 0x9e3779b9U + (result << 6U) +
                              (result >> 2U);
                    result ^= weight_hash + 0x9e3779b9U + (result << 6U) +
                              (result >> 2U);
                }
                return result;
            }
        };

        EndpointKey endpoint_key(
            const sasktran2::raytracing::GridWeightStencilView& weights) {
            EndpointKey result;
            result.reserve(weights.size());
            for (std::size_t index = 0; index < weights.size(); ++index) {
                const auto [atmosphere_index, weight] =
                    interpolation_pair(weights[index]);
                std::uint64_t bits;
                static_assert(sizeof(bits) == sizeof(weight));
                std::memcpy(&bits, &weight, sizeof(bits));
                result.emplace_back(atmosphere_index, bits);
            }
            return result;
        }

        struct ScalarLayerTransfer {
            double attenuation;
            double source_factor;
        };

        ScalarLayerTransfer scalar_layer_transfer(double optical_depth) {
            const double absolute_depth = std::abs(optical_depth);
            if (absolute_depth < 1.0e-12) {
                return {std::exp(-optical_depth), 1.0};
            }
            if (absolute_depth < 0.02) {
                const double source_factor =
                    1.0 +
                    optical_depth *
                        (-0.5 +
                         optical_depth *
                             (1.0 / 6.0 +
                              optical_depth * (-1.0 / 24.0 +
                                               optical_depth * (1.0 / 120.0))));
                return {1.0 - optical_depth * source_factor, source_factor};
            }
            const double transmission_loss = -std::expm1(-optical_depth);
            return {1.0 - transmission_loss, transmission_loss / optical_depth};
        }

        double constant_source_factor_derivative(double optical_depth,
                                                 double factor) {
            const double absolute_depth = std::abs(optical_depth);
            if (absolute_depth < 1.0e-12) {
                return -0.5;
            }
            if (absolute_depth < 0.02) {
                return -0.5 +
                       optical_depth *
                           (1.0 / 3.0 +
                            optical_depth *
                                (-1.0 / 8.0 + optical_depth * (1.0 / 30.0)));
            }
            return 1.0 / optical_depth - factor * (1.0 + 1.0 / optical_depth);
        }

        template <int N>
        EIGEN_STRONG_INLINE double fixed_dot(const double* left,
                                             const double* right) {
            double result = 0.0;
            for (int index = 0; index < N; ++index) {
                result += left[index] * right[index];
            }
            return result;
        }

        EIGEN_STRONG_INLINE double phase_dot(const double* coefficients,
                                             const double* basis, int order) {
            if (order == 3) {
                return fixed_dot<3>(coefficients, basis);
            }
            if (order == 16) {
                return fixed_dot<16>(coefficients, basis);
            }
            double result = 0.0;
            for (int degree = 0; degree < order; ++degree) {
                result += coefficients[degree] * basis[degree];
            }
            return result;
        }

        template <int N>
        EIGEN_STRONG_INLINE void fixed_axpy(double scale, const double* input,
                                            double* output) {
            for (int index = 0; index < N; ++index) {
                output[index] += scale * input[index];
            }
        }

        EIGEN_STRONG_INLINE void phase_axpy(double scale, const double* basis,
                                            int order, double* gradient) {
            if (order == 16) {
                fixed_axpy<16>(scale, basis, gradient);
                return;
            }
            for (int degree = 0; degree < order; ++degree) {
                gradient[degree] += scale * basis[degree];
            }
        }

        void append_legendre_basis(double cosine, int count,
                                   std::vector<double>& output) {
            output.push_back(1.0);
            if (count == 1) {
                return;
            }
            output.push_back(cosine);
            for (int degree = 2; degree < count; ++degree) {
                const double value =
                    ((2.0 * degree - 1.0) * cosine * output[output.size() - 1] -
                     (degree - 1.0) * output[output.size() - 2]) /
                    degree;
                output.push_back(value);
            }
        }
    } // namespace

    template <int NSTOKES>
    FirstOrderProvider<NSTOKES>::FirstOrderProvider(
        const sasktran2::Geometry1D& geometry,
        const sasktran2::raytracing::RayTracerBase& raytracer)
        : m_geometry(geometry), m_geometry_1d(&geometry),
          m_source(geometry, raytracer),
          m_solar_table(std::make_shared<
                        sasktran2::solartransmission::SolarTransmissionTable>(
              geometry, raytracer)),
          m_integrator(false), m_source_terms{&m_source} {}

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    FirstOrderProvider<NSTOKES>::FirstOrderProvider(
        const sasktran2::Geometry2D& geometry,
        const sasktran2::raytracing::RustRayTracer2D& raytracer,
        std::shared_ptr<sasktran2::solartransmission::SolarTransmissionTable2D>
            shared_solar_table)
        : m_geometry(geometry),
          m_source(geometry, raytracer, shared_solar_table),
          m_solar_table(
              shared_solar_table != nullptr
                  ? std::move(shared_solar_table)
                  : std::make_shared<
                        sasktran2::solartransmission::SolarTransmissionTable2D>(
                        geometry, raytracer)),
          m_integrator(false), m_source_terms{&m_source} {}
#endif

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::initialize_config(
        const sasktran2::Config& config) {
        m_geometry_initialized = false;
        m_atmosphere = nullptr;
        m_num_rays = 0;
        m_primal_scratch.clear();
        m_gradient_scratch.clear();
        m_vjp_radiance_scratch.clear();
        m_vjp_cotangent_scratch.clear();
        m_phase_product_scratch.clear();
        m_phase_order_scratch.clear();
        m_scalar_phase_orders.clear();
        m_uniform_phase_active.clear();
        m_uniform_phase_values.clear();
        m_uniform_albedo_active.clear();
        m_uniform_albedo_values.clear();
        m_cached_solar_transmission.clear();
        m_cached_solar_active.clear();
        m_scalar_layer_cache.clear();
        m_endpoint_medium_cache.clear();
        m_scalar_volume_cache.clear();
        m_scalar_cache_wavelength.clear();
        m_scalar_cache_index.clear();
        m_scalar_vjp_scratch.clear();
        m_source_geometry = nullptr;
        m_solar_offsets.clear();
        m_phase_basis.clear();
        m_endpoint_slots.clear();
        m_endpoint_stencils.clear();
        m_scalar_packed_rays.clear();
        m_scalar_packed_layers.clear();
        m_scalar_ground_geometry.clear();
        m_solar_interpolation.clear();
        m_solar_ground_hit.clear();
        m_solar_propagation_directions.clear();
        m_ground_horizontal_weights.clear();
        m_num_threads = config.num_threads();
        m_num_source_threads = config.num_source_threads();
        m_num_wavelength_threads = config.num_wavelength_threads();
        m_num_phase_moments = config.num_singlescatter_moments();
        m_solar_refraction = config.solar_refraction();
        m_endpoint_phase_basis = false;
        m_compact_scalar_requested =
            NSTOKES == 1 && !config.multiple_scatter_refraction() &&
            config.singlescatter_phasemode() ==
                sasktran2::Config::SingleScatterPhaseMode::from_legendre;
        m_use_compact_scalar = false;
        if (m_num_threads < 1 || m_num_source_threads < 1 ||
            m_num_wavelength_threads < 1) {
            throw std::invalid_argument(
                "Successive-orders first-order provider requires positive "
                "thread counts");
        }
        if (m_compact_scalar_requested) {
            m_solar_table->initialize_config(config);
        } else {
            m_source.initialize_config(config);
            m_source.set_wavelength_block_capacity(1);
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::initialize_geometry(
        const SourceGeometry1D& source_geometry) {
        m_geometry_initialized = false;
        m_atmosphere = nullptr;
        m_num_rays = 0;
        m_primal_scratch.clear();
        m_gradient_scratch.clear();
        m_vjp_radiance_scratch.clear();
        m_vjp_cotangent_scratch.clear();
        m_phase_product_scratch.clear();
        m_phase_order_scratch.clear();
        m_scalar_phase_orders.clear();
        m_uniform_phase_active.clear();
        m_uniform_phase_values.clear();
        m_uniform_albedo_active.clear();
        m_uniform_albedo_values.clear();
        m_cached_solar_transmission.clear();
        m_cached_solar_active.clear();
        m_scalar_layer_cache.clear();
        m_endpoint_medium_cache.clear();
        m_scalar_volume_cache.clear();
        m_scalar_cache_wavelength.clear();
        m_scalar_cache_index.clear();
        m_scalar_vjp_scratch.clear();
        m_source_geometry = nullptr;
        m_solar_offsets.clear();
        m_phase_basis.clear();
        m_endpoint_slots.clear();
        m_endpoint_stencils.clear();
        m_scalar_packed_rays.clear();
        m_scalar_packed_layers.clear();
        m_scalar_ground_geometry.clear();
        m_solar_interpolation.clear();
        m_solar_ground_hit.clear();
        m_solar_propagation_directions.clear();
        m_ground_horizontal_weights.clear();
        const auto& viewing = source_geometry.incoming_viewing_geometry();
        m_num_rays = static_cast<int>(viewing.traced_rays.size());
        m_source_geometry = &source_geometry;
        m_ground_horizontal_weights.resize(viewing.traced_rays.size());
        if (m_geometry_1d == nullptr) {
            const auto& geometry_2d =
                static_cast<const sasktran2::Geometry2D&>(m_geometry);
            for (std::size_t ray_index = 0;
                 ray_index < viewing.traced_rays.size(); ++ray_index) {
                const auto& ray = viewing.traced_rays[ray_index];
                if (ray.ground_is_hit && !ray.layers.empty()) {
                    geometry_2d.assign_horizontal_interpolation_weights(
                        ray.layers.front().exit,
                        m_ground_horizontal_weights[ray_index]);
                }
            }
        }
        m_use_compact_scalar =
            m_compact_scalar_requested &&
            std::all_of(viewing.traced_rays.begin(), viewing.traced_rays.end(),
                        [](const auto& ray) { return ray.is_straight; });
        if (m_compact_scalar_requested && !m_use_compact_scalar) {
            throw std::logic_error(
                "Compact scalar successive-orders geometry unexpectedly "
                "contains a refracted incoming ray");
        }
        m_use_lower_interpolation =
            m_geometry_1d != nullptr &&
            m_geometry_1d->altitude_grid().interpolation_method() ==
                sasktran2::grids::interpolation::lower;
        if (m_use_compact_scalar) {
            m_solar_table->initialize_geometry(viewing.traced_rays);
            m_solar_table->generate_interpolation(
                viewing.traced_rays, m_solar_interpolation, m_solar_ground_hit,
                m_solar_refraction ? &m_solar_propagation_directions : nullptr);
        } else {
            m_source.initialize_geometry(viewing);
        }
        if (!m_use_compact_scalar) {
            m_integrator.initialize_geometry(viewing.traced_rays, m_geometry);
            m_integrator.initialize_thread_storage(m_num_threads, 1);
            m_integrator.initialize_vjp_storage(m_num_threads, 1);
        }
        m_primal_scratch.resize(m_num_threads);
        m_gradient_scratch.resize(m_num_threads);
        m_vjp_radiance_scratch.resize(m_num_threads);
        m_vjp_cotangent_scratch.resize(m_num_threads);
        m_phase_product_scratch.resize(m_num_threads);
        m_phase_order_scratch.resize(m_num_threads);
        m_scalar_vjp_scratch.resize(m_num_threads);
        for (auto& scratch : m_primal_scratch) {
            scratch.resize(1, 0, true);
        }
        for (auto& scratch : m_vjp_radiance_scratch) {
            scratch.resize(NSTOKES, 1);
        }
        for (auto& scratch : m_vjp_cotangent_scratch) {
            scratch.resize(NSTOKES, 1);
        }
        if (m_use_compact_scalar) {
            m_solar_offsets.resize(static_cast<std::size_t>(m_num_rays) + 1);
            int solar_offset = 0;
            int maximum_layers = 0;
            for (int ray = 0; ray < m_num_rays; ++ray) {
                const auto& traced_ray = viewing.traced_rays[ray];
                m_solar_offsets[ray] = solar_offset;
                solar_offset += static_cast<int>(traced_ray.layers.size()) + 1;
                maximum_layers = std::max(
                    maximum_layers, static_cast<int>(traced_ray.layers.size()));
            }
            m_solar_offsets[m_num_rays] = solar_offset;
            m_endpoint_phase_basis = m_solar_refraction;
            if (m_endpoint_phase_basis &&
                m_solar_propagation_directions.size() !=
                    static_cast<std::size_t>(solar_offset)) {
                throw std::logic_error(
                    "Refracted successive-orders solar directions do not "
                    "match the source endpoints");
            }
            m_phase_basis.reserve(
                static_cast<std::size_t>(num_phase_basis_slots()) *
                m_num_phase_moments);
            const auto& sun = m_geometry.coordinates().sun_unit();
            for (int ray = 0; ray < m_num_rays; ++ray) {
                const auto& traced_ray = viewing.traced_rays[ray];
                if (!m_endpoint_phase_basis) {
                    const double cosine =
                        traced_ray.layers.empty()
                            ? 1.0
                            : std::clamp(sun.dot(traced_ray.layers.front()
                                                     .average_look_away),
                                         -1.0, 1.0);
                    append_legendre_basis(cosine, m_num_phase_moments,
                                          m_phase_basis);
                    continue;
                }
                if (traced_ray.layers.empty()) {
                    append_legendre_basis(1.0, m_num_phase_moments,
                                          m_phase_basis);
                    continue;
                }
                const int offset = m_solar_offsets[ray];
                const auto append_endpoint_basis =
                    [&](int solar_index, const Eigen::Vector3d& look_away) {
                        const double cosine = std::clamp(
                            -m_solar_propagation_directions[solar_index].dot(
                                look_away),
                            -1.0, 1.0);
                        append_legendre_basis(cosine, m_num_phase_moments,
                                              m_phase_basis);
                    };
                append_endpoint_basis(
                    offset, traced_ray.layers.front().average_look_away);
                for (int layer = 0; layer < traced_ray.layers.size(); ++layer) {
                    append_endpoint_basis(
                        offset + layer + 1,
                        traced_ray.layers[layer].average_look_away);
                }
            }
            std::vector<int> unique_endpoint_offsets;
            std::vector<InterpolationWeight> unique_endpoint_weights;
            if (!m_use_lower_interpolation) {
                m_endpoint_slots.assign(solar_offset, -1);
                unique_endpoint_offsets.push_back(0);
                std::unordered_map<EndpointKey, int, EndpointKeyHash> slots;
                const auto assign_endpoint =
                    [&](int solar_index,
                        const sasktran2::raytracing::GridWeightStencilView&
                            weights) {
                        auto key = endpoint_key(weights);
                        const auto [position, inserted] = slots.emplace(
                            std::move(key), static_cast<int>(slots.size()));
                        const int slot = position->second;
                        m_endpoint_slots[solar_index] = slot;
                        if (!inserted) {
                            return;
                        }
                        for (std::size_t index = 0; index < weights.size();
                             ++index) {
                            const auto [atmosphere_index, weight] =
                                weights[index];
                            unique_endpoint_weights.push_back(
                                {atmosphere_index, weight});
                        }
                        unique_endpoint_offsets.push_back(
                            static_cast<int>(unique_endpoint_weights.size()));
                    };
                for (int ray = 0; ray < m_num_rays; ++ray) {
                    const auto& traced_ray = viewing.traced_rays[ray];
                    if (traced_ray.layers.empty()) {
                        continue;
                    }
                    assign_endpoint(m_solar_offsets[ray],
                                    traced_ray.exit_weights(0));
                    for (int endpoint = 1;
                         endpoint <= static_cast<int>(traced_ray.layers.size());
                         ++endpoint) {
                        assign_endpoint(
                            m_solar_offsets[ray] + endpoint,
                            traced_ray.entrance_weights(endpoint - 1));
                    }
                }
            }
            if (!m_use_lower_interpolation) {
                m_scalar_packed_rays.resize(m_num_rays);
                m_scalar_packed_layers.reserve(
                    static_cast<std::size_t>(solar_offset - m_num_rays));
                for (int ray = 0; ray < m_num_rays; ++ray) {
                    const auto& traced_ray = viewing.traced_rays[ray];
                    auto& packed_ray = m_scalar_packed_rays[ray];
                    packed_ray.layer_begin = static_cast<std::uint32_t>(
                        m_scalar_packed_layers.size());
                    for (int layer = 0;
                         layer < static_cast<int>(traced_ray.layers.size());
                         ++layer) {
                        const auto& traced_layer = traced_ray.layers[layer];
                        m_scalar_packed_layers.push_back(
                            {traced_layer.od_quad_start,
                             traced_layer.od_quad_end});
                    }
                    packed_ray.layer_end = static_cast<std::uint32_t>(
                        m_scalar_packed_layers.size());
                    if (traced_ray.ground_is_hit &&
                        !traced_ray.layers.empty()) {
                        const auto& ground_layer = traced_ray.layers.front();
                        packed_ray.ground_geometry = static_cast<std::int32_t>(
                            m_scalar_ground_geometry.size());
                        m_scalar_ground_geometry.push_back(
                            {ground_layer.exit.position.normalized(),
                             ground_layer.average_look_away});
                    }
                }
            }
            // Construction can grow these vectors beyond their retained
            // endpoint data. Runtime views are acquired after setup, so
            // compaction preserves their order without invalidating a view.
            m_endpoint_slots.shrink_to_fit();
            unique_endpoint_offsets.shrink_to_fit();
            unique_endpoint_weights.shrink_to_fit();
            if (std::getenv("SASKTRAN2_PROFILE_MEMORY") != nullptr) {
                std::fprintf(stderr,
                             "SASKTRAN2_MEMORY "
                             "{\"kind\":\"endpoint_interpolation_weights\","
                             "\"count\":%zu,\"bytes\":%zu,"
                             "\"alignment_padding_bytes_saved\":%zu}\n",
                             unique_endpoint_weights.size(),
                             unique_endpoint_weights.size() *
                                 sizeof(InterpolationWeight),
                             unique_endpoint_weights.size() * 4);
            }

            if (std::getenv("SASKTRAN2_PROFILE_MEMORY") != nullptr) {
                std::size_t sizes[6] = {};
                std::size_t index16_stencils = 0;
                std::size_t index16_weights = 0;
                std::size_t corner_stencils = 0;
                std::size_t corner16_stencils = 0;
                int maximum_index = -1;
                int altitude_stride = 0;
#ifdef SKTRAN_RUST_SUPPORT
                if (const auto* geometry2d =
                        dynamic_cast<const sasktran2::Geometry2D*>(
                            &m_geometry)) {
                    altitude_stride = geometry2d->num_altitudes();
                }
#endif
                const std::size_t stencils =
                    unique_endpoint_offsets.empty()
                        ? 0
                        : unique_endpoint_offsets.size() - 1;
                for (std::size_t slot = 0; slot < stencils; ++slot) {
                    const int begin = unique_endpoint_offsets[slot];
                    const int end = unique_endpoint_offsets[slot + 1];
                    const int count = end - begin;
                    ++sizes[std::min(count, 5)];
                    bool index16 = true;
                    for (int index = begin; index < end; ++index) {
                        const int atmosphere_index =
                            unique_endpoint_weights[index].index;
                        maximum_index =
                            std::max(maximum_index, atmosphere_index);
                        index16 = index16 && atmosphere_index >= 0 &&
                                  atmosphere_index <=
                                      std::numeric_limits<std::uint16_t>::max();
                    }
                    if (index16) {
                        ++index16_stencils;
                        index16_weights += static_cast<std::size_t>(count);
                    }
                    if (count == 4 && altitude_stride > 0) {
                        const int base = unique_endpoint_weights[begin].index;
                        const bool corner =
                            base >= 0 &&
                            base / altitude_stride + 1 <
                                m_geometry.size() / altitude_stride &&
                            base % altitude_stride + 1 < altitude_stride &&
                            unique_endpoint_weights[begin + 1].index ==
                                base + 1 &&
                            unique_endpoint_weights[begin + 2].index ==
                                base + altitude_stride &&
                            unique_endpoint_weights[begin + 3].index ==
                                base + altitude_stride + 1;
                        corner_stencils += corner;
                        corner16_stencils +=
                            corner &&
                            base <= std::numeric_limits<std::uint16_t>::max();
                    }
                }
                std::fprintf(
                    stderr,
                    "SASKTRAN2_MEMORY "
                    "{\"kind\":\"endpoint_stencil_eligibility\","
                    "\"provider\":\"%p\",\"rays\":%d,\"stencils\":%zu,"
                    "\"weights\":%zu,\"altitude_stride\":%d,\"maximum_index\":%"
                    "d,"
                    "\"size0\":%zu,\"size1\":%zu,\"size2\":%zu,\"size3\":%zu,"
                    "\"size4\":%zu,\"size_gt4\":%zu,\"index16_stencils\":%zu,"
                    "\"index16_weights\":%zu,\"corner_stencils\":%zu,"
                    "\"corner16_stencils\":%zu,\"weights_bytes\":%zu,"
                    "\"offset_bytes\":%zu,\"slot_bytes\":%zu}\n",
                    static_cast<const void*>(this), m_num_rays, stencils,
                    unique_endpoint_weights.size(), altitude_stride,
                    maximum_index, sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4], sizes[5], index16_stencils, index16_weights,
                    corner_stencils, corner16_stencils,
                    unique_endpoint_weights.capacity() *
                        sizeof(InterpolationWeight),
                    unique_endpoint_offsets.capacity() * sizeof(int),
                    m_endpoint_slots.capacity() * sizeof(int));
            }
            int endpoint_altitude_stride = 0;
#ifdef SKTRAN_RUST_SUPPORT
            if (const auto* geometry2d =
                    dynamic_cast<const sasktran2::Geometry2D*>(&m_geometry)) {
                endpoint_altitude_stride = geometry2d->num_altitudes();
            }
#endif
            m_endpoint_stencils.assign(std::move(unique_endpoint_offsets),
                                       std::move(unique_endpoint_weights),
                                       endpoint_altitude_stride,
                                       m_geometry.size());
            if (std::getenv("SASKTRAN2_PROFILE_MEMORY") != nullptr) {
                std::fprintf(
                    stderr,
                    "SASKTRAN2_MEMORY {\"kind\":\"endpoint_stencil_storage\","
                    "\"provider\":\"%p\",\"stencils\":%zu,\"weights\":%zu,"
                    "\"base_bits\":%d,\"storage_bytes\":%zu}\n",
                    static_cast<const void*>(this), m_endpoint_stencils.size(),
                    m_endpoint_stencils.weight_count(),
                    m_endpoint_stencils.base_bits(),
                    m_endpoint_stencils.storage_bytes());
            }
            for (int thread = 0; thread < m_num_threads; ++thread) {
                auto& vjp = m_scalar_vjp_scratch[thread];
                vjp.optical_depth.resize(maximum_layers);
                vjp.attenuation.resize(maximum_layers);
                vjp.source_factor.resize(maximum_layers);
                vjp.prefix_attenuation.resize(maximum_layers);
                vjp.albedo.resize(maximum_layers);
                vjp.endpoints.resize(static_cast<std::size_t>(maximum_layers) +
                                     1);
                vjp.endpoint_cotangent.resize(maximum_layers + 1);
            }
        }
        m_geometry_initialized = true;
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
        bool volume_changed) {
        if (!m_geometry_initialized) {
            throw std::logic_error(
                "Successive-orders first-order geometry is not initialized");
        }
        m_atmosphere = nullptr;
        if (!m_use_compact_scalar) {
            // Successive orders evaluates this private exact source only via
            // its primal, JVP, and VJP hooks. Avoid the large per-layer active
            // derivative lists used exclusively by materialized Jacobians.
            m_source.initialize_atmosphere_native(atmosphere);
        }
        if (!m_use_compact_scalar) {
            m_integrator.initialize_atmosphere(atmosphere);
        }
        for (auto& gradient : m_gradient_scratch) {
            gradient.resize(atmosphere.num_deriv(), 1);
        }
        if (m_use_compact_scalar && volume_changed) {
            const int coefficient_count =
                atmosphere.num_scattering_deriv_groups() == 0
                    ? 0
                    : atmosphere.storage().total_extinction.rows() *
                          m_num_phase_moments;
            for (auto& product : m_phase_product_scratch) {
                product.resize(coefficient_count);
            }
            const int locations = atmosphere.storage().total_extinction.rows();
            for (auto& orders : m_phase_order_scratch) {
                if (coefficient_count == 0) {
                    std::vector<int>().swap(orders);
                } else {
                    orders.resize(locations);
                }
            }
            m_scalar_phase_orders.resize(static_cast<std::size_t>(locations) *
                                         atmosphere.num_wavel());
            m_uniform_phase_active.assign(atmosphere.num_wavel(), 0);
            m_uniform_phase_values.resize(atmosphere.num_wavel());
            m_uniform_albedo_active.assign(atmosphere.num_wavel(), 0);
            m_uniform_albedo_values.resize(atmosphere.num_wavel());
            const auto& storage = atmosphere.storage();
            for (int wavelength = 0; wavelength < atmosphere.num_wavel();
                 ++wavelength) {
                bool uniform_albedo = locations > 0;
                for (int location = 1; location < locations && uniform_albedo;
                     ++location) {
                    uniform_albedo = storage.ssa(location, wavelength) ==
                                     storage.ssa(0, wavelength);
                }
                if (uniform_albedo) {
                    m_uniform_albedo_active[wavelength] = 1;
                    m_uniform_albedo_values[wavelength] =
                        storage.ssa(0, wavelength);
                }
                for (int location = 0; location < locations; ++location) {
                    int order =
                        std::min(storage.max_order(location, wavelength),
                                 m_num_phase_moments);
                    while (order > 0 && storage.leg_coeff(order - 1, location,
                                                          wavelength) == 0.0) {
                        --order;
                    }
                    m_scalar_phase_orders[static_cast<std::size_t>(wavelength) *
                                              locations +
                                          location] = order;
                }
                bool uniform = locations > 0;
                for (int location = 1; location < locations && uniform;
                     ++location) {
                    for (int degree = 0; degree < m_num_phase_moments;
                         ++degree) {
                        if (storage.leg_coeff(degree, location, wavelength) !=
                            storage.leg_coeff(degree, 0, wavelength)) {
                            uniform = false;
                            break;
                        }
                    }
                }
                if (uniform) {
                    m_uniform_phase_active[wavelength] = 1;
                    const int order = m_scalar_phase_orders
                        [static_cast<std::size_t>(wavelength) * locations];
                    const double* coefficients =
                        &storage.leg_coeff(0, 0, wavelength);
                    const int phase_slots = num_phase_basis_slots();
                    auto& phase_values = m_uniform_phase_values[wavelength];
                    phase_values.resize(phase_slots);
                    for (int slot = 0; slot < phase_slots; ++slot) {
                        const double* basis = m_phase_basis.data() +
                                              static_cast<std::size_t>(slot) *
                                                  m_num_phase_moments;
                        phase_values[slot] =
                            phase_dot(coefficients, basis, order);
                    }
                } else {
                    // Nonuniform wavelengths do not use this fast-path data.
                    std::vector<double>().swap(
                        m_uniform_phase_values[wavelength]);
                }
            }
        } else if (!m_use_compact_scalar) {
            m_scalar_phase_orders.clear();
            m_uniform_phase_active.clear();
            m_uniform_phase_values.clear();
            m_uniform_albedo_active.clear();
            m_uniform_albedo_values.clear();
        }
        const bool cache_shape_changed =
            m_scalar_cache_index.size() !=
                static_cast<std::size_t>(atmosphere.num_wavel()) ||
            m_scalar_cache_wavelength.size() !=
                static_cast<std::size_t>(m_num_wavelength_threads);
        if (m_use_compact_scalar && (volume_changed || cache_shape_changed)) {
            if (cache_shape_changed) {
                m_cached_solar_transmission.clear();
                m_endpoint_medium_cache.clear();
                m_scalar_layer_cache.clear();
                m_scalar_volume_cache.clear();
                m_scalar_cache_index.assign(atmosphere.num_wavel(), -1);
                m_scalar_cache_wavelength.assign(m_num_wavelength_threads, -1);
            }
            m_cached_solar_transmission.resize(m_num_wavelength_threads);
            m_cached_solar_active.assign(m_num_wavelength_threads, 0);
            if (!m_endpoint_stencils.empty()) {
                m_endpoint_medium_cache.resize(m_num_wavelength_threads);
                for (auto& cache : m_endpoint_medium_cache) {
                    cache.active = false;
                }
            } else {
                m_endpoint_medium_cache.clear();
            }
            if (atmosphere.num_deriv() > 0) {
                m_scalar_layer_cache.resize(m_num_wavelength_threads);
                for (auto& cache : m_scalar_layer_cache) {
                    cache.active = false;
                }
            } else {
                m_scalar_layer_cache.clear();
            }
            m_scalar_volume_cache.resize(m_num_wavelength_threads);
            for (auto& cache : m_scalar_volume_cache) {
                cache.active = false;
            }
        }

        if (m_use_compact_scalar &&
            std::getenv("SASKTRAN2_PROFILE_MEMORY") != nullptr) {
            const std::size_t active_uniform = static_cast<std::size_t>(
                std::count(m_uniform_phase_active.begin(),
                           m_uniform_phase_active.end(),
                           static_cast<unsigned char>(1)));
            const std::size_t required_uniform_values =
                active_uniform *
                static_cast<std::size_t>(num_phase_basis_slots());
            std::size_t uniform_bytes = 0;
            for (const auto& values : m_uniform_phase_values) {
                uniform_bytes += values.capacity() * sizeof(double);
            }
            std::fprintf(
                stderr,
                "SASKTRAN2_MEMORY {\"kind\":\"first_order_spectral_metadata\","
                "\"provider\":\"%p\",\"instance_id\":%llu,\"revision\":%llu,"
                "\"volume_revision\":%llu,\"wavelengths\":%d,\"rays\":%d,"
                "\"phase_basis_slots\":%d,\"endpoint_phase_basis\":%s,"
                "\"uniform_wavelengths\":%zu,\"uniform_phase_allocated_bytes\":"
                "%zu,"
                "\"uniform_phase_used_bytes\":%zu,\"uniform_phase_unused_"
                "bytes\":%zu,"
                "\"phase_order_bytes\":%zu,\"phase_basis_geometry_bytes\":%zu,"
                "\"solar_interpolation_rows\":%lld,\"solar_interpolation_"
                "columns\":%lld,\"solar_interpolation_nonzeros\":%lld,"
                "\"solar_interpolation_geometry_bytes\":%zu,"
                "\"solar_table_shared_bytes\":%zu}\n",
                static_cast<const void*>(this),
                static_cast<unsigned long long>(atmosphere.instance_id()),
                static_cast<unsigned long long>(atmosphere.revision()),
                static_cast<unsigned long long>(atmosphere.volume_revision()),
                atmosphere.num_wavel(), m_num_rays, num_phase_basis_slots(),
                m_endpoint_phase_basis ? "true" : "false", active_uniform,
                uniform_bytes, required_uniform_values * sizeof(double),
                uniform_bytes - required_uniform_values * sizeof(double),
                m_scalar_phase_orders.capacity() * sizeof(int),
                m_phase_basis.capacity() * sizeof(double),
                static_cast<long long>(m_solar_interpolation.rows()),
                static_cast<long long>(m_solar_interpolation.cols()),
                static_cast<long long>(m_solar_interpolation.non_zeros()),
                m_solar_interpolation.storage_bytes(),
                m_solar_table == nullptr ? 0 : m_solar_table->storage_bytes());
        }
        m_atmosphere = &atmosphere;
    }

    template <int NSTOKES>
    void
    FirstOrderProvider<NSTOKES>::prepare_wavelength(int wavelength,
                                                    int wavelength_thread) {
        validate_ready(wavelength, wavelength_thread);
        if (!m_use_compact_scalar) {
            return;
        }
        // Only the current wavelength's map entry is written. Clearing the
        // old entry would race a different worker newly evaluating it.
        m_scalar_cache_index[wavelength] = wavelength_thread;
        if (m_scalar_cache_wavelength[wavelength_thread] != wavelength) {
            m_cached_solar_active[wavelength_thread] = 0;
            if (!m_endpoint_medium_cache.empty()) {
                m_endpoint_medium_cache[wavelength_thread].active = false;
            }
            if (!m_scalar_layer_cache.empty()) {
                m_scalar_layer_cache[wavelength_thread].active = false;
            }
            m_scalar_volume_cache[wavelength_thread].active = false;
            m_scalar_cache_wavelength[wavelength_thread] = wavelength;
        }
        const int solar_size = m_solar_offsets.back();
        m_cached_solar_transmission[wavelength_thread].resize(solar_size);
        if (!m_endpoint_medium_cache.empty()) {
            const int endpoint_count =
                static_cast<int>(m_endpoint_stencils.size());
            auto& cache = m_endpoint_medium_cache[wavelength_thread];
            cache.extinction.resize(endpoint_count);
            cache.albedo.resize(endpoint_count);
        }
        if (!m_scalar_layer_cache.empty()) {
            const int layer_count = solar_size - m_num_rays;
            auto& cache = m_scalar_layer_cache[wavelength_thread];
            cache.optical_depth.resize(layer_count);
            cache.source_factor.resize(layer_count);
        }

        if (std::getenv("SASKTRAN2_PROFILE_MEMORY") != nullptr) {
            const auto eigen_bytes = [](const auto& values) {
                return static_cast<std::size_t>(values.size()) * sizeof(double);
            };
            const std::size_t endpoint_bytes =
                m_endpoint_medium_cache.empty()
                    ? 0
                    : eigen_bytes(m_endpoint_medium_cache[wavelength_thread]
                                      .extinction) +
                          eigen_bytes(m_endpoint_medium_cache[wavelength_thread]
                                          .albedo);
            const std::size_t layer_bytes =
                m_scalar_layer_cache.empty()
                    ? 0
                    : eigen_bytes(m_scalar_layer_cache[wavelength_thread]
                                      .optical_depth) +
                          eigen_bytes(m_scalar_layer_cache[wavelength_thread]
                                          .source_factor);
            std::fprintf(
                stderr,
                "SASKTRAN2_MEMORY {\"kind\":\"first_order_worker_allocations\","
                "\"provider\":\"%p\",\"instance_id\":%llu,\"revision\":%llu,"
                "\"wavelength\":%d,\"worker\":%d,\"worker_count\":%d,"
                "\"solar_bytes\":%zu,\"endpoint_bytes\":%zu,\"layer_bytes\":%"
                "zu,"
                "\"ground_prefix_bytes\":%zu,\"endpoint_jvp_bytes\":%zu,"
                "\"solar_product_bytes\":%zu,\"solar_table_product_bytes\":%zu,"
                "\"phase_product_bytes\":%zu}\n",
                static_cast<const void*>(this),
                static_cast<unsigned long long>(m_atmosphere->instance_id()),
                static_cast<unsigned long long>(m_atmosphere->revision()),
                wavelength, wavelength_thread, m_num_wavelength_threads,
                eigen_bytes(m_cached_solar_transmission[wavelength_thread]),
                endpoint_bytes, layer_bytes,
                eigen_bytes(
                    m_scalar_volume_cache[wavelength_thread].ground_prefix),
                std::size_t{0}, std::size_t{0}, std::size_t{0},
                eigen_bytes(m_phase_product_scratch[wavelength_thread]));
        }
    }

    template <int NSTOKES>
    int FirstOrderProvider<NSTOKES>::scalar_cache_index(int wavelength) const {
        if (wavelength < 0 ||
            wavelength >= static_cast<int>(m_scalar_cache_index.size())) {
            throw std::logic_error(
                "Successive-orders scalar wavelength cache is not prepared");
        }
        const int slot = m_scalar_cache_index[wavelength];
        if (slot < 0 ||
            slot >= static_cast<int>(m_scalar_cache_wavelength.size()) ||
            m_scalar_cache_wavelength[slot] != wavelength) {
            throw std::logic_error(
                "Successive-orders scalar wavelength cache is not prepared");
        }
        return slot;
    }

    template <int NSTOKES>
    void
    FirstOrderProvider<NSTOKES>::validate_ready(int wavelength,
                                                int wavelength_thread) const {
        if (!m_geometry_initialized || m_atmosphere == nullptr ||
            wavelength < 0 || wavelength >= m_atmosphere->num_wavel() ||
            wavelength_thread < 0 ||
            wavelength_thread >= m_num_wavelength_threads) {
            throw std::invalid_argument(
                "Invalid successive-orders first-order calculation");
        }
    }

    template <int NSTOKES>
    int
    FirstOrderProvider<NSTOKES>::ray_thread_index(int wavelength_thread) const {
#ifdef SKTRAN_OPENMP_SUPPORT
        if (m_num_source_threads > 1) {
            return omp_get_thread_num() + wavelength_thread;
        }
#endif
        return wavelength_thread;
    }

    template <int NSTOKES>
    const Eigen::VectorXd&
    FirstOrderProvider<NSTOKES>::ensure_solar_transmission(
        int wavelength, int wavelength_thread) {
        const int cache_slot = scalar_cache_index(wavelength);
        auto& transmission = m_cached_solar_transmission.at(cache_slot);
        if (m_cached_solar_active.at(cache_slot) == 0) {
            FirstOrderScratchLease scratch(0, 0, m_solar_table->table_size());
            auto table = scratch.table();
            m_solar_table->apply(
                m_atmosphere->storage().total_extinction.col(wavelength),
                table);
            const double irradiance =
                m_atmosphere->storage().solar_irradiance(wavelength);
            m_solar_interpolation.apply(table, transmission);
            transmission.array() = (-transmission.array()).exp() * irradiance;
            for (int row = 0; row < m_solar_interpolation.rows(); ++row) {
                if (m_solar_ground_hit[row]) {
                    transmission(row) = 0.0;
                }
            }
            m_cached_solar_active[cache_slot] = 1;
            profile_product_scratch(this, scratch, "solar", wavelength,
                                    wavelength_thread);
        }
        return transmission;
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::ensure_endpoint_medium(int wavelength) {
        const int cache_slot = scalar_cache_index(wavelength);
        if (m_endpoint_medium_cache.empty()) {
            return;
        }
        auto& cache = m_endpoint_medium_cache[cache_slot];
        if (cache.active) {
            return;
        }
        const auto& storage = m_atmosphere->storage();
        m_endpoint_stencils.visit([&](const auto& stencils) {
            for (int slot = 0; slot < cache.extinction.size(); ++slot) {
                double extinction = 0.0;
                double albedo = 0.0;
                const auto weights = stencils[slot];
                for (std::size_t index = 0; index < weights.size(); ++index) {
                    const auto [atmosphere_index, weight] =
                        interpolation_pair(weights[index]);
                    extinction += weight * storage.total_extinction(
                                               atmosphere_index, wavelength);
                    albedo +=
                        weight * storage.ssa(atmosphere_index, wavelength);
                }
                cache.extinction(slot) = extinction;
                cache.albedo(slot) = albedo;
            }
        });
        cache.active = true;
    }

    template <int NSTOKES>
    typename FirstOrderProvider<NSTOKES>::ScalarEndpointContext
    FirstOrderProvider<NSTOKES>::prepare_scalar_endpoint_context(
        int wavelength, int cache_slot) const {
        ScalarEndpointContext context;
        const auto& storage = m_atmosphere->storage();
        const int locations = storage.total_extinction.rows();
        context.extinction = storage.total_extinction.col(wavelength).data();
        context.albedo = storage.ssa.col(wavelength).data();
        context.coefficient_stride = storage.leg_coeff.dimension(0);
        if (locations > 0) {
            context.coefficients = storage.leg_coeff.data() +
                                   static_cast<Eigen::Index>(wavelength) *
                                       storage.leg_coeff.dimension(1) *
                                       context.coefficient_stride;
            context.phase_orders =
                m_scalar_phase_orders.data() +
                static_cast<std::size_t>(wavelength) * locations;
        }
        context.maximum_orders = storage.max_order.col(wavelength).data();
        context.basis = m_phase_basis.data();
        context.endpoint_basis_stride =
            m_endpoint_phase_basis ? m_num_phase_moments : 0;
        if (m_uniform_phase_active[wavelength] != 0) {
            context.uniform_phase = m_uniform_phase_values[wavelength].data();
        }
        if (!m_endpoint_slots.empty()) {
            context.endpoint_slots = m_endpoint_slots.data();
        }
        if (!m_endpoint_medium_cache.empty()) {
            const auto& medium = m_endpoint_medium_cache[cache_slot];
            if (medium.active) {
                context.endpoint_extinction = medium.extinction.data();
                context.endpoint_albedo = medium.albedo.data();
            }
        }
        context.solar = m_cached_solar_transmission[cache_slot].data();
        context.ssa_deriv_start = m_atmosphere->ssa_deriv_start_index();
        return context;
    }

    template <int NSTOKES>
    typename FirstOrderProvider<NSTOKES>::ScalarEndpointContext
    FirstOrderProvider<NSTOKES>::scalar_endpoint_context_for_ray(
        const ScalarEndpointContext& context, int ray) const {
        auto result = context;
        if (!m_endpoint_phase_basis) {
            result.basis += static_cast<std::size_t>(ray) * m_num_phase_moments;
            if (result.uniform_phase != nullptr) {
                result.uniform_phase += ray;
            }
        }
        return result;
    }

    template <int NSTOKES>
    template <typename Weights>
    typename FirstOrderProvider<NSTOKES>::ScalarEndpoint
    FirstOrderProvider<NSTOKES>::scalar_endpoint(
        const ScalarEndpointContext& context, int layer, bool entrance,
        int solar_index, const Weights& weights) const {
        ScalarEndpoint result;
        const int endpoint_slot = context.endpoint_slots == nullptr
                                      ? -1
                                      : context.endpoint_slots[solar_index];
        const bool use_endpoint_medium =
            endpoint_slot >= 0 && context.endpoint_extinction != nullptr;
        if (use_endpoint_medium) {
            result.extinction = context.endpoint_extinction[endpoint_slot];
            result.albedo = context.endpoint_albedo[endpoint_slot];
        }
        const double* basis =
            context.basis + static_cast<std::size_t>(solar_index) *
                                context.endpoint_basis_stride;
        const bool uniform_phase = context.uniform_phase != nullptr;
        if (uniform_phase) {
            result.phase =
                context.uniform_phase[context.endpoint_basis_stride == 0
                                          ? 0
                                          : solar_index];
        }
        if (use_endpoint_medium && uniform_phase) {
            return result;
        }
        for (std::size_t index = 0; index < weights.size(); ++index) {
            const auto [atmosphere_index, weight] =
                interpolation_pair(weights[index]);
            if (weight == 0.0) {
                continue;
            }
            if (!use_endpoint_medium) {
                result.extinction +=
                    weight * context.extinction[atmosphere_index];
                result.albedo += weight * context.albedo[atmosphere_index];
            }
            if (!uniform_phase) {
                const int order = context.phase_orders[atmosphere_index];
                const double* coefficients =
                    context.coefficients +
                    atmosphere_index * context.coefficient_stride;
                const double phase = phase_dot(coefficients, basis, order);
                result.phase += weight * phase;
            }
        }
        (void)layer;
        (void)entrance;
        (void)solar_index;
        return result;
    }

    template <int NSTOKES>
    template <bool USE_ENDPOINT_MEDIUM, typename Weights>
    typename FirstOrderProvider<NSTOKES>::ScalarValueTangent
    FirstOrderProvider<NSTOKES>::scalar_endpoint_jvp(
        const ScalarEndpointContext& context, int layer, bool entrance,
        int solar_index, const Weights& weights,
        const double* __restrict extinction_direction,
        const double* __restrict albedo_direction,
        const double* __restrict solar_tangent, bool uniform_albedo_direction,
        double uniform_albedo_tangent, bool phase_tangent_active) const {
        ScalarValueTangent result;
        const double* basis =
            context.basis + static_cast<std::size_t>(solar_index) *
                                context.endpoint_basis_stride;
        double extinction = 0.0;
        double extinction_tangent = 0.0;
        double albedo = 0.0;
        double albedo_tangent =
            uniform_albedo_direction ? uniform_albedo_tangent : 0.0;
        double phase = 0.0;
        double phase_tangent = 0.0;
        const bool uniform_phase = context.uniform_phase != nullptr;
        if (uniform_phase) {
            phase = context.uniform_phase[context.endpoint_basis_stride == 0
                                              ? 0
                                              : solar_index];
        }
        if constexpr (USE_ENDPOINT_MEDIUM) {
            const int endpoint_slot = context.endpoint_slots[solar_index];
            extinction = context.endpoint_extinction[endpoint_slot];
            albedo = context.endpoint_albedo[endpoint_slot];
            extinction_tangent =
                context.endpoint_extinction_tangent[endpoint_slot];
            albedo_tangent = context.endpoint_albedo_tangent[endpoint_slot];
        }
        if constexpr (USE_ENDPOINT_MEDIUM) {
            if (uniform_phase && !phase_tangent_active) {
                const double solar = context.solar[solar_index];
                const double solar_jvp = solar_tangent[solar_index];
                const double phase_scale = phase * inverse_four_pi;
                const double extinction_albedo = extinction * albedo;
                result.value = extinction_albedo * solar * phase_scale;
                result.tangent = ((extinction_tangent * albedo +
                                   extinction * albedo_tangent) *
                                      solar +
                                  extinction_albedo * solar_jvp) *
                                 phase_scale;
                return result;
            }
        }
        const bool require_weight_interpolation =
            !USE_ENDPOINT_MEDIUM || !uniform_phase || phase_tangent_active;
        for (std::size_t index = 0;
             require_weight_interpolation && index < weights.size(); ++index) {
            const auto [atmosphere_index, weight] =
                interpolation_pair(weights[index]);
            if (weight == 0.0) {
                continue;
            }
            if constexpr (!USE_ENDPOINT_MEDIUM) {
                extinction += weight * context.extinction[atmosphere_index];
                extinction_tangent +=
                    weight * extinction_direction[atmosphere_index];
            }
            if constexpr (!USE_ENDPOINT_MEDIUM) {
                albedo += weight * context.albedo[atmosphere_index];
                if (!uniform_albedo_direction) {
                    albedo_tangent +=
                        weight * albedo_direction[atmosphere_index];
                }
            }
            if (phase_tangent_active) {
                const double* tangent_coefficients =
                    context.coefficient_tangent +
                    atmosphere_index * m_num_phase_moments;
                if (!uniform_phase) {
                    const int order = context.phase_orders[atmosphere_index];
                    const double* coefficients =
                        context.coefficients +
                        atmosphere_index * context.coefficient_stride;
                    phase += weight * phase_dot(coefficients, basis, order);
                }
                phase_tangent +=
                    weight *
                    phase_dot(tangent_coefficients, basis,
                              context.tangent_orders[atmosphere_index]);
            } else if (!uniform_phase) {
                const int order = context.phase_orders[atmosphere_index];
                const double* coefficients =
                    context.coefficients +
                    atmosphere_index * context.coefficient_stride;
                phase += weight * phase_dot(coefficients, basis, order);
            }
        }
        const double solar = context.solar[solar_index];
        const double solar_jvp = solar_tangent[solar_index];
        result.value = extinction * albedo * solar * phase * inverse_four_pi;
        result.tangent =
            (((extinction_tangent * albedo + extinction * albedo_tangent) *
                  solar +
              extinction * albedo * solar_jvp) *
                 phase +
             extinction * albedo * solar * phase_tangent) *
            inverse_four_pi;
        (void)layer;
        (void)entrance;
        return result;
    }

    template <int NSTOKES>
    template <bool WITH_PHASE_GRADIENT, typename Weights>
    void FirstOrderProvider<NSTOKES>::accumulate_scalar_endpoint_vjp(
        const ScalarEndpointContext& context, int layer, bool entrance,
        int solar_index, const Weights& weights, const ScalarEndpoint& endpoint,
        double source_cotangent, Eigen::Ref<Eigen::VectorXd> native_gradient,
        Eigen::Ref<Eigen::VectorXd> solar_gradient,
        Eigen::Ref<Eigen::VectorXd> coefficient_gradient) const {
        if (source_cotangent == 0.0) {
            return;
        }
        const double scale = source_cotangent * inverse_four_pi;
        const double extinction_cotangent = scale * endpoint.albedo *
                                            endpoint.solar_transmission *
                                            endpoint.phase;
        const double albedo_cotangent = scale * endpoint.extinction *
                                        endpoint.solar_transmission *
                                        endpoint.phase;
        const double solar_cotangent =
            scale * endpoint.extinction * endpoint.albedo * endpoint.phase;
        const double phase_cotangent = WITH_PHASE_GRADIENT
                                           ? scale * endpoint.extinction *
                                                 endpoint.albedo *
                                                 endpoint.solar_transmission
                                           : 0.0;
        const double* basis =
            context.basis + static_cast<std::size_t>(solar_index) *
                                context.endpoint_basis_stride;
        for (std::size_t index = 0; index < weights.size(); ++index) {
            const auto [atmosphere_index, weight] =
                interpolation_pair(weights[index]);
            if (weight == 0.0) {
                continue;
            }
            native_gradient(atmosphere_index) += weight * extinction_cotangent;
            native_gradient(context.ssa_deriv_start + atmosphere_index) +=
                weight * albedo_cotangent;
            if constexpr (WITH_PHASE_GRADIENT) {
                const int order =
                    std::min(context.maximum_orders[atmosphere_index],
                             m_num_phase_moments);
                const int coefficient_offset =
                    atmosphere_index * m_num_phase_moments;
                const double coefficient_cotangent = weight * phase_cotangent;
                phase_axpy(coefficient_cotangent, basis, order,
                           coefficient_gradient.data() + coefficient_offset);
            }
        }
        solar_gradient(solar_index) += solar_cotangent;
        (void)layer;
        (void)entrance;
    }

    template <int NSTOKES>
    bool FirstOrderProvider<NSTOKES>::ground_scattering_geometry(
        int solar_index, const ScalarGroundGeometry& ground, double& mu_in,
        double& mu_out, double& phi) const {
        Eigen::Vector3d direction_to_sun = m_geometry.coordinates().sun_unit();
        if (!m_solar_propagation_directions.empty()) {
            direction_to_sun = -m_solar_propagation_directions[solar_index];
        }
        sasktran2::Location ground_location;
        ground_location.position = ground.up;
        sasktran2::raytracing::calculate_csz_saz(
            direction_to_sun.normalized(), ground_location, ground.look_away,
            mu_in, phi, m_geometry.coordinates().geometry_type());
        mu_out = -ground.up.dot(ground.look_away.normalized());
        return mu_in > 0.0;
    }

    template <int NSTOKES>
    double FirstOrderProvider<NSTOKES>::ground_transport_albedo(int wavelength,
                                                                int ray) const {
        if (!m_atmosphere->surface().has_spatial_lambertian_albedo()) {
            return 1.0;
        }
        return m_atmosphere->surface().spatial_lambertian_albedo(
            wavelength, m_ground_horizontal_weights[ray]);
    }

    template <int NSTOKES>
    double FirstOrderProvider<NSTOKES>::ground_transport_albedo_tangent(
        int ray, Eigen::Ref<const Eigen::VectorXd> native_tangent) const {
        if (!m_atmosphere->surface().has_spatial_lambertian_albedo()) {
            return 0.0;
        }
        double result = 0.0;
        for (const auto& [horizontal_index, weight] :
             m_ground_horizontal_weights[ray]) {
            result += weight *
                      native_tangent(m_atmosphere->surface_deriv_start_index() +
                                     horizontal_index);
        }
        return result;
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::accumulate_ground_transport_albedo_vjp(
        int ray, double albedo_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        if (!m_atmosphere->surface().has_spatial_lambertian_albedo() ||
            albedo_cotangent == 0.0) {
            return;
        }
        for (const auto& [horizontal_index, weight] :
             m_ground_horizontal_weights[ray]) {
            native_gradient(m_atmosphere->surface_deriv_start_index() +
                            horizontal_index) += weight * albedo_cotangent;
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::calculate_scalar(
        int wavelength, int wavelength_thread,
        Eigen::Ref<Eigen::VectorXd> forcing, TransportOperator* transport) {
        prepare_wavelength(wavelength, wavelength_thread);
        if (transport != nullptr) {
            if (m_use_lower_interpolation) {
                calculate_scalar_impl<true, true>(wavelength, wavelength_thread,
                                                  forcing, transport);
            } else {
                calculate_scalar_impl<true, false>(
                    wavelength, wavelength_thread, forcing, transport);
            }
        } else if (m_use_lower_interpolation) {
            calculate_scalar_impl<false, true>(wavelength, wavelength_thread,
                                               forcing, nullptr);
        } else {
            calculate_scalar_impl<false, false>(wavelength, wavelength_thread,
                                                forcing, nullptr);
        }
    }

    template <int NSTOKES>
    template <bool WITH_TRANSPORT, bool WITH_FORCING>
    void FirstOrderProvider<NSTOKES>::calculate_scalar_uniform_impl(
        int wavelength, int wavelength_thread,
        Eigen::Ref<Eigen::VectorXd> forcing, TransportOperator* transport) {
        const int cache_slot = scalar_cache_index(wavelength);
        const auto& interpolation = m_source_geometry->incoming_interpolation();
        const auto& extinction = m_atmosphere->storage().total_extinction;
        const auto& ssa = m_atmosphere->storage().ssa;
        if constexpr (WITH_FORCING) {
            ensure_solar_transmission(wavelength, wavelength_thread);
            ensure_endpoint_medium(wavelength);
        }
        const auto& solar = m_cached_solar_transmission[cache_slot];
        const bool uniform_albedo = m_uniform_albedo_active[wavelength] != 0;
        const double uniform_albedo_value =
            uniform_albedo ? m_uniform_albedo_values[wavelength] : 0.0;
        if constexpr (WITH_TRANSPORT) {
            transport->values().setZero();
        }
        ScalarLayerCache* layer_cache = m_scalar_layer_cache.empty()
                                            ? nullptr
                                            : &m_scalar_layer_cache[cache_slot];
        if (layer_cache != nullptr) {
            layer_cache->active = false;
        }
#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const auto& packed_ray = m_scalar_packed_rays[ray];
            const auto& ray_interpolation = interpolation[ray];
            const double phase_scale =
                WITH_FORCING
                    ? m_uniform_phase_values[wavelength][ray] * inverse_four_pi
                    : 0.0;
            const auto endpoint_source = [&](int solar_index) {
                const int slot = m_endpoint_slots[solar_index];
                const auto& medium = m_endpoint_medium_cache[cache_slot];
                return medium.extinction(slot) * medium.albedo(slot) *
                       solar(solar_index) * phase_scale;
            };
            double prefix = 1.0;
            double radiance = 0.0;
            double shared_endpoint = 0.0;
            bool have_shared_endpoint = false;
            for (std::uint32_t flat_layer = packed_ray.layer_end;
                 flat_layer-- > packed_ray.layer_begin;) {
                const auto& layer = m_scalar_packed_layers[flat_layer];
                const auto local_layer = static_cast<std::size_t>(
                    flat_layer - packed_ray.layer_begin);
                const auto optical_depth_weights =
                    ray_interpolation.optical_depth_for_layer(local_layer);
                const auto atmosphere_weights =
                    ray_interpolation.atmosphere_for_layer(local_layer);
                double optical_depth = 0.0;
                double albedo = uniform_albedo_value;
                for (std::size_t index = 0;
                     index < optical_depth_weights.size(); ++index) {
                    const auto [atmosphere_index, weight] =
                        optical_depth_weights[index];
                    optical_depth +=
                        weight * extinction(atmosphere_index, wavelength);
                }
                if (!uniform_albedo) {
                    for (const auto& weight : atmosphere_weights) {
                        albedo +=
                            weight.weight() * ssa(weight.index, wavelength);
                    }
                }
                const auto transfer = scalar_layer_transfer(optical_depth);
                if (layer_cache != nullptr) {
                    layer_cache->optical_depth(flat_layer) = optical_depth;
                    layer_cache->source_factor(flat_layer) =
                        transfer.source_factor;
                }
                if constexpr (WITH_TRANSPORT) {
                    const double source_factor =
                        prefix * albedo * (1.0 - transfer.attenuation);
                    ray_interpolation.source_for_layer(local_layer)
                        .visit([&](const auto sources) {
                            for (const auto& source : sources) {
                                transport->values()(
                                    ray_interpolation.transport_value_offset +
                                    source.row_inner_index()) +=
                                    source.weight() * source_factor;
                            }
                        });
                }
                if constexpr (WITH_FORCING) {
                    const double layer_distance =
                        layer.source_quad_start + layer.source_quad_end;
                    if (layer_distance < minimum_layer_distance_m) {
                        have_shared_endpoint = false;
                    } else {
                        const int exit_solar = m_solar_offsets[ray] +
                                               static_cast<int>(local_layer);
                        const double start =
                            have_shared_endpoint
                                ? shared_endpoint
                                : endpoint_source(exit_solar + 1);
                        const double end = endpoint_source(exit_solar);
                        radiance += prefix * transfer.source_factor *
                                    (start * layer.source_quad_start +
                                     end * layer.source_quad_end);
                        shared_endpoint = end;
                        have_shared_endpoint = true;
                    }
                }
                prefix *= transfer.attenuation;
            }
            if constexpr (WITH_TRANSPORT) {
                if (ray_interpolation.ground_is_hit()) {
                    const double ground_albedo =
                        ground_transport_albedo(wavelength, ray);
                    ray_interpolation.ground_weights.visit(
                        [&](const auto sources) {
                            for (const auto& source : sources) {
                                transport->values()(
                                    ray_interpolation.transport_value_offset +
                                    source.row_inner_index()) +=
                                    source.weight() * prefix * ground_albedo;
                            }
                        });
                }
            }
            if constexpr (WITH_FORCING) {
                if (packed_ray.ground_geometry >= 0) {
                    const auto& ground =
                        m_scalar_ground_geometry[packed_ray.ground_geometry];
                    double mu_in;
                    double mu_out;
                    double phi;
                    if (ground_scattering_geometry(m_solar_offsets[ray], ground,
                                                   mu_in, mu_out, phi)) {
                        const auto brdf = m_atmosphere->surface().brdf(
                            wavelength, mu_in, mu_out, phi,
                            m_ground_horizontal_weights[ray]);
                        radiance += prefix * solar(m_solar_offsets[ray]) *
                                    mu_in * brdf(0, 0);
                    }
                }
                forcing(ray) = radiance;
            }
        }
        if (layer_cache != nullptr) {
            layer_cache->active = true;
        }
    }

    template <int NSTOKES>
    template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION, bool WITH_FORCING>
    void FirstOrderProvider<NSTOKES>::calculate_scalar_impl(
        int wavelength, int wavelength_thread,
        Eigen::Ref<Eigen::VectorXd> forcing, TransportOperator* transport) {
        const int cache_slot = scalar_cache_index(wavelength);
        if constexpr (!LOWER_INTERPOLATION) {
            if (m_uniform_phase_active[wavelength] != 0 &&
                !m_endpoint_phase_basis) {
                calculate_scalar_uniform_impl<WITH_TRANSPORT, WITH_FORCING>(
                    wavelength, wavelength_thread, forcing, transport);
                return;
            }
        }
        const auto& rays = m_source_geometry->incoming_rays();
        const auto& interpolation = m_source_geometry->incoming_interpolation();
        const auto& extinction = m_atmosphere->storage().total_extinction;
        const auto& ssa = m_atmosphere->storage().ssa;
        ScalarEndpointContext endpoint_context;
        if constexpr (WITH_FORCING) {
            ensure_solar_transmission(wavelength, wavelength_thread);
            ensure_endpoint_medium(wavelength);
            endpoint_context =
                prepare_scalar_endpoint_context(wavelength, cache_slot);
        }
        const auto& solar = m_cached_solar_transmission[cache_slot];
        const bool uniform_albedo = m_uniform_albedo_active[wavelength] != 0;
        const double uniform_albedo_value =
            uniform_albedo ? m_uniform_albedo_values[wavelength] : 0.0;
        ScalarVolumeCache* volume_cache = nullptr;
        if constexpr (WITH_TRANSPORT && !LOWER_INTERPOLATION) {
            if (!m_scalar_volume_cache.empty()) {
                volume_cache = &m_scalar_volume_cache[cache_slot];
            }
        }
        if constexpr (WITH_TRANSPORT) {
            transport->values().setZero();
        }
        if (volume_cache != nullptr) {
            volume_cache->ground_prefix.setZero(m_num_rays);
            volume_cache->active = false;
        }
        ScalarLayerCache* layer_cache = m_scalar_layer_cache.empty()
                                            ? nullptr
                                            : &m_scalar_layer_cache[cache_slot];
        if (layer_cache != nullptr) {
            layer_cache->active = false;
        }
#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            ScalarEndpointContext ray_endpoint_context;
            if constexpr (WITH_FORCING) {
                ray_endpoint_context =
                    scalar_endpoint_context_for_ray(endpoint_context, ray);
            }
            const auto& ray_interpolation = interpolation[ray];
            const ScalarPackedRay* packed_ray = nullptr;
            if constexpr (!LOWER_INTERPOLATION) {
                packed_ray = &m_scalar_packed_rays[ray];
            }
            const int flat_layer_offset = m_solar_offsets[ray] - ray;
            double prefix = 1.0;
            double radiance = 0.0;
            ScalarEndpoint shared_endpoint;
            bool have_shared_endpoint = false;
            for (int layer =
                     static_cast<int>(ray_interpolation.layers.size()) - 1;
                 layer >= 0; --layer) {
                double source_quad_start = 0.0;
                double source_quad_end = 0.0;
                double layer_distance = 0.0;
                if constexpr (WITH_FORCING) {
                    if constexpr (LOWER_INTERPOLATION) {
                        const auto& traced_layer = rays[ray].layers[layer];
                        layer_distance = traced_layer.layer_distance;
                        source_quad_start = layer_distance *
                                            traced_layer.od_quad_start_fraction;
                        source_quad_end =
                            layer_distance * traced_layer.od_quad_end_fraction;
                    } else {
                        const auto& packed_layer =
                            m_scalar_packed_layers[packed_ray->layer_begin +
                                                   layer];
                        source_quad_start = packed_layer.source_quad_start;
                        source_quad_end = packed_layer.source_quad_end;
                        layer_distance = source_quad_start + source_quad_end;
                    }
                }
                double optical_depth = 0.0;
                double albedo = uniform_albedo_value;
                const auto od_weights =
                    ray_interpolation.optical_depth_for_layer(layer);
                for (std::size_t index = 0; index < od_weights.size();
                     ++index) {
                    const auto [atmosphere_index, weight] = od_weights[index];
                    optical_depth +=
                        weight * extinction(atmosphere_index, wavelength);
                }
                if constexpr (WITH_TRANSPORT) {
                    if (!uniform_albedo) {
                        for (const auto& weight :
                             ray_interpolation.atmosphere_for_layer(layer)) {
                            albedo +=
                                weight.weight() * ssa(weight.index, wavelength);
                        }
                    }
                }
                const auto layer_transfer =
                    scalar_layer_transfer(optical_depth);
                const double attenuation = layer_transfer.attenuation;
                const double factor = layer_transfer.source_factor;
                if (layer_cache != nullptr) {
                    const int flat_layer = flat_layer_offset + layer;
                    layer_cache->optical_depth(flat_layer) = optical_depth;
                    layer_cache->source_factor(flat_layer) = factor;
                }
                if constexpr (WITH_TRANSPORT) {
                    const double source_factor =
                        prefix * albedo * (1.0 - attenuation);
                    ray_interpolation.source_for_layer(layer).visit(
                        [&](const auto sources) {
                            for (const auto& source : sources) {
                                const int transport_index =
                                    ray_interpolation.transport_value_offset +
                                    source.row_inner_index();
                                const double contribution =
                                    source.weight() * source_factor;
                                transport->values()(transport_index) +=
                                    contribution;
                            }
                        });
                }
                if constexpr (WITH_FORCING) {
                    if (layer_distance < minimum_layer_distance_m) {
                        have_shared_endpoint = false;
                    }
                    if (layer_distance >= minimum_layer_distance_m) {
                        if constexpr (LOWER_INTERPOLATION) {
                            const auto& traced_layer = rays[ray].layers[layer];
                            auto entrance_weights =
                                rays[ray].entrance_weights(layer);
                            auto exit_weights = rays[ray].exit_weights(layer);
                            const auto* start_weights = &entrance_weights;
                            const auto* end_weights = &exit_weights;
                            bool start_entrance = true;
                            bool end_entrance = false;
                            if (traced_layer.r_exit > traced_layer.r_entrance) {
                                end_weights = &entrance_weights;
                                end_entrance = true;
                            } else {
                                start_weights = &exit_weights;
                                start_entrance = false;
                            }
                            const int exit_solar = m_solar_offsets[ray] + layer;
                            auto start = scalar_endpoint(
                                ray_endpoint_context, layer, start_entrance,
                                exit_solar + 1, *start_weights);
                            auto end = scalar_endpoint(
                                ray_endpoint_context, layer, end_entrance,
                                exit_solar, *end_weights);
                            start.solar_transmission = solar(exit_solar + 1);
                            start.source = start.extinction * start.albedo *
                                           start.solar_transmission *
                                           start.phase * inverse_four_pi;
                            end.solar_transmission = solar(exit_solar);
                            end.source = end.extinction * end.albedo *
                                         end.solar_transmission * end.phase *
                                         inverse_four_pi;
                            radiance += prefix * factor *
                                        (start.source * source_quad_start +
                                         end.source * source_quad_end);
                        } else {
                            const int exit_solar = m_solar_offsets[ray] + layer;
                            auto start =
                                have_shared_endpoint
                                    ? shared_endpoint
                                    : scalar_endpoint(
                                          ray_endpoint_context, layer, true,
                                          exit_solar + 1,
                                          endpoint_weights(exit_solar + 1));
                            auto end = scalar_endpoint(
                                ray_endpoint_context, layer, false, exit_solar,
                                endpoint_weights(exit_solar));
                            if (!have_shared_endpoint) {
                                start.solar_transmission =
                                    solar(exit_solar + 1);
                                start.source = start.extinction * start.albedo *
                                               start.solar_transmission *
                                               start.phase * inverse_four_pi;
                            }
                            end.solar_transmission = solar(exit_solar);
                            end.source = end.extinction * end.albedo *
                                         end.solar_transmission * end.phase *
                                         inverse_four_pi;
                            radiance += prefix * factor *
                                        (start.source * source_quad_start +
                                         end.source * source_quad_end);
                            shared_endpoint = end;
                        }
                        have_shared_endpoint = true;
                    }
                }
                prefix *= attenuation;
            }
            if (volume_cache != nullptr) {
                volume_cache->ground_prefix(ray) = prefix;
            }
            if constexpr (WITH_TRANSPORT) {
                if (ray_interpolation.ground_is_hit()) {
                    const double ground_albedo =
                        ground_transport_albedo(wavelength, ray);
                    ray_interpolation.ground_weights.visit(
                        [&](const auto sources) {
                            for (const auto& source : sources) {
                                transport->values()(
                                    ray_interpolation.transport_value_offset +
                                    source.row_inner_index()) +=
                                    source.weight() * prefix * ground_albedo;
                            }
                        });
                }
            }
            if constexpr (WITH_FORCING) {
                if constexpr (LOWER_INTERPOLATION) {
                    const auto& traced_ray = rays[ray];
                    if (traced_ray.ground_is_hit &&
                        !traced_ray.layers.empty()) {
                        const auto& ground_layer = traced_ray.layers.front();
                        ScalarGroundGeometry ground{
                            ground_layer.exit.position.normalized(),
                            ground_layer.average_look_away};
                        double mu_in;
                        double mu_out;
                        double phi;
                        if (ground_scattering_geometry(m_solar_offsets[ray],
                                                       ground, mu_in, mu_out,
                                                       phi)) {
                            const auto brdf = m_atmosphere->surface().brdf(
                                wavelength, mu_in, mu_out, phi);
                            radiance += prefix * solar(m_solar_offsets[ray]) *
                                        mu_in * brdf(0, 0);
                        }
                    }
                } else if (packed_ray->ground_geometry >= 0) {
                    const auto& ground =
                        m_scalar_ground_geometry[packed_ray->ground_geometry];
                    double mu_in;
                    double mu_out;
                    double phi;
                    if (ground_scattering_geometry(m_solar_offsets[ray], ground,
                                                   mu_in, mu_out, phi)) {
                        const auto brdf = m_atmosphere->surface().brdf(
                            wavelength, mu_in, mu_out, phi,
                            m_ground_horizontal_weights[ray]);
                        radiance += prefix * solar(m_solar_offsets[ray]) *
                                    mu_in * brdf(0, 0);
                    }
                }
                forcing(ray) = radiance;
            }
        }
        if (volume_cache != nullptr) {
            volume_cache->active = true;
        }
        if (layer_cache != nullptr) {
            layer_cache->active = true;
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::calculate_scalar_jvp(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::Ref<Eigen::VectorXd> forcing,
        Eigen::Ref<Eigen::VectorXd> forcing_tangent,
        const Eigen::VectorXd* layer_state_projection,
        const Eigen::VectorXd* ground_state_projection,
        Eigen::VectorXd* direct_transport_tangent) {
        if (scalar_cache_index(wavelength) != wavelength_thread) {
            throw std::logic_error(
                "Successive-orders scalar JVP requires its prepared worker");
        }
        if (direct_transport_tangent != nullptr) {
            if (m_use_lower_interpolation) {
                calculate_scalar_jvp_impl<true, true>(
                    wavelength, wavelength_thread, native_tangent, forcing,
                    forcing_tangent, layer_state_projection,
                    ground_state_projection, direct_transport_tangent);
            } else {
                calculate_scalar_jvp_impl<true, false>(
                    wavelength, wavelength_thread, native_tangent, forcing,
                    forcing_tangent, layer_state_projection,
                    ground_state_projection, direct_transport_tangent);
            }
        } else if (m_use_lower_interpolation) {
            calculate_scalar_jvp_impl<false, true>(
                wavelength, wavelength_thread, native_tangent, forcing,
                forcing_tangent, nullptr, nullptr, nullptr);
        } else {
            calculate_scalar_jvp_impl<false, false>(
                wavelength, wavelength_thread, native_tangent, forcing,
                forcing_tangent, nullptr, nullptr, nullptr);
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::calculate_scalar_jvp_uniform_proportional(
        int wavelength, Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::Ref<const Eigen::VectorXd> solar_tangent,
        double extinction_direction_scale, double albedo, double albedo_tangent,
        const Eigen::VectorXd& layer_state_projection,
        const Eigen::VectorXd& ground_state_projection,
        Eigen::VectorXd& direct_transport_tangent,
        Eigen::Ref<Eigen::VectorXd> forcing_tangent) {
        const int cache_slot = scalar_cache_index(wavelength);
        const auto& interpolation = m_source_geometry->incoming_interpolation();
        const auto& solar = m_cached_solar_transmission[cache_slot];
        const auto& endpoint_medium = m_endpoint_medium_cache[cache_slot];
        const auto& layer_cache = m_scalar_layer_cache[cache_slot];
        if (!layer_cache.active ||
            layer_state_projection.size() != layer_cache.optical_depth.size() ||
            ground_state_projection.size() != m_num_rays ||
            direct_transport_tangent.size() != m_num_rays) {
            throw std::invalid_argument(
                "Invalid uniform scalar transport JVP storage");
        }
        direct_transport_tangent.setZero();
#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const auto& packed_ray = m_scalar_packed_rays[ray];
            const double phase_scale =
                m_uniform_phase_values[wavelength][ray] * inverse_four_pi;
            const auto endpoint = [&](int solar_index) {
                ScalarValueTangent result;
                const int slot = m_endpoint_slots[solar_index];
                const double extinction = endpoint_medium.extinction(slot);
                const double endpoint_albedo = endpoint_medium.albedo(slot);
                const double extinction_tangent =
                    extinction_direction_scale * extinction;
                const double extinction_albedo = extinction * endpoint_albedo;
                result.value =
                    extinction_albedo * solar(solar_index) * phase_scale;
                result.tangent =
                    ((extinction_tangent * endpoint_albedo +
                      extinction * albedo_tangent) *
                         solar(solar_index) +
                     extinction_albedo * solar_tangent(solar_index)) *
                    phase_scale;
                return result;
            };
            double prefix = 1.0;
            double prefix_tangent = 0.0;
            double radiance_tangent = 0.0;
            ScalarValueTangent shared_endpoint;
            bool have_shared_endpoint = false;
            for (std::uint32_t flat_layer = packed_ray.layer_end;
                 flat_layer-- > packed_ray.layer_begin;) {
                const auto& layer = m_scalar_packed_layers[flat_layer];
                const double optical_depth =
                    layer_cache.optical_depth(flat_layer);
                const double optical_depth_tangent =
                    extinction_direction_scale * optical_depth;
                const double factor = layer_cache.source_factor(flat_layer);
                const double attenuation = 1.0 - factor * optical_depth;
                const double attenuation_tangent =
                    -attenuation * optical_depth_tangent;
                const double source_factor_tangent =
                    prefix_tangent * albedo * (1.0 - attenuation) +
                    prefix * albedo_tangent * (1.0 - attenuation) -
                    prefix * albedo * attenuation_tangent;
                direct_transport_tangent(ray) +=
                    layer_state_projection(flat_layer) * source_factor_tangent;

                const double layer_distance =
                    layer.source_quad_start + layer.source_quad_end;
                if (layer_distance < minimum_layer_distance_m) {
                    have_shared_endpoint = false;
                } else {
                    const int local_layer =
                        static_cast<int>(flat_layer - packed_ray.layer_begin);
                    const int exit_solar = m_solar_offsets[ray] + local_layer;
                    const auto start = have_shared_endpoint
                                           ? shared_endpoint
                                           : endpoint(exit_solar + 1);
                    const auto end = endpoint(exit_solar);
                    const double factor_tangent =
                        constant_source_factor_derivative(optical_depth,
                                                          factor) *
                        optical_depth_tangent;
                    const double endpoint_value =
                        start.value * layer.source_quad_start +
                        end.value * layer.source_quad_end;
                    const double endpoint_tangent =
                        start.tangent * layer.source_quad_start +
                        end.tangent * layer.source_quad_end;
                    const double optical_endpoint_value =
                        start.value * layer.source_quad_start +
                        end.value * layer.source_quad_end;
                    radiance_tangent +=
                        prefix_tangent * factor * endpoint_value +
                        prefix * (factor * endpoint_tangent +
                                  factor_tangent * optical_endpoint_value);
                    shared_endpoint = end;
                    have_shared_endpoint = true;
                }
                prefix_tangent =
                    prefix_tangent * attenuation + prefix * attenuation_tangent;
                prefix *= attenuation;
            }
            if (interpolation[ray].ground_is_hit()) {
                const double ground_albedo =
                    ground_transport_albedo(wavelength, ray);
                const double ground_albedo_tangent =
                    ground_transport_albedo_tangent(ray, native_tangent);
                direct_transport_tangent(ray) +=
                    ground_state_projection(ray) *
                    (ground_albedo * prefix_tangent +
                     prefix * ground_albedo_tangent);
            }
            if (packed_ray.ground_geometry >= 0) {
                const auto& ground =
                    m_scalar_ground_geometry[packed_ray.ground_geometry];
                double mu_in;
                double mu_out;
                double phi;
                if (ground_scattering_geometry(m_solar_offsets[ray], ground,
                                               mu_in, mu_out, phi)) {
                    const auto brdf = m_atmosphere->surface().brdf(
                        wavelength, mu_in, mu_out, phi,
                        m_ground_horizontal_weights[ray]);
                    double brdf_tangent = 0.0;
                    for (int derivative = 0;
                         derivative < m_atmosphere->surface().num_deriv();
                         ++derivative) {
                        brdf_tangent +=
                            native_tangent(
                                m_atmosphere->surface_deriv_start_index() +
                                derivative) *
                            m_atmosphere->surface().d_brdf(
                                wavelength, mu_in, mu_out, phi, derivative,
                                m_ground_horizontal_weights[ray])(0, 0);
                    }
                    const double ground_source =
                        solar(m_solar_offsets[ray]) * mu_in * brdf(0, 0);
                    const double ground_tangent =
                        mu_in *
                        (solar_tangent(m_solar_offsets[ray]) * brdf(0, 0) +
                         solar(m_solar_offsets[ray]) * brdf_tangent);
                    radiance_tangent += prefix_tangent * ground_source +
                                        prefix * ground_tangent;
                }
            }
            forcing_tangent(ray) = radiance_tangent;
        }
    }

    template <int NSTOKES>
    template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION>
    void FirstOrderProvider<NSTOKES>::calculate_scalar_jvp_impl(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::Ref<Eigen::VectorXd> forcing,
        Eigen::Ref<Eigen::VectorXd> forcing_tangent,
        const Eigen::VectorXd* layer_state_projection,
        const Eigen::VectorXd* ground_state_projection,
        Eigen::VectorXd* direct_transport_tangent) {
        const int cache_slot = scalar_cache_index(wavelength);
        if constexpr (!LOWER_INTERPOLATION) {
            ensure_endpoint_medium(wavelength);
        }
        const auto& solar =
            ensure_solar_transmission(wavelength, wavelength_thread);
        FirstOrderScratchLease scratch(1, solar.size(),
                                       m_solar_table->table_size());
        auto solar_tangent = scratch.solar(0);
        auto table_tangent = scratch.table();
        const Eigen::Index solar_columns = m_solar_table->atmosphere_size();
        m_solar_table->apply(native_tangent.head(solar_columns), table_tangent);
        m_solar_interpolation.apply(table_tangent, solar_tangent);
        solar_tangent.array() *= -solar.array();
        auto& coefficient_tangent = m_phase_product_scratch[wavelength_thread];
        auto& tangent_orders = m_phase_order_scratch[wavelength_thread];
        const auto& storage = m_atmosphere->storage();
        const int locations = storage.total_extinction.rows();
        bool phase_tangent_active = false;
        if (m_atmosphere->num_scattering_deriv_groups() > 0) {
            coefficient_tangent.setZero();
            std::fill(tangent_orders.begin(), tangent_orders.end(), 0);
            for (int group = 0;
                 group < m_atmosphere->num_scattering_deriv_groups(); ++group) {
                const int native_offset =
                    m_atmosphere->scat_deriv_start_index() + group * locations;
                for (int location = 0; location < locations; ++location) {
                    const double direction =
                        native_tangent(native_offset + location);
                    if (direction == 0.0) {
                        continue;
                    }
                    const int order = std::min(
                        storage.d_max_order[group](location, wavelength),
                        m_num_phase_moments);
                    const int coefficient_offset =
                        location * m_num_phase_moments;
                    for (int degree = 0; degree < order; ++degree) {
                        coefficient_tangent(coefficient_offset + degree) +=
                            direction * storage.d_leg_coeff(degree, location,
                                                            wavelength, group);
                    }
                }
            }
            for (int location = 0; location < locations; ++location) {
                int order = m_num_phase_moments;
                const int coefficient_offset = location * m_num_phase_moments;
                while (order > 0 && coefficient_tangent(coefficient_offset +
                                                        order - 1) == 0.0) {
                    --order;
                }
                tangent_orders[location] = order;
                phase_tangent_active = phase_tangent_active || order > 0;
            }
        }
        const auto& rays = m_source_geometry->incoming_rays();
        const auto& interpolation = m_source_geometry->incoming_interpolation();
        const auto& ssa = m_atmosphere->storage().ssa;
        const double* __restrict albedo_values = ssa.col(wavelength).data();
        const double* __restrict extinction_direction = native_tangent.data();
        const double* __restrict albedo_direction =
            native_tangent.data() + m_atmosphere->ssa_deriv_start_index();
        const double* __restrict solar_direction = solar_tangent.data();
        const double* __restrict extinction_values =
            storage.total_extinction.col(wavelength).data();
        bool proportional_extinction_direction = locations > 0;
        double extinction_direction_scale = 0.0;
        int proportional_reference = 0;
        while (proportional_reference < locations &&
               extinction_values[proportional_reference] == 0.0) {
            proportional_extinction_direction =
                proportional_extinction_direction &&
                extinction_direction[proportional_reference] == 0.0;
            ++proportional_reference;
        }
        if (proportional_reference < locations &&
            proportional_extinction_direction) {
            extinction_direction_scale =
                extinction_direction[proportional_reference] /
                extinction_values[proportional_reference];
        }
        constexpr double proportional_tolerance =
            32.0 * std::numeric_limits<double>::epsilon();
        for (int location = proportional_reference + 1;
             location < locations && proportional_extinction_direction;
             ++location) {
            const double expected =
                extinction_direction_scale * extinction_values[location];
            const double scale = std::max(
                {std::abs(extinction_direction[location]), std::abs(expected),
                 std::numeric_limits<double>::min()});
            proportional_extinction_direction =
                std::abs(extinction_direction[location] - expected) <=
                proportional_tolerance * scale;
        }
        const bool uniform_albedo = m_uniform_albedo_active[wavelength] != 0;
        const double uniform_albedo_value =
            uniform_albedo ? m_uniform_albedo_values[wavelength] : 0.0;
        bool uniform_albedo_direction = locations > 0;
        const double uniform_albedo_tangent =
            uniform_albedo_direction ? albedo_direction[0] : 0.0;
        for (int location = 1; location < locations && uniform_albedo_direction;
             ++location) {
            uniform_albedo_direction =
                albedo_direction[location] == uniform_albedo_tangent;
        }
        if constexpr (WITH_TRANSPORT && !LOWER_INTERPOLATION) {
            if (m_uniform_phase_active[wavelength] != 0 &&
                !m_endpoint_phase_basis && !phase_tangent_active &&
                proportional_extinction_direction && uniform_albedo &&
                uniform_albedo_direction) {
                if (layer_state_projection == nullptr ||
                    ground_state_projection == nullptr ||
                    direct_transport_tangent == nullptr) {
                    throw std::invalid_argument(
                        "Missing uniform scalar transport JVP storage");
                }
                calculate_scalar_jvp_uniform_proportional(
                    wavelength, native_tangent, solar_tangent,
                    extinction_direction_scale, uniform_albedo_value,
                    uniform_albedo_tangent, *layer_state_projection,
                    *ground_state_projection, *direct_transport_tangent,
                    forcing_tangent);
                profile_product_scratch(this, scratch, "jvp", wavelength,
                                        wavelength_thread);
                return;
            }
        }
        const int endpoint_count = static_cast<int>(m_endpoint_stencils.size());
        scratch.prepare_endpoints(endpoint_count);
        auto endpoint_extinction_tangent = scratch.endpoint_extinction();
        auto endpoint_albedo_tangent = scratch.endpoint_albedo();
        profile_product_scratch(this, scratch, "jvp", wavelength,
                                wavelength_thread);
        m_endpoint_stencils.visit([&](const auto& stencils) {
            for (int slot = 0; slot < endpoint_extinction_tangent.size();
                 ++slot) {
                double extinction_value =
                    proportional_extinction_direction
                        ? extinction_direction_scale *
                              m_endpoint_medium_cache[cache_slot].extinction(
                                  slot)
                        : 0.0;
                double albedo_value = 0.0;
                const auto weights = stencils[slot];
                for (std::size_t index = 0; index < weights.size(); ++index) {
                    const auto [atmosphere_index, weight] =
                        interpolation_pair(weights[index]);
                    if (!proportional_extinction_direction) {
                        extinction_value +=
                            weight * extinction_direction[atmosphere_index];
                    }
                    if (!uniform_albedo_direction) {
                        albedo_value +=
                            weight * albedo_direction[atmosphere_index];
                    }
                }
                endpoint_extinction_tangent(slot) = extinction_value;
                endpoint_albedo_tangent(slot) = uniform_albedo_direction
                                                    ? uniform_albedo_tangent
                                                    : albedo_value;
            }
        });
        const double* __restrict layer_state =
            layer_state_projection == nullptr ? nullptr
                                              : layer_state_projection->data();
        const double* __restrict ground_state =
            ground_state_projection == nullptr
                ? nullptr
                : ground_state_projection->data();
        double* __restrict direct_transport_direction =
            direct_transport_tangent == nullptr
                ? nullptr
                : direct_transport_tangent->data();
        if (m_scalar_layer_cache.empty() ||
            !m_scalar_layer_cache[cache_slot].active) {
            throw std::logic_error(
                "Successive-orders scalar JVP requires a prepared primal");
        }
        const auto& layer_cache = m_scalar_layer_cache[cache_slot];
        if constexpr (WITH_TRANSPORT) {
            if (layer_state_projection == nullptr ||
                layer_state_projection->size() !=
                    layer_cache.optical_depth.size() ||
                ground_state_projection == nullptr ||
                ground_state_projection->size() != m_num_rays ||
                direct_transport_tangent == nullptr ||
                direct_transport_tangent->size() != m_num_rays) {
                throw std::invalid_argument(
                    "Invalid direct scalar transport JVP storage");
            }
            direct_transport_tangent->setZero();
        }

        auto endpoint_context =
            prepare_scalar_endpoint_context(wavelength, cache_slot);
        endpoint_context.coefficient_tangent = coefficient_tangent.data();
        endpoint_context.tangent_orders = tangent_orders.data();
        endpoint_context.endpoint_extinction_tangent =
            endpoint_extinction_tangent.data();
        endpoint_context.endpoint_albedo_tangent =
            endpoint_albedo_tangent.data();

#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const auto ray_endpoint_context =
                scalar_endpoint_context_for_ray(endpoint_context, ray);
            const auto& ray_interpolation = interpolation[ray];
            const ScalarPackedRay* packed_ray = nullptr;
            if constexpr (!LOWER_INTERPOLATION) {
                packed_ray = &m_scalar_packed_rays[ray];
            }
            const int flat_layer_offset = m_solar_offsets[ray] - ray;
            double prefix = 1.0;
            double prefix_tangent = 0.0;
            double radiance_tangent = 0.0;
            ScalarValueTangent shared_endpoint;
            bool have_shared_endpoint = false;
            for (int layer =
                     static_cast<int>(ray_interpolation.layers.size()) - 1;
                 layer >= 0; --layer) {
                double source_quad_start;
                double source_quad_end;
                double layer_distance;
                if constexpr (LOWER_INTERPOLATION) {
                    const auto& traced_layer = rays[ray].layers[layer];
                    layer_distance = traced_layer.layer_distance;
                    source_quad_start =
                        layer_distance * traced_layer.od_quad_start_fraction;
                    source_quad_end =
                        layer_distance * traced_layer.od_quad_end_fraction;
                } else {
                    const auto& packed_layer =
                        m_scalar_packed_layers[packed_ray->layer_begin + layer];
                    source_quad_start = packed_layer.source_quad_start;
                    source_quad_end = packed_layer.source_quad_end;
                    layer_distance = source_quad_start + source_quad_end;
                }
                const int flat_layer = flat_layer_offset + layer;
                const double optical_depth =
                    layer_cache.optical_depth(flat_layer);
                double optical_depth_tangent =
                    proportional_extinction_direction
                        ? extinction_direction_scale * optical_depth
                        : 0.0;
                double albedo = uniform_albedo_value;
                double albedo_tangent =
                    uniform_albedo_direction ? uniform_albedo_tangent : 0.0;
                if (!proportional_extinction_direction) {
                    const auto od_weights =
                        ray_interpolation.optical_depth_for_layer(layer);
                    for (std::size_t index = 0; index < od_weights.size();
                         ++index) {
                        const auto [atmosphere_index, weight] =
                            od_weights[index];
                        optical_depth_tangent +=
                            weight * extinction_direction[atmosphere_index];
                    }
                }
                if constexpr (WITH_TRANSPORT) {
                    if (!uniform_albedo || !uniform_albedo_direction) {
                        for (const auto& weight :
                             ray_interpolation.atmosphere_for_layer(layer)) {
                            if (!uniform_albedo) {
                                albedo += weight.weight() *
                                          albedo_values[weight.index];
                            }
                            if (!uniform_albedo_direction) {
                                albedo_tangent +=
                                    weight.weight() *
                                    albedo_direction[weight.index];
                            }
                        }
                    }
                }
                const double factor = layer_cache.source_factor(flat_layer);
                const double attenuation = 1.0 - factor * optical_depth;
                const double attenuation_tangent =
                    -attenuation * optical_depth_tangent;
                if constexpr (WITH_TRANSPORT) {
                    const double source_factor_tangent =
                        prefix_tangent * albedo * (1.0 - attenuation) +
                        prefix * albedo_tangent * (1.0 - attenuation) -
                        prefix * albedo * attenuation_tangent;
                    direct_transport_direction[ray] +=
                        layer_state[flat_layer] * source_factor_tangent;
                }
                if (layer_distance < minimum_layer_distance_m) {
                    have_shared_endpoint = false;
                }
                if (layer_distance >= minimum_layer_distance_m) {
                    ScalarValueTangent start;
                    ScalarValueTangent end;
                    if constexpr (LOWER_INTERPOLATION) {
                        const auto& traced_layer = rays[ray].layers[layer];
                        auto entrance_weights =
                            rays[ray].entrance_weights(layer);
                        auto exit_weights = rays[ray].exit_weights(layer);
                        const auto* start_weights = &entrance_weights;
                        const auto* end_weights = &exit_weights;
                        bool start_entrance = true;
                        bool end_entrance = false;
                        if (traced_layer.r_exit > traced_layer.r_entrance) {
                            end_weights = &entrance_weights;
                            end_entrance = true;
                        } else {
                            start_weights = &exit_weights;
                            start_entrance = false;
                        }
                        const int exit_solar = m_solar_offsets[ray] + layer;
                        start = scalar_endpoint_jvp<false>(
                            ray_endpoint_context, layer, start_entrance,
                            exit_solar + 1, *start_weights,
                            extinction_direction, albedo_direction,
                            solar_direction, uniform_albedo_direction,
                            uniform_albedo_tangent, phase_tangent_active);
                        end = scalar_endpoint_jvp<false>(
                            ray_endpoint_context, layer, end_entrance,
                            exit_solar, *end_weights, extinction_direction,
                            albedo_direction, solar_direction,
                            uniform_albedo_direction, uniform_albedo_tangent,
                            phase_tangent_active);
                    } else {
                        const int exit_solar = m_solar_offsets[ray] + layer;
                        start =
                            have_shared_endpoint
                                ? shared_endpoint
                                : scalar_endpoint_jvp<true>(
                                      ray_endpoint_context, layer, true,
                                      exit_solar + 1,
                                      endpoint_weights(exit_solar + 1),
                                      extinction_direction, albedo_direction,
                                      solar_direction, uniform_albedo_direction,
                                      uniform_albedo_tangent,
                                      phase_tangent_active);
                        end = scalar_endpoint_jvp<true>(
                            ray_endpoint_context, layer, false, exit_solar,
                            endpoint_weights(exit_solar), extinction_direction,
                            albedo_direction, solar_direction,
                            uniform_albedo_direction, uniform_albedo_tangent,
                            phase_tangent_active);
                    }
                    const double factor_tangent =
                        constant_source_factor_derivative(optical_depth,
                                                          factor) *
                        optical_depth_tangent;
                    const double endpoint_value =
                        start.value * source_quad_start +
                        end.value * source_quad_end;
                    const double endpoint_tangent =
                        start.tangent * source_quad_start +
                        end.tangent * source_quad_end;
                    const double optical_endpoint_value =
                        start.value * source_quad_start +
                        end.value * source_quad_end;
                    radiance_tangent +=
                        prefix_tangent * factor * endpoint_value +
                        prefix * (factor * endpoint_tangent +
                                  factor_tangent * optical_endpoint_value);
                    shared_endpoint = end;
                    have_shared_endpoint = true;
                }
                prefix_tangent =
                    prefix_tangent * attenuation + prefix * attenuation_tangent;
                prefix *= attenuation;
            }
            if constexpr (WITH_TRANSPORT) {
                if (ray_interpolation.ground_is_hit()) {
                    const double ground_albedo =
                        ground_transport_albedo(wavelength, ray);
                    const double ground_albedo_tangent =
                        ground_transport_albedo_tangent(ray, native_tangent);
                    direct_transport_direction[ray] +=
                        ground_state[ray] * (ground_albedo * prefix_tangent +
                                             prefix * ground_albedo_tangent);
                }
            }
            ScalarGroundGeometry local_ground;
            const ScalarGroundGeometry* ground = nullptr;
            if constexpr (LOWER_INTERPOLATION) {
                const auto& traced_ray = rays[ray];
                if (traced_ray.ground_is_hit && !traced_ray.layers.empty()) {
                    const auto& layer = traced_ray.layers.front();
                    local_ground = {layer.exit.position.normalized(),
                                    layer.average_look_away};
                    ground = &local_ground;
                }
            } else if (packed_ray->ground_geometry >= 0) {
                ground = &m_scalar_ground_geometry[packed_ray->ground_geometry];
            }
            if (ground != nullptr) {
                double mu_in;
                double mu_out;
                double phi;
                if (ground_scattering_geometry(m_solar_offsets[ray], *ground,
                                               mu_in, mu_out, phi)) {
                    const auto brdf = m_atmosphere->surface().brdf(
                        wavelength, mu_in, mu_out, phi,
                        m_ground_horizontal_weights[ray]);
                    double brdf_tangent = 0.0;
                    for (int derivative = 0;
                         derivative < m_atmosphere->surface().num_deriv();
                         ++derivative) {
                        brdf_tangent +=
                            native_tangent(
                                m_atmosphere->surface_deriv_start_index() +
                                derivative) *
                            m_atmosphere->surface().d_brdf(
                                wavelength, mu_in, mu_out, phi, derivative,
                                m_ground_horizontal_weights[ray])(0, 0);
                    }
                    const double ground_source =
                        solar(m_solar_offsets[ray]) * mu_in * brdf(0, 0);
                    const double ground_tangent =
                        mu_in *
                        (solar_tangent(m_solar_offsets[ray]) * brdf(0, 0) +
                         solar(m_solar_offsets[ray]) * brdf_tangent);
                    radiance_tangent += prefix_tangent * ground_source +
                                        prefix * ground_tangent;
                }
            }
            forcing_tangent(ray) = radiance_tangent;
        }
        (void)forcing;
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::accumulate_scalar_vjp(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient,
        const Eigen::VectorXd* transport_state,
        const Eigen::VectorXd* layer_state_projection,
        const Eigen::VectorXd* ground_state_projection) {
        if (scalar_cache_index(wavelength) != wavelength_thread) {
            throw std::logic_error(
                "Successive-orders scalar VJP requires its prepared worker");
        }
        const bool surface_only =
            m_atmosphere->storage().derivative_mappings_const().empty();
        if (transport_state != nullptr || layer_state_projection != nullptr) {
            if (m_use_lower_interpolation) {
                if (surface_only) {
                    accumulate_scalar_surface_vjp<true, true>(
                        wavelength, wavelength_thread, forcing_cotangent,
                        native_gradient, transport_state,
                        ground_state_projection);
                } else {
                    accumulate_scalar_vjp_dispatch<true, true>(
                        wavelength, wavelength_thread, forcing_cotangent,
                        native_gradient, transport_state,
                        layer_state_projection, ground_state_projection);
                }
            } else {
                if (surface_only) {
                    accumulate_scalar_surface_vjp<true, false>(
                        wavelength, wavelength_thread, forcing_cotangent,
                        native_gradient, transport_state,
                        ground_state_projection);
                } else {
                    accumulate_scalar_vjp_dispatch<true, false>(
                        wavelength, wavelength_thread, forcing_cotangent,
                        native_gradient, transport_state,
                        layer_state_projection, ground_state_projection);
                }
            }
        } else if (m_use_lower_interpolation) {
            if (surface_only) {
                accumulate_scalar_surface_vjp<false, true>(
                    wavelength, wavelength_thread, forcing_cotangent,
                    native_gradient, nullptr, nullptr);
            } else {
                accumulate_scalar_vjp_dispatch<false, true>(
                    wavelength, wavelength_thread, forcing_cotangent,
                    native_gradient, nullptr, nullptr, nullptr);
            }
        } else {
            if (surface_only) {
                accumulate_scalar_surface_vjp<false, false>(
                    wavelength, wavelength_thread, forcing_cotangent,
                    native_gradient, nullptr, nullptr);
            } else {
                accumulate_scalar_vjp_dispatch<false, false>(
                    wavelength, wavelength_thread, forcing_cotangent,
                    native_gradient, nullptr, nullptr, nullptr);
            }
        }
    }

    template <int NSTOKES>
    template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION>
    void FirstOrderProvider<NSTOKES>::accumulate_scalar_surface_vjp(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient,
        const Eigen::VectorXd* transport_state,
        const Eigen::VectorXd* ground_state_projection) {
        const int cache_slot = scalar_cache_index(wavelength);
        const int first_thread = wavelength_thread;
        const int last_thread = first_thread + m_num_source_threads;
        if (last_thread > static_cast<int>(m_gradient_scratch.size())) {
            throw std::logic_error(
                "Successive-orders compact scalar thread layout is invalid");
        }
        for (int thread = first_thread; thread < last_thread; ++thread) {
            m_gradient_scratch[thread].setZero();
        }
        const auto& rays = m_source_geometry->incoming_rays();
        const auto& interpolation = m_source_geometry->incoming_interpolation();
        const auto& extinction = m_atmosphere->storage().total_extinction;
        const auto& solar =
            ensure_solar_transmission(wavelength, wavelength_thread);
        const ScalarLayerCache* layer_cache =
            m_scalar_layer_cache.empty() ? nullptr
                                         : &m_scalar_layer_cache[cache_slot];
        const bool use_layer_cache =
            layer_cache != nullptr && layer_cache->active;
        const ScalarVolumeCache* volume_cache = nullptr;
        if constexpr (WITH_TRANSPORT && !LOWER_INTERPOLATION) {
            if (!m_scalar_volume_cache.empty() &&
                m_scalar_volume_cache[cache_slot].active) {
                volume_cache = &m_scalar_volume_cache[cache_slot];
            }
        }
#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const double forcing_gradient = forcing_cotangent(ray);
            if (forcing_gradient == 0.0) {
                continue;
            }
            const int thread = ray_thread_index(wavelength_thread);
            auto thread_gradient = m_gradient_scratch[thread].col(0);
            const auto& ray_interpolation = interpolation[ray];
            const ScalarPackedRay* packed_ray = nullptr;
            if constexpr (!LOWER_INTERPOLATION) {
                packed_ray = &m_scalar_packed_rays[ray];
            }
            double prefix = 1.0;
            if (volume_cache != nullptr) {
                prefix = volume_cache->ground_prefix(ray);
            } else {
                const int flat_layer_offset =
                    LOWER_INTERPOLATION
                        ? m_solar_offsets[ray] - ray
                        : static_cast<int>(packed_ray->layer_begin);
                for (int layer = 0;
                     layer < static_cast<int>(ray_interpolation.layers.size());
                     ++layer) {
                    double optical_depth = 0.0;
                    double attenuation;
                    if (use_layer_cache) {
                        const int flat_layer = flat_layer_offset + layer;
                        optical_depth = layer_cache->optical_depth(flat_layer);
                        attenuation =
                            1.0 - layer_cache->source_factor(flat_layer) *
                                      optical_depth;
                    } else {
                        const auto weights =
                            ray_interpolation.optical_depth_for_layer(layer);
                        for (std::size_t index = 0; index < weights.size();
                             ++index) {
                            const auto [atmosphere_index, weight] =
                                weights[index];
                            optical_depth +=
                                weight *
                                extinction(atmosphere_index, wavelength);
                        }
                        attenuation = std::exp(-optical_depth);
                    }
                    prefix *= attenuation;
                }
            }

            ScalarGroundGeometry local_ground;
            const ScalarGroundGeometry* ground = nullptr;
            if constexpr (LOWER_INTERPOLATION) {
                const auto& traced_ray = rays[ray];
                if (traced_ray.ground_is_hit && !traced_ray.layers.empty()) {
                    const auto& layer = traced_ray.layers.front();
                    local_ground = {layer.exit.position.normalized(),
                                    layer.average_look_away};
                    ground = &local_ground;
                }
            } else if (packed_ray->ground_geometry >= 0) {
                ground = &m_scalar_ground_geometry[packed_ray->ground_geometry];
            }
            if (ground != nullptr) {
                double mu_in;
                double mu_out;
                double phi;
                if (ground_scattering_geometry(m_solar_offsets[ray], *ground,
                                               mu_in, mu_out, phi)) {
                    for (int derivative = 0;
                         derivative < m_atmosphere->surface().num_deriv();
                         ++derivative) {
                        thread_gradient(
                            m_atmosphere->surface_deriv_start_index() +
                            derivative) +=
                            prefix * forcing_gradient *
                            solar(m_solar_offsets[ray]) * mu_in *
                            m_atmosphere->surface().d_brdf(
                                wavelength, mu_in, mu_out, phi, derivative,
                                m_ground_horizontal_weights[ray])(0, 0);
                    }
                }
            }
            if constexpr (WITH_TRANSPORT) {
                if (ray_interpolation.ground_is_hit()) {
                    double ground_value = 0.0;
                    if (ground_state_projection != nullptr) {
                        ground_value = (*ground_state_projection)(ray);
                    } else {
                        const auto transport_columns =
                            m_source_geometry->transport_columns_for_ray(ray);
                        transport_columns.visit([&](const auto& columns) {
                            ray_interpolation.ground_weights.visit(
                                [&](const auto sources) {
                                    for (const auto& source : sources) {
                                        ground_value +=
                                            source.weight() *
                                            (*transport_state)(
                                                columns
                                                    [source.row_inner_index()]);
                                    }
                                });
                        });
                    }
                    accumulate_ground_transport_albedo_vjp(
                        ray, prefix * ground_value * forcing_gradient,
                        thread_gradient);
                }
            }
        }
        for (int thread = first_thread; thread < last_thread; ++thread) {
            native_gradient += m_gradient_scratch[thread].col(0);
        }
    }

    template <int NSTOKES>
    template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION>
    void FirstOrderProvider<NSTOKES>::accumulate_scalar_vjp_dispatch(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient,
        const Eigen::VectorXd* transport_state,
        const Eigen::VectorXd* layer_state_projection,
        const Eigen::VectorXd* ground_state_projection) {
        if (m_atmosphere->num_scattering_deriv_groups() > 0) {
            accumulate_scalar_vjp_impl<WITH_TRANSPORT, LOWER_INTERPOLATION,
                                       true>(
                wavelength, wavelength_thread, forcing_cotangent,
                native_gradient, transport_state, layer_state_projection,
                ground_state_projection);
        } else {
            accumulate_scalar_vjp_impl<WITH_TRANSPORT, LOWER_INTERPOLATION,
                                       false>(
                wavelength, wavelength_thread, forcing_cotangent,
                native_gradient, transport_state, layer_state_projection,
                ground_state_projection);
        }
    }

    template <int NSTOKES>
    template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION,
              bool WITH_PHASE_GRADIENT>
    void FirstOrderProvider<NSTOKES>::accumulate_scalar_vjp_impl(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient,
        const Eigen::VectorXd* transport_state,
        const Eigen::VectorXd* layer_state_projection,
        const Eigen::VectorXd* ground_state_projection) {
        const int cache_slot = scalar_cache_index(wavelength);
        const int first_thread = wavelength_thread;
        const int last_thread = first_thread + m_num_source_threads;
        if (last_thread > static_cast<int>(m_gradient_scratch.size())) {
            throw std::logic_error(
                "Successive-orders compact scalar thread layout is invalid");
        }
        for (int thread = first_thread; thread < last_thread; ++thread) {
            m_gradient_scratch[thread].setZero();
            if constexpr (WITH_PHASE_GRADIENT) {
                m_phase_product_scratch[thread].setZero();
            }
        }
        const auto& rays = m_source_geometry->incoming_rays();
        const auto& interpolation = m_source_geometry->incoming_interpolation();
        const auto& extinction = m_atmosphere->storage().total_extinction;
        const auto& ssa = m_atmosphere->storage().ssa;
        const auto& solar =
            ensure_solar_transmission(wavelength, wavelength_thread);
        FirstOrderScratchLease product_scratch(
            m_num_source_threads, solar.size(), m_solar_table->table_size());
        for (int thread = 0; thread < m_num_source_threads; ++thread) {
            product_scratch.solar(thread).setZero();
        }
        profile_product_scratch(this, product_scratch, "vjp", wavelength,
                                wavelength_thread);
        ensure_endpoint_medium(wavelength);
        const auto endpoint_context =
            prepare_scalar_endpoint_context(wavelength, cache_slot);
        const bool uniform_phase = m_uniform_phase_active[wavelength] != 0;
        const bool uniform_albedo = m_uniform_albedo_active[wavelength] != 0;
        const double uniform_albedo_value =
            uniform_albedo ? m_uniform_albedo_values[wavelength] : 0.0;
        const ScalarLayerCache* layer_cache =
            m_scalar_layer_cache.empty() ? nullptr
                                         : &m_scalar_layer_cache[cache_slot];
        const bool use_layer_cache =
            layer_cache != nullptr && layer_cache->active;
        const double* __restrict layer_state =
            layer_state_projection == nullptr ? nullptr
                                              : layer_state_projection->data();
        const double* __restrict ground_state =
            ground_state_projection == nullptr
                ? nullptr
                : ground_state_projection->data();
#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const auto ray_endpoint_context =
                scalar_endpoint_context_for_ray(endpoint_context, ray);
            const int thread = ray_thread_index(wavelength_thread);
            auto thread_gradient = m_gradient_scratch[thread].col(0);
            auto solar_gradient = product_scratch.solar(thread - first_thread);
            auto& coefficient_gradient = m_phase_product_scratch[thread];
            auto& scratch = m_scalar_vjp_scratch[thread];
            const auto& ray_interpolation = interpolation[ray];
            const ScalarPackedRay* packed_ray = nullptr;
            if constexpr (!LOWER_INTERPOLATION) {
                packed_ray = &m_scalar_packed_rays[ray];
            }
            const int flat_layer_offset =
                LOWER_INTERPOLATION ? m_solar_offsets[ray] - ray
                                    : static_cast<int>(packed_ray->layer_begin);
            const auto transport_columns =
                transport_state == nullptr
                    ? TransportColumnView{}
                    : m_source_geometry->transport_columns_for_ray(ray);
            transport_columns.visit([&](const auto& columns) {
                const auto* transport_column_data = columns.data();
                const int num_layers =
                    static_cast<int>(ray_interpolation.layers.size());
                double prefix = 1.0;
                for (int layer = num_layers - 1; layer >= 0; --layer) {
                    double optical_depth = 0.0;
                    double attenuation = 0.0;
                    double factor = 1.0;
                    if (use_layer_cache) {
                        const int flat_layer = flat_layer_offset + layer;
                        optical_depth = layer_cache->optical_depth(flat_layer);
                        factor = layer_cache->source_factor(flat_layer);
                        attenuation = 1.0 - factor * optical_depth;
                    } else {
                        const auto od_weights =
                            ray_interpolation.optical_depth_for_layer(layer);
                        for (std::size_t index = 0; index < od_weights.size();
                             ++index) {
                            const auto [atmosphere_index, weight] =
                                od_weights[index];
                            optical_depth +=
                                weight *
                                extinction(atmosphere_index, wavelength);
                        }
                        const auto layer_transfer =
                            scalar_layer_transfer(optical_depth);
                        attenuation = layer_transfer.attenuation;
                        factor = layer_transfer.source_factor;
                    }
                    if (!use_layer_cache) {
                        scratch.optical_depth(layer) = optical_depth;
                        scratch.attenuation(layer) = attenuation;
                        scratch.source_factor(layer) = factor;
                    }
                    scratch.prefix_attenuation(layer) = prefix;
                    prefix *= attenuation;
                }
                if (!uniform_albedo) {
                    for (int layer = 0; layer < num_layers; ++layer) {
                        double albedo = 0.0;
                        for (const auto& weight :
                             ray_interpolation.atmosphere_for_layer(layer)) {
                            albedo +=
                                weight.weight() * ssa(weight.index, wavelength);
                        }
                        scratch.albedo(layer) = albedo;
                    }
                }
                if constexpr (!LOWER_INTERPOLATION) {
                    if (num_layers > 0) {
                        scratch.endpoint_cotangent.head(num_layers + 1)
                            .setZero();
                        for (int endpoint = 0; endpoint <= num_layers;
                             ++endpoint) {
                            const int solar_index =
                                m_solar_offsets[ray] + endpoint;
                            auto& value = scratch.endpoints[endpoint];
                            if (uniform_phase) {
                                const int slot =
                                    ray_endpoint_context
                                        .endpoint_slots[solar_index];
                                value.extinction =
                                    ray_endpoint_context
                                        .endpoint_extinction[slot];
                                value.albedo =
                                    ray_endpoint_context.endpoint_albedo[slot];
                                value.phase =
                                    ray_endpoint_context.uniform_phase
                                        [ray_endpoint_context
                                                     .endpoint_basis_stride == 0
                                             ? 0
                                             : solar_index];
                            } else {
                                const auto weights =
                                    endpoint_weights(solar_index);
                                value = scalar_endpoint(
                                    ray_endpoint_context,
                                    endpoint == 0 ? 0 : endpoint - 1,
                                    endpoint != 0, solar_index, weights);
                            }
                            value.solar_transmission = solar(solar_index);
                            value.source = value.extinction * value.albedo *
                                           value.solar_transmission *
                                           value.phase * inverse_four_pi;
                        }
                    }
                }
                const double forcing_gradient = forcing_cotangent(ray);
                double prefix_cotangent = 0.0;
                ScalarGroundGeometry local_ground;
                const ScalarGroundGeometry* ground = nullptr;
                if constexpr (LOWER_INTERPOLATION) {
                    const auto& traced_ray = rays[ray];
                    if (traced_ray.ground_is_hit &&
                        !traced_ray.layers.empty()) {
                        const auto& layer = traced_ray.layers.front();
                        local_ground = {layer.exit.position.normalized(),
                                        layer.average_look_away};
                        ground = &local_ground;
                    }
                } else if (packed_ray->ground_geometry >= 0) {
                    ground =
                        &m_scalar_ground_geometry[packed_ray->ground_geometry];
                }
                if (ground != nullptr) {
                    double mu_in;
                    double mu_out;
                    double phi;
                    if (ground_scattering_geometry(m_solar_offsets[ray],
                                                   *ground, mu_in, mu_out,
                                                   phi)) {
                        const auto brdf = m_atmosphere->surface().brdf(
                            wavelength, mu_in, mu_out, phi,
                            m_ground_horizontal_weights[ray]);
                        const double ground_source =
                            solar(m_solar_offsets[ray]) * mu_in * brdf(0, 0);
                        prefix_cotangent = forcing_gradient * ground_source;
                        solar_gradient(m_solar_offsets[ray]) +=
                            prefix * forcing_gradient * mu_in * brdf(0, 0);
                        for (int derivative = 0;
                             derivative < m_atmosphere->surface().num_deriv();
                             ++derivative) {
                            thread_gradient(
                                m_atmosphere->surface_deriv_start_index() +
                                derivative) +=
                                prefix * forcing_gradient *
                                solar(m_solar_offsets[ray]) * mu_in *
                                m_atmosphere->surface().d_brdf(
                                    wavelength, mu_in, mu_out, phi, derivative,
                                    m_ground_horizontal_weights[ray])(0, 0);
                        }
                    }
                }
                if constexpr (WITH_TRANSPORT) {
                    if (ray_interpolation.ground_is_hit()) {
                        double ground_value = 0.0;
                        if (ground_state != nullptr) {
                            ground_value = ground_state[ray];
                        } else {
                            ray_interpolation.ground_weights.visit(
                                [&](const auto sources) {
                                    for (const auto& source : sources) {
                                        ground_value +=
                                            source.weight() *
                                            (*transport_state)(
                                                transport_column_data
                                                    [source.row_inner_index()]);
                                    }
                                });
                        }
                        const double ground_albedo =
                            ground_transport_albedo(wavelength, ray);
                        prefix_cotangent +=
                            ground_albedo * ground_value * forcing_gradient;
                        accumulate_ground_transport_albedo_vjp(
                            ray, prefix * ground_value * forcing_gradient,
                            thread_gradient);
                    }
                }

                for (int layer = 0; layer < num_layers; ++layer) {
                    const sasktran2::raytracing::TracedLayer* traced_layer =
                        nullptr;
                    double layer_distance;
                    double source_quad_start;
                    double source_quad_end;
                    if constexpr (LOWER_INTERPOLATION) {
                        traced_layer = &rays[ray].layers[layer];
                        layer_distance = traced_layer->layer_distance;
                        source_quad_start = traced_layer->od_quad_start;
                        source_quad_end = traced_layer->od_quad_end;
                    } else {
                        const auto& packed_layer =
                            m_scalar_packed_layers[packed_ray->layer_begin +
                                                   layer];
                        source_quad_start = packed_layer.source_quad_start;
                        source_quad_end = packed_layer.source_quad_end;
                        layer_distance = source_quad_start + source_quad_end;
                    }
                    const int flat_layer = flat_layer_offset + layer;
                    const double optical_depth =
                        use_layer_cache ? layer_cache->optical_depth(flat_layer)
                                        : scratch.optical_depth(layer);
                    const double factor =
                        use_layer_cache ? layer_cache->source_factor(flat_layer)
                                        : scratch.source_factor(layer);
                    const double attenuation =
                        use_layer_cache ? 1.0 - factor * optical_depth
                                        : scratch.attenuation(layer);
                    const double layer_prefix =
                        scratch.prefix_attenuation(layer);
                    double layer_prefix_cotangent = 0.0;
                    double optical_depth_cotangent = 0.0;
                    double attenuation_cotangent = 0.0;
                    if constexpr (WITH_TRANSPORT) {
                        double source_factor_cotangent = 0.0;
                        if (layer_state != nullptr) {
                            source_factor_cotangent =
                                layer_state[flat_layer] * forcing_gradient;
                        } else {
                            ray_interpolation.source_for_layer(layer).visit(
                                [&](const auto sources) {
                                    for (const auto& source : sources) {
                                        source_factor_cotangent +=
                                            source.weight() * forcing_gradient *
                                            (*transport_state)(
                                                transport_column_data
                                                    [source.row_inner_index()]);
                                    }
                                });
                        }
                        const double albedo = uniform_albedo
                                                  ? uniform_albedo_value
                                                  : scratch.albedo(layer);
                        const double albedo_cotangent =
                            source_factor_cotangent * (1.0 - attenuation) *
                            layer_prefix;
                        attenuation_cotangent -=
                            source_factor_cotangent * albedo * layer_prefix;
                        layer_prefix_cotangent += source_factor_cotangent *
                                                  albedo * (1.0 - attenuation);
                        for (const auto& weight :
                             ray_interpolation.atmosphere_for_layer(layer)) {
                            thread_gradient(
                                m_atmosphere->ssa_deriv_start_index() +
                                weight.index) +=
                                weight.weight() * albedo_cotangent;
                        }
                    }
                    if (layer_distance >= minimum_layer_distance_m &&
                        forcing_gradient != 0.0) {
                        sasktran2::raytracing::GridWeightStencilView
                            entrance_weights;
                        sasktran2::raytracing::GridWeightStencilView
                            exit_weights;
                        const auto* start_weights = &entrance_weights;
                        const auto* end_weights = &exit_weights;
                        bool start_entrance = true;
                        bool end_entrance = false;
                        if constexpr (LOWER_INTERPOLATION) {
                            entrance_weights =
                                rays[ray].entrance_weights(layer);
                            exit_weights = rays[ray].exit_weights(layer);
                            if (traced_layer->r_exit >
                                traced_layer->r_entrance) {
                                end_weights = &entrance_weights;
                                end_entrance = true;
                            } else {
                                start_weights = &exit_weights;
                                start_entrance = false;
                            }
                        }
                        const int exit_solar = m_solar_offsets[ray] + layer;
                        auto start = !LOWER_INTERPOLATION
                                         ? scratch.endpoints[layer + 1]
                                         : scalar_endpoint(
                                               ray_endpoint_context, layer,
                                               start_entrance, exit_solar + 1,
                                               *start_weights);
                        auto end =
                            !LOWER_INTERPOLATION
                                ? scratch.endpoints[layer]
                                : scalar_endpoint(ray_endpoint_context, layer,
                                                  end_entrance, exit_solar,
                                                  *end_weights);
                        if constexpr (LOWER_INTERPOLATION) {
                            start.solar_transmission = solar(exit_solar + 1);
                            start.source = start.extinction * start.albedo *
                                           start.solar_transmission *
                                           start.phase * inverse_four_pi;
                            end.solar_transmission = solar(exit_solar);
                            end.source = end.extinction * end.albedo *
                                         end.solar_transmission * end.phase *
                                         inverse_four_pi;
                        }
                        const double factor_derivative =
                            constant_source_factor_derivative(optical_depth,
                                                              factor);
                        const double endpoint_value = [&]() {
                            if constexpr (LOWER_INTERPOLATION) {
                                return (start.source *
                                            traced_layer
                                                ->od_quad_start_fraction +
                                        end.source *
                                            traced_layer
                                                ->od_quad_end_fraction) *
                                       layer_distance;
                            } else {
                                return start.source * source_quad_start +
                                       end.source * source_quad_end;
                            }
                        }();
                        layer_prefix_cotangent +=
                            forcing_gradient * factor * endpoint_value;
                        optical_depth_cotangent +=
                            forcing_gradient * layer_prefix *
                            factor_derivative *
                            (start.source * source_quad_start +
                             end.source * source_quad_end);
                        const double start_cotangent = forcing_gradient *
                                                       layer_prefix * factor *
                                                       source_quad_start;
                        const double end_cotangent = forcing_gradient *
                                                     layer_prefix * factor *
                                                     source_quad_end;
                        if constexpr (!LOWER_INTERPOLATION) {
                            scratch.endpoint_cotangent(layer + 1) +=
                                start_cotangent;
                            scratch.endpoint_cotangent(layer) += end_cotangent;
                        } else {
                            accumulate_scalar_endpoint_vjp<WITH_PHASE_GRADIENT>(
                                ray_endpoint_context, layer, start_entrance,
                                exit_solar + 1, *start_weights, start,
                                start_cotangent, thread_gradient,
                                solar_gradient, coefficient_gradient);
                            accumulate_scalar_endpoint_vjp<WITH_PHASE_GRADIENT>(
                                ray_endpoint_context, layer, end_entrance,
                                exit_solar, *end_weights, end, end_cotangent,
                                thread_gradient, solar_gradient,
                                coefficient_gradient);
                        }
                    }
                    layer_prefix_cotangent += prefix_cotangent * attenuation;
                    attenuation_cotangent += prefix_cotangent * layer_prefix;
                    optical_depth_cotangent -=
                        attenuation * attenuation_cotangent;
                    const auto od_weights =
                        ray_interpolation.optical_depth_for_layer(layer);
                    for (std::size_t index = 0; index < od_weights.size();
                         ++index) {
                        const auto [atmosphere_index, weight] =
                            od_weights[index];
                        thread_gradient(atmosphere_index) +=
                            weight * optical_depth_cotangent;
                    }
                    prefix_cotangent = layer_prefix_cotangent;
                }
                if constexpr (!LOWER_INTERPOLATION) {
                    if (num_layers > 0) {
                        accumulate_scalar_endpoint_vjp<WITH_PHASE_GRADIENT>(
                            ray_endpoint_context, 0, false,
                            m_solar_offsets[ray],
                            endpoint_weights(m_solar_offsets[ray]),
                            scratch.endpoints[0], scratch.endpoint_cotangent(0),
                            thread_gradient, solar_gradient,
                            coefficient_gradient);
                        for (int endpoint = 1; endpoint <= num_layers;
                             ++endpoint) {
                            const int solar_index =
                                m_solar_offsets[ray] + endpoint;
                            accumulate_scalar_endpoint_vjp<WITH_PHASE_GRADIENT>(
                                ray_endpoint_context, endpoint - 1, true,
                                solar_index, endpoint_weights(solar_index),
                                scratch.endpoints[endpoint],
                                scratch.endpoint_cotangent(endpoint),
                                thread_gradient, solar_gradient,
                                coefficient_gradient);
                        }
                    }
                }
            });
        }

        for (int thread = first_thread; thread < last_thread; ++thread) {
            auto thread_gradient = m_gradient_scratch[thread].col(0);
            const auto& coefficient_gradient = m_phase_product_scratch[thread];
            const int locations =
                m_atmosphere->storage().total_extinction.rows();
            if constexpr (WITH_PHASE_GRADIENT) {
                for (int group = 0;
                     group < m_atmosphere->num_scattering_deriv_groups();
                     ++group) {
                    const int native_offset =
                        m_atmosphere->scat_deriv_start_index() +
                        group * locations;
                    for (int location = 0; location < locations; ++location) {
                        const int order =
                            std::min(m_atmosphere->storage().d_max_order[group](
                                         location, wavelength),
                                     m_num_phase_moments);
                        double value = 0.0;
                        const int coefficient_offset =
                            location * m_num_phase_moments;
                        for (int degree = 0; degree < order; ++degree) {
                            value += coefficient_gradient(coefficient_offset +
                                                          degree) *
                                     m_atmosphere->storage().d_leg_coeff(
                                         degree, location, wavelength, group);
                        }
                        thread_gradient(native_offset + location) += value;
                    }
                }
            }
            auto solar_gradient = product_scratch.solar(thread - first_thread);
            auto table_gradient = product_scratch.table();
            solar_gradient.array() *= solar.array();
            m_solar_interpolation.apply_transpose(solar_gradient,
                                                  table_gradient);
            const Eigen::Index solar_columns = m_solar_table->atmosphere_size();
            m_solar_table->accumulate_transpose(
                table_gradient, thread_gradient.head(solar_columns), -1.0);
            native_gradient += thread_gradient;
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::calculate_with_transport(
        int wavelength, int wavelength_thread, TransportOperator& transport,
        Eigen::Ref<Eigen::VectorXd> forcing) {
        validate_ready(wavelength, wavelength_thread);
        if (!m_use_compact_scalar || forcing.size() != size()) {
            throw std::invalid_argument(
                "Fused successive-orders assembly requires the compact "
                "scalar kernel");
        }
        calculate_scalar(wavelength, wavelength_thread, forcing, &transport);
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::calculate_transport_only(
        int wavelength, int wavelength_thread, TransportOperator& transport) {
        validate_ready(wavelength, wavelength_thread);
        if (!m_use_compact_scalar ||
            transport.sparsity().rows() != m_num_rays ||
            transport.sparsity().columns() !=
                m_source_geometry->total_num_outgoing() ||
            transport.values().size() !=
                m_source_geometry->transport_sparsity().nonzeros()) {
            throw std::invalid_argument(
                "Transport-only successive-orders assembly requires matching "
                "compact scalar transport");
        }
        prepare_wavelength(wavelength, wavelength_thread);
        // The false specialization never reads or writes this empty vector.
        // JVP/VJP prepare solar and endpoint-medium values when needed.
        Eigen::VectorXd unused_forcing;
        if (m_use_lower_interpolation) {
            calculate_scalar_impl<true, true, false>(
                wavelength, wavelength_thread, unused_forcing, &transport);
        } else {
            calculate_scalar_impl<true, false, false>(
                wavelength, wavelength_thread, unused_forcing, &transport);
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::project_transport_state(
        Eigen::Ref<const Eigen::VectorXd> transport_state,
        Eigen::VectorXd& layer_state_projection,
        Eigen::VectorXd& ground_state_projection) const {
        if (!m_use_compact_scalar ||
            transport_state.size() != m_source_geometry->total_num_outgoing()) {
            throw std::invalid_argument(
                "Invalid compact scalar transport state projection");
        }
        const double* state = transport_state.data();
        const int layer_count = m_solar_offsets.back() - m_num_rays;
        layer_state_projection.setZero(layer_count);
        ground_state_projection.setZero(m_num_rays);
        const auto& interpolation = m_source_geometry->incoming_interpolation();
        if (!m_use_lower_interpolation) {
#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
            for (int ray = 0; ray < m_num_rays; ++ray) {
                const auto& packed_ray = m_scalar_packed_rays[ray];
                const auto& ray_interpolation = interpolation[ray];
                const auto transport_columns =
                    m_source_geometry->transport_columns_for_ray(ray);
                transport_columns.visit([&](const auto& typed_columns) {
                    const auto* columns = typed_columns.data();
                    for (std::uint32_t flat_layer = packed_ray.layer_begin;
                         flat_layer < packed_ray.layer_end; ++flat_layer) {
                        const auto local_layer = static_cast<std::size_t>(
                            flat_layer - packed_ray.layer_begin);
                        double value = 0.0;
                        ray_interpolation.source_for_layer(local_layer)
                            .visit([&](const auto sources) {
                                for (const auto& source : sources) {
                                    value +=
                                        source.weight() *
                                        state
                                            [columns[source.row_inner_index()]];
                                }
                            });
                        layer_state_projection(flat_layer) = value;
                    }
                    double ground_value = 0.0;
                    ray_interpolation.ground_weights.visit(
                        [&](const auto sources) {
                            for (const auto& source : sources) {
                                ground_value +=
                                    source.weight() *
                                    state[columns[source.row_inner_index()]];
                            }
                        });
                    ground_state_projection(ray) = ground_value;
                });
            }
            return;
        }
#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const auto& ray_interpolation = interpolation[ray];
            const auto transport_columns =
                m_source_geometry->transport_columns_for_ray(ray);
            transport_columns.visit([&](const auto& typed_columns) {
                const auto* columns = typed_columns.data();
                const int flat_layer_offset = m_solar_offsets[ray] - ray;
                for (int layer = 0;
                     layer < static_cast<int>(ray_interpolation.layers.size());
                     ++layer) {
                    double value = 0.0;
                    ray_interpolation.source_for_layer(layer).visit(
                        [&](const auto sources) {
                            for (const auto& source : sources) {
                                value +=
                                    source.weight() *
                                    state[columns[source.row_inner_index()]];
                            }
                        });
                    layer_state_projection(flat_layer_offset + layer) = value;
                }
                double ground_value = 0.0;
                ray_interpolation.ground_weights.visit([&](const auto sources) {
                    for (const auto& source : sources) {
                        ground_value +=
                            source.weight() *
                            state[columns[source.row_inner_index()]];
                    }
                });
                ground_state_projection(ray) = ground_value;
            });
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::calculate_jvp_with_transport(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        const Eigen::VectorXd& layer_state_projection,
        const Eigen::VectorXd& ground_state_projection,
        Eigen::VectorXd& direct_transport_tangent,
        Eigen::Ref<Eigen::VectorXd> forcing,
        Eigen::Ref<Eigen::VectorXd> forcing_tangent) {
        validate_ready(wavelength, wavelength_thread);
        if (!m_use_compact_scalar ||
            native_tangent.size() != m_atmosphere->num_deriv() ||
            layer_state_projection.size() !=
                m_solar_offsets.back() - m_num_rays ||
            ground_state_projection.size() != m_num_rays ||
            direct_transport_tangent.size() != m_num_rays ||
            forcing.size() != size() || forcing_tangent.size() != size()) {
            throw std::invalid_argument(
                "Fused successive-orders JVP requires the compact scalar "
                "kernel");
        }
        calculate_scalar_jvp(wavelength, wavelength_thread, native_tangent,
                             forcing, forcing_tangent, &layer_state_projection,
                             &ground_state_projection,
                             &direct_transport_tangent);
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::accumulate_vjp_with_transport(
        int wavelength, int wavelength_thread,
        const Eigen::VectorXd& transport_state,
        Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) {
        validate_ready(wavelength, wavelength_thread);
        if (!m_use_compact_scalar || forcing_cotangent.size() != size() ||
            native_gradient.size() != m_atmosphere->num_deriv() ||
            transport_state.size() != m_source_geometry->total_num_outgoing()) {
            throw std::invalid_argument(
                "Fused successive-orders VJP requires the compact scalar "
                "kernel");
        }
        accumulate_scalar_vjp(wavelength, wavelength_thread, forcing_cotangent,
                              native_gradient, &transport_state, nullptr,
                              nullptr);
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::accumulate_vjp_with_projected_transport(
        int wavelength, int wavelength_thread,
        const Eigen::VectorXd& layer_state_projection,
        const Eigen::VectorXd& ground_state_projection,
        Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) {
        validate_ready(wavelength, wavelength_thread);
        if (!m_use_compact_scalar || forcing_cotangent.size() != size() ||
            native_gradient.size() != m_atmosphere->num_deriv() ||
            layer_state_projection.size() !=
                m_solar_offsets.back() - m_num_rays ||
            ground_state_projection.size() != m_num_rays) {
            throw std::invalid_argument(
                "Fused successive-orders VJP requires the compact scalar "
                "kernel");
        }
        accumulate_scalar_vjp(wavelength, wavelength_thread, forcing_cotangent,
                              native_gradient, nullptr, &layer_state_projection,
                              &ground_state_projection);
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::calculate(
        int wavelength, int wavelength_thread,
        Eigen::Ref<Eigen::VectorXd> forcing) {
        validate_ready(wavelength, wavelength_thread);
        if (forcing.size() != size()) {
            throw std::invalid_argument(
                "Invalid successive-orders first-order forcing size");
        }
        if (m_use_compact_scalar) {
            calculate_scalar(wavelength, wavelength_thread, forcing);
            return;
        }
        const sasktran2::WavelengthBlock<> block{wavelength, 1};
        m_source.calculate(block, wavelength_thread);

#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const int thread = ray_thread_index(wavelength_thread);
            auto& result = m_primal_scratch[thread];
            result.set_zero(1);
            m_integrator.integrate(result, m_source_terms, block, ray,
                                   wavelength_thread, thread);
            forcing.template segment<NSTOKES>(ray * NSTOKES) =
                result.value.col(0);
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::calculate_jvp(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::Ref<Eigen::VectorXd> forcing,
        Eigen::Ref<Eigen::VectorXd> forcing_tangent) {
        validate_ready(wavelength, wavelength_thread);
        if (native_tangent.size() != m_atmosphere->num_deriv() ||
            forcing.size() != size() || forcing_tangent.size() != size()) {
            throw std::invalid_argument(
                "Invalid successive-orders first-order JVP dimensions");
        }
        if (m_use_compact_scalar) {
            calculate_scalar_jvp(wavelength, wavelength_thread, native_tangent,
                                 forcing, forcing_tangent);
            return;
        }
        const sasktran2::WavelengthBlock<> block{wavelength, 1};
        m_source.calculate_jvp(block, wavelength_thread, native_tangent);

#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const int thread = ray_thread_index(wavelength_thread);
            sasktran2::RadianceJVP<NSTOKES> result;
            m_integrator.integrate_jvp(result, m_source_terms, wavelength, ray,
                                       wavelength_thread, thread,
                                       native_tangent);
            forcing.template segment<NSTOKES>(ray * NSTOKES) = result.value;
            forcing_tangent.template segment<NSTOKES>(ray * NSTOKES) =
                result.jvp;
        }
    }

    template <int NSTOKES>
    void FirstOrderProvider<NSTOKES>::accumulate_vjp(
        int wavelength, int wavelength_thread,
        Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) {
        validate_ready(wavelength, wavelength_thread);
        if (forcing_cotangent.size() != size() ||
            native_gradient.size() != m_atmosphere->num_deriv()) {
            throw std::invalid_argument(
                "Invalid successive-orders first-order VJP dimensions");
        }
        // `prepare_primal` may return a cached wavelength state while the
        // exact source's wavelength-thread slot has since been reused. Refresh
        // its solar transmission before every reverse product.
        if (m_use_compact_scalar) {
            accumulate_scalar_vjp(wavelength, wavelength_thread,
                                  forcing_cotangent, native_gradient, nullptr,
                                  nullptr, nullptr);
            return;
        }
        m_source.calculate_vjp(sasktran2::WavelengthBlock<>{wavelength, 1},
                               wavelength_thread);
        // The engine uses disjoint scratch slots for simultaneous wavelength
        // workers.  Clearing the whole vector here would race another
        // wavelength's reverse pass and could erase a partially accumulated
        // gradient.  A call owns only its wavelength slot when source
        // threading is disabled, or the contiguous source-thread slots that
        // start at that wavelength slot otherwise.
        const int first_thread = wavelength_thread;
        const int last_thread = first_thread + m_num_source_threads;
        if (last_thread > static_cast<int>(m_gradient_scratch.size())) {
            throw std::logic_error(
                "Successive-orders first-order thread layout is invalid");
        }
        for (int thread = first_thread; thread < last_thread; ++thread) {
            m_gradient_scratch[thread].setZero();
        }

#pragma omp parallel for if (m_num_source_threads > 1)                         \
    num_threads(m_num_source_threads) schedule(dynamic)
        for (int ray = 0; ray < m_num_rays; ++ray) {
            const int thread = ray_thread_index(wavelength_thread);
            auto& radiance = m_vjp_radiance_scratch[thread];
            auto& cotangent = m_vjp_cotangent_scratch[thread];
            cotangent.col(0) =
                forcing_cotangent.template segment<NSTOKES>(ray * NSTOKES);
            m_integrator.integrate_vjp_block(
                radiance, m_source_terms,
                sasktran2::WavelengthBlock<>{wavelength, 1}, ray,
                wavelength_thread, thread, cotangent,
                m_gradient_scratch[thread]);
        }

        for (int thread = first_thread; thread < last_thread; ++thread) {
            native_gradient += m_gradient_scratch[thread].col(0);
        }
    }

    template <int NSTOKES>
    std::size_t FirstOrderProvider<NSTOKES>::workspace_bytes() const {
        std::size_t result = m_solar_interpolation.storage_bytes() +
                             m_endpoint_stencils.storage_bytes();
        if (m_solar_table != nullptr) {
            result += m_solar_table->storage_bytes();
        }
        for (const auto& scratch : m_primal_scratch) {
            result += static_cast<std::size_t>(scratch.value.size() +
                                               scratch.deriv.size()) *
                      sizeof(double);
        }
        for (const auto& gradient : m_gradient_scratch) {
            result +=
                static_cast<std::size_t>(gradient.size()) * sizeof(double);
        }
        for (const auto& scratch : m_vjp_radiance_scratch) {
            result += static_cast<std::size_t>(scratch.size()) * sizeof(double);
        }
        for (const auto& scratch : m_vjp_cotangent_scratch) {
            result += static_cast<std::size_t>(scratch.size()) * sizeof(double);
        }
        for (const auto& scratch : m_phase_product_scratch) {
            result += static_cast<std::size_t>(scratch.size()) * sizeof(double);
        }
        for (const auto& scratch : m_phase_order_scratch) {
            result += scratch.capacity() * sizeof(int);
        }
        result += m_scalar_phase_orders.capacity() * sizeof(int);
        result += m_uniform_phase_active.capacity() * sizeof(unsigned char);
        result +=
            m_uniform_phase_values.capacity() * sizeof(std::vector<double>);
        for (const auto& values : m_uniform_phase_values) {
            result += values.capacity() * sizeof(double);
        }
        result += m_uniform_albedo_active.capacity() * sizeof(unsigned char);
        result += m_uniform_albedo_values.capacity() * sizeof(double);
        result += (m_scalar_cache_wavelength.capacity() +
                   m_scalar_cache_index.capacity()) *
                  sizeof(int);
        result +=
            m_scalar_packed_rays.capacity() * sizeof(ScalarPackedRay) +
            m_scalar_packed_layers.capacity() * sizeof(ScalarPackedLayer) +
            m_scalar_ground_geometry.capacity() * sizeof(ScalarGroundGeometry);
        for (const auto& transmission : m_cached_solar_transmission) {
            result +=
                static_cast<std::size_t>(transmission.size()) * sizeof(double);
        }
        result += m_cached_solar_active.size() * sizeof(unsigned char);
        for (const auto& cache : m_scalar_layer_cache) {
            result += static_cast<std::size_t>(cache.optical_depth.size() +
                                               cache.source_factor.size()) *
                      sizeof(double);
        }
        for (const auto& cache : m_endpoint_medium_cache) {
            result += static_cast<std::size_t>(cache.extinction.size() +
                                               cache.albedo.size()) *
                      sizeof(double);
        }
        for (const auto& cache : m_scalar_volume_cache) {
            result += static_cast<std::size_t>(cache.ground_prefix.size()) *
                      sizeof(double);
        }
        for (const auto& scratch : m_scalar_vjp_scratch) {
            result +=
                static_cast<std::size_t>(
                    scratch.optical_depth.size() + scratch.attenuation.size() +
                    scratch.source_factor.size() +
                    scratch.prefix_attenuation.size() + scratch.albedo.size() +
                    scratch.endpoint_cotangent.size()) *
                sizeof(double);
            result += scratch.endpoints.size() * sizeof(ScalarEndpoint);
        }
        return result;
    }

    template class FirstOrderProvider<1>;
    template class FirstOrderProvider<3>;

} // namespace sasktran2::successive_orders
