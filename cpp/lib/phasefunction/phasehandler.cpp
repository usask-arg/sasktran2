#include "sasktran2/raytracing.h"
#include <sasktran2/solartransmission.h>
#include <spdlog/spdlog.h>

#include <limits>

namespace sasktran2::solartransmission {
    namespace {
        bool same_nonzero_indices(
            const sasktran2::raytracing::GridWeightStencilView& lhs,
            const sasktran2::raytracing::GridWeightStencilView& rhs) {
            std::size_t lhs_index = 0;
            std::size_t rhs_index = 0;
            while (true) {
                while (lhs_index < lhs.size() && lhs[lhs_index].second == 0.0) {
                    ++lhs_index;
                }
                while (rhs_index < rhs.size() && rhs[rhs_index].second == 0.0) {
                    ++rhs_index;
                }
                if (lhs_index == lhs.size() || rhs_index == rhs.size()) {
                    return lhs_index == lhs.size() && rhs_index == rhs.size();
                }
                if (lhs[lhs_index].first != rhs[rhs_index].first) {
                    return false;
                }
                ++lhs_index;
                ++rhs_index;
            }
        }
    } // namespace

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
        const std::vector<std::vector<int>>& index_map,
        const std::vector<Eigen::Vector3d>* solar_propagation_directions) {
        initialize_geometry_impl(los_rays, index_map,
                                 solar_propagation_directions);
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_geometry_impl(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
        const std::vector<std::vector<int>>& index_map,
        const std::vector<Eigen::Vector3d>* solar_propagation_directions) {

        ZoneScopedN("Phase Handler Initialize Geometry");

        m_scatter_angles.clear();
        m_internal_to_geometry.clear();
        m_internal_to_cos_scatter.clear();
        m_geometry_layer_offsets.assign(los_rays.size() + 1, 0);
        std::size_t total_layers = 0;
        std::size_t num_endpoint_references = 0;
        std::size_t num_phase_entries = 0;
        const bool endpoint_solar_geometry =
            solar_propagation_directions != nullptr;
        if (endpoint_solar_geometry) {
            std::size_t endpoint_count = 0;
            for (const auto& ray : los_rays) {
                endpoint_count += ray.layers.size() + 1;
            }
            if (solar_propagation_directions->size() != endpoint_count) {
                throw std::invalid_argument(
                    "Single-scatter solar directions do not match the ray "
                    "endpoint geometry");
            }
        }
        int counting_scatter_index = 0;
        std::vector<int> last_scatter_by_geometry(m_geometry.size(), -1);
        const auto count_endpoint =
            [&](const raytracing::GridWeightStencilView& weights) {
                for (std::size_t index = 0; index < weights.size(); ++index) {
                    const auto weight = weights[index];
                    if (weight.second == 0.0) {
                        continue;
                    }
                    ++num_endpoint_references;
                    if (last_scatter_by_geometry[weight.first] !=
                        counting_scatter_index) {
                        last_scatter_by_geometry[weight.first] =
                            counting_scatter_index;
                        ++num_phase_entries;
                    }
                }
            };
        for (std::size_t ray_index = 0; ray_index < los_rays.size();
             ++ray_index) {
            const auto& ray = los_rays[ray_index];
            m_geometry_layer_offsets[ray_index] =
                static_cast<std::uint32_t>(total_layers);
            total_layers += ray.layers.size();
            if (total_layers > std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(
                    "Single-scatter phase geometry exceeds the packed layer "
                    "index range");
            }

            if (ray.layers.empty()) {
                continue;
            }
            if (endpoint_solar_geometry) {
                for (std::size_t layer_index = 0;
                     layer_index < ray.layers.size(); ++layer_index) {
                    count_endpoint(ray.entrance_weights(layer_index));
                    ++counting_scatter_index;
                    count_endpoint(ray.exit_weights(layer_index));
                    ++counting_scatter_index;
                }
            } else {
                for (std::size_t layer_index = 0;
                     layer_index < ray.layers.size(); ++layer_index) {
                    count_endpoint(ray.entrance_weights(layer_index));
                    const auto exit_weights = ray.exit_weights(layer_index);
                    const bool can_share_exit =
                        layer_index > 0 &&
                        same_nonzero_indices(
                            exit_weights,
                            ray.entrance_weights(layer_index - 1));
                    if (!can_share_exit) {
                        count_endpoint(exit_weights);
                    }
                    if (!ray.is_straight) {
                        ++counting_scatter_index;
                    }
                }
                if (ray.is_straight) {
                    ++counting_scatter_index;
                }
            }
            if (num_endpoint_references >
                std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(
                    "Single-scatter phase geometry exceeds the packed "
                    "endpoint index range");
            }
            if (num_phase_entries >
                static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::length_error(
                    "Single-scatter phase geometry exceeds the internal "
                    "index range");
            }
        }
        m_geometry_layer_offsets.back() =
            static_cast<std::uint32_t>(total_layers);
        m_geometry_entrance_offsets.assign(total_layers, 0);
        m_geometry_exit_offsets.assign(total_layers, 0);
        m_geometry_to_internal.clear();
        m_geometry_to_internal.reserve(num_endpoint_references);
        m_internal_to_geometry.reserve(num_phase_entries);
        m_internal_to_cos_scatter.reserve(num_phase_entries);

        int num_internal = 0;
        int num_scatter = 0;
        std::fill(last_scatter_by_geometry.begin(),
                  last_scatter_by_geometry.end(), -1);
        std::vector<int> internal_by_geometry(m_geometry.size(), -1);
        const auto append_endpoint =
            [&](const raytracing::GridWeightStencilView& weights) {
                for (std::size_t index = 0; index < weights.size(); ++index) {
                    const auto weight = weights[index];
                    if (weight.second == 0.0) {
                        continue;
                    }
                    if (last_scatter_by_geometry[weight.first] != num_scatter) {
                        last_scatter_by_geometry[weight.first] = num_scatter;
                        internal_by_geometry[weight.first] = num_internal;
                        m_internal_to_geometry.push_back(weight.first);
                        m_internal_to_cos_scatter.push_back(num_scatter);
                        ++num_internal;
                    }
                    m_geometry_to_internal.push_back(
                        internal_by_geometry[weight.first]);
                }
            };
        double theta, C1, C2, S1, S2;
        int negation;
        const auto append_scatter_angle = [&](const Eigen::Vector3d&
                                                  solar_propagation,
                                              const sasktran2::raytracing::
                                                  TracedLayer& scatter_layer,
                                              const sasktran2::raytracing::
                                                  TracedRay& ray) {
            const auto result =
                m_geometry.coordinates().stokes_standard_to_observer_z(
                    scatter_layer.average_look_away,
                    ray.observer_and_look.observer.position);
            math::stokes_scattering_factors(solar_propagation.normalized(),
                                            -scatter_layer.average_look_away,
                                            theta, C1, C2, S1, S2, negation);
            if constexpr (NSTOKES == 3) {
                const double adjusted_C2 =
                    C2 * result.first - S2 * result.second;
                const double adjusted_S2 =
                    C2 * result.second + S2 * result.first;
                m_scatter_angles.push_back({theta, adjusted_C2, adjusted_S2});
            } else {
                m_scatter_angles.push_back({theta});
            }
        };
        // First we need to iterate through and figure out how many internal
        // indices we will end up with and how many scatter angles we will need

        for (int i = 0; i < los_rays.size(); ++i) {
            const auto& ray = los_rays[i];

            // Empty rays don't need to be considered
            if (ray.layers.size() == 0) {
                continue;
            }

            if (endpoint_solar_geometry) {
                for (int j = 0; j < ray.layers.size(); ++j) {
                    const auto flat_layer = m_geometry_layer_offsets[i] +
                                            static_cast<std::uint32_t>(j);
                    const int exit_solar_index = index_map[i][j];
                    const int entrance_solar_index = exit_solar_index + 1;

                    append_scatter_angle(
                        (*solar_propagation_directions)[entrance_solar_index],
                        ray.layers[j], ray);
                    m_geometry_entrance_offsets[flat_layer] =
                        static_cast<std::uint32_t>(
                            m_geometry_to_internal.size());
                    append_endpoint(ray.entrance_weights(j));
                    ++num_scatter;

                    append_scatter_angle(
                        (*solar_propagation_directions)[exit_solar_index],
                        ray.layers[j], ray);
                    m_geometry_exit_offsets[flat_layer] =
                        static_cast<std::uint32_t>(
                            m_geometry_to_internal.size());
                    append_endpoint(ray.exit_weights(j));
                    ++num_scatter;
                }
                continue;
            }

            // Straight rays we can just use the first layer to get the
            // scattering angle, scattering angle does not change along the ray
            if (ray.is_straight) {
                auto result =
                    m_geometry.coordinates().stokes_standard_to_observer_z(
                        ray.layers[0].average_look_away,
                        ray.observer_and_look.observer.position);

                math::stokes_scattering_factors(
                    -1 * m_geometry.coordinates().sun_unit(),
                    -1 * ray.layers[0].average_look_away, theta, C1, C2, S1, S2,
                    negation);
                if constexpr (NSTOKES == 3) {
                    double adjusted_C2 = C2 * result.first - S2 * result.second;
                    double adjusted_S2 = C2 * result.second + S2 * result.first;

                    m_scatter_angles.push_back(
                        {theta, adjusted_C2, adjusted_S2});
                } else {
                    m_scatter_angles.push_back({theta});
                }
            }

            for (int j = 0; j < ray.layers.size(); ++j) {
                const auto flat_layer =
                    m_geometry_layer_offsets[i] + static_cast<std::uint32_t>(j);

                // If the ray isn't straight every layer has a scattering angle
                if (!ray.is_straight) {
                    const auto& scatter_layer = ray.layers[0];
                    auto result =
                        m_geometry.coordinates().stokes_standard_to_observer_z(
                            scatter_layer.average_look_away,
                            ray.observer_and_look.observer.position);

                    math::stokes_scattering_factors(
                        -1 * m_geometry.coordinates().sun_unit(),
                        -1 * scatter_layer.average_look_away, theta, C1, C2, S1,
                        S2, negation);
                    if constexpr (NSTOKES == 3) {
                        double adjusted_C2 =
                            C2 * result.first - S2 * result.second;
                        double adjusted_S2 =
                            C2 * result.second + S2 * result.first;

                        m_scatter_angles.push_back(
                            {theta, adjusted_C2, adjusted_S2});
                    } else {
                        m_scatter_angles.push_back({theta});
                    }
                }

                const auto entrance_weights = ray.entrance_weights(j);
                const auto exit_weights = ray.exit_weights(j);
                auto& entrance_offset = m_geometry_entrance_offsets[flat_layer];
                entrance_offset =
                    static_cast<std::uint32_t>(m_geometry_to_internal.size());

                append_endpoint(entrance_weights);

                const bool can_share_exit =
                    j > 0 && same_nonzero_indices(exit_weights,
                                                  ray.entrance_weights(j - 1));

                auto& exit_offset = m_geometry_exit_offsets[flat_layer];
                if (!can_share_exit) {
                    // End layer at TOA, need to use layer exit
                    exit_offset = static_cast<std::uint32_t>(
                        m_geometry_to_internal.size());
                    append_endpoint(exit_weights);
                } else {
                    // Assign to previous layers entrance
                    exit_offset = m_geometry_entrance_offsets[flat_layer - 1];
                }

                if (!ray.is_straight) {
                    ++num_scatter;
                }
            }
            if (ray.is_straight) {
                ++num_scatter;
            }
        }

        spdlog::debug(
            "Single-scatter phase geometry: {} layers, {} endpoint references "
            "(capacity {}), {} phase entries, {} scatter angles",
            total_layers, m_geometry_to_internal.size(),
            m_geometry_to_internal.capacity(), m_internal_to_geometry.size(),
            m_scatter_angles.size());

        if constexpr (NSTOKES == 3) {
            m_phase.resize(2, (int)m_internal_to_geometry.size(),
                           m_config->num_wavelength_threads());
        } else {
            m_phase.resize(1, (int)m_internal_to_geometry.size(),
                           m_config->num_wavelength_threads());
        }

        m_wigner_d00.resize(m_config->num_singlescatter_moments(),
                            m_scatter_angles.size());
        auto d00 = sasktran2::math::WignerDCalculator(0, 0);
        auto d02 = sasktran2::math::WignerDCalculator(0, 2);

        if constexpr (NSTOKES == 3) {
            m_wigner_d02.resize(m_config->num_singlescatter_moments(),
                                m_scatter_angles.size());
        }

        for (int i = 0; i < m_scatter_angles.size(); ++i) {
            ZoneScopedN("Phase Handler SS Wigner Calculation");
            d00.vec_d_emplace(m_scatter_angles[i][0],
                              m_config->num_singlescatter_moments(),
                              &m_wigner_d00(0, i));

            if constexpr (NSTOKES == 3) {
                d02.vec_d_emplace(m_scatter_angles[i][0],
                                  m_config->num_singlescatter_moments(),
                                  &m_wigner_d02(0, i));
            }
        }
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_atmosphere = &atmosphere;

        int numderiv = m_atmosphere->num_scattering_deriv_groups();

        if constexpr (NSTOKES == 1) {
            m_d_phase.resize(1, (int)m_internal_to_geometry.size(), numderiv,
                             m_config->num_wavelength_threads());
        } else {
            m_d_phase.resize(2, (int)m_internal_to_geometry.size(), numderiv,
                             m_config->num_wavelength_threads());
        }
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::calculate(int wavelidx, int threadidx) {

        // If we need to calculate the phase function, we do it
        if (m_config->singlescatter_phasemode() ==
            sasktran2::Config::SingleScatterPhaseMode::from_legendre) {
            // Set all elements with this threadidx to zero
            m_phase.chip(threadidx, 2).setZero();
            m_d_phase.chip(threadidx, 3).setZero();

            Eigen::array<Eigen::Index, 3> dims = m_phase.dimensions();
            Eigen::array<Eigen::Index, 3> offsets = {0, 1, 0};
            Eigen::array<Eigen::Index, 3> extents = {dims[0], 1, 1};

            for (int i = 0; i < m_internal_to_geometry.size(); ++i) {
                int scat_index = m_internal_to_cos_scatter[i];

                int atmo_index = m_internal_to_geometry[i];

                int max_order =
                    m_atmosphere->storage().max_order(atmo_index, wavelidx);

                if constexpr (NSTOKES == 1) {
                    m_phase(0, i, threadidx) =
                        Eigen::Map<const Eigen::VectorXd>(
                            &m_atmosphere->storage().leg_coeff(0, atmo_index,
                                                               wavelidx),
                            max_order)
                            .dot(m_wigner_d00(Eigen::seq(0, max_order - 1),
                                              scat_index));
                } else {
                    m_phase(0, i, threadidx) =
                        Eigen::Map<const Eigen::VectorXd, 0,
                                   Eigen::InnerStride<4>>(
                            &m_atmosphere->storage().leg_coeff(0, atmo_index,
                                                               wavelidx),
                            max_order)
                            .dot(m_wigner_d00(Eigen::seq(0, max_order - 1),
                                              scat_index));

                    m_phase(1, i, threadidx) =
                        Eigen::Map<const Eigen::VectorXd, 0,
                                   Eigen::InnerStride<4>>(
                            &m_atmosphere->storage().leg_coeff(3, atmo_index,
                                                               wavelidx),
                            max_order)
                            .dot(m_wigner_d02(Eigen::seq(0, max_order - 1),
                                              scat_index));
                }
                for (int d = 0; d < m_atmosphere->num_scattering_deriv_groups();
                     ++d) {

                    int d_max_order = m_atmosphere->storage().d_max_order[d](
                        atmo_index, wavelidx);
                    if constexpr (NSTOKES == 1) {
                        m_d_phase(0, i, d, threadidx) =
                            Eigen::Map<const Eigen::VectorXd>(
                                &m_atmosphere->storage().d_leg_coeff(
                                    0, atmo_index, wavelidx, d),
                                d_max_order)
                                .dot(
                                    m_wigner_d00(Eigen::seq(0, d_max_order - 1),
                                                 scat_index));
                    } else {
                        m_d_phase(0, i, d, threadidx) =
                            Eigen::Map<const Eigen::VectorXd, 0,
                                       Eigen::InnerStride<4>>(
                                &m_atmosphere->storage().d_leg_coeff(
                                    0, atmo_index, wavelidx, d),
                                d_max_order)
                                .dot(
                                    m_wigner_d00(Eigen::seq(0, d_max_order - 1),
                                                 scat_index));
                        m_d_phase(1, i, d, threadidx) =
                            Eigen::Map<const Eigen::VectorXd, 0,
                                       Eigen::InnerStride<4>>(
                                &m_atmosphere->storage().d_leg_coeff(
                                    3, atmo_index, wavelidx, d),
                                d_max_order)
                                .dot(
                                    m_wigner_d02(Eigen::seq(0, d_max_order - 1),
                                                 scat_index));
                    }
                }
            }
        } else {
            spdlog::error("Phase mode not implemented");
            throw std::runtime_error("Phase mode not implemented");
        }
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_wavelength_blocks(int batch_size) {
        m_wavelength_batch_capacity = batch_size;
        const int num_phase_components = NSTOKES == 1 ? 1 : 2;
        const int num_internal =
            static_cast<int>(m_internal_to_geometry.size());
        const int num_derivatives =
            m_atmosphere == nullptr
                ? 0
                : m_atmosphere->num_scattering_deriv_groups();

        m_phase_batch.resize(m_config->num_wavelength_threads());
        m_d_phase_batch.resize(m_config->num_wavelength_threads());
        m_active_wavelength_block_start.assign(
            m_config->num_wavelength_threads(), 0);
        m_active_wavelength_block_count.assign(
            m_config->num_wavelength_threads(), 0);
        for (int threadidx = 0; threadidx < m_config->num_wavelength_threads();
             ++threadidx) {
            m_phase_batch[threadidx].resize(num_phase_components * num_internal,
                                            batch_size);
            m_d_phase_batch[threadidx].resize(
                num_derivatives * num_phase_components * num_internal,
                batch_size);
        }
    }

    template <int NSTOKES>
    template <int N>
    void PhaseHandler<NSTOKES>::calculate_block(
        const sasktran2::WavelengthBlock<N>& batch, int threadidx) {
        if (m_config->singlescatter_phasemode() !=
            sasktran2::Config::SingleScatterPhaseMode::from_legendre) {
            throw std::runtime_error("Phase mode not implemented");
        }
        if (batch.count > m_wavelength_batch_capacity) {
            throw std::invalid_argument(
                "Wavelength batch exceeds phase storage capacity");
        }
        m_active_wavelength_block_start[threadidx] = batch.start;
        m_active_wavelength_block_count[threadidx] = batch.count;

        auto& phase = m_phase_batch[threadidx];
        auto& d_phase = m_d_phase_batch[threadidx];
        wavelength_left_cols(phase, batch).setZero();
        wavelength_left_cols(d_phase, batch).setZero();

        const int num_internal =
            static_cast<int>(m_internal_to_geometry.size());
        const int num_phase_components = NSTOKES == 1 ? 1 : 2;
        const auto phase_row = [num_internal](int component,
                                              int internal_index) {
            return component * num_internal + internal_index;
        };
        const auto derivative_row =
            [num_internal, num_phase_components](int derivative, int component,
                                                 int internal_index) {
                return (derivative * num_phase_components + component) *
                           num_internal +
                       internal_index;
            };

        for (int internal_index = 0; internal_index < num_internal;
             ++internal_index) {
            const int scatter_index = m_internal_to_cos_scatter[internal_index];
            const int atmosphere_index = m_internal_to_geometry[internal_index];

            for (int lane = 0; lane < batch.count; ++lane) {
                const int wavelength = batch.wavelength(lane);
                const int max_order = m_atmosphere->storage().max_order(
                    atmosphere_index, wavelength);

                if constexpr (NSTOKES == 1) {
                    phase(phase_row(0, internal_index), lane) =
                        Eigen::Map<const Eigen::VectorXd>(
                            &m_atmosphere->storage().leg_coeff(
                                0, atmosphere_index, wavelength),
                            max_order)
                            .dot(m_wigner_d00(Eigen::seq(0, max_order - 1),
                                              scatter_index));
                } else {
                    phase(phase_row(0, internal_index), lane) =
                        Eigen::Map<const Eigen::VectorXd, 0,
                                   Eigen::InnerStride<4>>(
                            &m_atmosphere->storage().leg_coeff(
                                0, atmosphere_index, wavelength),
                            max_order)
                            .dot(m_wigner_d00(Eigen::seq(0, max_order - 1),
                                              scatter_index));
                    phase(phase_row(1, internal_index), lane) =
                        Eigen::Map<const Eigen::VectorXd, 0,
                                   Eigen::InnerStride<4>>(
                            &m_atmosphere->storage().leg_coeff(
                                3, atmosphere_index, wavelength),
                            max_order)
                            .dot(m_wigner_d02(Eigen::seq(0, max_order - 1),
                                              scatter_index));
                }

                for (int derivative = 0;
                     derivative < m_atmosphere->num_scattering_deriv_groups();
                     ++derivative) {
                    const int derivative_max_order =
                        m_atmosphere->storage().d_max_order[derivative](
                            atmosphere_index, wavelength);
                    if constexpr (NSTOKES == 1) {
                        d_phase(derivative_row(derivative, 0, internal_index),
                                lane) =
                            Eigen::Map<const Eigen::VectorXd>(
                                &m_atmosphere->storage().d_leg_coeff(
                                    0, atmosphere_index, wavelength,
                                    derivative),
                                derivative_max_order)
                                .dot(m_wigner_d00(
                                    Eigen::seq(0, derivative_max_order - 1),
                                    scatter_index));
                    } else {
                        d_phase(derivative_row(derivative, 0, internal_index),
                                lane) =
                            Eigen::Map<const Eigen::VectorXd, 0,
                                       Eigen::InnerStride<4>>(
                                &m_atmosphere->storage().d_leg_coeff(
                                    0, atmosphere_index, wavelength,
                                    derivative),
                                derivative_max_order)
                                .dot(m_wigner_d00(
                                    Eigen::seq(0, derivative_max_order - 1),
                                    scatter_index));
                        d_phase(derivative_row(derivative, 1, internal_index),
                                lane) =
                            Eigen::Map<const Eigen::VectorXd, 0,
                                       Eigen::InnerStride<4>>(
                                &m_atmosphere->storage().d_leg_coeff(
                                    3, atmosphere_index, wavelength,
                                    derivative),
                                derivative_max_order)
                                .dot(m_wigner_d02(
                                    Eigen::seq(0, derivative_max_order - 1),
                                    scatter_index));
                    }
                }
            }
        }
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::scatter(
        int threadidx, int losidx, int layeridx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        scatter_impl(threadidx, losidx, layeridx, index_weights, is_entrance,
                     source);
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::scatter_impl(
        int threadidx, int losidx, int layeridx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        const double source_amplitude = source.value(0);
        const auto phase = scatter_and_accumulate_derivative(
            threadidx, losidx, layeridx, index_weights, is_entrance,
            source_amplitude, 1.0, source);
        source.value = source_amplitude * phase;
    }

    template <int NSTOKES>
    Eigen::Vector<double, NSTOKES>
    PhaseHandler<NSTOKES>::scatter_and_accumulate_derivative(
        int threadidx, int losidx, int layeridx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, double source_amplitude, double derivative_scale,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& target)
        const {
        return scatter_and_accumulate_derivative_impl(
            threadidx, losidx, layeridx, index_weights, is_entrance,
            source_amplitude, derivative_scale, target);
    }

    template <int NSTOKES>
    Eigen::Vector<double, NSTOKES>
    PhaseHandler<NSTOKES>::scatter_and_accumulate_derivative(
        int threadidx, int losidx, int layeridx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, double source_amplitude, double derivative_scale,
        sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& target) const {
        return scatter_and_accumulate_derivative_impl(
            threadidx, losidx, layeridx, index_weights, is_entrance,
            source_amplitude, derivative_scale, target);
    }

    template <int NSTOKES>
    template <typename Target>
    Eigen::Vector<double, NSTOKES>
    PhaseHandler<NSTOKES>::scatter_and_accumulate_derivative_impl(
        int threadidx, int losidx, int layeridx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, double source_amplitude, double derivative_scale,
        Target& target) const {
        const int* internal_indices =
            geometry_internal_indices(losidx, layeridx, is_entrance);

        if constexpr (NSTOKES == 1) {
            double phase_result = 0.0;

            if (index_weights.size() == 2) {
                const auto lower = index_weights[0];
                const auto upper = index_weights[1];
                if (lower.second == 0.0 || upper.second == 0.0) {
                    phase_result = m_phase(0, internal_indices[0], threadidx);
                } else {
                    phase_result = m_phase(0, internal_indices[0], threadidx) *
                                       lower.second +
                                   m_phase(0, internal_indices[1], threadidx) *
                                       upper.second;
                }
            } else {
                int internal_offset = 0;
                for (std::size_t index = 0; index < index_weights.size();
                     ++index) {
                    const auto index_weight = index_weights[index];
                    if (index_weight.second == 0.0) {
                        continue;
                    }
                    const int internal_index =
                        internal_indices[internal_offset++];
                    phase_result += m_phase(0, internal_index, threadidx) *
                                    index_weight.second;
                }
            }

            for (int d = 0; d < m_atmosphere->num_scattering_deriv_groups();
                 ++d) {
                const int derivative_start =
                    m_atmosphere->scat_deriv_start_index() +
                    d * m_geometry.size();
                int derivative_internal_offset = 0;
                for (std::size_t index = 0; index < index_weights.size();
                     ++index) {
                    const auto index_weight = index_weights[index];
                    if (index_weight.second == 0.0) {
                        continue;
                    }
                    const int internal_index =
                        internal_indices[derivative_internal_offset++];
                    target.deriv(0, derivative_start + index_weight.first) +=
                        derivative_scale * source_amplitude *
                        m_d_phase(0, internal_index, d, threadidx) *
                        index_weight.second;
                }
            }
            return Eigen::Vector<double, 1>::Constant(phase_result);

        } else if constexpr (NSTOKES == 3) {
            Eigen::Vector3d phase_result = Eigen::Vector3d::Zero();

            const auto accumulate_phase = [&](int internal_index,
                                              double weight) {
                const auto& scatter_angle =
                    m_scatter_angles[m_internal_to_cos_scatter[internal_index]];
                const double C2 = scatter_angle[1];
                const double S2 = scatter_angle[2];

                const double off_diag =
                    m_phase(1, internal_index, threadidx) * weight;

                phase_result(0) +=
                    m_phase(0, internal_index, threadidx) * weight;
                phase_result(1) += -C2 * off_diag;
                phase_result(2) += -S2 * off_diag;
            };

            if (index_weights.size() == 2) {
                const auto lower = index_weights[0];
                const auto upper = index_weights[1];
                if (lower.second == 0.0 || upper.second == 0.0) {
                    accumulate_phase(internal_indices[0], 1.0);
                } else {
                    accumulate_phase(internal_indices[0], lower.second);
                    accumulate_phase(internal_indices[1], upper.second);
                }
            } else {
                int phase_internal_offset = 0;
                for (std::size_t index = 0; index < index_weights.size();
                     ++index) {
                    const auto index_weight = index_weights[index];
                    if (index_weight.second == 0.0) {
                        continue;
                    }
                    accumulate_phase(internal_indices[phase_internal_offset++],
                                     index_weight.second);
                }
            }

            for (int d = 0; d < m_atmosphere->num_scattering_deriv_groups();
                 ++d) {
                const auto accumulate_derivative = [&](int internal_index,
                                                       int geometry_index,
                                                       double weight) {
                    const auto& scatter_angle = m_scatter_angles
                        [m_internal_to_cos_scatter[internal_index]];
                    const double C2 = scatter_angle[1];
                    const double S2 = scatter_angle[2];

                    int deriv_index = m_atmosphere->scat_deriv_start_index() +
                                      d * m_geometry.size() + geometry_index;
                    target.deriv(0, deriv_index) +=
                        derivative_scale * source_amplitude *
                        m_d_phase(0, internal_index, d, threadidx) * weight;

                    target.deriv(1, deriv_index) +=
                        derivative_scale * source_amplitude *
                        m_d_phase(1, internal_index, d, threadidx) * weight *
                        (-C2);

                    target.deriv(2, deriv_index) +=
                        derivative_scale * source_amplitude *
                        m_d_phase(1, internal_index, d, threadidx) * weight *
                        (-S2);
                };

                int derivative_internal_offset = 0;
                for (std::size_t index = 0; index < index_weights.size();
                     ++index) {
                    const auto index_weight = index_weights[index];
                    if (index_weight.second == 0.0) {
                        continue;
                    }
                    accumulate_derivative(
                        internal_indices[derivative_internal_offset++],
                        index_weight.first, index_weight.second);
                }
            }

            return phase_result;
        } else {
            return Eigen::Vector<double, NSTOKES>::Zero();
        }
    }

    template <int NSTOKES>
    Eigen::Vector<double, NSTOKES> PhaseHandler<NSTOKES>::scatter_value(
        int threadidx, int losidx, int layeridx, int wavelidx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance) const {
        Eigen::Vector<double, NSTOKES> phase =
            Eigen::Vector<double, NSTOKES>::Zero();
        const int* internal_indices =
            geometry_internal_indices(losidx, layeridx, is_entrance);

        const int batch_lane =
            m_wavelength_batch_capacity > 1
                ? wavelidx - m_active_wavelength_block_start.at(threadidx)
                : 0;
        if (m_wavelength_batch_capacity > 1 &&
            (batch_lane < 0 ||
             batch_lane >= m_active_wavelength_block_count.at(threadidx))) {
            throw std::out_of_range(
                "Phase wavelength is outside the active block");
        }
        const int num_internal =
            static_cast<int>(m_internal_to_geometry.size());
        const auto phase_value = [&](int component, int internal_index) {
            if (m_wavelength_batch_capacity > 1) {
                return m_phase_batch[threadidx](
                    component * num_internal + internal_index, batch_lane);
            }
            return m_phase(component, internal_index, threadidx);
        };
        const auto accumulate_phase = [&](int internal_index, double weight) {
            phase(0) += phase_value(0, internal_index) * weight;
            if constexpr (NSTOKES == 3) {
                const auto& scatter_angle =
                    m_scatter_angles[m_internal_to_cos_scatter[internal_index]];
                const double polarized =
                    phase_value(1, internal_index) * weight;
                phase(1) -= scatter_angle[1] * polarized;
                phase(2) -= scatter_angle[2] * polarized;
            }
        };

        if (index_weights.size() == 2 && (index_weights[0].second == 0.0 ||
                                          index_weights[1].second == 0.0)) {
            accumulate_phase(internal_indices[0], 1.0);
        } else {
            int internal_offset = 0;
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto weight = index_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                accumulate_phase(internal_indices[internal_offset++],
                                 weight.second);
            }
        }

        return phase;
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::scatter_jvp(
        int threadidx, int losidx, int layeridx, int wavelidx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::Vector<double, NSTOKES>& phase,
        Eigen::Vector<double, NSTOKES>& phase_jvp) const {
        phase = scatter_value(threadidx, losidx, layeridx, wavelidx,
                              index_weights, is_entrance);
        phase_jvp.setZero();
        const int* internal_indices =
            geometry_internal_indices(losidx, layeridx, is_entrance);

        const int batch_lane =
            m_wavelength_batch_capacity > 1
                ? wavelidx - m_active_wavelength_block_start.at(threadidx)
                : 0;
        if (m_wavelength_batch_capacity > 1 &&
            (batch_lane < 0 ||
             batch_lane >= m_active_wavelength_block_count.at(threadidx))) {
            throw std::out_of_range(
                "Phase wavelength is outside the active block");
        }
        const int num_internal =
            static_cast<int>(m_internal_to_geometry.size());
        const int num_phase_components = NSTOKES == 1 ? 1 : 2;
        const auto derivative_value = [&](int derivative, int component,
                                          int internal_index) {
            if (m_wavelength_batch_capacity > 1) {
                const int row =
                    (derivative * num_phase_components + component) *
                        num_internal +
                    internal_index;
                return m_d_phase_batch[threadidx](row, batch_lane);
            }
            return m_d_phase(component, internal_index, derivative, threadidx);
        };

        for (int derivative = 0;
             derivative < m_atmosphere->num_scattering_deriv_groups();
             ++derivative) {
            int internal_offset = 0;
            const int derivative_start =
                m_atmosphere->scat_deriv_start_index() +
                derivative * m_geometry.size();
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto weight = index_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                const int internal_index = internal_indices[internal_offset++];
                const double direction =
                    native_tangent(derivative_start + weight.first) *
                    weight.second;
                phase_jvp(0) +=
                    direction * derivative_value(derivative, 0, internal_index);
                if constexpr (NSTOKES == 3) {
                    const auto& scatter_angle = m_scatter_angles
                        [m_internal_to_cos_scatter[internal_index]];
                    const double polarized =
                        direction *
                        derivative_value(derivative, 1, internal_index);
                    phase_jvp(1) -= scatter_angle[1] * polarized;
                    phase_jvp(2) -= scatter_angle[2] * polarized;
                }
            }
        }
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::scatter_vjp(
        int threadidx, int losidx, int layeridx, int wavelidx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, const Eigen::Vector<double, NSTOKES>& phase_cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        const int* internal_indices =
            geometry_internal_indices(losidx, layeridx, is_entrance);
        const int batch_lane =
            m_wavelength_batch_capacity > 1
                ? wavelidx - m_active_wavelength_block_start.at(threadidx)
                : 0;
        if (m_wavelength_batch_capacity > 1 &&
            (batch_lane < 0 ||
             batch_lane >= m_active_wavelength_block_count.at(threadidx))) {
            throw std::out_of_range(
                "Phase wavelength is outside the active block");
        }
        const int num_internal =
            static_cast<int>(m_internal_to_geometry.size());
        const int num_phase_components = NSTOKES == 1 ? 1 : 2;
        const auto derivative_phase = [&](int derivative, int component,
                                          int internal_index) {
            if (m_wavelength_batch_capacity > 1) {
                const int row =
                    (derivative * num_phase_components + component) *
                        num_internal +
                    internal_index;
                return m_d_phase_batch[threadidx](row, batch_lane);
            }
            return m_d_phase(component, internal_index, derivative, threadidx);
        };
        for (int derivative = 0;
             derivative < m_atmosphere->num_scattering_deriv_groups();
             ++derivative) {
            int internal_offset = 0;
            const int derivative_start =
                m_atmosphere->scat_deriv_start_index() +
                derivative * m_geometry.size();
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto weight = index_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                const int internal_index = internal_indices[internal_offset++];
                double derivative_value =
                    phase_cotangent(0) *
                    derivative_phase(derivative, 0, internal_index);
                if constexpr (NSTOKES == 3) {
                    const auto& scatter_angle = m_scatter_angles
                        [m_internal_to_cos_scatter[internal_index]];
                    derivative_value +=
                        (-scatter_angle[1] * phase_cotangent(1) -
                         scatter_angle[2] * phase_cotangent(2)) *
                        derivative_phase(derivative, 1, internal_index);
                }
                native_gradient(derivative_start + weight.first) +=
                    weight.second * derivative_value;
            }
        }
    }

    template <int NSTOKES>
    template <int N>
    void PhaseHandler<NSTOKES>::scatter_and_accumulate_derivative_block(
        int threadidx, int losidx, int layeridx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, const sasktran2::WavelengthBlock<N>& batch,
        const Eigen::Ref<const Eigen::Matrix<double, 1, N, Eigen::RowMajor>>&
            source_amplitude,
        const Eigen::Ref<const Eigen::Matrix<double, 1, N, Eigen::RowMajor>>&
            derivative_scale,
        sasktran2::WavelengthBlockDual<NSTOKES>& target,
        Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>&
            phase_result) const {
        wavelength_left_cols(phase_result, batch).setZero();

        const auto& phase = m_phase_batch[threadidx];
        const auto& d_phase = m_d_phase_batch[threadidx];
        const int num_internal =
            static_cast<int>(m_internal_to_geometry.size());
        const int num_phase_components = NSTOKES == 1 ? 1 : 2;
        const auto phase_row = [num_internal](int component,
                                              int internal_index) {
            return component * num_internal + internal_index;
        };
        const auto derivative_row =
            [num_internal, num_phase_components](int derivative, int component,
                                                 int internal_index) {
                return (derivative * num_phase_components + component) *
                           num_internal +
                       internal_index;
            };
        const int* internal_indices =
            geometry_internal_indices(losidx, layeridx, is_entrance);

        const auto accumulate_phase = [&](int internal_index, double weight) {
            wavelength_head(phase_result.row(0), batch).array() +=
                weight *
                wavelength_head(phase.row(phase_row(0, internal_index)), batch)
                    .array();
            if constexpr (NSTOKES == 3) {
                const auto& scatter_angle =
                    m_scatter_angles[m_internal_to_cos_scatter[internal_index]];
                const auto off_diagonal =
                    weight * wavelength_head(
                                 phase.row(phase_row(1, internal_index)), batch)
                                 .array();
                wavelength_head(phase_result.row(1), batch).array() -=
                    scatter_angle[1] * off_diagonal;
                wavelength_head(phase_result.row(2), batch).array() -=
                    scatter_angle[2] * off_diagonal;
            }
        };

        if (index_weights.size() == 2 && (index_weights[0].second == 0.0 ||
                                          index_weights[1].second == 0.0)) {
            accumulate_phase(internal_indices[0], 1.0);
        } else {
            int internal_offset = 0;
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto weight = index_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                accumulate_phase(internal_indices[internal_offset++],
                                 weight.second);
            }
        }

        for (int derivative = 0;
             target.derivative_size() > 0 &&
             derivative < m_atmosphere->num_scattering_deriv_groups();
             ++derivative) {
            int internal_offset = 0;
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto weight = index_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                const int internal_index = internal_indices[internal_offset++];
                const int derivative_index =
                    m_atmosphere->scat_deriv_start_index() +
                    derivative * m_geometry.size() + weight.first;
                auto target_derivative =
                    target.derivative(derivative_index, batch);
                const auto common_factor = weight.second *
                                           derivative_scale.array() *
                                           source_amplitude.array();
                target_derivative.row(0).array() +=
                    common_factor *
                    wavelength_head(d_phase.row(derivative_row(derivative, 0,
                                                               internal_index)),
                                    batch)
                        .array();
                if constexpr (NSTOKES == 3) {
                    const auto& scatter_angle = m_scatter_angles
                        [m_internal_to_cos_scatter[internal_index]];
                    const auto polarized =
                        common_factor *
                        wavelength_head(d_phase.row(derivative_row(
                                            derivative, 1, internal_index)),
                                        batch)
                            .array();
                    target_derivative.row(1).array() -=
                        scatter_angle[1] * polarized;
                    target_derivative.row(2).array() -=
                        scatter_angle[2] * polarized;
                }
            }
        }
    }

#define SASKTRAN2_INSTANTIATE_PHASE_BLOCK(NSTOKES, N)                          \
    template void PhaseHandler<NSTOKES>::calculate_block<N>(                   \
        const sasktran2::WavelengthBlock<N>&, int);                            \
    template void                                                              \
    PhaseHandler<NSTOKES>::scatter_and_accumulate_derivative_block<N>(         \
        int, int, int, const raytracing::GridWeightStencilView&, bool,         \
        const sasktran2::WavelengthBlock<N>&,                                  \
        const Eigen::Ref<const Eigen::Matrix<double, 1, N, Eigen::RowMajor>>&, \
        const Eigen::Ref<const Eigen::Matrix<double, 1, N, Eigen::RowMajor>>&, \
        sasktran2::WavelengthBlockDual<NSTOKES>&,                              \
        Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>&)      \
        const;

    SASKTRAN2_INSTANTIATE_PHASE_BLOCK(1, Eigen::Dynamic)
    SASKTRAN2_INSTANTIATE_PHASE_BLOCK(1, 1)
    SASKTRAN2_INSTANTIATE_PHASE_BLOCK(1, 4)
    SASKTRAN2_INSTANTIATE_PHASE_BLOCK(3, Eigen::Dynamic)
    SASKTRAN2_INSTANTIATE_PHASE_BLOCK(3, 1)
    SASKTRAN2_INSTANTIATE_PHASE_BLOCK(3, 4)

#undef SASKTRAN2_INSTANTIATE_PHASE_BLOCK

    template class PhaseHandler<1>;
    template class PhaseHandler<3>;

} // namespace sasktran2::solartransmission
