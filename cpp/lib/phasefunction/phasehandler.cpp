#include "sasktran2/raytracing.h"
#include <sasktran2/solartransmission.h>
#include <spdlog/spdlog.h>

namespace sasktran2::solartransmission {
    namespace {
        std::vector<std::pair<int, double>>
        endpoint_weights(const sasktran2::raytracing::TracedRay& ray,
                         std::size_t layer_index, bool entrance) {
            const auto stencil = entrance ? ray.entrance_weights(layer_index)
                                          : ray.exit_weights(layer_index);
            std::vector<std::pair<int, double>> result;
            result.reserve(stencil.size());
            for (std::size_t index = 0; index < stencil.size(); ++index) {
                const auto weight = stencil[index];
                if (weight.second != 0.0) {
                    result.push_back(weight);
                }
            }
            return result;
        }
    } // namespace

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
        const std::vector<std::vector<int>>& index_map) {
        initialize_geometry_impl(los_rays, index_map);
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_geometry_impl(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
        const std::vector<std::vector<int>>& index_map) {

        ZoneScopedN("Phase Handler Initialize Geometry");

        m_scatter_angles.clear();
        m_internal_to_geometry.clear();
        m_internal_to_cos_scatter.clear();
        m_geometry_entrance_to_internal.clear();
        m_geometry_exit_to_internal.clear();

        int num_internal = 0;
        int num_scatter = 0;
        double theta, C1, C2, S1, S2;
        int negation;
        // First we need to iterate through and figure out how many internal
        // indices we will end up with and how many scatter angles we will need

        // Keep track of the entrance and exits separately
        m_geometry_entrance_to_internal.resize(los_rays.size());
        m_geometry_exit_to_internal.resize(los_rays.size());
        for (int i = 0; i < los_rays.size(); ++i) {
            const auto& ray = los_rays[i];

            // Empty rays don't need to be considered
            if (ray.layers.size() == 0) {
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

            // We have to map every exit/internal layer point to an internal
            // scattering angle
            m_geometry_entrance_to_internal[i].resize(ray.layers.size());
            m_geometry_exit_to_internal[i].resize(ray.layers.size());

            for (int j = 0; j < ray.layers.size(); ++j) {

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

                const auto entrance_weights = endpoint_weights(ray, j, true);
                const auto exit_weights = endpoint_weights(ray, j, false);
                m_geometry_entrance_to_internal[i][j].resize(
                    entrance_weights.size());
                m_geometry_exit_to_internal[i][j].resize(exit_weights.size());

                for (int k = 0; k < entrance_weights.size(); ++k) {
                    m_geometry_entrance_to_internal[i][j][k] = num_internal;
                    ++num_internal;

                    m_internal_to_geometry.push_back(entrance_weights[k].first);
                    m_internal_to_cos_scatter.push_back(num_scatter);
                }

                bool can_share_exit = false;
                if (j > 0) {
                    const auto previous_entrance_weights =
                        endpoint_weights(ray, j - 1, true);
                    can_share_exit =
                        exit_weights.size() == previous_entrance_weights.size();
                    for (int k = 0; can_share_exit && k < exit_weights.size();
                         ++k) {
                        can_share_exit = exit_weights[k].first ==
                                         previous_entrance_weights[k].first;
                    }
                }

                if (!can_share_exit) {
                    // End layer at TOA, need to use layer exit
                    for (int k = 0; k < exit_weights.size(); ++k) {
                        m_geometry_exit_to_internal[i][j][k] = num_internal;
                        ++num_internal;

                        m_internal_to_geometry.push_back(exit_weights[k].first);
                        m_internal_to_cos_scatter.push_back(num_scatter);
                    }
                } else {
                    // Assign to previous layers entrance
                    m_geometry_exit_to_internal[i][j] =
                        m_geometry_entrance_to_internal[i][j - 1];
                }

                if (!ray.is_straight) {
                    ++num_scatter;
                }
            }
            if (ray.is_straight) {
                ++num_scatter;
            }
        }

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
        const auto& internal_indices =
            is_entrance ? m_geometry_entrance_to_internal[losidx][layeridx]
                        : m_geometry_exit_to_internal[losidx][layeridx];

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

    template class PhaseHandler<1>;
    template class PhaseHandler<3>;

} // namespace sasktran2::solartransmission
