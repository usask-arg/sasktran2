#include "sasktran2/raytracing.h"
#include <sasktran2/solartransmission.h>
#include <spdlog/spdlog.h>

namespace sasktran2::solartransmission {
    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
        const std::vector<std::vector<int>>& index_map) {

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
            const sasktran2::raytracing::TracedRay& ray = los_rays[i];

            // Empty rays don't need to be considered
            if (ray.layers.size() == 0) {
                continue;
            }

            // Straight rays we can just use the first layer to get the
            // scattering angle, scattering angle does not change along the ray
            if (ray.is_straight) {
                math::stokes_scattering_factors(
                    -1 * m_geometry.coordinates().sun_unit(),
                    -1 * ray.layers[0].average_look_away, theta, C1, C2, S1, S2,
                    negation);
                if constexpr (NSTOKES == 3) {
                    m_scatter_angles.push_back({theta, C2, S2});
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
                    math::stokes_scattering_factors(
                        -1 * m_geometry.coordinates().sun_unit(),
                        -1 * ray.layers[0].average_look_away, theta, C1, C2, S1,
                        S2, negation);
                    if constexpr (NSTOKES == 3) {
                        m_scatter_angles.push_back({theta, C2, S2});
                    } else {
                        m_scatter_angles.push_back({theta});
                    }
                }

                const sasktran2::raytracing::SphericalLayer& layer =
                    ray.layers[j];
                m_geometry_entrance_to_internal[i][j].resize(
                    layer.entrance.interpolation_weights.size());
                m_geometry_exit_to_internal[i][j].resize(
                    layer.exit.interpolation_weights.size());

                for (int k = 0; k < layer.entrance.interpolation_weights.size();
                     ++k) {
                    m_geometry_entrance_to_internal[i][j][k] = num_internal;
                    ++num_internal;

                    m_internal_to_geometry.push_back(
                        layer.entrance.interpolation_weights[k].first);
                    m_internal_to_cos_scatter.push_back(num_scatter);
                }

                if (j == 0) {
                    // End layer at TOA, need to use layer exit
                    for (int k = 0; k < layer.exit.interpolation_weights.size();
                         ++k) {
                        m_geometry_exit_to_internal[i][j][k] = num_internal;
                        ++num_internal;

                        m_internal_to_geometry.push_back(
                            layer.exit.interpolation_weights[k].first);
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
            for (int j = 0; j < m_config->num_singlescatter_moments(); ++j) {
                m_wigner_d00(j, i) = d00.d(m_scatter_angles[i][0], j);

                if constexpr (NSTOKES == 3) {
                    m_wigner_d02(j, i) = d02.d(m_scatter_angles[i][0], j);
                }
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
        const std::vector<std::pair<int, double>>& index_weights,
        bool is_entrance,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        const auto& internal_indices =
            is_entrance ? m_geometry_entrance_to_internal[losidx][layeridx]
                        : m_geometry_exit_to_internal[losidx][layeridx];

        if constexpr (NSTOKES == 1) {
            double phase_result = 0.0;

            for (int i = 0; i < index_weights.size(); ++i) {
                const auto& index_weight = index_weights[i];
                int internal_index = internal_indices[i];
                phase_result +=
                    m_phase(0, internal_index, threadidx) * index_weight.second;
            }

            for (int d = 0; d < m_atmosphere->num_scattering_deriv_groups();
                 ++d) {
                for (int i = 0; i < index_weights.size(); ++i) {
                    const auto& index_weight = index_weights[i];
                    int internal_index = internal_indices[i];

                    source.deriv(0, m_atmosphere->scat_deriv_start_index() +
                                        d * m_geometry.size() +
                                        index_weight.first) +=
                        source.value(0) *
                        m_d_phase(0, internal_index, d, threadidx) *
                        index_weight.second;
                }
            }
            source.value(0) *= phase_result;

        } else if constexpr (NSTOKES == 3) {
            Eigen::Vector3d phase_result = Eigen::Vector3d::Zero();

            for (int i = 0; i < index_weights.size(); ++i) {
                const auto& index_weight = index_weights[i];
                int internal_index = internal_indices[i];
                double C2 = m_scatter_angles[m_internal_to_cos_scatter.at(
                    internal_index)][1];
                double S2 = m_scatter_angles[m_internal_to_cos_scatter.at(
                    internal_index)][2];

                double off_diag =
                    m_phase(1, internal_index, threadidx) * index_weight.second;

                phase_result(0) +=
                    m_phase(0, internal_index, threadidx) * index_weight.second;
                phase_result(1) += -C2 * off_diag;
                phase_result(2) += -S2 * off_diag;
            }
            for (int d = 0; d < m_atmosphere->num_scattering_deriv_groups();
                 ++d) {
                for (int i = 0; i < index_weights.size(); ++i) {
                    const auto& index_weight = index_weights[i];

                    int internal_index = internal_indices[i];
                    double C2 = m_scatter_angles[m_internal_to_cos_scatter.at(
                        internal_index)][1];
                    double S2 = m_scatter_angles[m_internal_to_cos_scatter.at(
                        internal_index)][2];

                    int deriv_index = m_atmosphere->scat_deriv_start_index() +
                                      d * m_geometry.size() +
                                      index_weight.first;
                    source.deriv(0, deriv_index) +=
                        source.value(0) *
                        m_d_phase(0, internal_index, d, threadidx) *
                        index_weight.second;

                    source.deriv(1, deriv_index) +=
                        source.value(0) *
                        m_d_phase(1, internal_index, d, threadidx) *
                        index_weight.second * (-C2);

                    source.deriv(2, deriv_index) +=
                        source.value(0) *
                        m_d_phase(1, internal_index, d, threadidx) *
                        index_weight.second * (-S2);
                }
            }

            source.value.array() = source.value(0) * phase_result.array();
        } else {
        }
    }

    template class PhaseHandler<1>;
    template class PhaseHandler<3>;

} // namespace sasktran2::solartransmission
