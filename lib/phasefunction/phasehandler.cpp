#include "sasktran2/raytracing.h"
#include <sasktran2/solartransmission.h>

namespace sasktran2::solartransmission {
    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays) {

        // Go through each ray and construct the mappings between internal and
        // atmosphere indicies
        int internal_index = 0;
        int cos_scatter_index = 0;
        double theta, C1, C2, S1, S2;
        int negation;
        for (const sasktran2::raytracing::TracedRay& ray : los_rays) {
            if (ray.layers.size() == 0) {
                continue;
            }
            // For each layer in the ray
            for (const sasktran2::raytracing::SphericalLayer& layer :
                 ray.layers) {
                const auto& entrance_weights =
                    layer.entrance.interpolation_weights;

                for (const auto& weight : entrance_weights) {
                    if (m_internal_to_geometry.count(internal_index) == 0) {
                        m_internal_to_geometry[internal_index] = weight.first;
                        m_geometry_to_internal[weight.first] = internal_index;
                        m_internal_to_cos_scatter[internal_index] =
                            cos_scatter_index;

                        ++internal_index;
                    }
                }

                const auto& exit_weights = layer.exit.interpolation_weights;
                for (const auto& weight : layer.exit.interpolation_weights) {
                    if (m_internal_to_geometry.count(internal_index) == 0) {
                        m_internal_to_geometry[internal_index] = weight.first;
                        m_geometry_to_internal[weight.first] = internal_index;
                        m_internal_to_cos_scatter[internal_index] =
                            cos_scatter_index;

                        ++internal_index;
                    }
                }

                if (!ray.is_straight) {
                    ++cos_scatter_index;

                    math::stokes_scattering_factors(
                        -1 * m_geometry.coordinates().sun_unit(),
                        -1 * layer.average_look_away, theta, C1, C2, S1, S2,
                        negation);
                    if constexpr (NSTOKES == 3) {
                        m_scatter_angles.push_back({theta, C1, C2});
                    } else {
                        m_scatter_angles.push_back({theta});
                    }
                }
            }
            if (ray.is_straight) {
                ++cos_scatter_index;

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
        }

        if constexpr (NSTOKES == 3) {
            m_phase.resize(2, m_internal_to_geometry.size(),
                           m_config->num_wavelength_threads());
        } else {
            m_phase.resize(1, m_internal_to_geometry.size(),
                           m_config->num_wavelength_threads());
        }
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_atmosphere = &atmosphere;

        int numderiv = m_atmosphere->num_scattering_deriv_groups();

        if constexpr (NSTOKES == 1) {
            m_d_phase.resize(1, m_internal_to_geometry.size(), numderiv,
                             m_config->num_wavelength_threads());
        } else {
            m_d_phase.resize(2, m_internal_to_geometry.size(), numderiv,
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

            sasktran2::math::WignerDCalculator d00(0, 0);
            sasktran2::math::WignerDCalculator d02(0, 2);
            for (int i = 0; i < m_internal_to_geometry.size(); ++i) {
                double theta =
                    m_scatter_angles[m_internal_to_cos_scatter[i]][0];

                int atmo_index = m_internal_to_geometry[i];

                for (int k = 0; k < m_atmosphere->storage().max_order(
                                        atmo_index, wavelidx);
                     ++k) {
                    if constexpr (NSTOKES == 1) {
                        m_phase(0, i, threadidx) +=
                            m_atmosphere->storage().leg_coeff(k, atmo_index,
                                                              wavelidx) *
                            d00.d(theta, k);
                    } else {
                        m_phase(0, i, threadidx) +=
                            m_atmosphere->storage().leg_coeff(k * 4, atmo_index,
                                                              wavelidx) *
                            d00.d(theta, k);
                        m_phase(1, i, threadidx) +=
                            m_atmosphere->storage().leg_coeff(
                                k * 4 + 3, atmo_index, wavelidx) *
                            d02.d(theta, k);
                    }

                    for (int d = 0;
                         d < m_atmosphere->num_scattering_deriv_groups(); ++d) {
                        if constexpr (NSTOKES == 1) {
                            m_d_phase(0, i, d, threadidx) +=
                                m_atmosphere->storage().d_leg_coeff(
                                    k, atmo_index, wavelidx, d) *
                                d00.d(theta, k);
                        } else {
                            m_d_phase(0, i, d, threadidx) +=
                                m_atmosphere->storage().d_leg_coeff(
                                    k * 4, atmo_index, wavelidx, d) *
                                d00.d(theta, k);
                            m_d_phase(1, i, d, threadidx) +=
                                m_atmosphere->storage().d_leg_coeff(
                                    k * 4 + 3, atmo_index, wavelidx, d) *
                                d02.d(theta, k);
                        }
                    }
                }
            }
        } else {
            // How to do this?
        }
    }

    template <int NSTOKES>
    void PhaseHandler<NSTOKES>::scatter(
        int threadidx, const std::vector<std::pair<int, double>>& index_weights,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        if constexpr (NSTOKES == 1) {
            double phase_result = 0.0;

            for (const auto& index_weight : index_weights) {
                int internal_index =
                    m_geometry_to_internal.at(index_weight.first);
                phase_result +=
                    m_phase(0, internal_index, threadidx) * index_weight.second;
            }

            for (int d = 0; d < m_atmosphere->num_scattering_deriv_groups();
                 ++d) {
                for (const auto& index_weight : index_weights) {
                    int internal_index =
                        m_geometry_to_internal.at(index_weight.first);

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

            for (const auto& index_weight : index_weights) {
                int internal_index =
                    m_geometry_to_internal.at(index_weight.first);
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
                for (const auto& index_weight : index_weights) {
                    int internal_index =
                        m_geometry_to_internal.at(index_weight.first);
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
