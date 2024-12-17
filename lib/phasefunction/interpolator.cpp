#include <sasktran2/atmosphere/grid_storage.h>
#include <sasktran2/math/wigner.h>
#include <sasktran2/math/scattering.h>

namespace sasktran2::atmosphere {

    template <int NSTOKES, bool ssonly>
    PhaseInterpolator<NSTOKES, ssonly>::PhaseInterpolator()
        : m_geometry_loaded(false) {}

    template <int NSTOKES, bool ssonly>
    void PhaseInterpolator<NSTOKES, ssonly>::load_scattering_angle(
        int num_legendre, const Eigen::Vector3d& incoming_ray,
        const Eigen::Vector3d& outgoing_ray, bool outgoing_facing_away) {
        if (m_geometry_loaded) {
            return;
        }
        m_geometry_loaded = true;

        if constexpr (NSTOKES == 1) {
            double cos_scatter = incoming_ray.dot(outgoing_ray);
            if (!outgoing_facing_away) {
                cos_scatter *= -1;
            }
            double theta = acos(cos_scatter);

            sasktran2::math::WignerDCalculator wigner(0, 0);

            m_scattering_weights.resize(1, num_legendre);

            for (int i = 0; i < num_legendre; ++i) {
                m_scattering_weights(0, i) = wigner.d(theta, i);
            }
        }

        if constexpr (NSTOKES == 3) {
            double C1, C2, S1, S2, theta;
            int negation;

            // INCOMING: This is look vector away from observer, i.e., opposite
            // the propagation direction. OUTGOING: This is direction towards
            // the sun, i.e., opposite the propagation direction

            math::stokes_scattering_factors(-1 * incoming_ray,
                                            -1 * outgoing_ray, theta, C1, C2,
                                            S1, S2, negation);

            // We have 4 greek coefficients (a1, a2, a3, b1), and we use the
            // expansion (sum over legendre poly assumed) P11 = a1 d_{0, 0} P12
            // = b1 d_{0, 2} P22 + P33 = (a2 + a3) d_{2, 2} P22 - P33 = (a2 -
            // a3) d_{2,-2}

            // Quick simplification
            // 2 P22 = (a2 + a3) d_{2, 2} + (a2 - a3) d_{2,-2}
            // 2 P33 = (a2 - a3) d_{2,-2} - (a2 - a3) d_{2,-2}

            // Create the wigner functions that we need
            sasktran2::math::WignerDCalculator d00(0, 0);
            sasktran2::math::WignerDCalculator d22(2, 2);
            sasktran2::math::WignerDCalculator d02(0, 2);
            sasktran2::math::WignerDCalculator d2m2(2, -2);

            // Now we have to create the (3x3, 4) scattering matrix multiplier
            // for each legendre expansion coefficient
            m_scattering_weights.resize(m_scattering_weights.rows(),
                                        4 * num_legendre);
            m_scattering_weights.setZero();

            for (int i = 0; i < num_legendre; ++i) {
                int a1_idx = i * 4;
                int a2_idx = a1_idx + 1;
                int a3_idx = a1_idx + 2;
                int b1_idx = a1_idx + 3;

                // 0,0 component
                m_scattering_weights(0, a1_idx) = d00.d(theta, i);

                // 1, 0 component
                m_scattering_weights(1, b1_idx) = -C2 * d02.d(theta, i);

                // 2, 0 component
                m_scattering_weights(2, b1_idx) = -S2 * d02.d(theta, i);

                if constexpr (!ssonly) {
                    // 0, 1 component
                    m_scattering_weights(3, b1_idx) = -C1 * d02.d(theta, i);

                    // 1, 1 component
                    m_scattering_weights(4, a2_idx) +=
                        0.5 * C1 * C2 * (d22.d(theta, i) + d2m2.d(theta, i));
                    m_scattering_weights(4, a3_idx) +=
                        0.5 * C1 * C2 * (d22.d(theta, i) - d2m2.d(theta, i));

                    m_scattering_weights(4, a2_idx) -=
                        0.5 * S1 * S2 * (d22.d(theta, i) - d2m2.d(theta, i));
                    m_scattering_weights(4, a3_idx) -=
                        0.5 * S1 * S2 * (d22.d(theta, i) + d2m2.d(theta, i));

                    // 2, 1 component
                    m_scattering_weights(5, a2_idx) +=
                        0.5 * C1 * S2 * (d22.d(theta, i) + d2m2.d(theta, i));
                    m_scattering_weights(5, a3_idx) +=
                        0.5 * C1 * S2 * (d22.d(theta, i) - d2m2.d(theta, i));

                    m_scattering_weights(5, a2_idx) +=
                        0.5 * S1 * C2 * (d22.d(theta, i) - d2m2.d(theta, i));
                    m_scattering_weights(5, a3_idx) +=
                        0.5 * S1 * C2 * (d22.d(theta, i) + d2m2.d(theta, i));

                    // 0, 2 component
                    m_scattering_weights(6, b1_idx) = S1 * d02.d(theta, i);

                    // 1, 2 component
                    m_scattering_weights(7, a2_idx) -=
                        0.5 * S1 * C2 * (d22.d(theta, i) + d2m2.d(theta, i));
                    m_scattering_weights(7, a3_idx) -=
                        0.5 * S1 * C2 * (d22.d(theta, i) - d2m2.d(theta, i));

                    m_scattering_weights(7, a2_idx) -=
                        0.5 * C1 * S2 * (d22.d(theta, i) - d2m2.d(theta, i));
                    m_scattering_weights(7, a3_idx) -=
                        0.5 * C1 * S2 * (d22.d(theta, i) + d2m2.d(theta, i));

                    // 2, 2 component
                    m_scattering_weights(8, a2_idx) -=
                        0.5 * S1 * S2 * (d22.d(theta, i) + d2m2.d(theta, i));
                    m_scattering_weights(8, a3_idx) -=
                        0.5 * S1 * S2 * (d22.d(theta, i) - d2m2.d(theta, i));

                    m_scattering_weights(8, a2_idx) +=
                        0.5 * C1 * C2 * (d22.d(theta, i) - d2m2.d(theta, i));
                    m_scattering_weights(8, a3_idx) +=
                        0.5 * C1 * C2 * (d22.d(theta, i) + d2m2.d(theta, i));
                }
            }
        }
    }

    template <int NSTOKES, bool ssonly>
    template <sasktran2::dualstorage S>
    void PhaseInterpolator<NSTOKES, ssonly>::scatter(
        const AtmosphereGridStorageFull<NSTOKES>& phase_storage, int wavelidx,
        const std::vector<std::pair<int, double>>& index_weights,
        sasktran2::Dual<double, S, NSTOKES>& source) const {

        using PhaseResult = typename std::conditional<
            ssonly, Eigen::Vector<sasktran2::types::leg_coeff, NSTOKES>,
            Eigen::Vector<sasktran2::types::leg_coeff,
                          NSTOKES * NSTOKES>>::type;

        int numscatderiv = phase_storage.numscatderiv;
        int derivstart = phase_storage.scatderivstart;

        const auto& d = phase_storage.leg_coeff.dimensions();

        Eigen::Map<const Eigen::MatrixXd> phase_matrix(
            &phase_storage.leg_coeff(0, 0, wavelidx), d[0], d[1]);

        if constexpr (NSTOKES == 1) {
            // Apply phase derivatives
            if constexpr (S == sasktran2::dualstorage::dense) {
                if (source.deriv.size() > 0) {
                    for (const auto& ele : index_weights) {
                        for (int i = 0; i < numscatderiv; ++i) {
                            Eigen::Map<const Eigen::MatrixXd> d_phase_matrix(
                                &phase_storage.d_leg_coeff(0, 0, wavelidx, i),
                                d[0], d[1]);

                            source.deriv(Eigen::all,
                                         derivstart + i * d[1] + ele.first) +=
                                (ele.second * source.value(0) *
                                 (m_scattering_weights *
                                  d_phase_matrix(Eigen::all, ele.first)))
                                    .template cast<double>();
                        }
                    }
                }
            } else {
                spdlog::warn("Phase derivatives with an input sparse source "
                             "derivative is not supported");
            }

            // Have to interpolate the phase
            sasktran2::types::leg_coeff phase_result = 0;

            for (const auto& ele : index_weights) {
                auto non_zero_leg =
                    Eigen::seq(0, phase_storage.max_order(ele.first, wavelidx));
                phase_result +=
                    ele.second *
                    (m_scattering_weights(non_zero_leg)
                         .dot(phase_matrix(non_zero_leg, ele.first)));
            }
            // Apply the phase function
            source.value(0) *= phase_result;
        } else {
            // Apply phase derivatives
            // Only apply phase derivatives in the single scatter case
            if constexpr (ssonly) {
                if constexpr (S == sasktran2::dualstorage::dense) {
                    if (source.deriv.size() > 0) {
                        for (const auto& ele : index_weights) {
                            for (int i = 0; i < numscatderiv; ++i) {
                                Eigen::Map<const Eigen::MatrixXd>
                                    d_phase_matrix(&phase_storage.d_leg_coeff(
                                                       0, 0, wavelidx, i),
                                                   d[0], d[1]);

                                source
                                    .deriv(Eigen::all,
                                           derivstart + i * d[1] + ele.first)
                                    .array() +=
                                    source.value(0) *
                                    (ele.second *
                                     (m_scattering_weights *
                                      d_phase_matrix(Eigen::all, ele.first))
                                         .array());
                            }
                        }
                    }
                }
            }
            PhaseResult phase_result;
            phase_result.setZero();

            // Interpolate the phase matrix to the geometric location we are at
            for (const auto& ele : index_weights) {
                phase_result.array() +=
                    ele.second *
                    (m_scattering_weights * phase_matrix(Eigen::all, ele.first))
                        .array();
            }

            // TODO: Derivatives here
            if constexpr (!ssonly) {
                source.value.applyOnTheLeft(
                    Eigen::Map<Eigen::Matrix<double, NSTOKES, NSTOKES>>(
                        phase_result.data(), NSTOKES, NSTOKES));
            } else {
                source.value = source.value(0) * phase_result;
            }

#ifdef SASKTRAN_DEBUG_ASSERTS
            if (source.value.hasNaN()) {
                spdlog::error("Error scattering");
            }
#endif
        }
    }

    template class PhaseInterpolator<1, true>;
    template class PhaseInterpolator<3, true>;

    template class PhaseInterpolator<1, false>;
    template class PhaseInterpolator<3, false>;

    template void PhaseInterpolator<1, false>::scatter(
        const AtmosphereGridStorageFull<1>&, int,
        const std::vector<std::pair<int, double>>&,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&) const;
    template void PhaseInterpolator<3, false>::scatter(
        const AtmosphereGridStorageFull<3>&, int,
        const std::vector<std::pair<int, double>>&,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>&) const;

    template void PhaseInterpolator<1, true>::scatter(
        const AtmosphereGridStorageFull<1>&, int,
        const std::vector<std::pair<int, double>>&,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&) const;
    template void PhaseInterpolator<3, true>::scatter(
        const AtmosphereGridStorageFull<3>&, int,
        const std::vector<std::pair<int, double>>&,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>&) const;

} // namespace sasktran2::atmosphere
