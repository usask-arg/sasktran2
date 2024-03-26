#pragma once
#include "sktran_disco/sktran_do.h"

namespace sasktran_disco {

    template <int NSTOKES, int CNSTR>
    inline void lp_triple_product(
        const Eigen::Block<Eigen::MatrixXd, NSTOKES, NSTOKES>& const_Splus,
        const Eigen::Block<Eigen::MatrixXd, NSTOKES, NSTOKES>& const_Sminus,
        const std::vector<LegendreCoefficient<NSTOKES>>& coeffs,
        const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1s,
        const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2s, int m,
        const sasktran_disco::LayerDual<double>& ssa, double weight, double mu,
        bool on_diagonal) {
        // This is what Eigen says to do in their documentation...
        auto& Splus =
            const_cast<Eigen::Block<Eigen::MatrixXd, NSTOKES, NSTOKES>&>(
                const_Splus);
        auto& Sminus =
            const_cast<Eigen::Block<Eigen::MatrixXd, NSTOKES, NSTOKES>&>(
                const_Sminus);

        Splus.setZero();
        Sminus.setZero();

        Eigen::Matrix<double, NSTOKES, NSTOKES> eta;
        eta.setZero();

        if constexpr (NSTOKES == 1) {
            // Start by setting Splus to zeta, and using temporary for eta
            for (int l = m; l < coeffs.size(); ++l) {
                const auto& coeff = coeffs[l];
                const auto& lp1 = lp1s[l];
                const auto& lp2 = lp2s[l];

                // Calculated product by hand
                Splus(0, 0) += lp1.P() * lp2.P() * coeff.a1;

                int negation_factor_upper = 1;
                if ((l - m) % 2 != 0) {
                    negation_factor_upper *= -1;
                }
                // Calculated product by hand
                eta(0, 0) +=
                    lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
            }
        }

        if constexpr (NSTOKES == 3) {
            /*
            for (int l = 0; l < m; ++l) {
                a1deriv(l) = 0.0;
                a2deriv(l, Eigen::all).setZero();
                a3deriv(l, Eigen::all).setZero();
                b1deriv(l, Eigen::all).setZero();
            }
            */

            /*
                    zeta = ((*this->M_WT)[j] * l_upwelling.value(s1, s2) -
                            kronDelta(ii, jj)) /
                           (*this->M_MU)[i];
                    eta = (*this->M_WT)[j] * l_downwelling.value(s1, s2) /
                          (*this->M_MU)[i];
                }

                layer.solution(m).cache.s_plus()(ii, jj) =
                    -1.0 * (zeta + eta);
                layer.solution(m).cache.s_minus()(ii, jj) =
                    -1.0 * (zeta - eta);
            }
            */

            // Start by setting Splus to zeta, and using temporary for eta
            for (int l = m; l < coeffs.size(); ++l) {
                const auto& coeff = coeffs[l];
                const auto& lp1 = lp1s[l];
                const auto& lp2 = lp2s[l];

                // Calculated product by hand
                Splus(0, 0) += lp1.P() * lp2.P() * coeff.a1;
                Splus(0, 1) += -1.0 * lp1.P() * lp2.R() * coeff.b1;
                Splus(0, 2) += 1.0 * lp1.P() * lp2.T() * coeff.b1;
                // 0, 3 is always 0

                Splus(1, 0) += -1.0 * lp1.R() * lp2.P() * coeff.b1;
                Splus(1, 1) +=
                    lp1.R() * lp2.R() * coeff.a2 + lp1.T() * lp2.T() * coeff.a3;
                Splus(1, 2) += -lp1.R() * lp2.T() * coeff.a2 -
                               lp1.T() * lp2.R() * coeff.a3;

                Splus(2, 0) += lp1.T() * lp2.P() * coeff.b1;
                Splus(2, 1) += -lp1.T() * lp2.R() * coeff.a2 -
                               lp1.R() * lp2.T() * coeff.a3;
                Splus(2, 2) +=
                    lp1.T() * lp2.T() * coeff.a2 + lp1.R() * lp2.R() * coeff.a3;

                int negation_factor_upper = 1;
                int negation_factor_lower = -1;
                if ((l - m) % 2 != 0) {
                    negation_factor_upper *= -1;
                    negation_factor_lower *= -1;
                }
                // Calculated product by hand
                eta(0, 0) +=
                    lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
                eta(0, 1) +=
                    -1.0 * lp1.P() * lp2.R() * coeff.b1 * negation_factor_upper;
                eta(0, 2) +=
                    1.0 * lp1.P() * lp2.T() * coeff.b1 * negation_factor_upper;
                // 0, 3 is always 0

                eta(1, 0) +=
                    -1.0 * lp1.R() * lp2.P() * coeff.b1 * negation_factor_upper;
                eta(1, 1) +=
                    lp1.R() * lp2.R() * coeff.a2 * negation_factor_upper +
                    lp1.T() * lp2.T() * coeff.a3 * negation_factor_lower;
                eta(1, 2) +=
                    -lp1.R() * lp2.T() * coeff.a2 * negation_factor_upper -
                    lp1.T() * lp2.R() * coeff.a3 * negation_factor_lower;

                eta(2, 0) +=
                    lp1.T() * lp2.P() * coeff.b1 * negation_factor_upper;
                eta(2, 1) +=
                    -lp1.T() * lp2.R() * coeff.a2 * negation_factor_upper -
                    lp1.R() * lp2.T() * coeff.a3 * negation_factor_lower;
                eta(2, 2) +=
                    lp1.T() * lp2.T() * coeff.a2 * negation_factor_upper +
                    lp1.R() * lp2.R() * coeff.a3 * negation_factor_lower;
            }
        }
        // Include ssa and weight
        Splus.array() *= -0.5 * ssa.value * weight / mu;
        eta.array() *= -0.5 * ssa.value * weight / mu;

        if (on_diagonal) {
            Splus.diagonal().array() += 1 / mu;
        }

        // Now copy Splus into Sminus
        Sminus = Splus;

        // and +/- eta
        Splus.array() += eta.array();
        Sminus.array() -= eta.array();
    }

} // namespace sasktran_disco
