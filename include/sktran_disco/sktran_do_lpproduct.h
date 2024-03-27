#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_types.h"

namespace sasktran_disco {
    template <int NSTOKES, int CNSTR> class DerivBlockIter {
        using Matrix =
            typename std::conditional<CNSTR == -1, Eigen::MatrixXd,
                                      Eigen::Matrix<double, CNSTR / 2 * NSTOKES,
                                                    CNSTR / 2 * NSTOKES>>::type;

      private:
        std::vector<sasktran_disco::RTESolutionCache<NSTOKES, CNSTR>>& m_cache;
        const sasktran_disco::InputDerivatives<NSTOKES>& m_input_derivs;
        LayerIndex m_layer_index;
        int m_i;
        int m_j;

        int m_num_deriv_layer;
        int m_layer_start;

      public:
        DerivBlockIter(
            std::vector<sasktran_disco::RTESolutionCache<NSTOKES, CNSTR>>&
                cache,
            const sasktran_disco::InputDerivatives<NSTOKES>& input_derivs,
            LayerIndex layer_index)
            : m_cache(cache), m_input_derivs(input_derivs),
              m_layer_index(layer_index) {
            m_num_deriv_layer = input_derivs.numDerivativeLayer(m_layer_index);
            m_layer_start = input_derivs.layerStartIndex(m_layer_index);
        }

        void set_block(int i, int j) {
            m_i = i;
            m_j = j;
        }

        size_t num_deriv() const { return m_num_deriv_layer; }

        const LayerInputDerivative<NSTOKES>& input_deriv(int deriv) {
            return m_input_derivs.layerDerivatives()[m_layer_start + deriv];
        }

        double& d_Splus(int deriv, int i, int j) {
            return m_cache[deriv].s_plus()(m_i * NSTOKES + i,
                                           m_j * NSTOKES + j);
        }

        double& d_Sminus(int deriv, int i, int j) {
            return m_cache[deriv].s_minus()(m_i * NSTOKES + i,
                                            m_j * NSTOKES + j);
        }

        Eigen::Block<Matrix, NSTOKES, NSTOKES> d_Splus(int deriv) {
            return m_cache[deriv].s_plus().template block<NSTOKES, NSTOKES>(
                m_i * NSTOKES, m_j * NSTOKES);
        }

        Eigen::Block<Matrix, NSTOKES, NSTOKES> d_Sminus(int deriv) {
            return m_cache[deriv].s_minus().template block<NSTOKES, NSTOKES>(
                m_i * NSTOKES, m_j * NSTOKES);
        }
    };

    template <int NSTOKES, int CNSTR>
    inline void lp_triple_product(
        const Eigen::Block<
            typename std::conditional<CNSTR == -1, Eigen::MatrixXd,
                                      Eigen::Matrix<double, CNSTR / 2 * NSTOKES,
                                                    CNSTR / 2 * NSTOKES>>::type,
            NSTOKES, NSTOKES>& const_Splus,
        const Eigen::Block<
            typename std::conditional<CNSTR == -1, Eigen::MatrixXd,
                                      Eigen::Matrix<double, CNSTR / 2 * NSTOKES,
                                                    CNSTR / 2 * NSTOKES>>::type,
            NSTOKES, NSTOKES>& const_Sminus,
        const std::vector<LegendreCoefficient<NSTOKES>>& coeffs,
        const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1s,
        const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2s, int m,
        const sasktran_disco::LayerDual<double>& ssa, double weight, double mu,
        DerivBlockIter<NSTOKES, CNSTR>& deriv_blocks, bool on_diagonal) {
        // This is what Eigen says to do in their documentation...
        auto& Splus = const_cast<Eigen::Block<
            typename std::conditional<CNSTR == -1, Eigen::MatrixXd,
                                      Eigen::Matrix<double, CNSTR / 2 * NSTOKES,
                                                    CNSTR / 2 * NSTOKES>>::type,
            NSTOKES, NSTOKES>&>(const_Splus);
        auto& Sminus = const_cast<Eigen::Block<
            typename std::conditional<CNSTR == -1, Eigen::MatrixXd,
                                      Eigen::Matrix<double, CNSTR / 2 * NSTOKES,
                                                    CNSTR / 2 * NSTOKES>>::type,
            NSTOKES, NSTOKES>&>(const_Sminus);

        Splus.setZero();
        Sminus.setZero();

        Eigen::Matrix<double, NSTOKES, NSTOKES> eta, d_eta;
        eta.setZero();

        for (int l = m; l < coeffs.size(); ++l) {
            if constexpr (NSTOKES == 1) {
                // Start by setting Splus to zeta, and using temporary for eta
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
            if constexpr (NSTOKES == 3) {
                // Start by setting Splus to zeta, and using temporary for eta
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
        // Include quadrature factors since we need these for the derivatives
        Splus.array() *= -0.5 * weight / mu;
        eta.array() *= -0.5 * weight / mu;

        // Now do the derivatives
        for (int k = 0; k < deriv_blocks.num_deriv(); ++k) {
            auto d_Splus = (deriv_blocks.d_Splus(k));
            auto d_Sminus = (deriv_blocks.d_Sminus(k));

            const sasktran_disco::LayerInputDerivative<NSTOKES>& deriv =
                deriv_blocks.input_deriv(k);

            d_Splus.setZero();
            d_Sminus.setZero();

            d_eta.setZero();

            for (int l = m; l < coeffs.size(); ++l) {
                const auto& coeff = deriv.d_legendre_coeff[l];
                const auto& lp1 = lp1s[l];
                const auto& lp2 = lp2s[l];
                if constexpr (NSTOKES == 1) {
                    d_Splus(0, 0) += lp1.P() * lp2.P() * coeff.a1;

                    int negation_factor_upper = 1;
                    if ((l - m) % 2 != 0) {
                        negation_factor_upper *= -1;
                    }
                    d_eta(0, 0) +=
                        lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
                }
                if constexpr (NSTOKES == 3) {
                    // Calculated product by hand
                    d_Splus(0, 0) += lp1.P() * lp2.P() * coeff.a1;
                    d_Splus(0, 1) += -1.0 * lp1.P() * lp2.R() * coeff.b1;
                    d_Splus(0, 2) += 1.0 * lp1.P() * lp2.T() * coeff.b1;
                    // 0, 3 is always 0

                    d_Splus(1, 0) += -1.0 * lp1.R() * lp2.P() * coeff.b1;
                    d_Splus(1, 1) += lp1.R() * lp2.R() * coeff.a2 +
                                     lp1.T() * lp2.T() * coeff.a3;
                    d_Splus(1, 2) += -lp1.R() * lp2.T() * coeff.a2 -
                                     lp1.T() * lp2.R() * coeff.a3;

                    d_Splus(2, 0) += lp1.T() * lp2.P() * coeff.b1;
                    d_Splus(2, 1) += -lp1.T() * lp2.R() * coeff.a2 -
                                     lp1.R() * lp2.T() * coeff.a3;
                    d_Splus(2, 2) += lp1.T() * lp2.T() * coeff.a2 +
                                     lp1.R() * lp2.R() * coeff.a3;

                    int negation_factor_upper = 1;
                    int negation_factor_lower = -1;
                    if ((l - m) % 2 != 0) {
                        negation_factor_upper *= -1;
                        negation_factor_lower *= -1;
                    }
                    // Calculated product by hand
                    d_eta(0, 0) +=
                        lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
                    d_eta(0, 1) += -1.0 * lp1.P() * lp2.R() * coeff.b1 *
                                   negation_factor_upper;
                    d_eta(0, 2) += 1.0 * lp1.P() * lp2.T() * coeff.b1 *
                                   negation_factor_upper;
                    // 0, 3 is always 0

                    d_eta(1, 0) += -1.0 * lp1.R() * lp2.P() * coeff.b1 *
                                   negation_factor_upper;
                    d_eta(1, 1) +=
                        lp1.R() * lp2.R() * coeff.a2 * negation_factor_upper +
                        lp1.T() * lp2.T() * coeff.a3 * negation_factor_lower;
                    d_eta(1, 2) +=
                        -lp1.R() * lp2.T() * coeff.a2 * negation_factor_upper -
                        lp1.T() * lp2.R() * coeff.a3 * negation_factor_lower;

                    d_eta(2, 0) +=
                        lp1.T() * lp2.P() * coeff.b1 * negation_factor_upper;
                    d_eta(2, 1) +=
                        -lp1.T() * lp2.R() * coeff.a2 * negation_factor_upper -
                        lp1.R() * lp2.T() * coeff.a3 * negation_factor_lower;
                    d_eta(2, 2) +=
                        lp1.T() * lp2.T() * coeff.a2 * negation_factor_upper +
                        lp1.R() * lp2.R() * coeff.a3 * negation_factor_lower;
                }
            }
            // Start by including the quadrature factors
            d_Splus.array() *= -0.5 * ssa.value * weight / mu;
            d_eta.array() *= -0.5 * ssa.value * weight / mu;

            // Add on the derivatives with respect to ssa
            d_Splus.array() += Splus.array() * deriv.d_SSA;
            d_eta.array() += eta.array() * deriv.d_SSA;

            // Now copy Splus into Sminus
            d_Sminus = d_Splus;

            // and +/- eta
            d_Splus.array() += d_eta.array();
            d_Sminus.array() -= d_eta.array();
        }

        // Now we can include ssa and finish the calculation
        Splus.array() *= ssa.value;
        eta.array() *= ssa.value;

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
