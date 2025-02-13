#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_polarization_types.h"
#include "sktran_disco/sktran_do_types.h"

namespace sasktran_disco {
    /**
     * Calculates the triple product of two Legendre polynomials and a Legendre
     * coefficient
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     * @tparam negation
     * @tparam T
     * @param coeff
     * @param lp1
     * @param lp2
     * @param l
     * @param m
     * @param const_result
     */
    template <int NSTOKES, int CNSTR, bool negation, typename T>
    inline void lp_triple_product(const LegendreCoefficient<NSTOKES>& coeff,
                                  const LegendrePhaseContainer<NSTOKES>& lp1,
                                  const LegendrePhaseContainer<NSTOKES>& lp2,
                                  int l, int m, const T& const_result) {
        auto& result = const_cast<T&>(const_result);

        if constexpr (NSTOKES == 1) {
            int negation_factor_upper = 1;
            if constexpr (negation) {
                if ((l - m) % 2 != 0 && negation) {
                    negation_factor_upper *= -1;
                }
            }

            result(0, 0) +=
                lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
        }
        if constexpr (NSTOKES == 3) {
            int negation_factor_upper = 1;
            int negation_factor_lower = 1;

            if constexpr (negation) {
                negation_factor_lower = -1;
                if ((l - m) % 2 != 0) {
                    negation_factor_upper *= -1;
                    negation_factor_lower *= -1;
                }
            }

            // Calculated product by hand
            result(0, 0) +=
                lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
            result(0, 1) +=
                -1.0 * lp1.P() * lp2.R() * coeff.b1 * negation_factor_upper;
            result(0, 2) +=
                1.0 * lp1.P() * lp2.T() * coeff.b1 * negation_factor_upper;
            // 0, 3 is always 0

            result(1, 0) +=
                -1.0 * lp1.R() * lp2.P() * coeff.b1 * negation_factor_upper;
            result(1, 1) +=
                lp1.R() * lp2.R() * coeff.a2 * negation_factor_upper +
                lp1.T() * lp2.T() * coeff.a3 * negation_factor_lower;
            result(1, 2) +=
                -lp1.R() * lp2.T() * coeff.a2 * negation_factor_upper -
                lp1.T() * lp2.R() * coeff.a3 * negation_factor_lower;

            result(2, 0) +=
                lp1.T() * lp2.P() * coeff.b1 * negation_factor_upper;
            result(2, 1) +=
                -lp1.T() * lp2.R() * coeff.a2 * negation_factor_upper -
                lp1.R() * lp2.T() * coeff.a3 * negation_factor_lower;
            result(2, 2) +=
                lp1.T() * lp2.T() * coeff.a2 * negation_factor_upper +
                lp1.R() * lp2.R() * coeff.a3 * negation_factor_lower;
        }
    }

    /**
     * @brief Handles derivatives of Splus/Sminus for the Legendre solution
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
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

    /**
     * @brief Calculates a block of the Splus/Sminus matrices in the homogeneous
     * solution and their derivatives
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     * @param const_Splus
     * @param const_Sminus
     * @param coeffs
     * @param lp1s
     * @param lp2s
     * @param m
     * @param ssa
     * @param weight
     * @param mu
     * @param deriv_blocks
     * @param on_diagonal
     */
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
            lp_triple_product<NSTOKES, CNSTR, false>(coeffs[l], lp1s[l],
                                                     lp2s[l], l, m, Splus);
            lp_triple_product<NSTOKES, CNSTR, true>(coeffs[l], lp1s[l], lp2s[l],
                                                    l, m, eta);
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

                lp_triple_product<NSTOKES, CNSTR, false>(coeff, lp1, lp2, l, m,
                                                         d_Splus);
                lp_triple_product<NSTOKES, CNSTR, true>(coeff, lp1, lp2, l, m,
                                                        d_eta);
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

    template <int NSTOKES, int CNSTR>
    inline void scat_phase_f(
        const std::vector<LegendreCoefficient<NSTOKES>>& coeffs,
        const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_mu,
        const sasktran_disco::VectorDim2<LegendrePhaseContainer<NSTOKES>>& lp,
        int m, LayerIndex p, const LayerDual<double>& ssa,
        const std::vector<double>& weights,
        const InputDerivatives<NSTOKES>& input_derivs,
        VectorLayerDual<double>& dual_lpsum_plus,
        VectorLayerDual<double>& dual_lpsum_minus) {
        ZoneScopedN("Scattering Phase Calculation");
        dual_lpsum_plus.value.setZero();
        dual_lpsum_minus.value.setZero();

        for (StreamIndex j = 0; j < coeffs.size() / 2; ++j) {
            Eigen::Map<Eigen::Matrix<double, NSTOKES, NSTOKES>> plus_matrix(
                dual_lpsum_plus.value.data() + j * NSTOKES * NSTOKES);
            Eigen::Map<Eigen::Matrix<double, NSTOKES, NSTOKES>> minus_matrix(
                dual_lpsum_minus.value.data() + j * NSTOKES * NSTOKES);

            for (int l = m; l < coeffs.size(); ++l) {
                lp_triple_product<NSTOKES, CNSTR, false>(
                    coeffs[l], lp_mu[l], lp[j][l], l, m, plus_matrix);
                lp_triple_product<NSTOKES, CNSTR, true>(
                    coeffs[l], lp_mu[l], lp[j][l], l, m, minus_matrix);
            }

            for (int k = 0; k < input_derivs.numDerivativeLayer(p); ++k) {
                const LayerInputDerivative<NSTOKES>& deriv =
                    input_derivs
                        .layerDerivatives()[input_derivs.layerStartIndex(p) +
                                            k];

                Eigen::Map<Eigen::Matrix<double, NSTOKES, NSTOKES>, 0,
                           Eigen::InnerStride<>>
                    d_plus_matrix(
                        &dual_lpsum_plus.deriv(k, j * NSTOKES * NSTOKES),
                        NSTOKES, NSTOKES,
                        Eigen::InnerStride<>(
                            input_derivs.numDerivativeLayer(p)));
                Eigen::Map<Eigen::Matrix<double, NSTOKES, NSTOKES>, 0,
                           Eigen::InnerStride<>>
                    d_minus_matrix(
                        &dual_lpsum_minus.deriv(k, j * NSTOKES * NSTOKES),
                        NSTOKES, NSTOKES,
                        Eigen::InnerStride<>(
                            input_derivs.numDerivativeLayer(p)));

                d_plus_matrix.setZero();
                d_minus_matrix.setZero();

                for (int l = m; l < coeffs.size(); ++l) {
                    lp_triple_product<NSTOKES, CNSTR, false>(
                        deriv.d_legendre_coeff[l], lp_mu[l], lp[j][l], l, m,
                        d_plus_matrix);
                    lp_triple_product<NSTOKES, CNSTR, true>(
                        deriv.d_legendre_coeff[l], lp_mu[l], lp[j][l], l, m,
                        d_minus_matrix);
                }

                d_plus_matrix.array() *= 0.5 * ssa.value * weights[j];
                d_minus_matrix.array() *= 0.5 * ssa.value * weights[j];

                d_plus_matrix.array() +=
                    plus_matrix.array() * deriv.d_SSA * 0.5 * weights[j];
                d_minus_matrix.array() +=
                    minus_matrix.array() * deriv.d_SSA * 0.5 * weights[j];
            }

            plus_matrix.array() *= 0.5 * ssa.value * weights[j];
            minus_matrix.array() *= 0.5 * ssa.value * weights[j];
        }
    }

    template <int NSTOKES, int CNSTR, bool negation>
    inline void
    single_scat_st(const std::vector<LegendreCoefficient<NSTOKES>>& coeffs,
                   const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_mu,
                   const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_csz,
                   int m, LayerIndex p, const LayerDual<double>& ssa,
                   double direct_intesity,
                   const InputDerivatives<NSTOKES>& input_derivs, double* value,
                   double* d_value, int nderiv) {
        Eigen::Matrix<double, NSTOKES, NSTOKES> result, d_result;
        result.setZero();

        double factor =
            (2 - kronDelta(m, 0)) * (1.0 / (4.0 * EIGEN_PI)) * direct_intesity;
        for (int l = m; l < coeffs.size(); ++l) {
            lp_triple_product<NSTOKES, CNSTR, negation>(
                coeffs[l], lp_mu[l], lp_csz[l], l, m, result);
        }
        result *= factor;

        for (int k = 0; k < input_derivs.numDerivativeLayer(p); ++k) {
            const LayerInputDerivative<NSTOKES>& deriv =
                input_derivs
                    .layerDerivatives()[input_derivs.layerStartIndex(p) + k];

            d_result.setZero();
            for (int l = m; l < coeffs.size(); ++l) {
                lp_triple_product<NSTOKES, CNSTR, negation>(
                    deriv.d_legendre_coeff[l], lp_mu[l], lp_csz[l], l, m,
                    d_result);
            }
            d_result *= factor * ssa.value;

            d_result.array() += result.array() * deriv.d_SSA;

            Eigen::Map<Eigen::Matrix<double, NSTOKES, 1>, 0,
                       Eigen::InnerStride<>>
                map_d_value(d_value + k, NSTOKES, 1,
                            Eigen::InnerStride<>(nderiv));
            map_d_value.array() = d_result(Eigen::all, 0).array();
        }

        Eigen::Map<Eigen::Matrix<double, NSTOKES, 1>> map_value(value);
        map_value.array() = result(Eigen::all, 0).array() * ssa.value;
    }

} // namespace sasktran_disco
