#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_linearization_types.h"
#include "sktran_disco/sktran_do_opticallayer.h"
#include "sktran_disco/sktran_do_lpproduct.h"

template <int NSTOKES, int CNSTR>
sasktran_disco::OpticalLayer<NSTOKES, CNSTR>::OpticalLayer(
    const PersistentConfiguration<NSTOKES, CNSTR>& config, LayerIndex index,
    double altitude_ceiling, double altitude_floor,
    const InputDerivatives<NSTOKES>& input_derivs,
    const ThreadData<NSTOKES, CNSTR>& thread_data)
    : OpticalLayerROP<NSTOKES>(config), M_ALT_CEILING(altitude_ceiling),
      M_ALT_FLOOR(altitude_floor), M_INDEX(index),
      m_solutions(thread_data.rte_solution(index)),
      m_input_derivs(input_derivs),
      m_layercache(thread_data.layer_cache(index)),
      m_postprocessing_cache(thread_data.postprocessing_cache(index)),
      m_dual_thickness(m_layercache.dual_thickness),
      m_dual_ssa(m_layercache.dual_ssa),
      m_average_secant(m_layercache.average_secant) {
    m_lephasef = std::make_unique<std::vector<LegendreCoefficient<NSTOKES>>>(
        config.nstr());
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayer<NSTOKES, CNSTR>::set_optical(
    double scat_ext, double tot_ext,
    const VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>&
        phasef_expansion,
    double optical_depth_ceiling, double optical_depth_floor) {
    M_SCAT_EXT = scat_ext;
    M_TOT_EXT = tot_ext;
    M_OPTICALDEPTH_CEILING = optical_depth_ceiling;
    M_OPTICALDEPTH_FLOOR = optical_depth_floor;

    M_OPTICAL_THICKNESS = optical_depth_floor - optical_depth_ceiling;

    for (int i = 0; i < phasef_expansion.size(); ++i) {
        (*m_lephasef)[i] = phasef_expansion[i];
    }

    M_SSA = scat_ext / tot_ext;

    double ssa_dither = this->m_userspec->getSSAEqual1Dither();
    if (1 - M_SSA < ssa_dither) {
        const_cast<double&>(M_SSA) = 1 - ssa_dither;
    }
}

template <int NSTOKES, int CNSTR>
sasktran_disco::OpticalLayer<NSTOKES, CNSTR>::OpticalLayer(
    const PersistentConfiguration<NSTOKES, CNSTR>& config, LayerIndex index,
    double scat_ext, double tot_ext,
    std::unique_ptr<VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>>
        phasef_expansion,
    double optical_depth_ceiling, double optical_depth_floor,
    double altitude_ceiling, double altitude_floor,
    const InputDerivatives<NSTOKES>& input_derivatives)
    : OpticalLayerROP<NSTOKES>(config), M_SSA(scat_ext / tot_ext),
      M_SCAT_EXT(scat_ext), M_TOT_EXT(tot_ext),
      M_OPTICALDEPTH_CEILING(optical_depth_ceiling),
      M_OPTICALDEPTH_FLOOR(optical_depth_floor),
      M_OPTICAL_THICKNESS(optical_depth_floor - optical_depth_ceiling),
      M_ALT_CEILING(altitude_ceiling), M_ALT_FLOOR(altitude_floor),
      m_lephasef(std::move(phasef_expansion)), M_INDEX(index),
      m_solutions(config.pool().thread_data().rte_solution(index)),
      m_input_derivs(input_derivatives),
      m_layercache(config.pool().thread_data().layer_cache(index)),
      m_postprocessing_cache(
          config.pool().thread_data().postprocessing_cache(index)),
      m_dual_thickness(m_layercache.dual_thickness),
      m_dual_ssa(m_layercache.dual_ssa),
      m_average_secant(m_layercache.average_secant) {
    // Check that SSA is not approximately zero. If so then emit warning once
    // and then dither SSA.
    double ssa_dither = this->m_userspec->getSSAEqual1Dither();
    if (1 - M_SSA < ssa_dither) {
        const_cast<double&>(M_SSA) = 1 - ssa_dither;
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayer<NSTOKES, CNSTR>::configureDerivative() {
    m_dual_thickness.resize(m_input_derivs.numDerivativeLayer(M_INDEX));
    m_dual_thickness.layer_index = M_INDEX;
    m_dual_thickness.layer_start =
        (uint)m_input_derivs.layerStartIndex(M_INDEX);

    m_dual_ssa.resize(m_input_derivs.numDerivativeLayer(M_INDEX));
    m_dual_ssa.layer_index = M_INDEX;
    m_dual_ssa.layer_start = (uint)m_input_derivs.layerStartIndex(M_INDEX);

    m_dual_thickness.value = M_OPTICAL_THICKNESS;
    m_dual_ssa.value = M_SSA;
    for (uint i = 0; i < m_input_derivs.numDerivativeLayer(M_INDEX); ++i) {
        m_dual_thickness.deriv(i) =
            m_input_derivs.layerDerivatives()[m_dual_thickness.layer_start + i]
                .d_optical_depth;
        m_dual_ssa.deriv(i) =
            m_input_derivs.layerDerivatives()[m_dual_ssa.layer_start + i].d_SSA;
    }

    LayerIndex p = index();
    uint numderiv = (uint)m_input_derivs.numDerivativeLayer(p);
    uint numtotalderiv = (uint)m_input_derivs.numDerivative();
    uint layerStart = (uint)m_input_derivs.layerStartIndex(p);
    m_postprocessing_cache.resize(this->M_NSTR, p, numderiv, layerStart,
                                  numtotalderiv);
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayer<NSTOKES, CNSTR>::integrate_source(
    AEOrder m, double mu, double obsod,
    const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_mu,
    Radiance<NSTOKES>& result,
    ReverseLinearizationTrace<NSTOKES>& reverse_trace,
    const sasktran_disco::Radiance<NSTOKES>* manual_ss_source,
    bool include_ss) const {
    LayerIndex p = index();
    uint numderiv = (uint)m_input_derivs.numDerivativeLayer(p);
    uint numtotalderiv = (uint)m_input_derivs.numDerivative();
    uint layerStart = (uint)m_input_derivs.layerStartIndex(p);

    double exit_optical_depth = mu > 0 ? opticalDepth(Location::CEILING)
                                       : opticalDepth(Location::FLOOR);

    double upper_bound;
    if (obsod < exit_optical_depth) {
        upper_bound = exit_optical_depth;
    } else {
        upper_bound = obsod;
    }

    double x = upper_bound - exit_optical_depth;
    const auto& transmission = ceiling_beam_transmittanc();

    const auto& average_secant = dual_average_secant();
    const LayerDual<double>& dual_thickness = this->dual_thickness();

    // Calculate legendre sums multiplied by stream weights
    VectorLayerDual<double>& dual_lpsum_plus =
        m_postprocessing_cache.dual_lpsum_plus;
    VectorLayerDual<double>& dual_lpsum_minus =
        m_postprocessing_cache.dual_lpsum_minus;

    scat_phase_f<NSTOKES, CNSTR>(*m_lephasef, lp_mu, (*this->M_LP_MU)[m], m, p,
                                 m_dual_ssa, (*this->M_WT), m_input_derivs,
                                 dual_lpsum_minus, dual_lpsum_plus);

    const auto& solution = m_solutions[m];

    const auto& dual_Aplus = solution.value.dual_green_A_plus();
    const auto& dual_Aminus = solution.value.dual_green_A_minus();

    // Eq (44)
    auto& V = m_postprocessing_cache.V;

    // Eq (43)
    auto& Q = m_postprocessing_cache.Q;

    if (include_ss) {
        Q.deriv.setZero();
        if constexpr (NSTOKES == 1) {
            single_scat_st<NSTOKES, CNSTR, true>(
                *m_lephasef, lp_mu, (*this->M_LP_CSZ)[m], m, p, m_dual_ssa, 1.0,
                m_input_derivs, &Q.value, &Q.deriv(layerStart, 0),
                Q.deriv.rows());
        } else {
            single_scat_st<NSTOKES, CNSTR, true>(
                *m_lephasef, lp_mu, (*this->M_LP_CSZ)[m], m, p, m_dual_ssa, 1.0,
                m_input_derivs, &Q.value(0), &Q.deriv(layerStart, 0),
                Q.deriv.rows());
        }
    } else {
        if (manual_ss_source == nullptr) {
            Q.setzero();
        } else {
            // Can just reassign Q to the manual source
            Q.value = (*manual_ss_source).value;
            Q.deriv = (*manual_ss_source).deriv;
        }
    }

    // could maybe use the full memory somehow? not sure
    auto& hp = m_postprocessing_cache.hp[0];
    auto& hm = m_postprocessing_cache.hm[0];
    auto& J = m_postprocessing_cache.J;

    auto& Y_plus = m_postprocessing_cache.Y_plus;
    auto& Y_minus = m_postprocessing_cache.Y_minus;

    using MatrixView = Eigen::Map<Eigen::MatrixXd>;
    using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;
    MatrixView Y_plus_matrix(Y_plus.value.data(), NSTOKES,
                             this->M_NSTR / 2 * NSTOKES);
    MatrixView Y_minus_matrix(Y_minus.value.data(), NSTOKES,
                              this->M_NSTR / 2 * NSTOKES);

    ConstMatrixView lpsum_plus_matrix(dual_lpsum_plus.value.data(), NSTOKES,
                                      this->M_NSTR / 2 * NSTOKES);
    ConstMatrixView lpsum_minus_matrix(dual_lpsum_minus.value.data(), NSTOKES,
                                       this->M_NSTR / 2 * NSTOKES);

    ConstMatrixView homog_plus_matrix(
        solution.value.dual_homog_plus().value.data(),
        this->M_NSTR / 2 * NSTOKES, this->M_NSTR / 2 * NSTOKES);
    ConstMatrixView homog_minus_matrix(
        solution.value.dual_homog_minus().value.data(),
        this->M_NSTR / 2 * NSTOKES, this->M_NSTR / 2 * NSTOKES);

    Y_plus_matrix = lpsum_plus_matrix * homog_plus_matrix +
                    lpsum_minus_matrix * homog_minus_matrix;
    Y_minus_matrix = lpsum_plus_matrix * homog_minus_matrix +
                     lpsum_minus_matrix * homog_plus_matrix;

    // Assign derivatives of Y_plus and Y_minus
    // Only have layer derivatives
    std::vector<
        Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>>
        Y_plus_deriv;
    std::vector<
        Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>>
        Y_minus_deriv;

    Y_plus_deriv.reserve(numderiv);
    Y_minus_deriv.reserve(numderiv);

    for (uint k = 0; k < numderiv; ++k) {
        Y_plus_deriv.emplace_back(&Y_plus.deriv(k, 0), NSTOKES,
                                  this->M_NSTR / 2 * NSTOKES,
                                  Eigen::InnerStride<>(numderiv));
        Y_minus_deriv.emplace_back(&Y_minus.deriv(k, 0), NSTOKES,
                                   this->M_NSTR / 2 * NSTOKES,
                                   Eigen::InnerStride<>(numderiv));

        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>
            lpsum_plus_deriv(&dual_lpsum_plus.deriv(k, 0), NSTOKES,
                             this->M_NSTR / 2 * NSTOKES,
                             Eigen::InnerStride<>(numderiv));

        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>
            lpsum_minus_deriv(&dual_lpsum_minus.deriv(k, 0), NSTOKES,
                              this->M_NSTR / 2 * NSTOKES,
                              Eigen::InnerStride<>(numderiv));

        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>
            homog_minus_deriv(&solution.value.dual_homog_minus().deriv(k, 0),
                              this->M_NSTR / 2 * NSTOKES,
                              this->M_NSTR / 2 * NSTOKES,
                              Eigen::InnerStride<>(numderiv));

        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>
            homog_plus_deriv(&solution.value.dual_homog_plus().deriv(k, 0),
                             this->M_NSTR / 2 * NSTOKES,
                             this->M_NSTR / 2 * NSTOKES,
                             Eigen::InnerStride<>(numderiv));

        Y_plus_deriv[k].noalias() = lpsum_plus_deriv * homog_plus_matrix +
                                    lpsum_plus_matrix * homog_plus_deriv +
                                    lpsum_minus_deriv * homog_minus_matrix +
                                    lpsum_minus_matrix * homog_minus_deriv;

        Y_minus_deriv[k].noalias() = lpsum_plus_deriv * homog_minus_matrix +
                                     lpsum_plus_matrix * homog_minus_deriv +
                                     lpsum_minus_deriv * homog_plus_matrix +
                                     lpsum_minus_matrix * homog_plus_deriv;
    }

    auto& Dm = m_postprocessing_cache.Dm[0];
    auto& Dp = m_postprocessing_cache.Dp[0];
    auto& Eform = m_postprocessing_cache.Eform[0];

    E(mu, x, M_OPTICAL_THICKNESS, transmission, Eform);

    const auto& dual_L = solution.boundary.L_coeffs;
    const auto& dual_M = solution.boundary.M_coeffs;

    J.setzero();
    V.setzero();

    int l_offset = p * this->M_NSTR * NSTOKES;
    int m_offset = l_offset + this->M_NSTR * NSTOKES / 2;

    for (SolutionIndex i = 0; i < this->M_NSTR / 2 * NSTOKES; ++i) {
        h_plus(m, mu, x, M_OPTICAL_THICKNESS, i, hp);
        h_minus(m, mu, x, M_OPTICAL_THICKNESS, i, hm);

        if constexpr (NSTOKES == 1) {
            J.value += Y_plus_matrix(0, i) * hp.value * dual_L.value(i);
            J.value += Y_minus_matrix(0, i) * hm.value * dual_M.value(i);
        } else {
            J.value +=
                Y_plus_matrix(Eigen::all, i) * hp.value * dual_L.value(i);
            J.value +=
                Y_minus_matrix(Eigen::all, i) * hm.value * dual_M.value(i);
        }

        // Y and Hp/Hm only have layer derivatives
        for (uint k = 0; k < numderiv; ++k) {
            J.deriv(layerStart + k, Eigen::all).noalias() +=
                Y_plus_deriv[k](Eigen::all, i) * hp.value * dual_L.value(i);
            J.deriv(layerStart + k, Eigen::all).noalias() +=
                Y_minus_deriv[k](Eigen::all, i) * hm.value * dual_M.value(i);

            J.deriv(layerStart + k, Eigen::all).noalias() +=
                Y_plus_matrix(Eigen::all, i) * hp.deriv(k) * dual_L.value(i);
            J.deriv(layerStart + k, Eigen::all).noalias() +=
                Y_minus_matrix(Eigen::all, i) * hm.deriv(k) * dual_M.value(i);
        }

        if (this->M_BACKPROP_BVP) {
            reverse_trace.bvp_coeff_weights()(l_offset + i, Eigen::all) +=
                Y_plus_matrix(Eigen::all, i) * hp.value;
            reverse_trace.bvp_coeff_weights()(m_offset + i, Eigen::all) +=
                Y_minus_matrix(Eigen::all, i) * hm.value;
        } else {
            // But L/M have full derivatives
            for (uint k = 0; k < numtotalderiv; ++k) {
                J.deriv(k, Eigen::all).noalias() +=
                    Y_plus_matrix(Eigen::all, i) * hp.value *
                    dual_L.deriv(k, i);
                J.deriv(k, Eigen::all).noalias() +=
                    Y_minus_matrix(Eigen::all, i) * hm.value *
                    dual_M.deriv(k, i);
            }
        }

        const auto& eigval = solution.value.dual_eigval();

        double expfactor =
            exp(-1 * (dual_thickness.value - x) * average_secant.value);

        Dp.value =
            (-1 * transmission.value * expfactor * hm.value + Eform.value) /
            (average_secant.value + solution.value.dual_eigval().value(i));
        Dm.value =
            (transmission.value * hp.value - Eform.value) /
            (average_secant.value - solution.value.dual_eigval().value(i));

        Dp.deriv.noalias() =
            average_secant.deriv *
            (transmission.value * hm.value * dual_thickness.value * expfactor -
             Dp.value) /
            (average_secant.value + eigval.value(i));
        Dp.deriv.noalias() +=
            Eform.deriv / (average_secant.value + eigval.value(i));
        Dp.deriv.noalias() += transmission.deriv * (-hm.value * expfactor) /
                              (average_secant.value + eigval.value(i));

        Dm.deriv.noalias() = average_secant.deriv * (-Dm.value) /
                             (average_secant.value - eigval.value(i));
        Dm.deriv.noalias() += transmission.deriv * (hp.value) /
                              (average_secant.value - eigval.value(i));
        Dm.deriv.noalias() +=
            Eform.deriv * -1.0 / (average_secant.value - eigval.value(i));

        for (uint k = 0; k < numderiv; ++k) {
            Dp.deriv(layerStart + k) +=
                eigval.deriv(k, i) * (-Dp.value) /
                (average_secant.value + eigval.value(i));
            Dp.deriv(layerStart + k) +=
                hm.deriv(k) * (-transmission.value * expfactor) /
                (average_secant.value + eigval.value(i));
            Dp.deriv(layerStart + k) +=
                (dual_thickness.value - x) / dual_thickness.value *
                dual_thickness.deriv(k) *
                (transmission.value * average_secant.value * hm.value *
                 expfactor) /
                (average_secant.value + eigval.value(i));

            Dm.deriv(layerStart + k) +=
                hp.deriv(k) * transmission.value /
                (average_secant.value - eigval.value(i));
            Dm.deriv(layerStart + k) +=
                eigval.deriv(k, i) * (Dm.value) /
                (average_secant.value - eigval.value(i));
        }

        if constexpr (NSTOKES == 1) {
            V.value += dual_Aplus.value(i) * Y_plus_matrix(0, i) * Dm.value +
                       dual_Aminus.value(i) * Y_minus_matrix(0, i) * Dp.value;
        } else {
            V.value +=
                dual_Aplus.value(i) * Y_plus_matrix(Eigen::all, i) * Dm.value +
                dual_Aminus.value(i) * Y_minus_matrix(Eigen::all, i) * Dp.value;
        }

        // Y and A only have layer derivatives
        for (uint k = 0; k < numderiv; ++k) {
            V.deriv(layerStart + k, Eigen::all).noalias() +=
                dual_Aplus.value(i) * Y_plus_deriv[k](Eigen::all, i) * Dm.value;
            V.deriv(layerStart + k, Eigen::all).noalias() +=
                dual_Aminus.value(i) * Y_minus_deriv[k](Eigen::all, i) *
                Dp.value;

            V.deriv(layerStart + k, Eigen::all).noalias() +=
                dual_Aplus.deriv(k, i) * Y_plus_matrix(Eigen::all, i) *
                Dm.value;
            V.deriv(layerStart + k, Eigen::all).noalias() +=
                dual_Aminus.deriv(k, i) * Y_minus_matrix(Eigen::all, i) *
                Dp.value;
        }
        // But D has full derivatives
        for (uint k = 0; k < numtotalderiv; ++k) {
            V.deriv(k, Eigen::all).noalias() += dual_Aplus.value(i) *
                                                Y_plus_matrix(Eigen::all, i) *
                                                Dm.deriv(k);
            V.deriv(k, Eigen::all).noalias() += dual_Aminus.value(i) *
                                                Y_minus_matrix(Eigen::all, i) *
                                                Dp.deriv(k);
        }
    }

    result.value = J.value + V.value + Q.value * Eform.value;

    result.deriv.noalias() = J.deriv + V.deriv + Q.deriv * Eform.value;
    for (uint k = 0; k < numtotalderiv; ++k) {
        if constexpr (NSTOKES == 1) {
            result.deriv(k) += Q.value * Eform.deriv(k);
        } else {
            result.deriv(k, Eigen::all).noalias() += Q.value * Eform.deriv(k);
        }
    }
}

#ifdef SASKTRAN_DISCO_FULL_COMPILE
template <>
void sasktran_disco::OpticalLayer<1, 2>::integrate_source(
    AEOrder m, double mu, double obsod,
    const std::vector<LegendrePhaseContainer<1>>& lp_mu, Radiance<1>& result,
    ReverseLinearizationTrace<1>& reverse_trace,
    const sasktran_disco::Radiance<1>* manual_ss_source,

    bool include_ss) const {
    // TODO: Harmonize this with the other integrate_source function, refactor a
    // lot of this calculation

    // Note that for the source function integration we switch convention for
    // upwelling/downwelling from the solution convention
    LayerIndex p = index();
    uint numderiv = (uint)m_input_derivs.numDerivativeLayer(p);
    uint numtotalderiv = (uint)m_input_derivs.numDerivative();
    uint layerStart = (uint)m_input_derivs.layerStartIndex(p);

    const auto& solution = m_solutions[m];
    auto& cache = m_postprocessing_cache;

    double exit_optical_depth = mu > 0 ? opticalDepth(Location::CEILING)
                                       : opticalDepth(Location::FLOOR);

    double upper_bound;
    if (obsod < exit_optical_depth) {
        upper_bound = exit_optical_depth;
    } else {
        upper_bound = obsod;
    }

    double x = upper_bound - exit_optical_depth;
    const auto& transmission = ceiling_beam_transmittanc();

    const auto& average_secant = dual_average_secant();
    const LayerDual<double>& dual_thickness = this->dual_thickness();

    // Calculate legendre sums multiplied by stream weights
    VectorLayerDual<double>& dual_lpsum_plus = cache.dual_lpsum_plus;
    VectorLayerDual<double>& dual_lpsum_minus = cache.dual_lpsum_minus;

    scat_phase_f<1, 2>(*m_lephasef, lp_mu, (*this->M_LP_MU)[m], m, p,
                       m_dual_ssa, (*this->M_WT), m_input_derivs,
                       dual_lpsum_minus, dual_lpsum_plus);

    const auto& dual_Aplus = solution.value.dual_green_A_plus();
    const auto& dual_Aminus = solution.value.dual_green_A_minus();

    // Eq (44)
    auto& V = cache.V;

    // Eq (43)
    auto& Q = cache.Q;

    if (include_ss) {
        Q.deriv.setZero();
        single_scat_st<1, 2, true>(*m_lephasef, lp_mu, (*this->M_LP_CSZ)[m], m,
                                   p, m_dual_ssa, 1.0, m_input_derivs, &Q.value,
                                   &Q.deriv(layerStart, 0), Q.deriv.rows());
    } else {
        if (manual_ss_source == nullptr) {
            Q.setzero();
        } else {
            // Can just reassign Q to the manual source
            Q.value = (*manual_ss_source).value;
            Q.deriv = (*manual_ss_source).deriv;
        }
    }

    auto& hp = cache.hp[0];
    auto& hm = cache.hm[0];
    auto& J = cache.J;

    auto& Y_plus = cache.Y_plus;
    auto& Y_minus = cache.Y_minus;

    Y_plus.value(0) =
        dual_lpsum_plus.value(0) * solution.value.dual_homog_plus().value(0) +
        dual_lpsum_minus.value(0) * solution.value.dual_homog_minus().value(0);
    Y_minus.value(0) =
        dual_lpsum_plus.value(0) * solution.value.dual_homog_minus().value(0) +
        dual_lpsum_minus.value(0) * solution.value.dual_homog_plus().value(0);

    // Assign derivatives of Y_plus and Y_minus
    // Only have layer derivatives

    for (uint k = 0; k < numderiv; ++k) {
        Y_plus.deriv(k, 0) = dual_lpsum_plus.deriv(k, 0) *
                                 solution.value.dual_homog_plus().value(0) +
                             dual_lpsum_minus.deriv(k, 0) *
                                 solution.value.dual_homog_minus().value(0) +
                             dual_lpsum_plus.value(0) *
                                 solution.value.dual_homog_plus().deriv(k, 0) +
                             dual_lpsum_minus.value(0) *
                                 solution.value.dual_homog_minus().deriv(k, 0);

        Y_minus.deriv(k, 0) =
            dual_lpsum_plus.deriv(k, 0) *
                solution.value.dual_homog_minus().value(0) +
            dual_lpsum_minus.deriv(k, 0) *
                solution.value.dual_homog_plus().value(0) +
            dual_lpsum_plus.value(0) *
                solution.value.dual_homog_minus().deriv(k, 0) +
            dual_lpsum_minus.value(0) *
                solution.value.dual_homog_plus().deriv(k, 0);
    }

    auto& Dm = cache.Dm[0];
    auto& Dp = cache.Dp[0];
    auto& Eform = cache.Eform[0];

    E(mu, x, M_OPTICAL_THICKNESS, transmission, Eform);

    const auto& dual_L = solution.boundary.L_coeffs;
    const auto& dual_M = solution.boundary.M_coeffs;

    J.setzero();
    V.setzero();

    h_plus(m, mu, x, M_OPTICAL_THICKNESS, 0, hp);
    h_minus(m, mu, x, M_OPTICAL_THICKNESS, 0, hm);

    J.value = Y_plus.value(0) * hp.value * dual_L.value(0);
    J.value += Y_minus.value(0) * hm.value * dual_M.value(0);

    // Y and Hp/Hm only have layer derivatives
    for (uint k = 0; k < numderiv; ++k) {
        J.deriv(layerStart + k, 0) +=
            Y_plus.deriv(k, 0) * hp.value * dual_L.value(0);
        J.deriv(layerStart + k, 0) +=
            Y_minus.deriv(k, 0) * hm.value * dual_M.value(0);

        J.deriv(layerStart + k, 0) +=
            Y_plus.value(0) * hp.deriv(k) * dual_L.value(0);
        J.deriv(layerStart + k, 0) +=
            Y_minus.value(0) * hm.deriv(k) * dual_M.value(0);
    }

    if (this->M_BACKPROP_BVP) {
        reverse_trace.bvp_coeff_weights()(p * 2, 0) +=
            Y_plus.value(0) * hp.value;
        reverse_trace.bvp_coeff_weights()(p * 2 + 1, 0) +=
            Y_minus.value(0) * hm.value;
    } else {
        // But L/M have full derivatives
        for (uint k = 0; k < numtotalderiv; ++k) {
            J.deriv(k, 0) += Y_plus.value(0) * hp.value * dual_L.deriv(k, 0);
            J.deriv(k, 0) += Y_minus.value(0) * hm.value * dual_M.deriv(k, 0);
        }
    }

    const auto& eigval = solution.value.dual_eigval();

    double expfactor =
        exp(-1 * (dual_thickness.value - x) * average_secant.value);

    Dp.value = (-1 * transmission.value * expfactor * hm.value + Eform.value) /
               (average_secant.value + solution.value.dual_eigval().value(0));
    Dm.value = (transmission.value * hp.value - Eform.value) /
               (average_secant.value - solution.value.dual_eigval().value(0));

    Dp.deriv.noalias() =
        average_secant.deriv *
        (transmission.value * hm.value * dual_thickness.value * expfactor -
         Dp.value) /
        (average_secant.value + eigval.value(0));
    Dp.deriv.noalias() +=
        Eform.deriv / (average_secant.value + eigval.value(0));
    Dp.deriv.noalias() += transmission.deriv * (-hm.value * expfactor) /
                          (average_secant.value + eigval.value(0));

    Dm.deriv.noalias() = average_secant.deriv * (-Dm.value) /
                         (average_secant.value - eigval.value(0));
    Dm.deriv.noalias() += transmission.deriv * (hp.value) /
                          (average_secant.value - eigval.value(0));
    Dm.deriv.noalias() +=
        Eform.deriv * -1.0 / (average_secant.value - eigval.value(0));

    for (uint k = 0; k < numderiv; ++k) {
        Dp.deriv(layerStart + k) += eigval.deriv(k, 0) * (-Dp.value) /
                                    (average_secant.value + eigval.value(0));
        Dp.deriv(layerStart + k) += hm.deriv(k) *
                                    (-transmission.value * expfactor) /
                                    (average_secant.value + eigval.value(0));
        Dp.deriv(layerStart + k) +=
            (dual_thickness.value - x) / dual_thickness.value *
            dual_thickness.deriv(k) *
            (transmission.value * average_secant.value * hm.value * expfactor) /
            (average_secant.value + eigval.value(0));

        Dm.deriv(layerStart + k) += hp.deriv(k) * transmission.value /
                                    (average_secant.value - eigval.value(0));
        Dm.deriv(layerStart + k) += eigval.deriv(k, 0) * (Dm.value) /
                                    (average_secant.value - eigval.value(0));
    }

    V.value += dual_Aplus.value(0) * Y_plus.value(0) * Dm.value +
               dual_Aminus.value(0) * Y_minus.value(0) * Dp.value;

    // Y and A only have layer derivatives
    for (uint k = 0; k < numderiv; ++k) {
        V.deriv(layerStart + k, 0) +=
            dual_Aplus.value(0) * Y_plus.deriv(k, 0) * Dm.value;
        V.deriv(layerStart + k, 0) +=
            dual_Aminus.value(0) * Y_minus.deriv(k, 0) * Dp.value;

        V.deriv(layerStart + k, 0) +=
            dual_Aplus.deriv(k, 0) * Y_plus.value(0) * Dm.value;
        V.deriv(layerStart + k, 0) +=
            dual_Aminus.deriv(k, 0) * Y_minus.value(0) * Dp.value;
    }
    // But D has full derivatives
    for (uint k = 0; k < numtotalderiv; ++k) {
        V.deriv(k, 0) += dual_Aplus.value(0) * Y_plus.value(0) * Dm.deriv(k);
        V.deriv(k, 0) += dual_Aminus.value(0) * Y_minus.value(0) * Dp.deriv(k);
    }

    result.value = J.value + V.value + Q.value * Eform.value;

    result.deriv.noalias() = J.deriv + V.deriv + Q.deriv * Eform.value;
    for (uint k = 0; k < numtotalderiv; ++k) {
        result.deriv(k) += Q.value * Eform.deriv(k);
    }
}
#endif

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayer<NSTOKES, CNSTR>::h_plus(
    AEOrder m, double mu, double x, double thickness, SolutionIndex j,
    LayerDual<HomogType>& xform) const {
    uint p = index();

    auto v_eigval = m_solutions[m].value.dual_eigval().value(j);
    const auto& d_eigval =
        m_solutions[m].value.dual_eigval().deriv(Eigen::all, j);

    const LayerDual<double>& dual_thickness = this->dual_thickness();

    // Exit depth is x, 0 is exit at the top of the layer
    // h_layer = h_low + (1 - x/thickness) * (h_upper - h_lower)
    // 0 = -dx * lh / thickness - x/(thickness*thickness) * dthickness * lh
    // dx = x/thickness * dthickness

    double layerfraction = 1 - x / thickness;
    HomogType den = (1.0 + mu * v_eigval);

    if (std::abs(den) > 0.0001) {
        HomogType e1 = std::exp(-x * v_eigval);
        HomogType e2 = std::exp(-1.0 * dual_thickness.value * v_eigval) *
                       std::exp(-(dual_thickness.value - x) / mu);

        xform.value = (e1 - e2) / (den);

        if (xform.deriv.size() > 0) {
            // Start with the numerator derivative
            xform.deriv.noalias() =
                (-1.0 * e1 *
                     (d_eigval * x +
                      dual_thickness.deriv * (1 - layerfraction) * v_eigval) +
                 e2 * ((v_eigval + layerfraction / mu) * dual_thickness.deriv +
                       dual_thickness.value * d_eigval)) /
                den;

            // Then add in the denominator derivative
            xform.deriv -= 1.0 / (den)*xform.value * mu * d_eigval;
        }
    } else { // limit as eigval -> -1/mu, take two terms in the taylor series
             // expansion
        HomogType e1 = std::exp(-x * v_eigval);
        xform.value = (dual_thickness.value - x) / mu * e1 *
                      (HomogType(1) - (dual_thickness.value - x) *
                                          (v_eigval + HomogType(1) / mu));
        if (xform.deriv.size() > 0) {
            xform.deriv.noalias() = -1.0 * x * d_eigval * xform.value;
            xform.deriv += 1 / mu * e1 * dual_thickness.deriv *
                           (HomogType(1) - (dual_thickness.value - x) *
                                               (v_eigval + HomogType(1) / mu));
            xform.deriv += 1 / mu * e1 * (dual_thickness.value - x) * -1.0 *
                           (v_eigval + 1 / mu) * dual_thickness.deriv;
        }

        // xform = dual_thickness / mu;
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayer<NSTOKES, CNSTR>::h_minus(
    AEOrder m, double mu, double x, double thickness, SolutionIndex j,
    LayerDual<HomogType>& xform) const {
    uint p = index();

    auto v_eigval = m_solutions[m].value.dual_eigval().value(j);
    const auto& d_eigval =
        m_solutions[m].value.dual_eigval().deriv(Eigen::all, j);

    const LayerDual<double>& dual_thickness = this->dual_thickness();
    double layerfraction = 1 - x / thickness;

    HomogType den = (1.0 - mu * v_eigval);

    if (std::abs(den) > 0.0001) {
        HomogType e1 = std::exp(-1.0 * v_eigval * (dual_thickness.value - x));
        HomogType e2 = std::exp(-1 / mu * (dual_thickness.value - x));

        xform.value = (e1 - e2) / den;

        if (xform.deriv.size() > 0) {
            xform.deriv.noalias() =
                (-1.0 * e1 *
                     (d_eigval * (dual_thickness.value - x) +
                      dual_thickness.deriv * layerfraction * v_eigval) +
                 e2 * dual_thickness.deriv * layerfraction / mu) /
                den;
            xform.deriv += 1.0 / den * xform.value * mu * d_eigval;
        }

    } else { // limit as eigval -> -1/mu, take two terms in the taylor series
             // expansion
        HomogType e1 = std::exp(-1.0 * v_eigval * (dual_thickness.value - x));

        xform.value = e1 * (dual_thickness.value - x) / mu *
                      (HomogType(1) - (dual_thickness.value - x) *
                                          (v_eigval - HomogType(1) / mu));
        if (xform.deriv.size() > 0) {
            xform.deriv.noalias() = -1.0 *
                                    (v_eigval * dual_thickness.deriv +
                                     dual_thickness.value * d_eigval) *
                                    xform.value;
            xform.deriv.noalias() +=
                e1 * dual_thickness.deriv / mu *
                (HomogType(1) -
                 (dual_thickness.value - x) * (v_eigval - HomogType(1) / mu));
            xform.deriv.noalias() += e1 * (dual_thickness.value - x) / mu *
                                     -1.0 * dual_thickness.deriv *
                                     (v_eigval - 1 / mu);
        }
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayer<NSTOKES, CNSTR>::E(
    double mu, double x, double thickness, const Dual<double>& transmission,
    Dual<double>& xform) const {
    uint p = index();
    uint layerStart = (uint)m_input_derivs.layerStartIndex(p);
    uint numDeriv = (uint)m_input_derivs.numDerivativeLayer(p);

    const LayerDual<double>& dual_thickness = this->dual_thickness();
    const Dual<double>& avg_secant = this->dual_average_secant();
    double layerfraction = 1 - x / thickness;

    // Start by assigning cross derivative terms of xform
    double e1 = exp(-x * avg_secant.value);
    double e2 = exp(-1.0 * dual_thickness.value * avg_secant.value) *
                exp(-1.0 * (dual_thickness.value - x) / mu);
    double den = (1.0 + mu * avg_secant.value);

    xform.value = transmission.value / den * (e1 - e2);

    if (xform.deriv.size() > 0) {
        xform.deriv.noalias() = transmission.deriv / den * (e1 - e2);

        xform.deriv.noalias() +=
            avg_secant.deriv *
            ((transmission.value / den * e2 * dual_thickness.value) -
             (xform.value / den * mu) +
             (transmission.value / den * e1 * (-1.0 * x)));
    }

    // Then add in the terms that do not contain cross derivatives
    if (numDeriv > 0) {
        xform
            .deriv(Eigen::seq(layerStart, layerStart + numDeriv - 1),
                   Eigen::all)
            .noalias() += transmission.value / den *
                          (e2 * dual_thickness.deriv *
                               (avg_secant.value + layerfraction / mu) -
                           e1 * dual_thickness.deriv * (1 - layerfraction) *
                               avg_secant.value);
    }
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::OpticalLayer);
