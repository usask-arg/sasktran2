#pragma once

#include <cstddef>
#include <sasktran2/internal_common.h>
#include "storage.h"

namespace sasktran2::twostream {

    /**
     * Solves the pentadiagonal system
     *
     * @param bvp
     * @return int
     */
    inline int pentadiagonal_solve(BVPCoeffs& bvp, Eigen::MatrixXd& rhs,
                                   bool transpose = false) {

        int N = bvp.a.size();
        // Diagonals are e,c,d,a,b
        auto& d = bvp.d;

        // auto& a = transpose? bvp.c : bvp.a;
        // auto& c = transpose? bvp.a : bvp.c;
        // auto& b = transpose? bvp.e : bvp.b;
        // auto& e = transpose? bvp.b : bvp.e;

        Eigen::Map<Eigen::VectorXd> a(
            transpose ? bvp.c.data() + 1 : bvp.a.data(), N);
        Eigen::Map<Eigen::VectorXd> c(
            transpose ? bvp.a.data() - 1 : bvp.c.data(), N);
        Eigen::Map<Eigen::VectorXd> b(
            transpose ? bvp.e.data() + 2 : bvp.b.data(), N);
        Eigen::Map<Eigen::VectorXd> e(
            transpose ? bvp.b.data() - 2 : bvp.e.data(), N);

        // auto& y = bvp.rhs;
        auto& y = rhs;
        auto& z = bvp.z;

        auto& mu = bvp.mu;
        auto& alpha = bvp.alpha;
        auto& beta = bvp.beta;
        auto& gamma = bvp.gamma;

        mu(0) = d(0);
        alpha(0) = a(0) / mu(0);

        beta(0) = b(0) / mu(0);
        z(0, Eigen::all) = y(0, Eigen::all) / mu(0);

        gamma(1) = c(1);
        mu(1) = d(1) - alpha(0) * gamma(1);
        alpha(1) = (a(1) - beta(0) * gamma(1)) / mu(1);
        beta(1) = b(1) / mu(1);
        z(1, Eigen::all) =
            (y(1, Eigen::all) - z(0, Eigen::all) * gamma(1)) / mu(1);

        for (int i = 2; i < N - 2; ++i) {
            gamma(i) = c(i) - alpha(i - 2) * e(i);
            mu(i) = d(i) - beta(i - 2) * e(i) - alpha(i - 1) * gamma(i);
            alpha(i) = (a(i) - beta(i - 1) * gamma(i)) / mu(i);
            beta(i) = b(i) / mu(i);
            z(i, Eigen::all) = (y(i, Eigen::all) - z(i - 2, Eigen::all) * e(i) -
                                z(i - 1, Eigen::all) * gamma(i)) /
                               mu(i);
        }

        if (N >= 4) {
            gamma(N - 2) = c(N - 2) - alpha(N - 4) * e(N - 2);
            mu(N - 2) =
                d(N - 2) - beta(N - 4) * e(N - 2) - alpha(N - 3) * gamma(N - 2);
            alpha(N - 2) = (a(N - 2) - beta(N - 3) * gamma(N - 2)) / mu(N - 2);

            gamma(N - 1) = c(N - 1) - alpha(N - 3) * e(N - 1);
            mu(N - 1) =
                d(N - 1) - beta(N - 3) * e(N - 1) - alpha(N - 2) * gamma(N - 1);
        }

        z(N - 2, Eigen::all) =
            (y(N - 2, Eigen::all) - z(N - 4, Eigen::all) * e(N - 2) -
             z(N - 3, Eigen::all) * gamma(N - 2)) /
            mu(N - 2);
        z(N - 1, Eigen::all) =
            (y(N - 1, Eigen::all) - z(N - 3, Eigen::all) * e(N - 1) -
             z(N - 2, Eigen::all) * gamma(N - 1)) /
            mu(N - 1);

        y(N - 1, Eigen::all) = z(N - 1, Eigen::all);
        y(N - 2, Eigen::all) =
            z(N - 2, Eigen::all) - alpha(N - 2) * y(N - 1, Eigen::all);

        for (int i = N - 3; i >= 0; --i) {
            y(i, Eigen::all) = z(i, Eigen::all) -
                               alpha(i) * y(i + 1, Eigen::all) -
                               beta(i) * y(i + 2, Eigen::all);
        }

        // There are situations where this linear system is underdetermined,
        // this is usually only the case when we are evaluating an azimuthal
        // moment m, that is greater than the largest phase expansion
        // coefficient.  We could check the determinant of the system in the
        // beginning, but it is more computationally efficient just to solve the
        // system, then if the resulting coefficients are nan just set
        // everything to 0 since there will be 0 contribution from this order
        // anyways
#ifdef SASKTRAN_DEBUG_ASSERTS
        if (y.hasNaN()) {
            spdlog::warn("Pentadiagonal solver failed");
            y.setZero();
        }
#endif

        return 0;
    }

    /**
     * @brief Constructs the homogeneous and particular solutions for each layer
     *
     * @param input
     * @param solution
     */
    inline void solve_layers(const Input& input, Solution& solution) {
        solution.homog[0].d.array() =
            (input.ssa.array() * input.b1.array() * input.mu - 1 / input.mu);
        solution.homog[0].s.array() = 1 / input.mu * (input.ssa.array() - 1);

        solution.homog[1].d.array() = -1 / input.mu;
        solution.homog[1].s.array() =
            1.0 / (2.0 * input.mu) *
            (input.ssa.array() * input.b1.array() * (1 - input.mu * input.mu) -
             2);

        // Gradients of d by ssa
        solution.homog[0].d_d_by_ssa.array() = input.b1.array() * input.mu;

        solution.homog[1].d_d_by_ssa.array() = 0;

        // Gradients of s by ssa
        solution.homog[0].d_s_by_ssa.array() = 1 / input.mu;
        solution.homog[1].d_s_by_ssa.array() = 1.0 / (2.0 * input.mu) *
                                               input.b1.array() *
                                               (1 - input.mu * input.mu);

        for (auto& homog : solution.homog) {
            // Eigenvalue and gradient
            homog.k.array() = (homog.s.array() * homog.d.array()).sqrt();
            homog.d_k_by_ssa.array() =
                (0.5 / homog.k.transpose().array()) *
                (homog.s.transpose().array() * homog.d_d_by_ssa.array() +
                 homog.d.transpose().array() * homog.d_s_by_ssa.array());

            // Eigenvectors, and gradients
            homog.X_plus.array() =
                0.5 * (1 - homog.s.array() / homog.k.array());
            homog.X_minus.array() =
                0.5 * (1 + homog.s.array() / homog.k.array());

            homog.d_X_plus_by_ssa.array() =
                (-0.5 / homog.k.transpose().array() *
                 homog.d_s_by_ssa.array()) +
                0.5 / homog.d.transpose().array() * homog.d_k_by_ssa.array();
            homog.d_X_minus_by_ssa.array() =
                (0.5 / homog.k.transpose().array() * homog.d_s_by_ssa.array()) -
                0.5 / homog.d.transpose().array() * homog.d_k_by_ssa.array();

            // Convenience quantity for later
            homog.omega.array() = (-homog.k.array() * input.od.array()).exp();
            homog.d_omega_by_ssa.array() = -homog.omega.transpose().array() *
                                           input.od.transpose().array() *
                                           homog.d_k_by_ssa.array();
        }

        // Q_plus and Q_minus
        solution.particular[0].Q_plus.array() =
            1.0 / (4 * EIGEN_PI) *
            (input.ssa.array() +
             input.b1.array() * input.ssa.array() * input.csz * input.mu);
        solution.particular[0].Q_minus.array() =
            1.0 / (4 * EIGEN_PI) *
            (input.ssa.array() -
             input.b1.array() * input.ssa.array() * input.csz * input.mu);

        solution.particular[1].Q_plus.array() =
            1.0 / (4 * EIGEN_PI) * input.ssa.array() * input.b1.array() *
            sqrt(1 - input.mu * input.mu) * sqrt(1 - input.csz * input.csz);
        solution.particular[1].Q_minus.array() =
            solution.particular[1].Q_plus.array();

        // And their gradients
        solution.particular[0].d_Q_plus_by_ssa.array() =
            1.0 / (4.0 * EIGEN_PI) *
            (1 + input.b1.transpose().array() * input.csz * input.mu);
        solution.particular[0].d_Q_minus_by_ssa.array() =
            1.0 / (4.0 * EIGEN_PI) *
            (1 - input.b1.transpose().array() * input.csz * input.mu);

        solution.particular[1].d_Q_plus_by_ssa.array() =
            1.0 / (4.0 * EIGEN_PI) * input.b1.transpose().array() *
            sqrt(1 - input.mu * input.mu) * sqrt(1 - input.csz * input.csz);
        solution.particular[1].d_Q_minus_by_ssa.array() =
            solution.particular[1].d_Q_plus_by_ssa.array();

        for (int i = 0; i < 2; ++i) {
            auto& homog = solution.homog[i];
            auto& particular = solution.particular[i];

            particular.norm.array() =
                input.mu * (homog.X_plus.array() * homog.X_plus.array() -
                            homog.X_minus.array() * homog.X_minus.array());

            particular.d_norm_by_ssa.array() =
                2 * input.mu *
                (homog.X_plus.transpose().array() *
                     homog.d_X_plus_by_ssa.array() -
                 homog.X_minus.transpose().array() *
                     homog.d_X_minus_by_ssa.array());

            // A_plus and A_minus and their gradients by ssa
            particular.A_plus.array() =
                (particular.Q_plus.array() * homog.X_plus.array() +
                 particular.Q_minus.array() * homog.X_minus.array()) /
                particular.norm.array();
            particular.A_minus.array() =
                (particular.Q_minus.array() * homog.X_plus.array() +
                 particular.Q_plus.array() * homog.X_minus.array()) /
                particular.norm.array();

            particular.d_A_plus_by_ssa.array() =
                (particular.d_Q_plus_by_ssa.array() *
                     homog.X_plus.transpose().array() +
                 particular.d_Q_minus_by_ssa.array() *
                     homog.X_minus.transpose().array() +
                 particular.Q_plus.transpose().array() *
                     homog.d_X_plus_by_ssa.array() +
                 particular.Q_minus.transpose().array() *
                     homog.d_X_minus_by_ssa.array() -
                 particular.d_norm_by_ssa.array() *
                     particular.A_plus.transpose().array()) /
                particular.norm.transpose().array();

            particular.d_A_minus_by_ssa.array() =
                (particular.d_Q_minus_by_ssa.array() *
                     homog.X_plus.transpose().array() +
                 particular.d_Q_plus_by_ssa.array() *
                     homog.X_minus.transpose().array() +
                 particular.Q_minus.transpose().array() *
                     homog.d_X_plus_by_ssa.array() +
                 particular.Q_plus.transpose().array() *
                     homog.d_X_minus_by_ssa.array() -
                 particular.d_norm_by_ssa.array() *
                     particular.A_minus.transpose().array()) /
                particular.norm.transpose().array();

            // Solution multipliers and their derivatives
            particular.C_plus.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (homog.omega.array() - input.expsec.array()) /
                (input.average_secant.array() - homog.k.array());
            particular.C_minus.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (1 - homog.omega.array() * input.expsec.array()) /
                (input.average_secant.array() + homog.k.array());

            particular.d_C_plus_by_ssa.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1))
                        .transpose()
                        .array() *
                    (homog.d_omega_by_ssa.array()) /
                    (input.average_secant.transpose().array() -
                     homog.k.transpose().array()) +
                homog.d_k_by_ssa.array() *
                    particular.C_plus.transpose().array() /
                    (input.average_secant.transpose().array() -
                     homog.k.transpose().array());

            particular.d_C_minus_by_ssa.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1))
                        .transpose()
                        .array() *
                    (-homog.d_omega_by_ssa.array() *
                     input.expsec.transpose().array()) /
                    (input.average_secant.transpose().array() +
                     homog.k.transpose().array()) -
                homog.d_k_by_ssa.array() *
                    particular.C_minus.transpose().array() /
                    (input.average_secant.transpose().array() +
                     homog.k.transpose().array());

            // Particular solutions and their derivatives
            particular.G_plus_top.array() = particular.A_minus.array() *
                                            particular.C_minus.array() *
                                            homog.X_minus.array();
            particular.G_plus_bottom.array() = particular.A_plus.array() *
                                               particular.C_plus.array() *
                                               homog.X_plus.array();
            particular.G_minus_top.array() = particular.A_minus.array() *
                                             particular.C_minus.array() *
                                             homog.X_plus.array();
            particular.G_minus_bottom.array() = particular.A_plus.array() *
                                                particular.C_plus.array() *
                                                homog.X_minus.array();

            particular.d_G_plus_top_by_ssa.array() =
                (particular.d_A_minus_by_ssa.array() *
                     particular.C_minus.transpose().array() *
                     homog.X_minus.transpose().array() +
                 particular.A_minus.transpose().array() *
                     particular.d_C_minus_by_ssa.array() *
                     homog.X_minus.transpose().array() +
                 particular.A_minus.transpose().array() *
                     particular.C_minus.transpose().array() *
                     homog.d_X_minus_by_ssa.array());

            particular.d_G_plus_bottom_by_ssa.array() =
                (particular.d_A_plus_by_ssa.array() *
                     particular.C_plus.transpose().array() *
                     homog.X_plus.transpose().array() +
                 particular.A_plus.transpose().array() *
                     particular.d_C_plus_by_ssa.array() *
                     homog.X_plus.transpose().array() +
                 particular.A_plus.transpose().array() *
                     particular.C_plus.transpose().array() *
                     homog.d_X_plus_by_ssa.array());

            particular.d_G_minus_top_by_ssa.array() =
                (particular.d_A_minus_by_ssa.array() *
                     particular.C_minus.transpose().array() *
                     homog.X_plus.transpose().array() +
                 particular.A_minus.transpose().array() *
                     particular.d_C_minus_by_ssa.array() *
                     homog.X_plus.transpose().array() +
                 particular.A_minus.transpose().array() *
                     particular.C_minus.transpose().array() *
                     homog.d_X_plus_by_ssa.array());

            particular.d_G_minus_bottom_by_ssa.array() =
                (particular.d_A_plus_by_ssa.array() *
                     particular.C_plus.transpose().array() *
                     homog.X_minus.transpose().array() +
                 particular.A_plus.transpose().array() *
                     particular.d_C_plus_by_ssa.array() *
                     homog.X_minus.transpose().array() +
                 particular.A_plus.transpose().array() *
                     particular.C_plus.transpose().array() *
                     homog.d_X_minus_by_ssa.array());
        }
    }

    /**
     * Solves the boundary value problem
     *
     * @param input
     * @param solution
     */
    inline void solve_bvp(const Input& input, Solution& solution) {
        // Solve the BVP

        for (int i = 0; i < 2; ++i) {
            auto& homog = solution.homog[i];
            auto& particular = solution.particular[i];
            auto& bvp = solution.bvp_coeffs[i];

            bvp.rhs(0, 0) = -1.0 * particular.G_plus_top(0);

            bvp.rhs(Eigen::seq(1, Eigen::last - 1, 2), 0) =
                particular.G_minus_top(Eigen::seq(1, Eigen::last)) -
                particular.G_minus_bottom(Eigen::seq(0, Eigen::last - 1));
            bvp.rhs(Eigen::seq(2, Eigen::last, 2), 0) =
                particular.G_plus_top(Eigen::seq(1, Eigen::last)) -
                particular.G_plus_bottom(Eigen::seq(0, Eigen::last - 1));

            bvp.rhs(2 * input.nlyr - 1, 0) =
                input.csz * input.albedo / EIGEN_PI *
                    input.transmission(Eigen::last) -
                (particular.G_minus_bottom(Eigen::last) -
                 2 * input.mu * input.albedo *
                     particular.G_minus_bottom(Eigen::last));

            // Now the LHS and gradients
            bvp.d(0) = homog.X_plus(0);
            bvp.a(0) = homog.X_minus(0) * homog.omega(0);

            bvp.d_d_by_ssa(0, 0) = homog.d_X_plus_by_ssa(0);
            bvp.d_a_by_ssa(0, 0) = homog.d_X_minus_by_ssa(0) * homog.omega(0) +
                                   homog.X_minus(0) * homog.d_omega_by_ssa(0);

            for (int j = 0; j < input.nlyr - 1; ++j) {
                int j2 = j * 2;

                bvp.c(j2 + 1) = homog.X_minus(j) * homog.omega(j);
                bvp.d(j2 + 1) = homog.X_plus(j);
                bvp.a(j2 + 1) = -homog.X_minus(j + 1);
                bvp.b(j2 + 1) = -homog.X_plus(j + 1) * homog.omega(j + 1);

                bvp.e(j2 + 2) = homog.X_plus(j) * homog.omega(j);
                bvp.c(j2 + 2) = homog.X_minus(j);
                bvp.d(j2 + 2) = -homog.X_plus(j + 1);
                bvp.a(j2 + 2) = -homog.X_minus(j + 1) * homog.omega(j + 1);

                bvp.d_c_by_ssa(j2 + 1, j) =
                    homog.d_X_minus_by_ssa(j) * homog.omega(j) +
                    homog.X_minus(j) * homog.d_omega_by_ssa(j);
                bvp.d_d_by_ssa(j2 + 1, j) = homog.d_X_plus_by_ssa(j);
                bvp.d_a_by_ssa(j2 + 1, j + 1) = -homog.d_X_minus_by_ssa(j + 1);
                bvp.d_b_by_ssa(j2 + 1, j + 1) =
                    -homog.d_X_plus_by_ssa(j + 1) * homog.omega(j + 1) -
                    homog.X_plus(j + 1) * homog.d_omega_by_ssa(j + 1);

                bvp.d_e_by_ssa(j2 + 2, j) =
                    homog.d_X_plus_by_ssa(j) * homog.omega(j) +
                    homog.X_plus(j) * homog.d_omega_by_ssa(j);
                bvp.d_c_by_ssa(j2 + 2, j) = homog.d_X_minus_by_ssa(j);
                bvp.d_d_by_ssa(j2 + 2, j + 1) = -homog.d_X_plus_by_ssa(j + 1);
                bvp.d_a_by_ssa(j2 + 2, j + 1) =
                    -homog.d_X_minus_by_ssa(j + 1) * homog.omega(j + 1) -
                    homog.X_minus(j + 1) * homog.d_omega_by_ssa(j + 1);
            }

            bvp.c(2 * input.nlyr - 1) =
                (homog.X_minus(Eigen::last) -
                 2 * input.mu * input.albedo * homog.X_plus(Eigen::last)) *
                homog.omega(Eigen::last);
            bvp.d(2 * input.nlyr - 1) =
                (homog.X_plus(Eigen::last) -
                 2 * input.mu * input.albedo * homog.X_minus(Eigen::last));

            bvp.d_d_by_ssa(2 * input.nlyr - 1, input.nlyr - 1) =
                homog.d_X_plus_by_ssa(Eigen::last) -
                2 * input.mu * input.albedo *
                    homog.d_X_minus_by_ssa(Eigen::last);

            bvp.d_c_by_ssa(2 * input.nlyr - 1, input.nlyr - 1) =
                (homog.d_X_minus_by_ssa(Eigen::last) -
                 2 * input.mu * input.albedo *
                     homog.d_X_plus_by_ssa(Eigen::last)) *
                    homog.omega(Eigen::last) +
                (homog.X_minus(Eigen::last) -
                 2 * input.mu * input.albedo * homog.X_plus(Eigen::last)) *
                    homog.d_omega_by_ssa(Eigen::last);

            pentadiagonal_solve(bvp, bvp.rhs);
        }
    }

    /**
     * @brief Calls solve_layers and then solve_bvp
     *
     * @param input
     * @param solution
     */
    inline void solve(const Input& input, Solution& solution) {
        // Assign the homogneous and particular solutions
        solve_layers(input, solution);

        // Solve the BVP
        solve_bvp(input, solution);
    }

    /**
     * @brief Calculates the post processed integrated sources for a given line
     * of sight
     *
     * @param input
     * @param viewing_zenith
     * @param azimuth
     * @param solution
     * @param sources
     */
    inline void post_process(const Input& input, double viewing_zenith,
                             double azimuth, const Solution& solution,
                             Sources& sources) {
        // Lpsums
        sources.lpsum_plus[0].array() =
            0.5 * (input.ssa.array() - input.ssa.array() * input.b1.array() *
                                           viewing_zenith * input.mu);
        sources.lpsum_minus[0].array() =
            0.5 * (input.ssa.array() + input.ssa.array() * input.b1.array() *
                                           viewing_zenith * input.mu);

        sources.lpsum_plus[1].array() =
            0.25 * input.ssa.array() * input.b1.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));
        sources.lpsum_minus[1].array() =
            0.25 * input.ssa.array() * input.b1.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));

        // And their derivatives
        sources.d_lpsum_plus_by_ssa[0].array() =
            0.5 * (1 - input.b1.array() * viewing_zenith * input.mu);
        sources.d_lpsum_minus_by_ssa[0].array() =
            0.5 * (1 + input.b1.array() * viewing_zenith * input.mu);

        sources.d_lpsum_plus_by_ssa[1].array() =
            0.25 * input.b1.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));
        sources.d_lpsum_minus_by_ssa[1].array() =
            0.25 * input.b1.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));

        sources.beamtrans.array() = (-input.od.array() / viewing_zenith).exp();

        sources.E_minus.array() =
            (1 - input.expsec.array() * sources.beamtrans.array()) /
            (1 + input.average_secant.array() * viewing_zenith);

        sources.source.array().setZero();
        sources.d_source_by_ssa.array().setZero();
        for (int i = 0; i < 2; ++i) {
            double azi_factor = cos(i * azimuth);
            const auto& homog = solution.homog[i];
            const auto& part = solution.particular[i];

            sources.Y_plus[i].array() =
                sources.lpsum_plus[i].array() * homog.X_plus.array() +
                sources.lpsum_minus[i].array() * homog.X_minus.array();
            sources.Y_minus[i].array() =
                sources.lpsum_plus[i].array() * homog.X_minus.array() +
                sources.lpsum_minus[i].array() * homog.X_plus.array();

            // and their gradients
            sources.d_Y_plus_by_ssa[i].array() =
                sources.d_lpsum_plus_by_ssa[i].array() *
                    homog.X_plus.transpose().array() +
                sources.d_lpsum_minus_by_ssa[i].array() *
                    homog.X_minus.transpose().array() +
                sources.lpsum_plus[i].transpose().array() *
                    homog.d_X_plus_by_ssa.array() +
                sources.lpsum_minus[i].transpose().array() *
                    homog.d_X_minus_by_ssa.array();

            sources.d_Y_minus_by_ssa[i].array() =
                sources.d_lpsum_plus_by_ssa[i].array() *
                    homog.X_minus.transpose().array() +
                sources.d_lpsum_minus_by_ssa[i].array() *
                    homog.X_plus.transpose().array() +
                sources.lpsum_plus[i].transpose().array() *
                    homog.d_X_minus_by_ssa.array() +
                sources.lpsum_minus[i].transpose().array() *
                    homog.d_X_plus_by_ssa.array();

            // Homogeneous multipliers
            sources.H_minus[i].array() =
                (homog.omega.array() - sources.beamtrans.array()) /
                (1 - homog.k.array() * viewing_zenith);
            sources.H_plus[i].array() =
                (1 - homog.omega.array() * sources.beamtrans.array()) /
                (1 + homog.k.array() * viewing_zenith);

            // And their gradients
            sources.d_H_minus_by_ssa[i].array() =
                homog.d_omega_by_ssa.array() /
                    (1 - homog.k.transpose().array() * viewing_zenith) +
                homog.d_k_by_ssa.array() * viewing_zenith *
                    sources.H_minus[i].transpose().array() /
                    (1 - homog.k.transpose().array() * viewing_zenith);

            sources.d_H_plus_by_ssa[i].array() =
                -homog.d_omega_by_ssa.array() *
                    sources.beamtrans.transpose().array() /
                    (1 + homog.k.transpose().array() * viewing_zenith) -
                homog.d_k_by_ssa.array() * viewing_zenith *
                    sources.H_plus[i].transpose().array() /
                    (1 + homog.k.transpose().array() * viewing_zenith);

            // Particular multipliers
            sources.D_plus[i].array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (sources.E_minus.array() -
                 input.expsec.array() * sources.H_minus[i].array()) /
                (input.average_secant.array() + homog.k.array());
            sources.D_minus[i].array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (sources.H_plus[i].array() - sources.E_minus.array()) /
                (input.average_secant.array() - homog.k.array());

            // And their gradients
            sources.d_D_plus_by_ssa[i].array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1))
                         .transpose()
                         .array() *
                     (-sources.d_H_minus_by_ssa[i].array() *
                      input.expsec.transpose().array()) -
                 homog.d_k_by_ssa.array() *
                     sources.D_plus[i].transpose().array()) /
                (input.average_secant.transpose().array() +
                 homog.k.transpose().array());

            sources.d_D_minus_by_ssa[i].array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1))
                         .transpose()
                         .array() *
                     sources.d_H_plus_by_ssa[i].array() +
                 homog.d_k_by_ssa.array() *
                     sources.D_minus[i].transpose().array()) /
                (input.average_secant.transpose().array() -
                 homog.k.transpose().array());

            // The particular source itself
            sources.V[i].array() =
                part.A_plus.array() * sources.Y_plus[i].array() *
                    sources.D_minus[i].array() +
                part.A_minus.array() * sources.Y_minus[i].array() *
                    sources.D_plus[i].array();

            // with it's derivative
            sources.d_V_by_ssa[i].array() =
                (part.d_A_plus_by_ssa.array() *
                     sources.Y_plus[i].transpose().array() *
                     sources.D_minus[i].transpose().array() +
                 part.A_plus.transpose().array() *
                     sources.d_Y_plus_by_ssa[i].array() *
                     sources.D_minus[i].transpose().array() +
                 part.A_plus.transpose().array() *
                     sources.Y_plus[i].transpose().array() *
                     sources.d_D_minus_by_ssa[i].array()) +
                (part.d_A_minus_by_ssa.array() *
                     sources.Y_minus[i].transpose().array() *
                     sources.D_plus[i].transpose().array() +
                 part.A_minus.transpose().array() *
                     sources.d_Y_minus_by_ssa[i].array() *
                     sources.D_plus[i].transpose().array() +
                 part.A_minus.transpose().array() *
                     sources.Y_minus[i].transpose().array() *
                     sources.d_D_plus_by_ssa[i].array());

            // Add in homogeneous sources
            sources.source.array() +=
                azi_factor *
                (solution.bvp_coeffs[i]
                     .rhs(Eigen::seq(0, Eigen::last, 2), 0)
                     .array() *
                 sources.Y_plus[i].array() * sources.H_plus[i].array());

            // Append the derivative
            sources.d_source_by_ssa.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(0, Eigen::last, 2), 0)
                                  .transpose()
                                  .array() *
                              (sources.d_Y_plus_by_ssa[i].array() *
                                   sources.H_plus[i].transpose().array() +
                               sources.Y_plus[i].transpose().array() *
                                   sources.d_H_plus_by_ssa[i].array()));

            sources.source.array() +=
                azi_factor *
                (solution.bvp_coeffs[i]
                     .rhs(Eigen::seq(1, Eigen::last, 2), 0)
                     .array() *
                 sources.Y_minus[i].array() * sources.H_minus[i].array());

            sources.d_source_by_ssa.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(1, Eigen::last, 2), 0)
                                  .transpose()
                                  .array() *
                              (sources.d_Y_minus_by_ssa[i].array() *
                                   sources.H_minus[i].transpose().array() +
                               sources.Y_minus[i].transpose().array() *
                                   sources.d_H_minus_by_ssa[i].array()));

            // And keep track of the gradients for backprop
            sources.d_bvp_coeff[i](Eigen::seq(0, 2 * input.nlyr - 1, 2))
                .array() = azi_factor * (sources.Y_plus[i].transpose().array() *
                                         sources.H_plus[i].transpose().array());
            sources.d_bvp_coeff[i](Eigen::seq(1, 2 * input.nlyr - 1, 2))
                .array() =
                azi_factor * (sources.Y_minus[i].transpose().array() *
                              sources.H_minus[i].transpose().array());

            sources.source.array() += azi_factor * (sources.V[i].array());
            sources.d_source_by_ssa.array() +=
                azi_factor * (sources.d_V_by_ssa[i].array());
        }
    }
} // namespace sasktran2::twostream
