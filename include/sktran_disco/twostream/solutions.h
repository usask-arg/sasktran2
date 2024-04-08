#pragma once

#include <sasktran2/internal_common.h>
#include "storage.h"

namespace sasktran2::twostream {

    /**
     * Solves the pentadiagonal system
     *
     * @param bvp
     * @return int
     */
    inline int pentadiagonal_solve(BVPCoeffs& bvp) {

        int N = bvp.a.size();

        auto& a = bvp.a;
        auto& b = bvp.b;
        auto& e = bvp.e;
        auto& c = bvp.c;
        auto& d = bvp.d;

        auto& y = bvp.rhs;
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

        for (auto& homog : solution.homog) {
            homog.k.array() = (homog.s.array() * homog.d.array()).sqrt();

            homog.X_plus.array() =
                0.5 * (1 - homog.s.array() / homog.k.array());
            homog.X_minus.array() =
                0.5 * (1 + homog.s.array() / homog.k.array());

            homog.omega.array() = (-homog.k.array() * input.od.array()).exp();
        }

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

        for (int i = 0; i < 2; ++i) {
            auto& homog = solution.homog[i];
            auto& particular = solution.particular[i];

            particular.norm.array() =
                input.mu * (homog.X_plus.array() * homog.X_plus.array() -
                            homog.X_minus.array() * homog.X_minus.array());

            particular.A_plus.array() =
                (particular.Q_plus.array() * homog.X_plus.array() +
                 particular.Q_minus.array() * homog.X_minus.array()) /
                particular.norm.array();
            particular.A_minus.array() =
                (particular.Q_minus.array() * homog.X_plus.array() +
                 particular.Q_plus.array() * homog.X_minus.array()) /
                particular.norm.array();

            particular.C_plus.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (homog.omega.array() - input.expsec.array()) /
                (input.average_secant.array() - homog.k.array());
            particular.C_minus.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (1 - homog.omega.array() * input.expsec.array()) /
                (input.average_secant.array() + homog.k.array());

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

            // Now the LHS

            bvp.d(0) = homog.X_plus(0);
            bvp.a(0) = homog.X_minus(0) * homog.omega(0);

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
            }

            bvp.c(2 * input.nlyr - 1) =
                (homog.X_minus(Eigen::last) -
                 2 * input.mu * input.albedo * homog.X_plus(Eigen::last)) *
                homog.omega(Eigen::last);
            bvp.d(2 * input.nlyr - 1) =
                (homog.X_plus(Eigen::last) -
                 2 * input.mu * input.albedo * homog.X_minus(Eigen::last));

            pentadiagonal_solve(bvp);
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

        sources.beamtrans.array() = (-input.od.array() / viewing_zenith).exp();

        sources.E_minus.array() =
            (1 - input.expsec.array() * sources.beamtrans.array()) /
            (1 + input.average_secant.array() * viewing_zenith);

        sources.source.array().setZero();
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

            sources.H_minus[i].array() =
                (homog.omega.array() - sources.beamtrans.array()) /
                (1 - homog.k.array() * viewing_zenith);
            sources.H_plus[i].array() =
                (1 - homog.omega.array() * sources.beamtrans.array()) /
                (1 + homog.k.array() * viewing_zenith);

            sources.D_plus[i].array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (sources.E_minus.array() -
                 input.expsec.array() * sources.H_minus[i].array()) /
                (input.average_secant.array() + homog.k.array());
            sources.D_minus[i].array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (sources.H_plus[i].array() - sources.E_minus.array()) /
                (input.average_secant.array() - homog.k.array());

            sources.V[i].array() =
                part.A_plus.array() * sources.Y_plus[i].array() *
                    sources.D_minus[i].array() +
                part.A_minus.array() * sources.Y_minus[i].array() *
                    sources.D_plus[i].array();

            sources.source.array() +=
                azi_factor *
                (solution.bvp_coeffs[i]
                     .rhs(Eigen::seq(0, Eigen::last, 2), 0)
                     .array() *
                 sources.Y_plus[i].array() * sources.H_plus[i].array());
            sources.source.array() +=
                azi_factor *
                (solution.bvp_coeffs[i]
                     .rhs(Eigen::seq(1, Eigen::last, 2), 0)
                     .array() *
                 sources.Y_minus[i].array() * sources.H_minus[i].array());
            sources.source.array() += azi_factor * (sources.V[i].array());
        }
    }
} // namespace sasktran2::twostream
