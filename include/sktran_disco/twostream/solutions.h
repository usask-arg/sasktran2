#pragma once

#include <cstddef>
#include <sasktran2/internal_common.h>
#include "sktran_disco/sktran_do_types.h"
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
        solution.homog[0].d.value.array() =
            (input.ssa.array() * input.b1.array() * input.mu - 1 / input.mu);
        solution.homog[0].s.value.array() =
            1 / input.mu * (input.ssa.array() - 1);

        solution.homog[1].d.value.array() = -1 / input.mu;
        solution.homog[1].s.value.array() =
            1.0 / (2.0 * input.mu) *
            (input.ssa.array() * input.b1.array() * (1 - input.mu * input.mu) -
             2);

        // Gradients of d by ssa
        solution.homog[0].d.d_ssa.array() = input.b1.array() * input.mu;

        solution.homog[1].d.d_ssa.array() = 0;

        // Gradients of s by ssa
        solution.homog[0].s.d_ssa.array() = 1 / input.mu;
        solution.homog[1].s.d_ssa.array() = 1.0 / (2.0 * input.mu) *
                                            input.b1.array() *
                                            (1 - input.mu * input.mu);

        // Gradients of d/s by b1
        solution.homog[0].d.d_b1.array() = input.mu * input.ssa.array();
        // solution.homog[1].d.d_b1.array().setZero(); // Should always be 0

        // solution.homog[0].s.d_b1.array().setZero(); // Should already get set
        // to 0 by default
        solution.homog[1].s.d_b1.array() = 1.0 / (2.0 * input.mu) *
                                           input.ssa.array() *
                                           (1 - input.mu * input.mu);

        for (auto& homog : solution.homog) {
            // Eigenvalue and gradient
            homog.k.value.array() =
                (homog.s.value.array() * homog.d.value.array()).sqrt();
            homog.k.d_ssa.array() =
                (0.5 / homog.k.value.array()) *
                (homog.s.value.array() * homog.d.d_ssa.array() +
                 homog.d.value.array() * homog.s.d_ssa.array());
            homog.k.d_b1.array() =
                (0.5 / homog.k.value.array()) *
                (homog.s.value.array() * homog.d.d_b1.array() +
                 homog.d.value.array() * homog.s.d_b1.array());

            // Eigenvectors, and gradients
            homog.X_plus.value.array() =
                0.5 * (1 - homog.s.value.array() / homog.k.value.array());
            homog.X_minus.value.array() =
                0.5 * (1 + homog.s.value.array() / homog.k.value.array());

            homog.X_plus.d_ssa.array() =
                (-0.5 / homog.k.value.array() * homog.s.d_ssa.array()) +
                0.5 / homog.d.value.array() * homog.k.d_ssa.array();
            homog.X_minus.d_ssa.array() =
                (0.5 / homog.k.value.array() * homog.s.d_ssa.array()) -
                0.5 / homog.d.value.array() * homog.k.d_ssa.array();

            homog.X_plus.d_b1.array() =
                (-0.5 / homog.k.value.array() * homog.s.d_b1.array()) +
                0.5 / homog.d.value.array() * homog.k.d_b1.array();
            homog.X_minus.d_b1.array() =
                (0.5 / homog.k.value.array() * homog.s.d_b1.array()) -
                0.5 / homog.d.value.array() * homog.k.d_b1.array();

            // Convenience quantity for later
            homog.omega.value.array() =
                (-homog.k.value.array() * input.od.array()).exp();
            homog.omega.d_ssa.array() = -homog.omega.value.array() *
                                        input.od.array() *
                                        homog.k.d_ssa.array();

            homog.omega.d_b1.array() = -homog.omega.value.array() *
                                       input.od.array() * homog.k.d_b1.array();

            homog.omega.d_od.array() =
                -homog.omega.value.array() * homog.k.value.array();
        }

        // Q_plus and Q_minus
        solution.particular[0].Q_plus.value.array() =
            1.0 / (4 * EIGEN_PI) *
            (input.ssa.array() +
             input.b1.array() * input.ssa.array() * input.csz * input.mu);
        solution.particular[0].Q_minus.value.array() =
            1.0 / (4 * EIGEN_PI) *
            (input.ssa.array() -
             input.b1.array() * input.ssa.array() * input.csz * input.mu);

        solution.particular[1].Q_plus.value.array() =
            1.0 / (4 * EIGEN_PI) * input.ssa.array() * input.b1.array() *
            sqrt(1 - input.mu * input.mu) * sqrt(1 - input.csz * input.csz);
        solution.particular[1].Q_minus.value.array() =
            solution.particular[1].Q_plus.value.array();

        // And their gradients
        solution.particular[0].Q_plus.d_ssa.array() =
            1.0 / (4.0 * EIGEN_PI) *
            (1 + input.b1.array() * input.csz * input.mu);
        solution.particular[0].Q_minus.d_ssa.array() =
            1.0 / (4.0 * EIGEN_PI) *
            (1 - input.b1.array() * input.csz * input.mu);

        solution.particular[1].Q_plus.d_ssa.array() =
            1.0 / (4.0 * EIGEN_PI) * input.b1.array() *
            sqrt(1 - input.mu * input.mu) * sqrt(1 - input.csz * input.csz);
        solution.particular[1].Q_minus.d_ssa.array() =
            solution.particular[1].Q_plus.d_ssa.array();

        solution.particular[0].Q_plus.d_b1.array() =
            1.0 / (4.0 * EIGEN_PI) * input.ssa.array() * input.csz * input.mu;
        solution.particular[0].Q_minus.d_b1.array() =
            -1.0 / (4.0 * EIGEN_PI) * input.ssa.array() * input.csz * input.mu;

        solution.particular[1].Q_plus.d_b1.array() =
            1.0 / (4.0 * EIGEN_PI) * input.ssa.array() *
            sqrt(1 - input.mu * input.mu) * sqrt(1 - input.csz * input.csz);
        solution.particular[1].Q_minus.d_b1.array() =
            solution.particular[1].Q_plus.d_b1.array();

        for (int i = 0; i < 2; ++i) {
            auto& homog = solution.homog[i];
            auto& particular = solution.particular[i];

            particular.norm.value.array() =
                input.mu *
                (homog.X_plus.value.array() * homog.X_plus.value.array() -
                 homog.X_minus.value.array() * homog.X_minus.value.array());

            particular.norm.d_ssa.array() =
                2 * input.mu *
                (homog.X_plus.value.array() * homog.X_plus.d_ssa.array() -
                 homog.X_minus.value.array() * homog.X_minus.d_ssa.array());

            particular.norm.d_b1.array() =
                2 * input.mu *
                (homog.X_plus.value.array() * homog.X_plus.d_b1.array() -
                 homog.X_minus.value.array() * homog.X_minus.d_b1.array());

            // A_plus and A_minus and their gradients by ssa
            particular.A_plus.value.array() =
                (particular.Q_plus.value.array() * homog.X_plus.value.array() +
                 particular.Q_minus.value.array() *
                     homog.X_minus.value.array()) /
                particular.norm.value.array();
            particular.A_minus.value.array() =
                (particular.Q_minus.value.array() * homog.X_plus.value.array() +
                 particular.Q_plus.value.array() *
                     homog.X_minus.value.array()) /
                particular.norm.value.array();

            particular.A_plus.d_ssa.array() =
                (particular.Q_plus.d_ssa.array() * homog.X_plus.value.array() +
                 particular.Q_minus.d_ssa.array() *
                     homog.X_minus.value.array() +
                 particular.Q_plus.value.array() * homog.X_plus.d_ssa.array() +
                 particular.Q_minus.value.array() *
                     homog.X_minus.d_ssa.array() -
                 particular.norm.d_ssa.array() *
                     particular.A_plus.value.array()) /
                particular.norm.value.array();

            particular.A_minus.d_ssa.array() =
                (particular.Q_minus.d_ssa.array() * homog.X_plus.value.array() +
                 particular.Q_plus.d_ssa.array() * homog.X_minus.value.array() +
                 particular.Q_minus.value.array() * homog.X_plus.d_ssa.array() +
                 particular.Q_plus.value.array() * homog.X_minus.d_ssa.array() -
                 particular.norm.d_ssa.array() *
                     particular.A_minus.value.array()) /
                particular.norm.value.array();

            particular.A_plus.d_b1.array() =
                (particular.Q_plus.d_b1.array() * homog.X_plus.value.array() +
                 particular.Q_minus.d_b1.array() * homog.X_minus.value.array() +
                 particular.Q_plus.value.array() * homog.X_plus.d_b1.array() +
                 particular.Q_minus.value.array() * homog.X_minus.d_b1.array() -
                 particular.norm.d_b1.array() *
                     particular.A_plus.value.array()) /
                particular.norm.value.array();

            particular.A_minus.d_b1.array() =
                (particular.Q_minus.d_b1.array() * homog.X_plus.value.array() +
                 particular.Q_plus.d_b1.array() * homog.X_minus.value.array() +
                 particular.Q_minus.value.array() * homog.X_plus.d_b1.array() +
                 particular.Q_plus.value.array() * homog.X_minus.d_b1.array() -
                 particular.norm.d_b1.array() *
                     particular.A_minus.value.array()) /
                particular.norm.value.array();

            // Solution multipliers and their derivatives
            particular.C_plus.value.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (homog.omega.value.array() - input.expsec.array()) /
                (input.average_secant.array() - homog.k.value.array());
            particular.C_minus.value.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (1 - homog.omega.value.array() * input.expsec.array()) /
                (input.average_secant.array() + homog.k.value.array());

            particular.C_plus.d_ssa.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                    (homog.omega.d_ssa.array()) /
                    (input.average_secant.array() - homog.k.value.array()) +
                homog.k.d_ssa.array() * particular.C_plus.value.array() /
                    (input.average_secant.array() - homog.k.value.array());

            particular.C_minus.d_ssa.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                    (-homog.omega.d_ssa.array() * input.expsec.array()) /
                    (input.average_secant.array() + homog.k.value.array()) -
                homog.k.d_ssa.array() * particular.C_minus.value.array() /
                    (input.average_secant.array() + homog.k.value.array());

            particular.C_plus.d_b1.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                    (homog.omega.d_b1.array()) /
                    (input.average_secant.array() - homog.k.value.array()) +
                homog.k.d_b1.array() * particular.C_plus.value.array() /
                    (input.average_secant.array() - homog.k.value.array());

            particular.C_minus.d_b1.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                    (-homog.omega.d_b1.array() * input.expsec.array()) /
                    (input.average_secant.array() + homog.k.value.array()) -
                homog.k.d_b1.array() * particular.C_minus.value.array() /
                    (input.average_secant.array() + homog.k.value.array());

            // OD Derivatives through omega and expsec.  expsec.d_od = -secant *
            // expsec
            particular.C_plus.d_od.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (homog.omega.d_od.array() +
                 input.average_secant.array() * input.expsec.array()) /
                (input.average_secant.array() - homog.k.value.array());

            particular.C_minus.d_od.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (-homog.omega.d_od.array() * input.expsec.array() +
                 input.average_secant.array() * homog.omega.value.array() *
                     input.expsec.array()) /
                (input.average_secant.array() + homog.k.value.array());

            // Transmission derivatives are straightfoward
            particular.C_plus.d_transmission(Eigen::seq(0, Eigen::last - 1))
                .array() =
                (homog.omega.value.array() - input.expsec.array()) /
                (input.average_secant.array() - homog.k.value.array());

            particular.C_minus.d_transmission(Eigen::seq(0, Eigen::last - 1))
                .array() =
                (1 - homog.omega.value.array() * input.expsec.array()) /
                (input.average_secant.array() + homog.k.value.array());

            // Secant derivatives are a little more complicated
            // expsec.d_sec = -od * expsec
            particular.C_plus.d_secant.array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                     input.od.array() * input.expsec.array() -
                 particular.C_plus.value.array()) /
                (input.average_secant.array() - homog.k.value.array());

            particular.C_minus.d_secant.array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                     input.od.array() * input.expsec.array() *
                     homog.omega.value.array() -
                 particular.C_minus.value.array()) /
                (input.average_secant.array() + homog.k.value.array());

            // Particular solutions and their derivatives
            particular.G_plus_top.value.array() =
                particular.A_minus.value.array() *
                particular.C_minus.value.array() * homog.X_minus.value.array();
            particular.G_plus_bottom.value.array() =
                particular.A_plus.value.array() *
                particular.C_plus.value.array() * homog.X_plus.value.array();
            particular.G_minus_top.value.array() =
                particular.A_minus.value.array() *
                particular.C_minus.value.array() * homog.X_plus.value.array();
            particular.G_minus_bottom.value.array() =
                particular.A_plus.value.array() *
                particular.C_plus.value.array() * homog.X_minus.value.array();

            particular.G_plus_top.d_ssa.array() =
                (particular.A_minus.d_ssa.array() *
                     particular.C_minus.value.array() *
                     homog.X_minus.value.array() +
                 particular.A_minus.value.array() *
                     particular.C_minus.d_ssa.array() *
                     homog.X_minus.value.array() +
                 particular.A_minus.value.array() *
                     particular.C_minus.value.array() *
                     homog.X_minus.d_ssa.array());

            particular.G_plus_bottom.d_ssa.array() =
                (particular.A_plus.d_ssa.array() *
                     particular.C_plus.value.array() *
                     homog.X_plus.value.array() +
                 particular.A_plus.value.array() *
                     particular.C_plus.d_ssa.array() *
                     homog.X_plus.value.array() +
                 particular.A_plus.value.array() *
                     particular.C_plus.value.array() *
                     homog.X_plus.d_ssa.array());

            particular.G_minus_top.d_ssa.array() =
                (particular.A_minus.d_ssa.array() *
                     particular.C_minus.value.array() *
                     homog.X_plus.value.array() +
                 particular.A_minus.value.array() *
                     particular.C_minus.d_ssa.array() *
                     homog.X_plus.value.array() +
                 particular.A_minus.value.array() *
                     particular.C_minus.value.array() *
                     homog.X_plus.d_ssa.array());

            particular.G_minus_bottom.d_ssa.array() =
                (particular.A_plus.d_ssa.array() *
                     particular.C_plus.value.array() *
                     homog.X_minus.value.array() +
                 particular.A_plus.value.array() *
                     particular.C_plus.d_ssa.array() *
                     homog.X_minus.value.array() +
                 particular.A_plus.value.array() *
                     particular.C_plus.value.array() *
                     homog.X_minus.d_ssa.array());

            // b1 derivatives
            particular.G_plus_top.d_b1.array() =
                (particular.A_minus.d_b1.array() *
                     particular.C_minus.value.array() *
                     homog.X_minus.value.array() +
                 particular.A_minus.value.array() *
                     particular.C_minus.d_b1.array() *
                     homog.X_minus.value.array() +
                 particular.A_minus.value.array() *
                     particular.C_minus.value.array() *
                     homog.X_minus.d_b1.array());

            particular.G_plus_bottom.d_b1.array() =
                (particular.A_plus.d_b1.array() *
                     particular.C_plus.value.array() *
                     homog.X_plus.value.array() +
                 particular.A_plus.value.array() *
                     particular.C_plus.d_b1.array() *
                     homog.X_plus.value.array() +
                 particular.A_plus.value.array() *
                     particular.C_plus.value.array() *
                     homog.X_plus.d_b1.array());

            particular.G_minus_top.d_b1.array() =
                (particular.A_minus.d_b1.array() *
                     particular.C_minus.value.array() *
                     homog.X_plus.value.array() +
                 particular.A_minus.value.array() *
                     particular.C_minus.d_b1.array() *
                     homog.X_plus.value.array() +
                 particular.A_minus.value.array() *
                     particular.C_minus.value.array() *
                     homog.X_plus.d_b1.array());

            particular.G_minus_bottom.d_b1.array() =
                (particular.A_plus.d_b1.array() *
                     particular.C_plus.value.array() *
                     homog.X_minus.value.array() +
                 particular.A_plus.value.array() *
                     particular.C_plus.d_b1.array() *
                     homog.X_minus.value.array() +
                 particular.A_plus.value.array() *
                     particular.C_plus.value.array() *
                     homog.X_minus.d_b1.array());

            // G's only take od, transmission, and secant derivatives from C's
            particular.G_plus_top.d_od.array() =
                particular.A_minus.value.array() *
                particular.C_minus.d_od.array() * homog.X_minus.value.array();

            particular.G_plus_bottom.d_od.array() =
                particular.A_plus.value.array() *
                particular.C_plus.d_od.array() * homog.X_plus.value.array();

            particular.G_minus_top.d_od.array() =
                particular.A_minus.value.array() *
                particular.C_minus.d_od.array() * homog.X_plus.value.array();

            particular.G_minus_bottom.d_od.array() =
                particular.A_plus.value.array() *
                particular.C_plus.d_od.array() * homog.X_minus.value.array();

            particular.G_plus_top.d_transmission(Eigen::seq(0, Eigen::last - 1))
                .array() = particular.A_minus.value.array() *
                           particular.C_minus
                               .d_transmission(Eigen::seq(0, Eigen::last - 1))
                               .array() *
                           homog.X_minus.value.array();

            particular.G_plus_bottom
                .d_transmission(Eigen::seq(0, Eigen::last - 1))
                .array() =
                particular.A_plus.value.array() *
                particular.C_plus.d_transmission(Eigen::seq(0, Eigen::last - 1))
                    .array() *
                homog.X_plus.value.array();

            particular.G_minus_top
                .d_transmission(Eigen::seq(0, Eigen::last - 1))
                .array() = particular.A_minus.value.array() *
                           particular.C_minus
                               .d_transmission(Eigen::seq(0, Eigen::last - 1))
                               .array() *
                           homog.X_plus.value.array();

            particular.G_minus_bottom
                .d_transmission(Eigen::seq(0, Eigen::last - 1))
                .array() =
                particular.A_plus.value.array() *
                particular.C_plus.d_transmission(Eigen::seq(0, Eigen::last - 1))
                    .array() *
                homog.X_minus.value.array();

            particular.G_plus_top.d_secant.array() =
                particular.A_minus.value.array() *
                particular.C_minus.d_secant.array() *
                homog.X_minus.value.array();

            particular.G_plus_bottom.d_secant.array() =
                particular.A_plus.value.array() *
                particular.C_plus.d_secant.array() * homog.X_plus.value.array();

            particular.G_minus_top.d_secant.array() =
                particular.A_minus.value.array() *
                particular.C_minus.d_secant.array() *
                homog.X_plus.value.array();

            particular.G_minus_bottom.d_secant.array() =
                particular.A_plus.value.array() *
                particular.C_plus.d_secant.array() *
                homog.X_minus.value.array();
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

            bvp.rhs(0, 0) = -1.0 * particular.G_plus_top.value(0);

            bvp.rhs(Eigen::seq(1, Eigen::last - 1, 2), 0) =
                particular.G_minus_top.value(Eigen::seq(1, Eigen::last)) -
                particular.G_minus_bottom.value(Eigen::seq(0, Eigen::last - 1));
            bvp.rhs(Eigen::seq(2, Eigen::last, 2), 0) =
                particular.G_plus_top.value(Eigen::seq(1, Eigen::last)) -
                particular.G_plus_bottom.value(Eigen::seq(0, Eigen::last - 1));

            bvp.rhs(2 * input.nlyr - 1, 0) =
                sasktran_disco::kronDelta(i, 0) * input.csz * input.albedo /
                    EIGEN_PI * input.transmission(Eigen::last) -
                (particular.G_minus_bottom.value(Eigen::last) -
                 2 * sasktran_disco::kronDelta(i, 0) * input.mu * input.albedo *
                     particular.G_plus_bottom.value(Eigen::last));

            // Now the LHS and gradients
            bvp.d(0) = homog.X_plus.value(0);
            bvp.a(0) = homog.X_minus.value(0) * homog.omega.value(0);

            bvp.d_d_by_ssa(0, 0) = homog.X_plus.d_ssa(0);
            bvp.d_a_by_ssa(0, 0) =
                homog.X_minus.d_ssa(0) * homog.omega.value(0) +
                homog.X_minus.value(0) * homog.omega.d_ssa(0);

            bvp.d_d_by_b1(0, 0) = homog.X_plus.d_b1(0);
            bvp.d_a_by_b1(0, 0) = homog.X_minus.d_b1(0) * homog.omega.value(0) +
                                  homog.X_minus.value(0) * homog.omega.d_b1(0);
            ;

            bvp.d_a_by_od(0, 0) = homog.X_minus.value(0) * homog.omega.d_od(0);

            for (int j = 0; j < input.nlyr - 1; ++j) {
                int j2 = j * 2;

                bvp.c(j2 + 1) = homog.X_minus.value(j) * homog.omega.value(j);
                bvp.d(j2 + 1) = homog.X_plus.value(j);
                bvp.a(j2 + 1) = -homog.X_minus.value(j + 1);
                bvp.b(j2 + 1) =
                    -homog.X_plus.value(j + 1) * homog.omega.value(j + 1);

                bvp.e(j2 + 2) = homog.X_plus.value(j) * homog.omega.value(j);
                bvp.c(j2 + 2) = homog.X_minus.value(j);
                bvp.d(j2 + 2) = -homog.X_plus.value(j + 1);
                bvp.a(j2 + 2) =
                    -homog.X_minus.value(j + 1) * homog.omega.value(j + 1);

                bvp.d_c_by_ssa(j2 + 1, j) =
                    homog.X_minus.d_ssa(j) * homog.omega.value(j) +
                    homog.X_minus.value(j) * homog.omega.d_ssa(j);
                bvp.d_d_by_ssa(j2 + 1, j) = homog.X_plus.d_ssa(j);
                bvp.d_a_by_ssa(j2 + 1, j + 1) = -homog.X_minus.d_ssa(j + 1);
                bvp.d_b_by_ssa(j2 + 1, j + 1) =
                    -homog.X_plus.d_ssa(j + 1) * homog.omega.value(j + 1) -
                    homog.X_plus.value(j + 1) * homog.omega.d_ssa(j + 1);

                bvp.d_e_by_ssa(j2 + 2, j) =
                    homog.X_plus.d_ssa(j) * homog.omega.value(j) +
                    homog.X_plus.value(j) * homog.omega.d_ssa(j);
                bvp.d_c_by_ssa(j2 + 2, j) = homog.X_minus.d_ssa(j);
                bvp.d_d_by_ssa(j2 + 2, j + 1) = -homog.X_plus.d_ssa(j + 1);
                bvp.d_a_by_ssa(j2 + 2, j + 1) =
                    -homog.X_minus.d_ssa(j + 1) * homog.omega.value(j + 1) -
                    homog.X_minus.value(j + 1) * homog.omega.d_ssa(j + 1);

                bvp.d_c_by_b1(j2 + 1, j) =
                    homog.X_minus.d_b1(j) * homog.omega.value(j) +
                    homog.X_minus.value(j) * homog.omega.d_b1(j);
                bvp.d_d_by_b1(j2 + 1, j) = homog.X_plus.d_b1(j);
                bvp.d_a_by_b1(j2 + 1, j + 1) = -homog.X_minus.d_b1(j + 1);
                bvp.d_b_by_b1(j2 + 1, j + 1) =
                    -homog.X_plus.d_b1(j + 1) * homog.omega.value(j + 1) -
                    homog.X_plus.value(j + 1) * homog.omega.d_b1(j + 1);

                bvp.d_e_by_b1(j2 + 2, j) =
                    homog.X_plus.d_b1(j) * homog.omega.value(j) +
                    homog.X_plus.value(j) * homog.omega.d_b1(j);
                bvp.d_c_by_b1(j2 + 2, j) = homog.X_minus.d_b1(j);
                bvp.d_d_by_b1(j2 + 2, j + 1) = -homog.X_plus.d_b1(j + 1);
                bvp.d_a_by_b1(j2 + 2, j + 1) =
                    -homog.X_minus.d_b1(j + 1) * homog.omega.value(j + 1) -
                    homog.X_minus.value(j + 1) * homog.omega.d_b1(j + 1);

                bvp.d_c_by_od(j2 + 1, j) =
                    homog.X_minus.value(j) * homog.omega.d_od(j);
                bvp.d_b_by_od(j2 + 1, j + 1) =
                    -homog.X_plus.value(j + 1) * homog.omega.d_od(j + 1);

                bvp.d_e_by_od(j2 + 2, j) =
                    homog.X_plus.value(j) * homog.omega.d_od(j);
                bvp.d_a_by_od(j2 + 2, j + 1) =
                    -homog.X_minus.value(j + 1) * homog.omega.d_od(j + 1);
            }

            bvp.c(2 * input.nlyr - 1) =
                (homog.X_minus.value(Eigen::last) -
                 2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                     homog.X_plus.value(Eigen::last)) *
                homog.omega.value(Eigen::last);
            bvp.d(2 * input.nlyr - 1) =
                (homog.X_plus.value(Eigen::last) -
                 2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                     homog.X_minus.value(Eigen::last));

            bvp.d_d_by_ssa(2 * input.nlyr - 1, input.nlyr - 1) =
                homog.X_plus.d_ssa(Eigen::last) -
                2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                    homog.X_minus.d_ssa(Eigen::last);

            bvp.d_c_by_ssa(2 * input.nlyr - 1, input.nlyr - 1) =
                (homog.X_minus.d_ssa(Eigen::last) -
                 2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                     homog.X_plus.d_ssa(Eigen::last)) *
                    homog.omega.value(Eigen::last) +
                (homog.X_minus.value(Eigen::last) -
                 2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                     homog.X_plus.value(Eigen::last)) *
                    homog.omega.d_ssa(Eigen::last);

            bvp.d_d_by_b1(2 * input.nlyr - 1, input.nlyr - 1) =
                homog.X_plus.d_b1(Eigen::last) -
                2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                    homog.X_minus.d_b1(Eigen::last);

            bvp.d_c_by_b1(2 * input.nlyr - 1, input.nlyr - 1) =
                (homog.X_minus.d_b1(Eigen::last) -
                 2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                     homog.X_plus.d_b1(Eigen::last)) *
                    homog.omega.value(Eigen::last) +
                (homog.X_minus.value(Eigen::last) -
                 2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                     homog.X_plus.value(Eigen::last)) *
                    homog.omega.d_b1(Eigen::last);

            bvp.d_c_by_od(2 * input.nlyr - 1, input.nlyr - 1) =
                (homog.X_minus.value(Eigen::last) -
                 2 * input.mu * input.albedo * sasktran_disco::kronDelta(i, 0) *
                     homog.X_plus.value(Eigen::last)) *
                homog.omega.d_od(Eigen::last);

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
        sources.lpsum_plus[0].value.array() =
            0.5 * (input.ssa.array() - input.ssa.array() * input.b1.array() *
                                           viewing_zenith * input.mu);
        sources.lpsum_minus[0].value.array() =
            0.5 * (input.ssa.array() + input.ssa.array() * input.b1.array() *
                                           viewing_zenith * input.mu);

        sources.lpsum_plus[1].value.array() =
            0.25 * input.ssa.array() * input.b1.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));
        sources.lpsum_minus[1].value.array() =
            0.25 * input.ssa.array() * input.b1.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));

        // And their derivatives
        sources.lpsum_plus[0].d_ssa.array() =
            0.5 * (1 - input.b1.array() * viewing_zenith * input.mu);
        sources.lpsum_minus[0].d_ssa.array() =
            0.5 * (1 + input.b1.array() * viewing_zenith * input.mu);

        sources.lpsum_plus[1].d_ssa.array() =
            0.25 * input.b1.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));
        sources.lpsum_minus[1].d_ssa.array() =
            0.25 * input.b1.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));

        sources.lpsum_plus[0].d_b1.array() =
            -0.5 * input.ssa.array() * viewing_zenith * input.mu;

        sources.lpsum_minus[0].d_b1.array() =
            0.5 * input.ssa.array() * viewing_zenith * input.mu;

        sources.lpsum_plus[1].d_b1.array() =
            0.25 * input.ssa.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));

        sources.lpsum_minus[1].d_b1.array() =
            0.25 * input.ssa.array() *
            sqrt((1 - viewing_zenith * viewing_zenith) *
                 (1 - input.mu * input.mu));

        sources.beamtrans.value.array() =
            (-input.od.array() / viewing_zenith).exp();

        sources.beamtrans.d_od.array() =
            -sources.beamtrans.value.array() / viewing_zenith;

        sources.E_minus.value.array() =
            (1 - input.expsec.array() * sources.beamtrans.value.array()) /
            (1 + input.average_secant.array() * viewing_zenith);

        sources.E_minus.d_od.array() =
            (input.average_secant.array() * input.expsec.array() *
                 sources.beamtrans.value.array() -
             input.expsec.array() * sources.beamtrans.d_od.array()) /
            (1 + input.average_secant.array() * viewing_zenith);

        sources.E_minus.d_secant.array() =
            (input.od.array() * input.expsec.array() *
                 sources.beamtrans.value.array() -
             viewing_zenith * sources.E_minus.value.array()) /
            (1 + input.average_secant.array() * viewing_zenith);

        sources.source.value.array().setZero();
        sources.source.d_ssa.array().setZero();
        sources.source.d_od.array().setZero();
        sources.source.d_b1.array().setZero();
        sources.source.d_transmission.array().setZero();
        sources.source.d_secant.array().setZero();
        for (int i = 0; i < 2; ++i) {
            double azi_factor = cos(i * azimuth);
            const auto& homog = solution.homog[i];
            const auto& part = solution.particular[i];

            sources.Y_plus[i].value.array() =
                sources.lpsum_plus[i].value.array() *
                    homog.X_plus.value.array() +
                sources.lpsum_minus[i].value.array() *
                    homog.X_minus.value.array();
            sources.Y_minus[i].value.array() =
                sources.lpsum_plus[i].value.array() *
                    homog.X_minus.value.array() +
                sources.lpsum_minus[i].value.array() *
                    homog.X_plus.value.array();

            // and their gradients
            sources.Y_plus[i].d_ssa.array() =
                sources.lpsum_plus[i].d_ssa.array() *
                    homog.X_plus.value.array() +
                sources.lpsum_minus[i].d_ssa.array() *
                    homog.X_minus.value.array() +
                sources.lpsum_plus[i].value.array() *
                    homog.X_plus.d_ssa.array() +
                sources.lpsum_minus[i].value.array() *
                    homog.X_minus.d_ssa.array();

            sources.Y_minus[i].d_ssa.array() =
                sources.lpsum_plus[i].d_ssa.array() *
                    homog.X_minus.value.array() +
                sources.lpsum_minus[i].d_ssa.array() *
                    homog.X_plus.value.array() +
                sources.lpsum_plus[i].value.array() *
                    homog.X_minus.d_ssa.array() +
                sources.lpsum_minus[i].value.array() *
                    homog.X_plus.d_ssa.array();

            // by b1
            sources.Y_plus[i].d_b1.array() =
                sources.lpsum_plus[i].d_b1.array() *
                    homog.X_plus.value.array() +
                sources.lpsum_minus[i].d_b1.array() *
                    homog.X_minus.value.array() +
                sources.lpsum_plus[i].value.array() *
                    homog.X_plus.d_b1.array() +
                sources.lpsum_minus[i].value.array() *
                    homog.X_minus.d_b1.array();

            sources.Y_minus[i].d_b1.array() =
                sources.lpsum_plus[i].d_b1.array() *
                    homog.X_minus.value.array() +
                sources.lpsum_minus[i].d_b1.array() *
                    homog.X_plus.value.array() +
                sources.lpsum_plus[i].value.array() *
                    homog.X_minus.d_b1.array() +
                sources.lpsum_minus[i].value.array() *
                    homog.X_plus.d_b1.array();

            // Homogeneous multipliers
            sources.H_minus[i].value.array() =
                (homog.omega.value.array() - sources.beamtrans.value.array()) /
                (1 - homog.k.value.array() * viewing_zenith);
            sources.H_plus[i].value.array() =
                (1 -
                 homog.omega.value.array() * sources.beamtrans.value.array()) /
                (1 + homog.k.value.array() * viewing_zenith);

            // And their gradients
            sources.H_minus[i].d_ssa.array() =
                homog.omega.d_ssa.array() /
                    (1 - homog.k.value.array() * viewing_zenith) +
                homog.k.d_ssa.array() * viewing_zenith *
                    sources.H_minus[i].value.array() /
                    (1 - homog.k.value.array() * viewing_zenith);

            sources.H_plus[i].d_ssa.array() =
                -homog.omega.d_ssa.array() * sources.beamtrans.value.array() /
                    (1 + homog.k.value.array() * viewing_zenith) -
                homog.k.d_ssa.array() * viewing_zenith *
                    sources.H_plus[i].value.array() /
                    (1 + homog.k.value.array() * viewing_zenith);

            // Also by b1
            sources.H_minus[i].d_b1.array() =
                homog.omega.d_b1.array() /
                    (1 - homog.k.value.array() * viewing_zenith) +
                homog.k.d_b1.array() * viewing_zenith *
                    sources.H_minus[i].value.array() /
                    (1 - homog.k.value.array() * viewing_zenith);

            sources.H_plus[i].d_b1.array() =
                -homog.omega.d_b1.array() * sources.beamtrans.value.array() /
                    (1 + homog.k.value.array() * viewing_zenith) -
                homog.k.d_b1.array() * viewing_zenith *
                    sources.H_plus[i].value.array() /
                    (1 + homog.k.value.array() * viewing_zenith);

            // omega and beamtrans have od derivatives
            sources.H_minus[i].d_od.array() =
                (homog.omega.d_od.array() - sources.beamtrans.d_od.array()) /
                (1 - homog.k.value.array() * viewing_zenith);

            sources.H_plus[i].d_od.array() =
                (-homog.omega.d_od.array() * sources.beamtrans.value.array() -
                 homog.omega.value.array() * sources.beamtrans.d_od.array()) /
                (1 + homog.k.value.array() * viewing_zenith);

            // Particular multipliers
            sources.D_plus[i].value.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (sources.E_minus.value.array() -
                 input.expsec.array() * sources.H_minus[i].value.array()) /
                (input.average_secant.array() + homog.k.value.array());
            sources.D_minus[i].value.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (sources.H_plus[i].value.array() -
                 sources.E_minus.value.array()) /
                (input.average_secant.array() - homog.k.value.array());

            // And their gradients with ssa
            sources.D_plus[i].d_ssa.array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                     (-sources.H_minus[i].d_ssa.array() *
                      input.expsec.array()) -
                 homog.k.d_ssa.array() * sources.D_plus[i].value.array()) /
                (input.average_secant.array() + homog.k.value.array());

            sources.D_minus[i].d_ssa.array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                     sources.H_plus[i].d_ssa.array() +
                 homog.k.d_ssa.array() * sources.D_minus[i].value.array()) /
                (input.average_secant.array() - homog.k.value.array());

            // and b1
            sources.D_plus[i].d_b1.array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                     (-sources.H_minus[i].d_b1.array() * input.expsec.array()) -
                 homog.k.d_b1.array() * sources.D_plus[i].value.array()) /
                (input.average_secant.array() + homog.k.value.array());

            sources.D_minus[i].d_b1.array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                     sources.H_plus[i].d_b1.array() +
                 homog.k.d_b1.array() * sources.D_minus[i].value.array()) /
                (input.average_secant.array() - homog.k.value.array());

            // D's also have od derivatives through E_minus, expsec, and H
            sources.D_plus[i].d_od.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (sources.E_minus.d_od.array() -
                 input.expsec.array() * sources.H_minus[i].d_od.array() +
                 input.average_secant.array() * input.expsec.array() *
                     sources.H_minus[i].value.array()) /
                (input.average_secant.array() + homog.k.value.array());

            sources.D_minus[i].d_od.array() =
                input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                (sources.H_plus[i].d_od.array() -
                 sources.E_minus.d_od.array()) /
                (input.average_secant.array() - homog.k.value.array());

            // Then we have straightforward transmission derivatives
            sources.D_plus[i].d_transmission(Eigen::seq(0, Eigen::last - 1)) =
                (sources.E_minus.value.array() -
                 input.expsec.array() * sources.H_minus[i].value.array()) /
                (input.average_secant.array() + homog.k.value.array());

            sources.D_minus[i].d_transmission(Eigen::seq(0, Eigen::last - 1)) =
                (sources.H_plus[i].value.array() -
                 sources.E_minus.value.array()) /
                (input.average_secant.array() - homog.k.value.array());

            // And more annoying secant derivatives
            // E and expsec have secant derivatives
            sources.D_plus[i].d_secant.array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                     (sources.E_minus.d_secant.array() +
                      input.od.array() * input.expsec.array() *
                          sources.H_minus[i].value.array()) -
                 sources.D_plus[i].value.array()) /
                (input.average_secant.array() + homog.k.value.array());

            sources.D_minus[i].d_secant.array() =
                (input.transmission(Eigen::seq(0, Eigen::last - 1)).array() *
                     (-sources.E_minus.d_secant.array()) -
                 sources.D_minus[i].value.array()) /
                (input.average_secant.array() - homog.k.value.array());

            // The particular source itself
            sources.V[i].value.array() =
                part.A_plus.value.array() * sources.Y_plus[i].value.array() *
                    sources.D_minus[i].value.array() +
                part.A_minus.value.array() * sources.Y_minus[i].value.array() *
                    sources.D_plus[i].value.array();

            // with it's derivative
            sources.V[i].d_ssa.array() =
                (part.A_plus.d_ssa.array() * sources.Y_plus[i].value.array() *
                     sources.D_minus[i].value.array() +
                 part.A_plus.value.array() * sources.Y_plus[i].d_ssa.array() *
                     sources.D_minus[i].value.array() +
                 part.A_plus.value.array() * sources.Y_plus[i].value.array() *
                     sources.D_minus[i].d_ssa.array()) +
                (part.A_minus.d_ssa.array() * sources.Y_minus[i].value.array() *
                     sources.D_plus[i].value.array() +
                 part.A_minus.value.array() * sources.Y_minus[i].d_ssa.array() *
                     sources.D_plus[i].value.array() +
                 part.A_minus.value.array() * sources.Y_minus[i].value.array() *
                     sources.D_plus[i].d_ssa.array());

            sources.V[i].d_b1.array() =
                (part.A_plus.d_b1.array() * sources.Y_plus[i].value.array() *
                     sources.D_minus[i].value.array() +
                 part.A_plus.value.array() * sources.Y_plus[i].d_b1.array() *
                     sources.D_minus[i].value.array() +
                 part.A_plus.value.array() * sources.Y_plus[i].value.array() *
                     sources.D_minus[i].d_b1.array()) +
                (part.A_minus.d_b1.array() * sources.Y_minus[i].value.array() *
                     sources.D_plus[i].value.array() +
                 part.A_minus.value.array() * sources.Y_minus[i].d_b1.array() *
                     sources.D_plus[i].value.array() +
                 part.A_minus.value.array() * sources.Y_minus[i].value.array() *
                     sources.D_plus[i].d_b1.array());

            // D Has extra derivatives wrt to od, transmission, and secant
            sources.V[i].d_od.array() =
                (part.A_plus.value.array() * sources.Y_plus[i].value.array() *
                     sources.D_minus[i].d_od.array() +
                 part.A_minus.value.array() * sources.Y_minus[i].value.array() *
                     sources.D_plus[i].d_od.array());

            sources.V[i]
                .d_transmission(Eigen::seq(0, Eigen::last - 1))
                .array() =
                (part.A_plus.value.array() * sources.Y_plus[i].value.array() *
                     sources.D_minus[i]
                         .d_transmission(Eigen::seq(0, Eigen::last - 1))
                         .array() +
                 part.A_minus.value.array() * sources.Y_minus[i].value.array() *
                     sources.D_plus[i]
                         .d_transmission(Eigen::seq(0, Eigen::last - 1))
                         .array());

            sources.V[i].d_secant.array() =
                (part.A_plus.value.array() * sources.Y_plus[i].value.array() *
                     sources.D_minus[i].d_secant.array() +
                 part.A_minus.value.array() * sources.Y_minus[i].value.array() *
                     sources.D_plus[i].d_secant.array());

            // Add in homogeneous sources
            sources.source.value.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(0, Eigen::last, 2), 0)
                                  .array() *
                              sources.Y_plus[i].value.array() *
                              sources.H_plus[i].value.array());

            // Append the derivative
            sources.source.d_ssa.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(0, Eigen::last, 2), 0)
                                  .array() *
                              (sources.Y_plus[i].d_ssa.array() *
                                   sources.H_plus[i].value.array() +
                               sources.Y_plus[i].value.array() *
                                   sources.H_plus[i].d_ssa.array()));

            sources.source.d_b1.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(0, Eigen::last, 2), 0)
                                  .array() *
                              (sources.Y_plus[i].d_b1.array() *
                                   sources.H_plus[i].value.array() +
                               sources.Y_plus[i].value.array() *
                                   sources.H_plus[i].d_b1.array()));

            // H_plus has extra derivatives wrt to od
            sources.source.d_od.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(0, Eigen::last, 2), 0)
                                  .array() *
                              (sources.Y_plus[i].value.array() *
                               sources.H_plus[i].d_od.array()));

            // Second term
            sources.source.value.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(1, Eigen::last, 2), 0)
                                  .array() *
                              sources.Y_minus[i].value.array() *
                              sources.H_minus[i].value.array());

            sources.source.d_ssa.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(1, Eigen::last, 2), 0)
                                  .array() *
                              (sources.Y_minus[i].d_ssa.array() *
                                   sources.H_minus[i].value.array() +
                               sources.Y_minus[i].value.array() *
                                   sources.H_minus[i].d_ssa.array()));

            sources.source.d_b1.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(1, Eigen::last, 2), 0)
                                  .array() *
                              (sources.Y_minus[i].d_b1.array() *
                                   sources.H_minus[i].value.array() +
                               sources.Y_minus[i].value.array() *
                                   sources.H_minus[i].d_b1.array()));

            // Extra derivatives wrt to od
            sources.source.d_od.array() +=
                azi_factor * (solution.bvp_coeffs[i]
                                  .rhs(Eigen::seq(1, Eigen::last, 2), 0)
                                  .array() *
                              (sources.Y_minus[i].value.array() *
                               sources.H_minus[i].d_od.array()));

            // And keep track of the gradients for backprop
            sources.d_bvp_coeff[i](Eigen::seq(0, 2 * input.nlyr - 1, 2))
                .array() =
                azi_factor * (sources.Y_plus[i].value.transpose().array() *
                              sources.H_plus[i].value.transpose().array());
            sources.d_bvp_coeff[i](Eigen::seq(1, 2 * input.nlyr - 1, 2))
                .array() =
                azi_factor * (sources.Y_minus[i].value.transpose().array() *
                              sources.H_minus[i].value.transpose().array());

            sources.source.value.array() +=
                azi_factor * (sources.V[i].value.array());
            sources.source.d_ssa.array() +=
                azi_factor * (sources.V[i].d_ssa.array());
            sources.source.d_b1.array() +=
                azi_factor * (sources.V[i].d_b1.array());
            sources.source.d_od.array() +=
                azi_factor * (sources.V[i].d_od.array());
            sources.source.d_transmission.array() +=
                azi_factor * (sources.V[i].d_transmission.array());
            sources.source.d_secant.array() +=
                azi_factor * (sources.V[i].d_secant.array());
        }
    }
} // namespace sasktran2::twostream
