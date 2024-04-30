#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2.h>
#include "sktran_disco/twostream/solutions.h"
#include "storage.h"

namespace sasktran2::twostream::backprop {

    struct GradientMap {

        Eigen::Map<Eigen::RowVectorXd> d_extinction;
        Eigen::Map<Eigen::RowVectorXd> d_ssa;
        Eigen::Map<Eigen::RowVectorXd> d_b1;
        double& d_albedo;

        GradientMap(const sasktran2::atmosphere::Atmosphere<1>& atmo,
                    double* deriv_start)
            : d_extinction(deriv_start, 1,
                           atmo.storage().total_extinction.rows() - 1),
              d_ssa(deriv_start + atmo.storage().total_extinction.rows(), 1,
                    atmo.storage().total_extinction.rows() - 1),
              d_b1(deriv_start + 2 * atmo.storage().total_extinction.rows(), 1,
                   atmo.storage().total_extinction.rows() - 1),
              d_albedo(
                  deriv_start[atmo.storage().total_extinction.rows() - 1]) {}
    };

    inline void map_to_atmosphere(
        const Input& input, const GradientMap& internal_grad,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>& source) {
        int n = input.nlyr + 1;
        const auto& M = input.geometry_layers->interpolating_matrix();

        // Go through each atmosphere grid element
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n - 1; ++j) {
                double a = M(j, i);

                if (a == 0) {
                    continue;
                }

                // Have a weight contribution

                // For extinction we just have to multiply by the weight and the
                // layer thickness
                source.deriv(i) += a * internal_grad.d_extinction(j) *
                                   (input.geometry_layers->layer_ceiling()(j) -
                                    input.geometry_layers->layer_floor()(j));

                // For SSA,
                source.deriv(i + n) +=
                    a * internal_grad.d_ssa(j) *
                    input.atmosphere->storage().total_extinction(
                        i, input.wavelidx) /
                    (input.od(j) / (input.geometry_layers->layer_ceiling()(j) -
                                    input.geometry_layers->layer_floor()(j)));

                source.deriv(i) +=
                    a * internal_grad.d_ssa(j) *
                    (input.atmosphere->storage().ssa(i, input.wavelidx) -
                     input.ssa(j)) /
                    (input.od(j) / (input.geometry_layers->layer_ceiling()(j) -
                                    input.geometry_layers->layer_floor()(j)));

                // for b1, find the contribution to atmo b1
                double denom = (input.ssa(j) * input.od(j) /
                                (input.geometry_layers->layer_ceiling()(j) -
                                 input.geometry_layers->layer_floor()(j)));
                double b1_deriv =
                    a * internal_grad.d_b1(j) *
                    (input.atmosphere->storage().ssa(i, input.wavelidx) *
                     input.atmosphere->storage().total_extinction(
                         i, input.wavelidx)) /
                    denom;

                // to extinction
                source.deriv(i) +=
                    a * internal_grad.d_b1(j) *
                    input.atmosphere->storage().ssa(i, input.wavelidx) *
                    (input.atmosphere->storage().leg_coeff(1, i,
                                                           input.wavelidx) -
                     input.b1(j)) /
                    denom;

                // to ssa
                source.deriv(i) += a * internal_grad.d_b1(j) *
                                   input.atmosphere->storage().total_extinction(
                                       i, input.wavelidx) *
                                   (input.atmosphere->storage().leg_coeff(
                                        1, i, input.wavelidx) -
                                    input.b1(j)) /
                                   denom;

                // To scattering derivatives
                for (int k = 0;
                     k < input.atmosphere->num_scattering_deriv_groups(); ++k) {
                    source.deriv(input.atmosphere->scat_deriv_start_index() +
                                 k * n + i) +=
                        input.atmosphere->storage().d_leg_coeff(
                            1, i, input.wavelidx, k) *
                        b1_deriv;
                }
            }
        }

        // Albedo derivatives
        for (int i = 0; i < input.atmosphere->surface().num_deriv(); ++i) {
            source.deriv(input.atmosphere->surface_deriv_start_index() + i) +=
                internal_grad.d_albedo;
        }
    }

    /**
     * Takes a gradient vector with respect to the input optical depth and maps
     * it back to the gradient vector with respect to the input atmosphere
     *
     * @param input
     * @param d_od
     * @param grad
     */
    template <typename Derived>
    inline void od(const Input& input, const Eigen::MatrixBase<Derived>& d_od,
                   GradientMap& grad) {
        // return;
        grad.d_extinction += d_od;
    }

    /**
     * Takes a gradient vector with respect to the input ssa and maps it back to
     * the gradient vector with respect to the input atmosphere
     *
     * @param input
     * @param d_ssa
     * @param grad
     */
    template <typename Derived>
    inline void ssa(const Input& input, const Eigen::MatrixBase<Derived>& d_ssa,
                    GradientMap& grad) {
        grad.d_ssa += d_ssa;
    }

    /**
     * Takes a gradient vector with respect to the input b1 and appends it to
     * the full b1 derivative
     *
     * @param input
     * @param d_ssa
     * @param grad
     */
    template <typename Derived>
    inline void b1(const Input& input, const Eigen::MatrixBase<Derived>& d_b1,
                   GradientMap& grad) {

        grad.d_b1 += d_b1;
    }

    template <typename Derived>
    inline void transmission(const Input& input,
                             const Eigen::MatrixBase<Derived>& d_transmission,
                             GradientMap& grad) {
        // Transmission is calculated through
        /*
            transmission(0) = 0;
            transmission(Eigen::seq(1, nlyr)) =
                -1 * (geometry_layers->chapman_factors() * od);

            transmission.array() = transmission.array().exp();
        */

        // d_transmission_i / by d_od = -chapman_factors(i, :) * transmission

        // First transmission is always 1
        for (int i = 1; i < d_transmission.size(); ++i) {
            od(input,
               d_transmission(i) *
                   input.geometry_layers->chapman_factors().row(i - 1) *
                   (-input.transmission(i)),
               grad);
        }
    }

    template <typename Derived>
    inline void secant(const Input& input,
                       const Eigen::MatrixBase<Derived>& d_secant,
                       GradientMap& grad) {
        // Let M = chapman factors
        const auto& M = input.geometry_layers->chapman_factors();
        // Then y = M * od

        // secant(i) = (y(i) - y(i-1)) / od(i)

        // With the special case of secant(0) = y(0) / od(0)

        for (int i = 0; i < d_secant.size(); ++i) {
            // directly add the dense derivative to extinction
            grad.d_extinction(i) +=
                -d_secant(i) * input.average_secant(i) / input.od(i);

            // always add in the term for y(0)
            grad.d_extinction += M.row(i) / input.od(i) * d_secant(i);

            if (i > 0) {
                // add in the term for y(i-1)
                grad.d_extinction += -M.row(i - 1) / input.od(i) * d_secant(i);
            }
        }
    }

    inline void homog_X(const Input& input,
                        const std::array<HomogSolution, 2>& solution,
                        std::array<Eigen::RowVectorXd, 2>& d_X_plus,
                        std::array<Eigen::RowVectorXd, 2>& d_X_minus,
                        GradientMap& grad) {
        // We also have
        // X_plus = 0.5 * (1 - s / k)
        // X_minus = 0.5 * (1 + s / k)

        // The gradient wrt to s is  -/+ 0.5 / k

        // The gradient wrt to k is +/- 0.5 s / k^2, but since k^2 = s*d, this
        // is +/- 0.5 / d

        for (int i = 0; i < 2; ++i) {
            // Now do X_plus and X_minus

            ssa(input,
                d_X_plus[i].cwiseProduct(solution[i].X_plus.d_ssa.transpose()),
                grad);
            ssa(input,
                d_X_minus[i].cwiseProduct(
                    solution[i].X_minus.d_ssa.transpose()),
                grad);
        }
    }
    inline void homog_k(const Input& input,
                        const std::array<HomogSolution, 2>& solution,
                        std::array<Eigen::RowVectorXd, 2>& d_k,
                        GradientMap& grad) {
        // To start, we have k = sqrt(s * d) and so then
        // The gradient wrt to s is 0.5 / k * d
        // The gradient wrt to d is 0.5 / k * s

        for (int i = 0; i < 2; ++i) {
            // This is now grad wrt to ssa
            ssa(input, d_k[i].cwiseProduct(solution[i].k.d_ssa.transpose()),
                grad);
        }
    }
    inline void homog_omega(const Input& input,
                            const std::array<HomogSolution, 2>& solution,
                            std::array<Eigen::RowVectorXd, 2>& d_omega,
                            GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_omega[i].cwiseProduct(solution[i].omega.d_ssa.transpose()),
                grad);
        }
    }

    inline void particular_Q_plus(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_Q_plus, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_Q_plus[i].cwiseProduct(solution[i].Q_plus.d_ssa.transpose()),
                grad);
        }
    }

    inline void particular_Q_minus(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_Q_minus, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_Q_minus[i].cwiseProduct(
                    solution[i].Q_minus.d_ssa.transpose()),
                grad);
        }
    }

    inline void particular_norm(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_norm, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_norm[i].cwiseProduct(solution[i].norm.d_ssa.transpose()),
                grad);
        }
    }

    inline void particular_Aplus(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_Aplus, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_Aplus[i].cwiseProduct(solution[i].A_plus.d_ssa.transpose()),
                grad);
        }
    }

    inline void particular_Aminus(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_Aminus, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_Aminus[i].cwiseProduct(solution[i].A_minus.d_ssa.transpose()),
                grad);
        }
    }

    inline void
    particular_G_plus_top(const Input& input,
                          const std::array<ParticularSolution, 2>& solution,
                          const std::array<Eigen::RowVectorXd, 2>& d_G_plus_top,
                          GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            ssa(input,
                d_G_plus_top[i].cwiseProduct(
                    solution[i].G_plus_top.d_ssa.transpose()),
                grad);
            b1(input,
               d_G_plus_top[i].cwiseProduct(
                   solution[i].G_plus_top.d_b1.transpose()),
               grad);
            od(input,
               d_G_plus_top[i].cwiseProduct(
                   solution[i].G_plus_top.d_od.transpose()),
               grad);
            transmission(input,
                         d_G_plus_top[i].cwiseProduct(
                             solution[i]
                                 .G_plus_top
                                 .d_transmission(Eigen::seq(0, Eigen::last - 1))
                                 .transpose()),
                         grad);
            secant(input,
                   d_G_plus_top[i].cwiseProduct(
                       solution[i].G_plus_top.d_secant.transpose()),
                   grad);
        }
    }

    inline void particular_G_plus_bottom(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        const std::array<Eigen::RowVectorXd, 2>& d_G_plus_bottom,
        GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            ssa(input,
                d_G_plus_bottom[i].cwiseProduct(
                    solution[i].G_plus_bottom.d_ssa.transpose()),
                grad);
            b1(input,
               d_G_plus_bottom[i].cwiseProduct(
                   solution[i].G_plus_bottom.d_b1.transpose()),
               grad);
            od(input,
               d_G_plus_bottom[i].cwiseProduct(
                   solution[i].G_plus_bottom.d_od.transpose()),
               grad);
            transmission(input,
                         d_G_plus_bottom[i].cwiseProduct(
                             solution[i]
                                 .G_plus_bottom
                                 .d_transmission(Eigen::seq(0, Eigen::last - 1))
                                 .transpose()),
                         grad);
            secant(input,
                   d_G_plus_bottom[i].cwiseProduct(
                       solution[i].G_plus_bottom.d_secant.transpose()),
                   grad);
        }
    }

    inline void particular_G_minus_top(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        const std::array<Eigen::RowVectorXd, 2>& d_G_minus_top,
        GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            ssa(input,
                d_G_minus_top[i].cwiseProduct(
                    solution[i].G_minus_top.d_ssa.transpose()),
                grad);
            b1(input,
               d_G_minus_top[i].cwiseProduct(
                   solution[i].G_minus_top.d_b1.transpose()),
               grad);
            od(input,
               d_G_minus_top[i].cwiseProduct(
                   solution[i].G_minus_top.d_od.transpose()),
               grad);
            transmission(input,
                         d_G_minus_top[i].cwiseProduct(
                             solution[i]
                                 .G_minus_top
                                 .d_transmission(Eigen::seq(0, Eigen::last - 1))
                                 .transpose()),
                         grad);
            secant(input,
                   d_G_minus_top[i].cwiseProduct(
                       solution[i].G_minus_top.d_secant.transpose()),
                   grad);
        }
    }

    inline void particular_G_minus_bottom(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        const std::array<Eigen::RowVectorXd, 2>& d_G_minus_bottom,
        GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_G_minus_bottom[i].cwiseProduct(
                    solution[i].G_minus_bottom.d_ssa.transpose()),
                grad);
            od(input,
               d_G_minus_bottom[i].cwiseProduct(
                   solution[i].G_minus_bottom.d_od.transpose()),
               grad);
            b1(input,
               d_G_minus_bottom[i].cwiseProduct(
                   solution[i].G_minus_bottom.d_b1.transpose()),
               grad);
            transmission(input,
                         d_G_minus_bottom[i].cwiseProduct(
                             solution[i]
                                 .G_minus_bottom
                                 .d_transmission(Eigen::seq(0, Eigen::last - 1))
                                 .transpose()),
                         grad);
            secant(input,
                   d_G_minus_bottom[i].cwiseProduct(
                       solution[i].G_minus_bottom.d_secant.transpose()),
                   grad);
        }
    }

    inline void bvp(const Input& input, std::array<BVPCoeffs, 2>& bvp_coeffs,
                    const std::array<HomogSolution, 2>& homog,
                    const std::array<ParticularSolution, 2>& solution,
                    std::array<Eigen::MatrixXd, 2>& d_coeffs,
                    std::vector<GradientMap>& grad) {
        for (int i = 0; i < 2; ++i) {
            // We have to solve the adjoint equation to get the gradient
            pentadiagonal_solve(bvp_coeffs[i], d_coeffs[i], true);

            for (int j = 0; j < grad.size(); ++j) {
                // Now we have to backprop these through the boundary value
                // problem We have 6 BVP vectors that we have to account for (e,
                // c, d, a, b) for the system and RHS for the solution

                // Let's start with the RHS, we have contributions to all four G
                // terms
                bvp_coeffs[i].d_G_plus_top(0) = -d_coeffs[i](0, j);

                bvp_coeffs[i]
                    .d_G_plus_top(Eigen::seq(1, input.nlyr - 1))
                    .array() = d_coeffs[i](Eigen::seq(2, Eigen::last, 2), j)
                                   .transpose()
                                   .array();

                bvp_coeffs[i]
                    .d_G_plus_bottom(Eigen::seq(0, input.nlyr - 2))
                    .array() = -d_coeffs[i](Eigen::seq(2, Eigen::last, 2), j)
                                    .transpose()
                                    .array();

                bvp_coeffs[i]
                    .d_G_minus_top(Eigen::seq(1, input.nlyr - 1))
                    .array() = d_coeffs[i](Eigen::seq(1, Eigen::last - 1, 2), j)
                                   .transpose()
                                   .array();
                bvp_coeffs[i]
                    .d_G_minus_bottom(Eigen::seq(0, input.nlyr - 2))
                    .array() =
                    -d_coeffs[i](Eigen::seq(1, Eigen::last - 1, 2), j)
                         .transpose()
                         .array();

                bvp_coeffs[i].d_G_minus_bottom(input.nlyr - 1) =
                    -d_coeffs[i](Eigen::last, j) * (1);

                bvp_coeffs[i].d_G_plus_bottom(input.nlyr - 1) =
                    -d_coeffs[i](Eigen::last, j) *
                    (-2 * input.mu * input.albedo *
                     sasktran_disco::kronDelta(i, 0));

                // Backprop the SSA factors
                ssa(input,
                    bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                        solution[i].G_plus_top.d_ssa.transpose()),
                    grad[j]);
                ssa(input,
                    bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                        solution[i].G_plus_bottom.d_ssa.transpose()),
                    grad[j]);
                ssa(input,
                    bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                        solution[i].G_minus_top.d_ssa.transpose()),
                    grad[j]);
                ssa(input,
                    bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                        solution[i].G_minus_bottom.d_ssa.transpose()),
                    grad[j]);

                // b1 factors
                b1(input,
                   bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                       solution[i].G_plus_top.d_b1.transpose()),
                   grad[j]);
                b1(input,
                   bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                       solution[i].G_plus_bottom.d_b1.transpose()),
                   grad[j]);
                b1(input,
                   bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                       solution[i].G_minus_top.d_b1.transpose()),
                   grad[j]);
                b1(input,
                   bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                       solution[i].G_minus_bottom.d_b1.transpose()),
                   grad[j]);

                // And the OD factors
                od(input,
                   bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                       solution[i].G_plus_top.d_od.transpose()),
                   grad[j]);
                od(input,
                   bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                       solution[i].G_plus_bottom.d_od.transpose()),
                   grad[j]);
                od(input,
                   bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                       solution[i].G_minus_top.d_od.transpose()),
                   grad[j]);
                od(input,
                   bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                       solution[i].G_minus_bottom.d_od.transpose()),
                   grad[j]);

                // And the transmission factors

                bvp_coeffs[0].d_temp_transmission(
                    Eigen::seq(0, input.nlyr - 1)) +=
                    bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                        solution[i]
                            .G_plus_top
                            .d_transmission(Eigen::seq(0, Eigen::last - 1))
                            .transpose()) +
                    bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                        solution[i]
                            .G_plus_bottom
                            .d_transmission(Eigen::seq(0, Eigen::last - 1))
                            .transpose()) +
                    bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                        solution[i]
                            .G_minus_top
                            .d_transmission(Eigen::seq(0, Eigen::last - 1))
                            .transpose()) +
                    bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                        solution[i]
                            .G_minus_bottom
                            .d_transmission(Eigen::seq(0, Eigen::last - 1))
                            .transpose());

                bvp_coeffs[0].d_temp_transmission(input.nlyr) +=
                    input.csz * input.albedo * sasktran_disco::kronDelta(i, 0) /
                    EIGEN_PI * d_coeffs[i](2 * input.nlyr - 1, 0);

                // And finally the secant factors
                bvp_coeffs[0].d_temp_secant +=
                    bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                        solution[i].G_plus_top.d_secant.transpose()) +
                    bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                        solution[i].G_plus_bottom.d_secant.transpose()) +
                    bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                        solution[i].G_minus_top.d_secant.transpose()) +
                    bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                        solution[i].G_minus_bottom.d_secant.transpose());

                // Directly assign the albedo derivative
                grad[j].d_albedo +=
                    d_coeffs[i](2 * input.nlyr - 1, 0) *
                    (input.csz / EIGEN_PI * input.transmission(Eigen::last) +
                     2 * input.mu *
                         solution[i].G_plus_bottom.value(Eigen::last));

                // Now for the hard part, the pentadiagonal system (e, c, d, a,
                // b) The derivative term is -(dA)z, where A is our
                // pentadiagonal matrix and z is our BVP solution The
                // derivatives of A are given by bvp_coeffs[i].d_d_by_ssa etc.

                /** SSA Derivatives */
                bvp_coeffs[i].d_temp_ssa.setZero();
                // TOA terms
                bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(0, 0) *
                                               d_coeffs[i](0, j) *
                                               bvp_coeffs[i].d_d_by_ssa(0, 0);
                bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(1, 0) *
                                               d_coeffs[i](0, j) *
                                               bvp_coeffs[i].d_a_by_ssa(0, 0);
                // Continuity terms
                for (int k = 0; k < input.nlyr - 1; ++k) {
                    // d factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs(2 * k + 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_d_by_ssa(2 * k + 1, k);

                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs(2 * k + 2, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_d_by_ssa(2 * k + 2, k + 1);

                    // a factors
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) + 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_a_by_ssa(2 * k + 1, k + 1);
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) + 1, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_a_by_ssa(2 * k + 2, k + 1);

                    // b factors
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) + 2, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_b_by_ssa(2 * k + 1, k + 1);

                    // c factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) - 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_c_by_ssa(2 * k + 1, k);
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) - 1, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_c_by_ssa(2 * k + 2, k);

                    // e factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) - 2, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_e_by_ssa(2 * k + 2, k);
                }

                // Ground term
                bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                    -bvp_coeffs[i].rhs(2 * input.nlyr - 1, 0) *
                    d_coeffs[i](2 * input.nlyr - 1, j) *
                    bvp_coeffs[i].d_d_by_ssa(2 * input.nlyr - 1,
                                             input.nlyr - 1);
                bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                    -bvp_coeffs[i].rhs(2 * input.nlyr - 1 - 1, 0) *
                    d_coeffs[i](2 * input.nlyr - 1, j) *
                    bvp_coeffs[i].d_c_by_ssa(2 * input.nlyr - 1,
                                             input.nlyr - 1);

                ssa(input, bvp_coeffs[i].d_temp_ssa, grad[j]);

                /** OD Derivatives */
                bvp_coeffs[i].d_temp_ssa.setZero();
                // TOA terms
                bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(0, 0) *
                                               d_coeffs[i](0, j) *
                                               bvp_coeffs[i].d_d_by_od(0, 0);
                bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(1, 0) *
                                               d_coeffs[i](0, j) *
                                               bvp_coeffs[i].d_a_by_od(0, 0);
                // Continuity terms
                for (int k = 0; k < input.nlyr - 1; ++k) {
                    // d factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs(2 * k + 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_d_by_od(2 * k + 1, k);

                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs(2 * k + 2, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_d_by_od(2 * k + 2, k + 1);

                    // a factors
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) + 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_a_by_od(2 * k + 1, k + 1);
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) + 1, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_a_by_od(2 * k + 2, k + 1);

                    // b factors
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) + 2, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_b_by_od(2 * k + 1, k + 1);

                    // c factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) - 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_c_by_od(2 * k + 1, k);
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) - 1, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_c_by_od(2 * k + 2, k);

                    // e factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) - 2, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_e_by_od(2 * k + 2, k);
                }

                // Ground term
                bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                    -bvp_coeffs[i].rhs(2 * input.nlyr - 1, 0) *
                    d_coeffs[i](2 * input.nlyr - 1, j) *
                    bvp_coeffs[i].d_d_by_od(2 * input.nlyr - 1, input.nlyr - 1);
                bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                    -bvp_coeffs[i].rhs(2 * input.nlyr - 1 - 1, 0) *
                    d_coeffs[i](2 * input.nlyr - 1, j) *
                    bvp_coeffs[i].d_c_by_od(2 * input.nlyr - 1, input.nlyr - 1);

                od(input, bvp_coeffs[i].d_temp_ssa, grad[j]);

                /** b1 derivatives */
                bvp_coeffs[i].d_temp_ssa.setZero();
                // TOA terms
                bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(0, 0) *
                                               d_coeffs[i](0, j) *
                                               bvp_coeffs[i].d_d_by_b1(0, 0);
                bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(1, 0) *
                                               d_coeffs[i](0, j) *
                                               bvp_coeffs[i].d_a_by_b1(0, 0);
                // Continuity terms
                for (int k = 0; k < input.nlyr - 1; ++k) {
                    // d factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs(2 * k + 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_d_by_b1(2 * k + 1, k);

                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs(2 * k + 2, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_d_by_b1(2 * k + 2, k + 1);

                    // a factors
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) + 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_a_by_b1(2 * k + 1, k + 1);
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) + 1, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_a_by_b1(2 * k + 2, k + 1);

                    // b factors
                    bvp_coeffs[i].d_temp_ssa(k + 1) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) + 2, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_b_by_b1(2 * k + 1, k + 1);

                    // c factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 1) - 1, 0) *
                        d_coeffs[i](2 * k + 1, j) *
                        bvp_coeffs[i].d_c_by_b1(2 * k + 1, k);
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) - 1, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_c_by_b1(2 * k + 2, k);

                    // e factors
                    bvp_coeffs[i].d_temp_ssa(k) +=
                        -bvp_coeffs[i].rhs((2 * k + 2) - 2, 0) *
                        d_coeffs[i](2 * k + 2, j) *
                        bvp_coeffs[i].d_e_by_b1(2 * k + 2, k);
                }

                // Ground term
                bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                    -bvp_coeffs[i].rhs(2 * input.nlyr - 1, 0) *
                    d_coeffs[i](2 * input.nlyr - 1, j) *
                    bvp_coeffs[i].d_d_by_b1(2 * input.nlyr - 1, input.nlyr - 1);
                bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                    -bvp_coeffs[i].rhs(2 * input.nlyr - 1 - 1, 0) *
                    d_coeffs[i](2 * input.nlyr - 1, j) *
                    bvp_coeffs[i].d_c_by_b1(2 * input.nlyr - 1, input.nlyr - 1);

                b1(input, bvp_coeffs[i].d_temp_ssa, grad[j]);

                // BVP matrix albedo derivatives
                // c and d have albedo derivatives

                // for c
                double d_albedo = -2 * input.mu *
                                  homog[i].X_plus.value(Eigen::last) *
                                  homog[i].omega.value(Eigen::last) *
                                  sasktran_disco::kronDelta(i, 0);
                grad[j].d_albedo +=
                    -bvp_coeffs[i].rhs(2 * input.nlyr - 1 - 1, 0) *
                    d_coeffs[i](2 * input.nlyr - 1, j) * d_albedo;
                // for d
                d_albedo = -2 * input.mu * homog[i].X_minus.value(Eigen::last) *
                           sasktran_disco::kronDelta(i, 0);
                grad[j].d_albedo += -bvp_coeffs[i].rhs(2 * input.nlyr - 1, 0) *
                                    d_coeffs[i](2 * input.nlyr - 1, j) *
                                    d_albedo;
            }
        }
        transmission(input, bvp_coeffs[0].d_temp_transmission, grad[0]);
        secant(input, bvp_coeffs[0].d_temp_secant, grad[0]);
    }

    inline void full(const Input& input, Solution& solution,
                     const Sources& sources,
                     const Eigen::RowVectorXd& source_weights,
                     std::array<Eigen::MatrixXd, 2>& d_coeffs,
                     GradientMap& grad, double ground_weight = 0.0) {
        // Direct source terms
        ssa(input,
            sources.source.d_ssa.transpose().cwiseProduct(source_weights),
            grad);
        od(input, sources.source.d_od.transpose().cwiseProduct(source_weights),
           grad);
        b1(input, sources.source.d_b1.transpose().cwiseProduct(source_weights),
           grad);

        solution.bvp_coeffs[0].d_temp_secant =
            sources.source.d_secant.transpose().cwiseProduct(source_weights);

        solution.bvp_coeffs[0].d_temp_transmission(
            Eigen::seq(0, input.nlyr - 1)) =
            sources.source.d_transmission(Eigen::seq(0, input.nlyr - 1))
                .transpose()
                .cwiseProduct(source_weights);
        solution.bvp_coeffs[0].d_temp_transmission(input.nlyr) = 0.0;

        // Transmission contributions to the ground multiple scatter source

        solution.bvp_coeffs[0].d_temp_transmission(input.nlyr - 1) +=
            solution.particular[0].G_plus_bottom.d_transmission(input.nlyr -
                                                                1) *
            ground_weight;

        // Secant contributions to the ground multiple scatter source
        solution.bvp_coeffs[0].d_temp_secant(input.nlyr - 1) +=
            solution.particular[0].G_plus_bottom.d_secant(input.nlyr - 1) *
            ground_weight;

        // Direct ground multiple scatter contributions

        grad.d_ssa(input.nlyr - 1) +=
            (solution.particular[0].G_plus_bottom.d_ssa(input.nlyr - 1) +
             solution.bvp_coeffs[0].rhs(Eigen::last - 1, 0) *
                 (solution.homog[0].X_plus.d_ssa(Eigen::last) *
                      solution.homog[0].omega.value(Eigen::last) +
                  solution.homog[0].X_plus.value(Eigen::last) *
                      solution.homog[0].omega.d_ssa(Eigen::last)) +
             solution.bvp_coeffs[0].rhs(Eigen::last, 0) *
                 (solution.homog[0].X_minus.d_ssa(Eigen::last))) *
            ground_weight;

        grad.d_b1(input.nlyr - 1) +=
            (solution.particular[0].G_plus_bottom.d_b1(input.nlyr - 1) +
             solution.bvp_coeffs[0].rhs(Eigen::last - 1, 0) *
                 (solution.homog[0].X_plus.d_b1(Eigen::last) *
                      solution.homog[0].omega.value(Eigen::last) +
                  solution.homog[0].X_plus.value(Eigen::last) *
                      solution.homog[0].omega.d_b1(Eigen::last)) +
             solution.bvp_coeffs[0].rhs(Eigen::last, 0) *
                 (solution.homog[0].X_minus.d_b1(Eigen::last))) *
            ground_weight;

        grad.d_extinction(input.nlyr - 1) +=
            (solution.particular[0].G_plus_bottom.d_od(input.nlyr - 1) +
             solution.bvp_coeffs[0].rhs(Eigen::last - 1, 0) *
                 solution.homog[0].X_plus.value(Eigen::last) *
                 solution.homog[0].omega.d_od(Eigen::last)) *
            ground_weight;

        // BVP adjoint solution
        d_coeffs[0].col(0) = sources.d_bvp_coeff[0];
        d_coeffs[1].col(0) = sources.d_bvp_coeff[1];

        d_coeffs[0](Eigen::seq(0, 2 * input.nlyr - 1, 2), 0).array() *=
            source_weights.transpose().array();
        d_coeffs[0](Eigen::seq(1, 2 * input.nlyr - 1, 2), 0).array() *=
            source_weights.transpose().array();

        d_coeffs[1](Eigen::seq(0, 2 * input.nlyr - 1, 2), 0).array() *=
            source_weights.transpose().array();
        d_coeffs[1](Eigen::seq(1, 2 * input.nlyr - 1, 2), 0).array() *=
            source_weights.transpose().array();

        // Albedo contributions to BVP coefficients

        d_coeffs[0](2 * input.nlyr - 1, 0) +=
            ground_weight * (solution.homog[0].X_minus.value(Eigen::last));

        d_coeffs[0](2 * input.nlyr - 2, 0) +=
            ground_weight * (solution.homog[0].X_plus.value(Eigen::last) *
                             solution.homog[0].omega.value(Eigen::last));

        std::vector<GradientMap> grad_vec;
        grad_vec.push_back(grad);

        bvp(input, solution.bvp_coeffs, solution.homog, solution.particular,
            d_coeffs, grad_vec);
    }
} // namespace sasktran2::twostream::backprop
