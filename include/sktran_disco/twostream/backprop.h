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
        // Eigen::Map<Eigen::RowVectorXd> d_albedo;

        GradientMap(const sasktran2::atmosphere::Atmosphere<1>& atmo,
                    double* deriv_start)
            : d_extinction(deriv_start, 1,
                           atmo.storage().total_extinction.rows()),
              d_ssa(deriv_start + atmo.storage().total_extinction.rows(), 1,
                    atmo.storage().total_extinction.rows()),
              d_b1(deriv_start + 2 * atmo.storage().total_extinction.rows(), 1,
                   atmo.storage().total_extinction.rows()) {}
    };

    /**
     * Takes a gradient vector with respect to the input optical depth and maps
     * it back to the gradient vector with respect to the input atmosphere
     *
     * @param input
     * @param d_od
     * @param grad
     */
    inline void od(const Input& input, const Eigen::RowVectorXd& d_od,
                   GradientMap& grad) {
        // Map back the derivative of OD to the derivative of extinction
        // OD = (interpolating_matrix) * (extinction) * (layer_delta_height)
        // jacobian = (interpolating_matrix) * (layer_delta_height)

        grad.d_extinction +=
            (d_od.cwiseProduct((input.geometry_layers->layer_ceiling() -
                                input.geometry_layers->layer_floor())
                                   .transpose()) *
             input.geometry_layers->interpolating_matrix());
    }

    /**
     * Takes a gradient vector with respect to the input ssa and maps it back to
     * the gradient vector with respect to the input atmosphere
     *
     * @param input
     * @param d_ssa
     * @param grad
     */
    inline void ssa(const Input& input, const Eigen::RowVectorXd& d_ssa,
                    GradientMap& grad) {
        // The equation for layer SSA is, let M = (interpolating matrix)  and .
        // be elementwise multiplication, and k, w be the sasktran2 atmosphere
        // quantities

        // SSA = M * (k . w) / (M * k) = f(w, k)

        // The jacobian with respect to w is

        // K_w = M with each row weighted by k/OD
        grad.d_ssa +=
            ((d_ssa.cwiseQuotient(input.od.transpose())) *
             input.geometry_layers->interpolating_matrix())
                .cwiseProduct(input.atmosphere->storage()
                                  .total_extinction.col(input.wavelidx)
                                  .transpose());

        // TODO: Add the jacobian with respect to k
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
            d_X_plus[i].array() *= solution[i].d_X_plus_by_ssa.array();
            d_X_minus[i].array() *= solution[i].d_X_minus_by_ssa.array();

            ssa(input, d_X_plus[i], grad);
            ssa(input, d_X_minus[i], grad);
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
            d_k[i].array() *= solution[i].d_k_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_k[i], grad);
        }
    }
    inline void homog_omega(const Input& input,
                            const std::array<HomogSolution, 2>& solution,
                            std::array<Eigen::RowVectorXd, 2>& d_omega,
                            GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_omega[i].array() *= solution[i].d_omega_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_omega[i], grad);
        }
    }

    inline void particular_Q_plus(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_Q_plus, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_Q_plus[i].array() *= solution[i].d_Q_plus_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_Q_plus[i], grad);
        }
    }

    inline void particular_Q_minus(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_Q_minus, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_Q_minus[i].array() *= solution[i].d_Q_minus_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_Q_minus[i], grad);
        }
    }

    inline void particular_norm(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_norm, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_norm[i].array() *= solution[i].d_norm_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_norm[i], grad);
        }
    }

    inline void particular_Aplus(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_Aplus, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_Aplus[i].array() *= solution[i].d_A_plus_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_Aplus[i], grad);
        }
    }

    inline void particular_Aminus(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_Aminus, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_Aminus[i].array() *= solution[i].d_A_minus_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_Aminus[i], grad);
        }
    }

    inline void particular_G_plus_top(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_G_plus_top, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_G_plus_top[i].array() *= solution[i].d_G_plus_top_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_G_plus_top[i], grad);
        }
    }

    inline void particular_G_plus_bottom(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_G_plus_bottom, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_G_plus_bottom[i].array() *=
                solution[i].d_G_plus_bottom_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_G_plus_bottom[i], grad);
        }
    }

    inline void particular_G_minus_top(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_G_minus_top, GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_G_minus_top[i].array() *=
                solution[i].d_G_minus_top_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_G_minus_top[i], grad);
        }
    }

    inline void particular_G_minus_bottom(
        const Input& input, const std::array<ParticularSolution, 2>& solution,
        std::array<Eigen::RowVectorXd, 2>& d_G_minus_bottom,
        GradientMap& grad) {
        for (int i = 0; i < 2; ++i) {
            d_G_minus_bottom[i].array() *=
                solution[i].d_G_minus_bottom_by_ssa.array();

            // This is now grad wrt to ssa
            ssa(input, d_G_minus_bottom[i], grad);
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

                // TODO: When albedo!=0 this changes
                bvp_coeffs[i].d_G_minus_bottom(input.nlyr - 1) =
                    -d_coeffs[i](Eigen::last, j);

                bvp_coeffs[i].d_G_plus_top.array() *=
                    solution[i].d_G_plus_top_by_ssa.array();
                ssa(input, bvp_coeffs[i].d_G_plus_top, grad[j]);

                bvp_coeffs[i].d_G_plus_bottom.array() *=
                    solution[i].d_G_plus_bottom_by_ssa.array();
                ssa(input, bvp_coeffs[i].d_G_plus_bottom, grad[j]);

                bvp_coeffs[i].d_G_minus_top.array() *=
                    solution[i].d_G_minus_top_by_ssa.array();
                ssa(input, bvp_coeffs[i].d_G_minus_top, grad[j]);

                bvp_coeffs[i].d_G_minus_bottom.array() *=
                    solution[i].d_G_minus_bottom_by_ssa.array();
                ssa(input, bvp_coeffs[i].d_G_minus_bottom, grad[j]);

                // Now for the hard part, the pentadiagonal system (e, c, d, a,
                // b) The derivative term is -(dA)z, where A is our
                // pentadiagonal matrix and z is our BVP solution The
                // derivatives of A are given by bvp_coeffs[i].d_d_by_ssa etc.

                bvp_coeffs[i].d_temp = (-bvp_coeffs[i].rhs.cwiseProduct(
                                            d_coeffs[i](Eigen::all, j)))
                                           .transpose() *
                                       bvp_coeffs[i].d_d_by_ssa;
                ssa(input, bvp_coeffs[i].d_temp, grad[j]);

                bvp_coeffs[i].d_temp =
                    (-bvp_coeffs[i]
                          .rhs(Eigen::seq(1, Eigen::last), 0)
                          .cwiseProduct(
                              d_coeffs[i](Eigen::seq(0, Eigen::last - 1), j)))
                        .transpose() *
                    bvp_coeffs[i].d_a_by_ssa(Eigen::seq(0, Eigen::last - 1),
                                             Eigen::all);
                ssa(input, bvp_coeffs[i].d_temp, grad[j]);

                bvp_coeffs[i].d_temp =
                    (-bvp_coeffs[i]
                          .rhs(Eigen::seq(2, Eigen::last), 0)
                          .cwiseProduct(
                              d_coeffs[i](Eigen::seq(0, Eigen::last - 2), j)))
                        .transpose() *
                    bvp_coeffs[i].d_b_by_ssa(Eigen::seq(0, Eigen::last - 2),
                                             Eigen::all);
                ssa(input, bvp_coeffs[i].d_temp, grad[j]);

                bvp_coeffs[i].d_temp =
                    (-bvp_coeffs[i]
                          .rhs(Eigen::seq(0, Eigen::last - 1), 0)
                          .cwiseProduct(
                              d_coeffs[i](Eigen::seq(1, Eigen::last), j)))
                        .transpose() *
                    bvp_coeffs[i].d_c_by_ssa(Eigen::seq(1, Eigen::last),
                                             Eigen::all);
                ssa(input, bvp_coeffs[i].d_temp, grad[j]);

                bvp_coeffs[i].d_temp =
                    (-bvp_coeffs[i]
                          .rhs(Eigen::seq(0, Eigen::last - 2), 0)
                          .cwiseProduct(
                              d_coeffs[i](Eigen::seq(2, Eigen::last), j)))
                        .transpose() *
                    bvp_coeffs[i].d_e_by_ssa(Eigen::seq(2, Eigen::last),
                                             Eigen::all);
                ssa(input, bvp_coeffs[i].d_temp, grad[j]);
            }
        }
    }

    inline void full(const Input& input, Solution& solution,
                     const Sources& sources,
                     const Eigen::RowVectorXd& source_weights,
                     GradientMap& grad) {
        Eigen::RowVectorXd temp =
            sources.d_source_by_ssa.cwiseProduct(source_weights);

        ssa(input, temp, grad);

        std::array<Eigen::MatrixXd, 2> d_coeffs;
        d_coeffs[0].resize(input.nlyr * 2, 1);
        d_coeffs[1].resize(input.nlyr * 2, 1);

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

        std::vector<GradientMap> grad_vec;
        grad_vec.push_back(grad);

        bvp(input, solution.bvp_coeffs, solution.homog, solution.particular,
            d_coeffs, grad_vec);
    }
} // namespace sasktran2::twostream::backprop
