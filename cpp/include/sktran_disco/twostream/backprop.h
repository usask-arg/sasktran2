#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2.h>
#include "sktran_disco/twostream/solutions.h"
#include "storage.h"
#include "meta.h"

namespace sasktran2::twostream::backprop {

    template <SourceType Source> struct GradientMap {

        static constexpr int scalar_offset(int n_layer) {
            if constexpr (has_thermal<Source>()) {
                return 5 * n_layer;
            }

            return 3 * n_layer;
        }

        Eigen::Map<Eigen::RowVectorXd> d_extinction;
        Eigen::Map<Eigen::RowVectorXd> d_ssa;
        Eigen::Map<Eigen::RowVectorXd> d_b1;
        Eigen::Map<Eigen::RowVectorXd> d_thermal_b0;
        Eigen::Map<Eigen::RowVectorXd> d_thermal_b1;
        double& d_albedo;
        double& d_thermal_surf;

        GradientMap(const sasktran2::atmosphere::Atmosphere<1>& atmo,
                    double* deriv_start)
            : d_extinction(nullptr, 0), d_ssa(nullptr, 0), d_b1(nullptr, 0),
              d_thermal_b0(nullptr, 0), d_thermal_b1(nullptr, 0),
              d_albedo(
                  *(deriv_start +
                    scalar_offset(atmo.storage().total_extinction.rows() - 1))),
              d_thermal_surf(
                  *(deriv_start +
                    scalar_offset(atmo.storage().total_extinction.rows() - 1) +
                    (has_thermal<Source>() ? 1 : 0))) {
            const int n_layer = atmo.storage().total_extinction.rows() - 1;

            new (&d_extinction)
                Eigen::Map<Eigen::RowVectorXd>(deriv_start, n_layer);

            new (&d_ssa)
                Eigen::Map<Eigen::RowVectorXd>(deriv_start + n_layer, n_layer);

            new (&d_b1) Eigen::Map<Eigen::RowVectorXd>(
                deriv_start + 2 * n_layer, n_layer);

            if constexpr (has_thermal<Source>()) {
                new (&d_thermal_b0) Eigen::Map<Eigen::RowVectorXd>(
                    deriv_start + 3 * n_layer, n_layer);

                new (&d_thermal_b1) Eigen::Map<Eigen::RowVectorXd>(
                    deriv_start + 4 * n_layer, n_layer);
            }
        }
    };

    template <SourceType Source>
    inline void map_to_atmosphere(
        const Input<Source>& input, const GradientMap<Source>& internal_grad,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>& source) {
        ZoneScopedN("Twostream Backprop Map To Atmosphere");
        int n = input.nlyr + 1;
        int nlyr = input.nlyr;

        // Go through each atmosphere grid element
        for (int i = 0; i < nlyr; ++i) {
            const int top_atmo_idx = nlyr - i; // Top atmo idx of the layer
            const int bottom_atmo_idx =
                nlyr - i - 1; // Bottom atmo idx of the layer

            const double layer_thickness =
                input.geometry_layers->layer_ceiling()(i) -
                input.geometry_layers->layer_floor()(i);

            for (const auto idx : {top_atmo_idx, bottom_atmo_idx}) {
                // For extinction contribution is weight (0.5) * layer_thickness
                // equally to top and bottom
                source.deriv(idx) +=
                    0.5 * internal_grad.d_extinction(i) * layer_thickness;

                // For SSA, mapping to ssa
                source.deriv(idx + n) +=
                    0.5 * internal_grad.d_ssa(i) *
                    input.atmosphere->storage().total_extinction(
                        idx, input.wavelidx) /
                    (input.od(i) / layer_thickness);

                // and extinction contribution
                source.deriv(idx) +=
                    0.5 * internal_grad.d_ssa(i) *
                    (input.atmosphere->storage().ssa(idx, input.wavelidx) -
                     input.ssa(i)) /
                    (input.od(i) / layer_thickness);

                // for b1, find the contribution to atmo b1
                double denom = (input.ssa(i) * input.od(i) / layer_thickness);
                double b1_deriv =
                    0.5 * internal_grad.d_b1(i) *
                    (input.atmosphere->storage().ssa(idx, input.wavelidx) *
                     input.atmosphere->storage().total_extinction(
                         idx, input.wavelidx)) /
                    denom;

                // to extinction
                source.deriv(idx) +=
                    0.5 * internal_grad.d_b1(i) *
                    input.atmosphere->storage().ssa(idx, input.wavelidx) *
                    (input.atmosphere->storage().leg_coeff(1, idx,
                                                           input.wavelidx) -
                     input.b1(i)) /
                    denom;

                // to ssa
                source.deriv(idx + n) +=
                    0.5 * internal_grad.d_b1(i) *
                    input.atmosphere->storage().total_extinction(
                        idx, input.wavelidx) *
                    (input.atmosphere->storage().leg_coeff(1, idx,
                                                           input.wavelidx) -
                     input.b1(i)) /
                    denom;

                // To scattering derivatives
                for (int k = 0;
                     k < input.atmosphere->num_scattering_deriv_groups(); ++k) {
                    source.deriv(input.atmosphere->scat_deriv_start_index() +
                                 k * n + idx) +=
                        input.atmosphere->storage().d_leg_coeff(
                            1, idx, input.wavelidx, k) *
                        b1_deriv;
                }
            }

            if constexpr (has_thermal<Source>()) {
                const double thermal_top =
                    input.atmosphere->storage().emission_source(top_atmo_idx,
                                                                input.wavelidx);
                const double thermal_bottom =
                    input.atmosphere->storage().emission_source(bottom_atmo_idx,
                                                                input.wavelidx);

                // This is where it gets a little trickier
                // b0 = thermal_top, b1 = -log(thermal_bottom / thermal_top) /
                // layer_thickness

                source.deriv(input.atmosphere->emission_deriv_start_index() +
                             top_atmo_idx) += internal_grad.d_thermal_b0(i);

                // But the b1 derivatives have contributions to top, bottom, and
                // layer thickness (which maps to extinction)
                double b1_deriv = internal_grad.d_thermal_b1(i);

                // derivatives are
                // db1 / d thermal_top = 1 / (layer_thickness * thermal_top)
                // db1 / d thermal_bottom = -1 / (layer_thickness *
                // thermal_bottom) db1 / d layer_thickness = -b1 /
                // layer_thickness

                source.deriv(input.atmosphere->emission_deriv_start_index() +
                             top_atmo_idx) +=
                    b1_deriv / (layer_thickness * thermal_top);
                source.deriv(input.atmosphere->emission_deriv_start_index() +
                             bottom_atmo_idx) +=
                    -b1_deriv / (layer_thickness * thermal_bottom);

                // extinction contributes equally to top and bottom
                for (const auto& idx : {top_atmo_idx, bottom_atmo_idx}) {
                    source.deriv(idx) +=
                        b1_deriv * input.b1_thermal(i) * 0.5 / layer_thickness;
                }
            }
        }

        // Albedo derivatives
        for (int i = 0; i < input.atmosphere->surface().num_deriv(); ++i) {
            source.deriv(input.atmosphere->surface_deriv_start_index() + i) +=
                internal_grad.d_albedo;
        }

        if constexpr (has_thermal<Source>()) {
            // Thermal surface derivatives
            source.deriv(
                input.atmosphere->surface_emission_deriv_start_index()) +=
                internal_grad.d_thermal_surf;
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
    template <typename Derived, SourceType Source>
    inline void od(const Input<Source>& input,
                   const Eigen::MatrixBase<Derived>& d_od,
                   GradientMap<Source>& grad) {
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
    template <typename Derived, SourceType Source>
    inline void ssa(const Input<Source>& input,
                    const Eigen::MatrixBase<Derived>& d_ssa,
                    GradientMap<Source>& grad) {
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
    template <typename Derived, SourceType Source>
    inline void b1(const Input<Source>& input,
                   const Eigen::MatrixBase<Derived>& d_b1,
                   GradientMap<Source>& grad) {

        grad.d_b1 += d_b1;
    }

    template <typename Derived, SourceType Source>
    inline void transmission(const Input<Source>& input,
                             const Eigen::MatrixBase<Derived>& d_transmission,
                             GradientMap<Source>& grad) {
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

    template <typename Derived, SourceType Source>
    inline void secant(const Input<Source>& input,
                       const Eigen::MatrixBase<Derived>& d_secant,
                       GradientMap<Source>& grad) {
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

    template <typename Derived, SourceType Source>
    inline void thermal_b0(const Input<Source>& input,
                           const Eigen::MatrixBase<Derived>& d_thermal_b0,
                           GradientMap<Source>& grad) {
        // b0 doesn't depend on OD, so just add the gradient to the thermal_b0
        // gradient
        grad.d_thermal_b0 += d_thermal_b0;
    }

    template <typename Derived, SourceType Source>
    inline void thermal_b1(const Input<Source>& input,
                           const Eigen::MatrixBase<Derived>& d_thermal_b1,
                           GradientMap<Source>& grad) {
        // In terms of layer quantities we don't have to map b1 either, when we
        // move to atmosphere mapping then we have to consider the contribution
        // from od?
        grad.d_thermal_b1 += d_thermal_b1;
    }

    template <SourceType Source>
    inline void
    homog_X(const Input<Source>& input,
            const std::array<HomogSolution<Source>, num_azimuth<Source>()>&
                solution,
            std::array<Eigen::RowVectorXd, num_azimuth<Source>()>& d_X_plus,
            std::array<Eigen::RowVectorXd, num_azimuth<Source>()>& d_X_minus,
            GradientMap<Source>& grad) {
        // We also have
        // X_plus = 0.5 * (1 - s / k)
        // X_minus = 0.5 * (1 + s / k)

        // The gradient wrt to s is  -/+ 0.5 / k

        // The gradient wrt to k is +/- 0.5 s / k^2, but since k^2 = s*d, this
        // is +/- 0.5 / d

        for (int i = 0; i < num_azimuth<Source>(); ++i) {
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
    template <SourceType Source>
    inline void
    homog_k(const Input<Source>& input,
            const std::array<HomogSolution<Source>, num_azimuth<Source>()>&
                solution,
            std::array<Eigen::RowVectorXd, num_azimuth<Source>()>& d_k,
            GradientMap<Source>& grad) {
        // To start, we have k = sqrt(s * d) and so then
        // The gradient wrt to s is 0.5 / k * d
        // The gradient wrt to d is 0.5 / k * s

        for (int i = 0; i < num_azimuth<Source>(); ++i) {
            // This is now grad wrt to ssa
            ssa(input, d_k[i].cwiseProduct(solution[i].k.d_ssa.transpose()),
                grad);
        }
    }

    template <SourceType Source>
    inline void
    homog_omega(const Input<Source>& input,
                const std::array<HomogSolution<Source>, num_azimuth<Source>()>&
                    solution,
                std::array<Eigen::RowVectorXd, num_azimuth<Source>()>& d_omega,
                GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_omega[i].cwiseProduct(solution[i].omega.d_ssa.transpose()),
                grad);
        }
    }

    template <SourceType Source>
    inline void particular_Q_plus(
        const Input<Source>& input,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        std::array<Eigen::RowVectorXd, num_azimuth<Source>()>& d_Q_plus,
        GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_Q_plus[i].cwiseProduct(solution[i].Q_plus.d_ssa.transpose()),
                grad);
        }
    }

    template <SourceType Source>
    inline void particular_Q_minus(
        const Input<Source>& input,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        std::array<Eigen::RowVectorXd, num_azimuth<Source>()>& d_Q_minus,
        GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_Q_minus[i].cwiseProduct(
                    solution[i].Q_minus.d_ssa.transpose()),
                grad);
        }
    }

    template <SourceType Source>
    inline void particular_norm(
        const Input<Source>& input,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        std::array<Eigen::RowVectorXd, num_azimuth<Source>()>& d_norm,
        GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_norm[i].cwiseProduct(solution[i].norm.d_ssa.transpose()),
                grad);
        }
    }

    template <SourceType Source>
    inline void particular_Aplus(
        const Input<Source>& input,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        std::array<Eigen::RowVectorXd, num_azimuth<Source>()>& d_Aplus,
        GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_Aplus[i].cwiseProduct(solution[i].A_plus.d_ssa.transpose()),
                grad);
        }
    }

    template <SourceType Source>
    inline void
    particular_Aminus(const Input<Source>& input,
                      const std::array<ParticularSolution<Source>,
                                       num_azimuth<Source>()>& solution,
                      std::array<Eigen::RowVectorXd, 2>& d_Aminus,
                      GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
            // This is now grad wrt to ssa
            ssa(input,
                d_Aminus[i].cwiseProduct(solution[i].A_minus.d_ssa.transpose()),
                grad);
        }
    }

    template <SourceType Source>
    inline void particular_G_plus_top(
        const Input<Source>& input,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        const std::array<Eigen::RowVectorXd, num_azimuth<Source>()>&
            d_G_plus_top,
        GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
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

            if constexpr (has_solar<Source>()) {
                transmission(input,
                             d_G_plus_top[i].cwiseProduct(
                                 solution[i]
                                     .G_plus_top
                                     .d_transmission(Eigen::seq(
                                         0, Eigen::placeholders::last - 1))
                                     .transpose()),
                             grad);
                secant(input,
                       d_G_plus_top[i].cwiseProduct(
                           solution[i].G_plus_top.d_secant.transpose()),
                       grad);
            }

            if constexpr (has_thermal<Source>()) {
                thermal_b0(input,
                           d_G_plus_top[i].cwiseProduct(
                               solution[i].G_plus_top.d_thermal_b0.transpose()),
                           grad);
                thermal_b1(input,
                           d_G_plus_top[i].cwiseProduct(
                               solution[i].G_plus_top.d_thermal_b1.transpose()),
                           grad);
            }
        }
    }

    template <SourceType Source>
    inline void particular_G_plus_bottom(
        const Input<Source>& input,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        const std::array<Eigen::RowVectorXd, num_azimuth<Source>()>&
            d_G_plus_bottom,
        GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
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

            if constexpr (has_solar<Source>()) {
                transmission(input,
                             d_G_plus_bottom[i].cwiseProduct(
                                 solution[i]
                                     .G_plus_bottom
                                     .d_transmission(Eigen::seq(
                                         0, Eigen::placeholders::last - 1))
                                     .transpose()),
                             grad);
                secant(input,
                       d_G_plus_bottom[i].cwiseProduct(
                           solution[i].G_plus_bottom.d_secant.transpose()),
                       grad);
            }

            if constexpr (has_thermal<Source>()) {
                thermal_b0(
                    input,
                    d_G_plus_bottom[i].cwiseProduct(
                        solution[i].G_plus_bottom.d_thermal_b0.transpose()),
                    grad);
                thermal_b1(
                    input,
                    d_G_plus_bottom[i].cwiseProduct(
                        solution[i].G_plus_bottom.d_thermal_b1.transpose()),
                    grad);
            }
        }
    }

    template <SourceType Source>
    inline void particular_G_minus_top(
        const Input<Source>& input,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        const std::array<Eigen::RowVectorXd, num_azimuth<Source>()>&
            d_G_minus_top,
        GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
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

            if constexpr (has_solar<Source>()) {
                transmission(input,
                             d_G_minus_top[i].cwiseProduct(
                                 solution[i]
                                     .G_minus_top
                                     .d_transmission(Eigen::seq(
                                         0, Eigen::placeholders::last - 1))
                                     .transpose()),
                             grad);
                secant(input,
                       d_G_minus_top[i].cwiseProduct(
                           solution[i].G_minus_top.d_secant.transpose()),
                       grad);
            }

            if constexpr (has_thermal<Source>()) {
                thermal_b0(
                    input,
                    d_G_minus_top[i].cwiseProduct(
                        solution[i].G_minus_top.d_thermal_b0.transpose()),
                    grad);
                thermal_b1(
                    input,
                    d_G_minus_top[i].cwiseProduct(
                        solution[i].G_minus_top.d_thermal_b1.transpose()),
                    grad);
            }
        }
    }

    template <SourceType Source>
    inline void particular_G_minus_bottom(
        const Input<Source>& input,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        const std::array<Eigen::RowVectorXd, num_azimuth<Source>()>&
            d_G_minus_bottom,
        GradientMap<Source>& grad) {
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
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

            if constexpr (has_solar<Source>()) {
                transmission(input,
                             d_G_minus_bottom[i].cwiseProduct(
                                 solution[i]
                                     .G_minus_bottom
                                     .d_transmission(Eigen::seq(
                                         0, Eigen::placeholders::last - 1))
                                     .transpose()),
                             grad);
                secant(input,
                       d_G_minus_bottom[i].cwiseProduct(
                           solution[i].G_minus_bottom.d_secant.transpose()),
                       grad);
            }

            if constexpr (has_thermal<Source>()) {
                thermal_b0(
                    input,
                    d_G_minus_bottom[i].cwiseProduct(
                        solution[i].G_minus_bottom.d_thermal_b0.transpose()),
                    grad);
                thermal_b1(
                    input,
                    d_G_minus_bottom[i].cwiseProduct(
                        solution[i].G_minus_bottom.d_thermal_b1.transpose()),
                    grad);
            }
        }
    }

    template <SourceType Source>
    inline void
    bvp(const Input<Source>& input,
        std::array<BVPCoeffs<Source>, num_azimuth<Source>()>& bvp_coeffs,
        const std::array<HomogSolution<Source>, num_azimuth<Source>()>& homog,
        const std::array<ParticularSolution<Source>, num_azimuth<Source>()>&
            solution,
        std::array<Eigen::MatrixXd, num_azimuth<Source>()>& d_coeffs,
        GradientMap<Source>& grad) {
        ZoneScopedN("Twostream Backprop BVP");
        for (int i = 0; i < num_azimuth<Source>(); ++i) {
            // We have to solve the adjoint equation to get the gradient
            pentadiagonal_solve<Source, true>(bvp_coeffs[i], d_coeffs[i]);
            // Now we have to backprop these through the boundary value
            // problem We have 6 BVP vectors that we have to account for (e,
            // c, d, a, b) for the system and RHS for the solution

            // Let's start with the RHS, we have contributions to all four G
            // terms
            bvp_coeffs[i].d_G_plus_top(0) = -d_coeffs[i](0, 0);

            bvp_coeffs[i].d_G_plus_top(Eigen::seq(1, input.nlyr - 1)).array() =
                d_coeffs[i](Eigen::seq(2, Eigen::placeholders::last, 2), 0)
                    .transpose()
                    .array();

            bvp_coeffs[i]
                .d_G_plus_bottom(Eigen::seq(0, input.nlyr - 2))
                .array() =
                -d_coeffs[i](Eigen::seq(2, Eigen::placeholders::last, 2), 0)
                     .transpose()
                     .array();

            bvp_coeffs[i].d_G_minus_top(Eigen::seq(1, input.nlyr - 1)).array() =
                d_coeffs[i](Eigen::seq(1, Eigen::placeholders::last - 1, 2), 0)
                    .transpose()
                    .array();
            bvp_coeffs[i]
                .d_G_minus_bottom(Eigen::seq(0, input.nlyr - 2))
                .array() =
                -d_coeffs[i](Eigen::seq(1, Eigen::placeholders::last - 1, 2), 0)
                     .transpose()
                     .array();

            if constexpr (has_solar<Source>()) {
                bvp_coeffs[i].d_G_minus_bottom(input.nlyr - 1) =
                    -d_coeffs[i](Eigen::placeholders::last, 0) * (1);

                bvp_coeffs[i].d_G_plus_bottom(input.nlyr - 1) =
                    -d_coeffs[i](Eigen::placeholders::last, 0) *
                    (-2 * input.mu * input.albedo *
                     sasktran_disco::kronDelta(i, 0));
            } else {
                bvp_coeffs[i].d_G_minus_bottom(input.nlyr - 1) = 0.0;
                bvp_coeffs[i].d_G_plus_bottom(input.nlyr - 1) = 0.0;
            }

            // Backprop the SSA factors
            ssa(input,
                bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                    solution[i].G_plus_top.d_ssa.transpose()),
                grad);
            ssa(input,
                bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                    solution[i].G_plus_bottom.d_ssa.transpose()),
                grad);
            ssa(input,
                bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                    solution[i].G_minus_top.d_ssa.transpose()),
                grad);
            ssa(input,
                bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                    solution[i].G_minus_bottom.d_ssa.transpose()),
                grad);

            // b1 factors
            b1(input,
               bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                   solution[i].G_plus_top.d_b1.transpose()),
               grad);
            b1(input,
               bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                   solution[i].G_plus_bottom.d_b1.transpose()),
               grad);
            b1(input,
               bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                   solution[i].G_minus_top.d_b1.transpose()),
               grad);
            b1(input,
               bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                   solution[i].G_minus_bottom.d_b1.transpose()),
               grad);

            // And the OD factors
            od(input,
               bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                   solution[i].G_plus_top.d_od.transpose()),
               grad);
            od(input,
               bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                   solution[i].G_plus_bottom.d_od.transpose()),
               grad);
            od(input,
               bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                   solution[i].G_minus_top.d_od.transpose()),
               grad);
            od(input,
               bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                   solution[i].G_minus_bottom.d_od.transpose()),
               grad);

            // And the transmission factors

            if constexpr (has_solar<Source>()) {
                bvp_coeffs[0].d_temp_transmission(
                    Eigen::seq(0, input.nlyr - 1)) +=
                    bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                        solution[i]
                            .G_plus_top
                            .d_transmission(
                                Eigen::seq(0, Eigen::placeholders::last - 1))
                            .transpose()) +
                    bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                        solution[i]
                            .G_plus_bottom
                            .d_transmission(
                                Eigen::seq(0, Eigen::placeholders::last - 1))
                            .transpose()) +
                    bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                        solution[i]
                            .G_minus_top
                            .d_transmission(
                                Eigen::seq(0, Eigen::placeholders::last - 1))
                            .transpose()) +
                    bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                        solution[i]
                            .G_minus_bottom
                            .d_transmission(
                                Eigen::seq(0, Eigen::placeholders::last - 1))
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
            }

            if constexpr (has_thermal<Source>()) {
                bvp_coeffs[0].d_temp_thermal_b0 +=
                    bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                        solution[i].G_plus_top.d_thermal_b0.transpose()) +
                    bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                        solution[i].G_plus_bottom.d_thermal_b0.transpose()) +
                    bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                        solution[i].G_minus_top.d_thermal_b0.transpose()) +
                    bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                        solution[i].G_minus_bottom.d_thermal_b0.transpose());

                // And finally the secant factors
                bvp_coeffs[0].d_temp_thermal_b1 +=
                    bvp_coeffs[i].d_G_plus_top.cwiseProduct(
                        solution[i].G_plus_top.d_thermal_b1.transpose()) +
                    bvp_coeffs[i].d_G_plus_bottom.cwiseProduct(
                        solution[i].G_plus_bottom.d_thermal_b1.transpose()) +
                    bvp_coeffs[i].d_G_minus_top.cwiseProduct(
                        solution[i].G_minus_top.d_thermal_b1.transpose()) +
                    bvp_coeffs[i].d_G_minus_bottom.cwiseProduct(
                        solution[i].G_minus_bottom.d_thermal_b1.transpose());
            }

            if constexpr (has_solar<Source>()) {
                // Directly assign the albedo derivative
                grad.d_albedo +=
                    d_coeffs[i](2 * input.nlyr - 1, 0) *
                    (input.csz / EIGEN_PI *
                         input.transmission(Eigen::placeholders::last) +
                     2 * input.mu *
                         solution[i].G_plus_bottom.value(
                             Eigen::placeholders::last));
            }
            // Now for the hard part, the pentadiagonal system (e, c, d, a,
            // b) The derivative term is -(dA)z, where A is our
            // pentadiagonal matrix and z is our BVP solution The
            // derivatives of A are given by bvp_coeffs[i].d_d_by_ssa etc.

            /** SSA Derivatives */
            bvp_coeffs[i].d_temp_ssa.setZero();
            // TOA terms
            bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(0, 0) *
                                           d_coeffs[i](0, 0) *
                                           bvp_coeffs[i].d_d_by_ssa(0, 0);
            bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(1, 0) *
                                           d_coeffs[i](0, 0) *
                                           bvp_coeffs[i].d_a_by_ssa(0, 0);
            // Continuity terms
            for (int k = 0; k < input.nlyr - 1; ++k) {
                // d factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs(2 * k + 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_d_by_ssa(2 * k + 1, k);

                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs(2 * k + 2, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_d_by_ssa(2 * k + 2, k + 1);

                // a factors
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) + 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_a_by_ssa(2 * k + 1, k + 1);
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) + 1, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_a_by_ssa(2 * k + 2, k + 1);

                // b factors
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) + 2, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_b_by_ssa(2 * k + 1, k + 1);

                // c factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) - 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_c_by_ssa(2 * k + 1, k);
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) - 1, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_c_by_ssa(2 * k + 2, k);

                // e factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) - 2, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_e_by_ssa(2 * k + 2, k);
            }

            // Ground term
            bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                -bvp_coeffs[i].rhs(2 * input.nlyr - 1, 0) *
                d_coeffs[i](2 * input.nlyr - 1, 0) *
                bvp_coeffs[i].d_d_by_ssa(2 * input.nlyr - 1, input.nlyr - 1);
            bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                -bvp_coeffs[i].rhs(2 * input.nlyr - 1 - 1, 0) *
                d_coeffs[i](2 * input.nlyr - 1, 0) *
                bvp_coeffs[i].d_c_by_ssa(2 * input.nlyr - 1, input.nlyr - 1);

            ssa(input, bvp_coeffs[i].d_temp_ssa, grad);

            /** OD Derivatives */
            bvp_coeffs[i].d_temp_ssa.setZero();
            // TOA terms
            bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(0, 0) *
                                           d_coeffs[i](0, 0) *
                                           bvp_coeffs[i].d_d_by_od(0, 0);
            bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(1, 0) *
                                           d_coeffs[i](0, 0) *
                                           bvp_coeffs[i].d_a_by_od(0, 0);
            // Continuity terms
            for (int k = 0; k < input.nlyr - 1; ++k) {
                // d factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs(2 * k + 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_d_by_od(2 * k + 1, k);

                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs(2 * k + 2, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_d_by_od(2 * k + 2, k + 1);

                // a factors
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) + 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_a_by_od(2 * k + 1, k + 1);
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) + 1, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_a_by_od(2 * k + 2, k + 1);

                // b factors
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) + 2, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_b_by_od(2 * k + 1, k + 1);

                // c factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) - 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_c_by_od(2 * k + 1, k);
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) - 1, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_c_by_od(2 * k + 2, k);

                // e factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) - 2, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_e_by_od(2 * k + 2, k);
            }

            // Ground term
            bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                -bvp_coeffs[i].rhs(2 * input.nlyr - 1, 0) *
                d_coeffs[i](2 * input.nlyr - 1, 0) *
                bvp_coeffs[i].d_d_by_od(2 * input.nlyr - 1, input.nlyr - 1);
            bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                -bvp_coeffs[i].rhs(2 * input.nlyr - 1 - 1, 0) *
                d_coeffs[i](2 * input.nlyr - 1, 0) *
                bvp_coeffs[i].d_c_by_od(2 * input.nlyr - 1, input.nlyr - 1);

            od(input, bvp_coeffs[i].d_temp_ssa, grad);

            /** b1 derivatives */
            bvp_coeffs[i].d_temp_ssa.setZero();
            // TOA terms
            bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(0, 0) *
                                           d_coeffs[i](0, 0) *
                                           bvp_coeffs[i].d_d_by_b1(0, 0);
            bvp_coeffs[i].d_temp_ssa(0) += -bvp_coeffs[i].rhs(1, 0) *
                                           d_coeffs[i](0, 0) *
                                           bvp_coeffs[i].d_a_by_b1(0, 0);
            // Continuity terms
            for (int k = 0; k < input.nlyr - 1; ++k) {
                // d factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs(2 * k + 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_d_by_b1(2 * k + 1, k);

                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs(2 * k + 2, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_d_by_b1(2 * k + 2, k + 1);

                // a factors
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) + 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_a_by_b1(2 * k + 1, k + 1);
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) + 1, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_a_by_b1(2 * k + 2, k + 1);

                // b factors
                bvp_coeffs[i].d_temp_ssa(k + 1) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) + 2, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_b_by_b1(2 * k + 1, k + 1);

                // c factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 1) - 1, 0) *
                    d_coeffs[i](2 * k + 1, 0) *
                    bvp_coeffs[i].d_c_by_b1(2 * k + 1, k);
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) - 1, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_c_by_b1(2 * k + 2, k);

                // e factors
                bvp_coeffs[i].d_temp_ssa(k) +=
                    -bvp_coeffs[i].rhs((2 * k + 2) - 2, 0) *
                    d_coeffs[i](2 * k + 2, 0) *
                    bvp_coeffs[i].d_e_by_b1(2 * k + 2, k);
            }

            // Ground term
            bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                -bvp_coeffs[i].rhs(2 * input.nlyr - 1, 0) *
                d_coeffs[i](2 * input.nlyr - 1, 0) *
                bvp_coeffs[i].d_d_by_b1(2 * input.nlyr - 1, input.nlyr - 1);
            bvp_coeffs[i].d_temp_ssa(input.nlyr - 1) +=
                -bvp_coeffs[i].rhs(2 * input.nlyr - 1 - 1, 0) *
                d_coeffs[i](2 * input.nlyr - 1, 0) *
                bvp_coeffs[i].d_c_by_b1(2 * input.nlyr - 1, input.nlyr - 1);

            b1(input, bvp_coeffs[i].d_temp_ssa, grad);

            // BVP matrix albedo derivatives
            // c and d have albedo derivatives

            // for c
            double d_albedo = -2 * input.mu *
                              homog[i].X_plus.value(Eigen::placeholders::last) *
                              homog[i].omega.value(Eigen::placeholders::last) *
                              sasktran_disco::kronDelta(i, 0);
            grad.d_albedo += -bvp_coeffs[i].rhs(2 * input.nlyr - 1 - 1, 0) *
                             d_coeffs[i](2 * input.nlyr - 1, 0) * d_albedo;
            // for d
            d_albedo = -2 * input.mu *
                       homog[i].X_minus.value(Eigen::placeholders::last) *
                       sasktran_disco::kronDelta(i, 0);
            grad.d_albedo += -bvp_coeffs[i].rhs(2 * input.nlyr - 1, 0) *
                             d_coeffs[i](2 * input.nlyr - 1, 0) * d_albedo;
        }
        if constexpr (has_solar<Source>()) {
            transmission(input, bvp_coeffs[0].d_temp_transmission, grad);
            secant(input, bvp_coeffs[0].d_temp_secant, grad);
        }

        if constexpr (has_thermal<Source>()) {
            thermal_b0(input, bvp_coeffs[0].d_temp_thermal_b0, grad);
            thermal_b1(input, bvp_coeffs[0].d_temp_thermal_b1, grad);
        }
    }

    template <SourceType Source>
    inline void
    full(const Input<Source>& input, Solution<Source>& solution,
         const Sources<Source>& sources,
         const Eigen::RowVectorXd& source_weights,
         std::array<Eigen::MatrixXd, num_azimuth<Source>()>& d_coeffs,
         GradientMap<Source>& grad, double ground_weight = 0.0) {
        ZoneScopedN("Twostream full backprop");
        // Direct source terms
        ssa(input,
            sources.source.d_ssa.transpose().cwiseProduct(source_weights),
            grad);
        od(input, sources.source.d_od.transpose().cwiseProduct(source_weights),
           grad);
        b1(input, sources.source.d_b1.transpose().cwiseProduct(source_weights),
           grad);

        if constexpr (has_solar<Source>()) {
            solution.bvp_coeffs[0].d_temp_secant =
                sources.source.d_secant.transpose().cwiseProduct(
                    source_weights);

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
        }

        if constexpr (has_thermal<Source>()) {
            solution.bvp_coeffs[0].d_temp_thermal_b0 =
                sources.source.d_thermal_b0.transpose().cwiseProduct(
                    source_weights);
            solution.bvp_coeffs[0].d_temp_thermal_b1 =
                sources.source.d_thermal_b1.transpose().cwiseProduct(
                    source_weights);
        }

        // Direct ground multiple scatter contributions

        grad.d_ssa(input.nlyr - 1) +=
            (solution.particular[0].G_plus_bottom.d_ssa(input.nlyr - 1) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                 (solution.homog[0].X_plus.d_ssa(Eigen::placeholders::last) *
                      solution.homog[0].omega.value(Eigen::placeholders::last) +
                  solution.homog[0].X_plus.value(Eigen::placeholders::last) *
                      solution.homog[0].omega.d_ssa(
                          Eigen::placeholders::last)) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
                 (solution.homog[0].X_minus.d_ssa(Eigen::placeholders::last))) *
            ground_weight;

        grad.d_b1(input.nlyr - 1) +=
            (solution.particular[0].G_plus_bottom.d_b1(input.nlyr - 1) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                 (solution.homog[0].X_plus.d_b1(Eigen::placeholders::last) *
                      solution.homog[0].omega.value(Eigen::placeholders::last) +
                  solution.homog[0].X_plus.value(Eigen::placeholders::last) *
                      solution.homog[0].omega.d_b1(Eigen::placeholders::last)) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
                 (solution.homog[0].X_minus.d_b1(Eigen::placeholders::last))) *
            ground_weight;

        grad.d_extinction(input.nlyr - 1) +=
            (solution.particular[0].G_plus_bottom.d_od(input.nlyr - 1) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                 solution.homog[0].X_plus.value(Eigen::placeholders::last) *
                 solution.homog[0].omega.d_od(Eigen::placeholders::last)) *
            ground_weight;

        // BVP adjoint solution
        d_coeffs[0].col(0) = sources.d_bvp_coeff[0];

        if constexpr (num_azimuth<Source>() > 1) {
            d_coeffs[1].col(0) = sources.d_bvp_coeff[1];
            d_coeffs[1](Eigen::seq(0, 2 * input.nlyr - 1, 2), 0).array() *=
                source_weights.transpose().array();
            d_coeffs[1](Eigen::seq(1, 2 * input.nlyr - 1, 2), 0).array() *=
                source_weights.transpose().array();
        }

        d_coeffs[0](Eigen::seq(0, 2 * input.nlyr - 1, 2), 0).array() *=
            source_weights.transpose().array();
        d_coeffs[0](Eigen::seq(1, 2 * input.nlyr - 1, 2), 0).array() *=
            source_weights.transpose().array();

        // Albedo contributions to BVP coefficients

        d_coeffs[0](2 * input.nlyr - 1, 0) +=
            ground_weight *
            (solution.homog[0].X_minus.value(Eigen::placeholders::last));

        d_coeffs[0](2 * input.nlyr - 2, 0) +=
            ground_weight *
            (solution.homog[0].X_plus.value(Eigen::placeholders::last) *
             solution.homog[0].omega.value(Eigen::placeholders::last));

        bvp(input, solution.bvp_coeffs, solution.homog, solution.particular,
            d_coeffs, grad);
    }
} // namespace sasktran2::twostream::backprop
