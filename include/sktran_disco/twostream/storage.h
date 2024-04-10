#pragma once
#include "sasktran2/atmosphere/atmosphere.h"
#include "sktran_disco/sktran_do_geometrylayerarray.h"
#include <sasktran2/internal_common.h>

namespace sasktran2::twostream {
    /**
     *  Stores all of the input parameters necessary for the two stream model,
     * as well as some derived cached values
     *
     */
    struct Input {
        const sasktran_disco::GeometryLayerArray<1>*
            geometry_layers; /** Pointer to the geometry layer information */
        const sasktran2::atmosphere::Atmosphere<1>*
            atmosphere; /** Pointer to the sasktran atmosphere */

        Eigen::VectorXd od;  /** Optical depth in each [layer]*/
        Eigen::VectorXd ssa; /** Single scatter albedo in each [layer]*/
        Eigen::VectorXd b1;  /** First legendre coefficient in each [layer]*/

        Eigen::VectorXd
            transmission; /** solar transmission to the top of each [layer+1],
                             last element is transmission to the surface*/
        Eigen::VectorXd average_secant; /** Average secant in each [layer] */
        Eigen::VectorXd expsec;         /** exp(-od * secant) in each [layer] */

        double csz;    /** Cosine of the solar zenith angle*/
        double mu;     /** Quadrature angle */
        double albedo; /** Surface albedo */

        int nlyr;     /** Number of layers */
        int wavelidx; /** Wavelength index */

        /** Initializes the storage for a given number of layers*/
        void init(int nlyr) {
            this->nlyr = nlyr;
            od.resize(nlyr);
            ssa.resize(nlyr);
            b1.resize(nlyr);
            transmission.resize(nlyr + 1);
            average_secant.resize(nlyr);

            albedo = 0.0;

            mu = 0.5;
        }

        /**
         * Sets the input parameters based on the wavelength index from the
         * sasktran2 atmosphere object
         *
         * @param wavelidx
         */
        void calculate(int wavelidx) {
            this->wavelidx = wavelidx;

            // Start by interpolating extinction to the layers
            od.array() =
                (geometry_layers->interpolating_matrix() *
                 (atmosphere->storage().total_extinction.col(wavelidx)))
                    .array();

            // And interpolate the scattering extinction
            ssa = (geometry_layers->interpolating_matrix() *
                   (atmosphere->storage().ssa.col(wavelidx).cwiseProduct(
                       atmosphere->storage().total_extinction.col(wavelidx))));

            // We also need to interpolate the first legendre coefficient
            // multiplied by the scattering extinction
            int stride = atmosphere->storage().leg_coeff.dimension(0);

            Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> grid_b1(
                &atmosphere->storage().leg_coeff(1, 0, wavelidx), nlyr,
                Eigen::InnerStride<>(stride));

            /*
            b1 = (geometry_layers->interpolating_matrix() *
                  (atmosphere->storage().ssa.col(wavelidx).cwiseProduct(
                       atmosphere->storage().total_extinction.col(wavelidx)))
                      .cwiseProduct(grid_b1));
            */

            // Then we can divide by the total scattering extinction to get b1
            b1.array() /= ssa.array();
            b1.setConstant(0.1);

            // Then ssa will be the scattering extinction divided by the total
            // extinction
            ssa.array() /= od.array();

            // And we multiply by the thickness to get the optical depth
            od.array() *= (geometry_layers->layer_ceiling().array() -
                           geometry_layers->layer_floor().array());

            transmission(0) = 0;
            transmission(Eigen::seq(1, nlyr)) =
                -1 * (geometry_layers->chapman_factors() * od);

            average_secant.array() =
                (transmission(Eigen::seq(0, nlyr - 1)).array() -
                 transmission(Eigen::seq(1, nlyr)).array()) /
                od.array();

            transmission.array() = transmission.array().exp();

            expsec = (-1.0 * average_secant.array() * od.array()).exp();
        }
    };

    /**
     *  Storage for the homogeneous solution and cached quantities
     *
     */
    struct HomogSolution {
        Eigen::VectorXd k;       /** Eigen value, [lyr] */
        Eigen::VectorXd X_plus;  /** Position eigen vector [lyr] */
        Eigen::VectorXd X_minus; /** Negative eigen vector [wlyr] */

        // Temporaries
        Eigen::VectorXd d;     /** Temporary [lyr] */
        Eigen::VectorXd s;     /** Temporary [lyr] */
        Eigen::VectorXd omega; /** exp(-k*od) [lyr]*/

        // Gradients
        Eigen::RowVectorXd
            d_d_by_ssa; /** Gradient of d with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_s_by_ssa; /** Gradient of s with respect to ssa [lyr] */

        Eigen::RowVectorXd
            d_k_by_ssa; /** Gradient of k with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_X_plus_by_ssa; /** Gradient of X_plus with respect to ssa [lyr] */
        Eigen::RowVectorXd d_X_minus_by_ssa; /** Gradient of X_minus with
                                                respect to ssa [lyr] */

        Eigen::RowVectorXd
            d_omega_by_ssa; /** Gradient of omega with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_omega_by_od; /** Gradient of omega with respect to od [lyr] */

        /**
         *  Initializes the storage for a given number of layers
         *
         * @param nlyr
         */
        void init(int nlyr) {
            k.resize(nlyr);
            X_plus.resize(nlyr);
            X_minus.resize(nlyr);

            d.resize(nlyr);
            s.resize(nlyr);

            omega.resize(nlyr);

            d_d_by_ssa.resize(nlyr);
            d_s_by_ssa.resize(nlyr);
            d_k_by_ssa.resize(nlyr);

            d_X_minus_by_ssa.resize(nlyr);
            d_X_plus_by_ssa.resize(nlyr);

            d_omega_by_od.resize(nlyr);
            d_omega_by_ssa.resize(nlyr);
        }
    };

    /**
     *  Storage for the particular solution and cached quantities
     *
     */
    struct ParticularSolution {
        Eigen::VectorXd G_plus_top;     /** G+ at top of the [lyr] */
        Eigen::VectorXd G_plus_bottom;  /** G+ at bottom of the [lyr]*/
        Eigen::VectorXd G_minus_top;    /** G- at the top of the [lyr] */
        Eigen::VectorXd G_minus_bottom; /** G- at the bottom of the [lyr] */
        Eigen::VectorXd A_plus;         /** A+ in the layer [lyr] */
        Eigen::VectorXd A_minus;        /** A- in the [lyr] */

        // Temporaries
        Eigen::VectorXd Q_plus;  /** Q+ in the [lyr] */
        Eigen::VectorXd Q_minus; /** Q- in the [lyr] */
        Eigen::VectorXd norm;    /** Normalization factor [lyr] */
        Eigen::VectorXd C_plus;  /** C+ [lyr]*/
        Eigen::VectorXd C_minus; /** C- [lyr]*/

        // Gradients
        Eigen::RowVectorXd
            d_Q_plus_by_ssa; /** Gradient of Q+ with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_Q_minus_by_ssa; /** Gradient of Q- with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_norm_by_ssa; /** Gradient of norm with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_A_plus_by_ssa; /** Gradient of A+ with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_A_minus_by_ssa; /** Gradient of A- with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_C_plus_by_ssa; /** Gradient of C+ with respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_C_minus_by_ssa; /** Gradient of C- with respect to ssa [lyr] */

        Eigen::RowVectorXd d_G_plus_top_by_ssa; /** Gradient of G+ at top with
                                                   respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_G_plus_bottom_by_ssa; /** Gradient of G+ at bottom with respect to
                                       ssa [lyr] */
        Eigen::RowVectorXd d_G_minus_top_by_ssa; /** Gradient of G- at top with
                                                    respect to ssa [lyr] */
        Eigen::RowVectorXd
            d_G_minus_bottom_by_ssa; /** Gradient of G- at bottom with respect
                                        to ssa [lyr] */

        /**
         *  Initializes the storage for a given number of layers
         *
         * @param nlyr
         */
        void init(int nlyr) {
            G_plus_top.resize(nlyr);
            G_plus_bottom.resize(nlyr);
            G_minus_top.resize(nlyr);
            G_minus_bottom.resize(nlyr);
            A_plus.resize(nlyr);
            A_minus.resize(nlyr);

            Q_plus.resize(nlyr);
            Q_minus.resize(nlyr);
            norm.resize(nlyr);
            C_plus.resize(nlyr);
            C_minus.resize(nlyr);

            d_Q_minus_by_ssa.resize(nlyr);
            d_Q_plus_by_ssa.resize(nlyr);
            d_norm_by_ssa.resize(nlyr);
            d_A_plus_by_ssa.resize(nlyr);
            d_A_minus_by_ssa.resize(nlyr);
            d_C_minus_by_ssa.resize(nlyr);
            d_C_plus_by_ssa.resize(nlyr);

            d_G_plus_bottom_by_ssa.resize(nlyr);
            d_G_plus_top_by_ssa.resize(nlyr);
            d_G_minus_bottom_by_ssa.resize(nlyr);
            d_G_minus_top_by_ssa.resize(nlyr);
        }
    };

    /**
     *  Storage for the boundary value problem coefficients and cached
     * quantities
     *
     */
    struct BVPCoeffs {
        Eigen::VectorXd e, c, d, b, a; /** Pentadiagonal matrix diagonals*/

        Eigen::VectorXd gamma, mu, alpha,
            beta; /** Temporary storage for the pentadiagonal solver */

        Eigen::MatrixXd z; /** Temporary storage for the pentadiagonal solver */
        Eigen::MatrixXd rhs; /** RHS of the pentadiagonal system */

        // Gradient mappings from the RHS to the BVP coefficients
        Eigen::RowVectorXd d_G_plus_top;
        Eigen::RowVectorXd d_G_plus_bottom;
        Eigen::RowVectorXd d_G_minus_top;
        Eigen::RowVectorXd d_G_minus_bottom;

        // Gradient mappings from the adjoint solution to the diagonals
        Eigen::RowVectorXd d_e, d_c, d_d, d_b, d_a;

        // Gradient mappings of the diagonals with respect to ssa
        Eigen::MatrixXd d_e_by_ssa, d_c_by_ssa, d_d_by_ssa, d_b_by_ssa,
            d_a_by_ssa;

        // Temporary storage
        Eigen::RowVectorXd d_temp;

        void init(int nlyr) {
            z.resize(nlyr * 2, 1);
            rhs.resize(nlyr * 2, 1);

            e.resize(nlyr * 2);
            c.resize(nlyr * 2);
            d.resize(nlyr * 2);
            b.resize(nlyr * 2);
            a.resize(nlyr * 2);

            e.setZero();
            c.setZero();
            d.setZero();
            b.setZero();
            a.setZero();

            gamma.resize(nlyr * 2);
            mu.resize(nlyr * 2);
            alpha.resize(nlyr * 2);
            beta.resize(nlyr * 2);

            d_G_plus_top.resize(nlyr);
            d_G_plus_bottom.resize(nlyr);
            d_G_minus_top.resize(nlyr);
            d_G_minus_bottom.resize(nlyr);

            d_G_minus_top.setZero();
            d_G_minus_bottom.setZero();
            d_G_plus_top.setZero();
            d_G_plus_bottom.setZero();

            d_e.resize(nlyr * 2);
            d_c.resize(nlyr * 2);
            d_d.resize(nlyr * 2);
            d_b.resize(nlyr * 2);
            d_a.resize(nlyr * 2);

            d_e.setZero();
            d_c.setZero();
            d_d.setZero();
            d_b.setZero();
            d_a.setZero();

            d_temp.resize(nlyr);
            d_temp.setZero();

            d_e_by_ssa.resize(nlyr * 2, nlyr);
            d_c_by_ssa.resize(nlyr * 2, nlyr);
            d_d_by_ssa.resize(nlyr * 2, nlyr);
            d_b_by_ssa.resize(nlyr * 2, nlyr);
            d_a_by_ssa.resize(nlyr * 2, nlyr);

            d_e_by_ssa.setZero();
            d_c_by_ssa.setZero();
            d_d_by_ssa.setZero();
            d_b_by_ssa.setZero();
            d_a_by_ssa.setZero();
        }
    };

    /**
     *  Storage for the entire discrete ordinates solution at both azimuthal
     * orders
     *
     */
    struct Solution {
        std::array<HomogSolution, 2>
            homog; /** Homogeneous solutino for both azimuthal orders */
        std::array<ParticularSolution, 2>
            particular; /** Particular solution for both azimuthal orders*/
        std::array<BVPCoeffs, 2>
            bvp_coeffs; /** BVP solution for both azimuthal orders */

        /**
         * Initializes the storage for a given number of lyaers
         *
         * @param nlyr
         */
        void init(int nlyr) {
            for (auto& h : homog) {
                h.init(nlyr);
            }
            for (auto& p : particular) {
                p.init(nlyr);
            }
            for (auto& b : bvp_coeffs) {
                b.init(nlyr);
            }
        }
    };

    /**
     * Storage for the post-processed sources for a given line of sight
     *
     */
    struct Sources {
        Eigen::VectorXd source; /** Integrated sources in each [lyr] */

        Eigen::VectorXd beamtrans; /** LOS transmission factors exp(-od /
                                      viewing_zenith) for each [layer] */

        Eigen::VectorXd final_weight_factors; /** Final weight factors for each
                                                [layer] */

        std::array<Eigen::VectorXd, 2> lpsum_plus,
            lpsum_minus; /** LP triple products for upwelling and downwelling
                            [lyr] */
        std::array<Eigen::VectorXd, 2> Y_plus,
            Y_minus; /** "Interpolated" homogenous solutions in each [lyr] */

        std::array<Eigen::VectorXd, 2> H_plus,
            H_minus; /** Homogenous solution multipliers [lyr] */
        std::array<Eigen::VectorXd, 2> D_plus,
            D_minus; /** Particular solution multipliers [lyr] */

        std::array<Eigen::VectorXd, 2> V; /** Particular source [lyr] */

        Eigen::VectorXd E_minus; /** Solar multiplier [lyr] */

        // Gradients
        Eigen::RowVectorXd d_source_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_lpsum_plus_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_lpsum_minus_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_Y_plus_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_Y_minus_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_H_plus_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_H_minus_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_D_plus_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_D_minus_by_ssa;
        std::array<Eigen::RowVectorXd, 2> d_V_by_ssa;

        // Backprop factors
        std::array<Eigen::RowVectorXd, 2> d_bvp_coeff;

        /**
         * Initializes the storage for a given number of layers
         *
         * @param nlyr
         */
        void init(int nlyr) {
            source.resize(nlyr);
            beamtrans.resize(nlyr);
            final_weight_factors.resize(nlyr);

            for (auto& l : lpsum_plus) {
                l.resize(nlyr);
            }
            for (auto& l : lpsum_minus) {
                l.resize(nlyr);
            }
            for (auto& l : Y_plus) {
                l.resize(nlyr);
            }
            for (auto& l : Y_minus) {
                l.resize(nlyr);
            }
            for (auto& l : H_plus) {
                l.resize(nlyr);
            }
            for (auto& l : H_minus) {
                l.resize(nlyr);
            }
            for (auto& l : D_plus) {
                l.resize(nlyr);
            }
            for (auto& l : D_minus) {
                l.resize(nlyr);
            }
            for (auto& l : V) {
                l.resize(nlyr);
            }

            E_minus.resize(nlyr);

            for (auto& d : d_bvp_coeff) {
                d.resize(nlyr * 2);
            }

            d_source_by_ssa.resize(nlyr);
            for (auto& d : d_lpsum_plus_by_ssa) {
                d.resize(nlyr);
            }
            for (auto& d : d_lpsum_minus_by_ssa) {
                d.resize(nlyr);
            }
            for (auto& d : d_Y_plus_by_ssa) {
                d.resize(nlyr);
            }
            for (auto& d : d_Y_minus_by_ssa) {
                d.resize(nlyr);
            }
            for (auto& d : d_H_plus_by_ssa) {
                d.resize(nlyr);
            }
            for (auto& d : d_H_minus_by_ssa) {
                d.resize(nlyr);
            }
            for (auto& d : d_D_plus_by_ssa) {
                d.resize(nlyr);
            }
            for (auto& d : d_D_minus_by_ssa) {
                d.resize(nlyr);
            }
            for (auto& d : d_V_by_ssa) {
                d.resize(nlyr);
            }
        }
    };
} // namespace sasktran2::twostream
