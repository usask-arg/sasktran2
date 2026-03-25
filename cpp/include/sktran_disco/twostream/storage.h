#pragma once
#include "sasktran2/atmosphere/atmosphere.h"
#include "sktran_disco/sktran_do_geometrylayerarray.h"
#include <sasktran2/internal_common.h>
#include "meta.h"

namespace sasktran2::twostream {

    /**
     *  Storage for a quantity in the twostream model that has layer derivatives
     * with respect to (optionally) the single scatter albedo, optical depth,
     * and first legendre coefficient
     *
     */
    template <bool ssa_deriv, bool od_deriv, bool b1_deriv, bool trans_deriv,
              bool secant_deriv, bool thermal_b0_deriv, bool thermal_b1_deriv,
              SourceType source_type>
    struct LayerVector {
        Eigen::VectorXd value;
        Eigen::VectorXd d_ssa;
        Eigen::VectorXd d_od;
        Eigen::VectorXd d_b1;
        Eigen::VectorXd d_thermal_b1;
        Eigen::VectorXd d_thermal_b0;
        Eigen::VectorXd d_transmission;
        Eigen::VectorXd d_secant;

        void resize(int n) {
            value.resize(n);
            value.setZero();

            if constexpr (ssa_deriv) {
                d_ssa.resize(n);
                d_ssa.setZero();
            }
            if constexpr (od_deriv) {
                d_od.resize(n);
                d_od.setZero();
            }
            if constexpr (b1_deriv) {
                d_b1.resize(n);
                d_b1.setZero();
            }

            if constexpr (trans_deriv && has_solar<source_type>()) {
                d_transmission.resize(n + 1);
                d_transmission.setZero();
            }

            if constexpr (secant_deriv && has_solar<source_type>()) {
                d_secant.resize(n);
                d_secant.setZero();
            }

            if constexpr (thermal_b0_deriv && has_thermal<source_type>()) {
                d_thermal_b0.resize(n);
                d_thermal_b0.setZero();
            }

            if constexpr (thermal_b1_deriv && has_thermal<source_type>()) {
                d_thermal_b1.resize(n);
                d_thermal_b1.setZero();
            }
        }
    };

    /**
     *  Stores all of the input parameters necessary for the two stream model,
     * as well as some derived cached values
     *
     */
    template <SourceType source_type> struct Input {
        const sasktran_disco::GeometryLayerArray<1>*
            geometry_layers; /** Pointer to the geometry layer information */
        const sasktran2::atmosphere::Atmosphere<1>*
            atmosphere; /** Pointer to the sasktran atmosphere */

        Eigen::VectorXd od;  /** Optical depth in each [layer]*/
        Eigen::VectorXd ssa; /** Single scatter albedo in each [layer]*/
        Eigen::VectorXd b1;  /** First legendre coefficient in each [layer]*/
        Eigen::VectorXd b0_thermal; /** thermal source = b0 exp(-b1 x) */
        Eigen::VectorXd b1_thermal; /** thermal source = b0 exp(-b1 x) */

        Eigen::VectorXd
            transmission; /** solar transmission to the top of each [layer+1],
                             last element is transmission to the surface*/
        Eigen::VectorXd average_secant; /** Average secant in each [layer] */
        Eigen::VectorXd expsec;         /** exp(-od * secant) in each [layer] */
        Eigen::VectorXd
            exp_thermal; /** exp(-b1_thermal * od) in each [layer] */

        Eigen::SparseMatrix<double>
            interpolating_matrix; /** Matrix to interpolate from the
                           atmosphere grid to the geometry
                           layers */

        double csz;          /** Cosine of the solar zenith angle*/
        double mu;           /** Quadrature angle */
        double albedo;       /** Surface albedo */
        double thermal_surf; /** Surface thermal emission */

        int nlyr;     /** Number of layers */
        int wavelidx; /** Wavelength index */

        /** Initializes the storage for a given number of layers*/
        void init(int nlyr) {
            this->nlyr = nlyr;
            od.resize(nlyr);
            ssa.resize(nlyr);
            b1.resize(nlyr);

            if constexpr (has_thermal<source_type>()) {
                b0_thermal.resize(nlyr);
                b1_thermal.resize(nlyr);
                exp_thermal.resize(nlyr);
            }

            if constexpr (has_solar<source_type>()) {
                transmission.resize(nlyr + 1);
                average_secant.resize(nlyr);
                expsec.resize(nlyr);
            }

            albedo = 0.0;

            mu = 0.5;

            interpolating_matrix =
                geometry_layers->interpolating_matrix().sparseView();
        }

        void calculate_base(int wavelidx) {
            this->wavelidx = wavelidx;

            const auto& k =
                atmosphere->storage().total_extinction.col(wavelidx);
            const auto& ssa_atmo = atmosphere->storage().ssa.col(wavelidx);
            const auto& thermal_source =
                atmosphere->storage().emission_source.col(wavelidx);

            int stride = atmosphere->storage().leg_coeff.dimension(0);

            Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> grid_b1(
                &atmosphere->storage().leg_coeff(1, 0, wavelidx), nlyr + 1,
                Eigen::InnerStride<>(stride));

            // Layers start at TOA and go down
            for (int i = 0; i < nlyr; ++i) {
                const int top_atmo_idx = nlyr - i; // Top atmo idx of the layer
                const int bottom_atmo_idx = nlyr - i - 1; // Bottom atmo idx

                const double layer_thickness =
                    geometry_layers->layer_ceiling()(i) -
                    geometry_layers->layer_floor()(i);

                const double avg_k =
                    (k(top_atmo_idx) + k(bottom_atmo_idx)) / 2.0;
                const double avg_k_scat =
                    (k(top_atmo_idx) * ssa_atmo(top_atmo_idx) +
                     k(bottom_atmo_idx) * ssa_atmo(bottom_atmo_idx)) /
                    2.0;

                const double avg_b1 =
                    (k(top_atmo_idx) * ssa_atmo(top_atmo_idx) *
                         grid_b1(top_atmo_idx) +
                     k(bottom_atmo_idx) * ssa_atmo(bottom_atmo_idx) *
                         grid_b1(bottom_atmo_idx)) /
                    2.0;

                od(i) = avg_k * layer_thickness;
                ssa(i) = avg_k_scat / avg_k;
                b1(i) = avg_b1 / avg_k_scat;

                if constexpr (has_thermal<source_type>()) {
                    const double thermal_top = thermal_source(top_atmo_idx);
                    const double thermal_bottom =
                        thermal_source(bottom_atmo_idx);

                    // from top of layer to bottom, thermal = b0 exp(-b1 x)
                    // at x = 0, thermal = b0, at x = layer thickness,
                    // thermal = b0 exp(-b1 * layer thickness)
                    // => b1 = -log(thermal_bottom / thermal_top) /
                    // layer_thickness
                    b0_thermal(i) = thermal_top;
                    b1_thermal(i) = -std::log(thermal_bottom / thermal_top) /
                                    layer_thickness;
                }
            }

            ssa = ssa.cwiseMin(1 - 1e-9);
            albedo =
                atmosphere->surface().brdf(wavelidx, 0, 0, 0)(0, 0) * EIGEN_PI;

            thermal_surf = atmosphere->surface().emission()[wavelidx];
        }

        void calculate_derived(int wavelidx) {
            if constexpr (has_solar<source_type>()) {
                transmission(0) = 0;
                transmission(Eigen::seq(1, nlyr)).noalias() =
                    -1 * (geometry_layers->chapman_factors() * od);

                average_secant.array() =
                    (transmission(Eigen::seq(0, nlyr - 1)).array() -
                     transmission(Eigen::seq(1, nlyr)).array()) /
                    od.array();

                transmission.array() =
                    transmission.array().exp() *
                    atmosphere->storage().solar_irradiance(wavelidx);

                expsec.array() =
                    (-1.0 * average_secant.array() * od.array()).exp();
            }
            if constexpr (has_thermal<source_type>()) {
                exp_thermal.array() = (-b1_thermal.array() * od.array()).exp();
            }
        }

        /**
         * Sets the input parameters based on the wavelength index from the
         * sasktran2 atmosphere object
         *
         * @param wavelidx
         */
        void calculate(int wavelidx) {
            calculate_base(wavelidx);
            calculate_derived(wavelidx);
        }
    };

    /**
     *  Storage for the homogeneous solution and cached quantities
     *
     */
    template <SourceType source_type> struct HomogSolution {
        LayerVector<true, false, true, false, false, false, false, source_type>
            k; /** Eigen value, [lyr] */
        LayerVector<true, false, true, false, false, false, false, source_type>
            X_plus; /** Position eigen vector [lyr] */
        LayerVector<true, false, true, false, false, false, false, source_type>
            X_minus; /** Negative eigen vector [wlyr] */

        // Temporaries
        LayerVector<true, false, true, false, false, false, false, source_type>
            d; /** Temporary [lyr] */
        LayerVector<true, false, true, false, false, false, false, source_type>
            s; /** Temporary [lyr] */
        LayerVector<true, true, true, false, false, false, false, source_type>
            omega; /** exp(-k*od) [lyr]*/

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
        }
    };

    /**
     *  Storage for the particular solution and cached quantities
     *
     */
    template <SourceType source_type> struct ParticularSolution {
        LayerVector<true, true, true, true, true, true, true, source_type>
            G_plus_top; /** G+ at top of the [lyr] */
        LayerVector<true, true, true, true, true, true, true, source_type>
            G_plus_bottom; /** G+ at bottom of the [lyr]*/
        LayerVector<true, true, true, true, true, true, true, source_type>
            G_minus_top; /** G- at the top of the [lyr] */
        LayerVector<true, true, true, true, true, true, true, source_type>
            G_minus_bottom; /** G- at the bottom of the [lyr] */
        LayerVector<true, false, true, false, false, false, false, source_type>
            A_plus; /** A+ in the layer [lyr] */
        LayerVector<true, false, true, false, false, false, false, source_type>
            A_minus; /** A- in the [lyr] */

        LayerVector<true, false, false, false, false, false, false, source_type>
            A_thermal; /** A+ in the layer [lyr] */

        // Temporaries
        LayerVector<true, false, true, false, false, false, false, source_type>
            Q_plus; /** Q+ in the [lyr] */
        LayerVector<true, false, true, false, false, false, false, source_type>
            Q_minus; /** Q- in the [lyr] */
        LayerVector<true, false, true, false, false, false, false, source_type>
            norm; /** Normalization factor [lyr] */
        LayerVector<true, true, true, true, true, false, false, source_type>
            C_plus; /** C+ [lyr]*/
        LayerVector<true, true, true, true, true, false, false, source_type>
            C_minus; /** C- [lyr]*/

        LayerVector<true, true, true, false, false, true, true, source_type>
            C_plus_thermal; /** C+ [lyr]*/
        LayerVector<true, true, true, false, false, true, true, source_type>
            C_minus_thermal; /** C- [lyr]*/

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

            if constexpr (has_thermal<source_type>()) {
                A_thermal.resize(nlyr);

                C_plus_thermal.resize(nlyr);
                C_minus_thermal.resize(nlyr);
            }
            if constexpr (has_solar<source_type>()) {
                A_plus.resize(nlyr);
                A_minus.resize(nlyr);

                Q_plus.resize(nlyr);
                Q_minus.resize(nlyr);

                C_plus.resize(nlyr);
                C_minus.resize(nlyr);
            }
            norm.resize(nlyr);
            ;
        }
    };

    struct BVPDerivMatrix {
        Eigen::VectorXd _storage;

        void init(int nlyr) {
            _storage.resize((nlyr + 1) * 5);

            _storage.setZero();
        }

        double& operator()(int a, int b) {
            int j = (a - 1) / 2;
            int a_offset = (a) % 3;
            int b_offset = b - j;

            int index = j * 5 + 3 * a_offset + b_offset;

            return _storage[index];
        }
    };
    /**
     *  Storage for the boundary value problem coefficients and cached
     * quantities
     *
     */
    template <SourceType source_type> struct BVPCoeffs {
        Eigen::VectorXd e, c, d, b, a; /** Pentadiagonal matrix diagonals*/

        Eigen::VectorXd gamma, mu, alpha,
            beta; /** Temporary storage for the pentadiagonal solver */

        Eigen::MatrixXd z; /** Temporary storage for the pentadiagonal solver */
        Eigen::MatrixXd rhs; /** RHS of the pentadiagonal system */

        // Gradient mappings from the RHS to the G terms
        Eigen::RowVectorXd d_G_plus_top;
        Eigen::RowVectorXd d_G_plus_bottom;
        Eigen::RowVectorXd d_G_minus_top;
        Eigen::RowVectorXd d_G_minus_bottom;

        // Gradient mappings from the adjoint solution to the diagonals
        Eigen::RowVectorXd d_e, d_c, d_d, d_b, d_a;

        // Gradient mappings of the diagonals with respect to ssa
        BVPDerivMatrix d_e_by_ssa, d_c_by_ssa, d_d_by_ssa, d_b_by_ssa,
            d_a_by_ssa;

        // Gradient mappings of the diagonals with respect to b1
        BVPDerivMatrix d_e_by_b1, d_c_by_b1, d_d_by_b1, d_b_by_b1, d_a_by_b1;

        // Gradient mappings of the diagonals with respect to od
        BVPDerivMatrix d_e_by_od, d_c_by_od, d_d_by_od, d_b_by_od, d_a_by_od;

        // Temporary storage
        Eigen::RowVectorXd d_temp_ssa;
        Eigen::RowVectorXd d_temp_od;
        Eigen::RowVectorXd d_temp_transmission;
        Eigen::RowVectorXd d_temp_secant;
        Eigen::RowVectorXd d_temp_thermal_b0;
        Eigen::RowVectorXd d_temp_thermal_b1;

        void init(int nlyr) {
            z.resize(nlyr * 2, 1);
            rhs.resize(nlyr * 2, 1);

            z.setZero();
            rhs.setZero();

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

            gamma.setZero();
            mu.setZero();
            alpha.setZero();
            beta.setZero();

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

            d_temp_ssa.resize(nlyr);
            d_temp_ssa.setZero();
            d_temp_od.resize(nlyr);
            d_temp_od.setZero();

            if constexpr (has_thermal<source_type>()) {
                d_temp_thermal_b0.resize(nlyr);
                d_temp_thermal_b0.setZero();
                d_temp_thermal_b1.resize(nlyr);
                d_temp_thermal_b1.setZero();
            }

            if constexpr (has_solar<source_type>()) {
                d_temp_transmission.resize(nlyr + 1);
                d_temp_transmission.setZero();
                d_temp_secant.resize(nlyr);
                d_temp_secant.setZero();
            }

            d_e_by_ssa.init(nlyr);
            d_c_by_ssa.init(nlyr);
            d_d_by_ssa.init(nlyr);
            d_b_by_ssa.init(nlyr);
            d_a_by_ssa.init(nlyr);

            d_e_by_b1.init(nlyr);
            d_c_by_b1.init(nlyr);
            d_d_by_b1.init(nlyr);
            d_b_by_b1.init(nlyr);
            d_a_by_b1.init(nlyr);

            d_e_by_od.init(nlyr);
            d_c_by_od.init(nlyr);
            d_d_by_od.init(nlyr);
            d_b_by_od.init(nlyr);
            d_a_by_od.init(nlyr);
        }
    };

    /**
     *  Storage for the entire discrete ordinates solution at both azimuthal
     * orders
     *
     */
    template <SourceType source_type> struct Solution {
        std::array<HomogSolution<source_type>, num_azimuth<source_type>()>
            homog; /** Homogeneous solutino for both azimuthal orders */
        std::array<ParticularSolution<source_type>, num_azimuth<source_type>()>
            particular; /** Particular solution for both azimuthal orders*/
        std::array<BVPCoeffs<source_type>, num_azimuth<source_type>()>
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
    template <SourceType source_type> struct Sources {
        LayerVector<true, true, true, true, true, true, true, source_type>
            source; /** Integrated sources in each [lyr] */

        LayerVector<false, true, false, false, false, false, false, source_type>
            beamtrans; /** LOS transmission factors exp(-od /
         viewing_zenith) for each [layer] */

        Eigen::VectorXd final_weight_factors; /** Final weight factors for each
                                                [layer] */

        std::array<LayerVector<true, false, true, false, false, false, false,
                               source_type>,
                   num_azimuth<source_type>()>
            lpsum_plus, lpsum_minus; /** LP triple products for upwelling and
                                        downwelling [lyr] */
        std::array<LayerVector<true, false, true, false, false, false, false,
                               source_type>,
                   num_azimuth<source_type>()>
            Y_plus,
            Y_minus; /** "Interpolated" homogenous solutions in each [lyr] */

        std::array<LayerVector<true, true, true, false, false, false, false,
                               source_type>,
                   num_azimuth<source_type>()>
            H_plus, H_minus; /** Homogenous solution multipliers [lyr] */
        std::array<LayerVector<true, true, true, true, true, false, false,
                               source_type>,
                   num_azimuth<source_type>()>
            D_plus, D_minus; /** Particular solution multipliers [lyr] */

        std::array<LayerVector<true, true, true, false, false, true, true,
                               source_type>,
                   num_azimuth<source_type>()>
            D_plus_thermal,
            D_minus_thermal; /** Particular solution multipliers [lyr] */

        std::array<LayerVector<true, true, true, true, true, false, false,
                               source_type>,
                   num_azimuth<source_type>()>
            V; /** Particular source [lyr] */

        LayerVector<false, true, false, false, true, false, false, source_type>
            E_minus; /** Solar multiplier [lyr] */

        LayerVector<false, true, false, false, false, false, true, source_type>
            E_thermal; /** thermal multiplier [lyr] */

        // Backprop factors
        std::array<Eigen::RowVectorXd, num_azimuth<source_type>()> d_bvp_coeff;

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

            if constexpr (has_thermal<source_type>()) {
                for (auto& l : D_plus_thermal) {
                    l.resize(nlyr);
                }
                for (auto& l : D_minus_thermal) {
                    l.resize(nlyr);
                }
                E_thermal.resize(nlyr);
            }
            E_minus.resize(nlyr);

            for (auto& d : d_bvp_coeff) {
                d.resize(nlyr * 2);
                d.setZero();
            }
        }
    };
} // namespace sasktran2::twostream
