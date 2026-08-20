#pragma once

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <complex>
#include <vector>

namespace sasktran2::math {
    class UnitSphere;
}

namespace sasktran2::successive_orders {

    /** Geometry-only real spherical-harmonic transform for scalar scattering.
     *
     * The normalization follows SASKTRAN's Legendre convention: the dot
     * product of all degree-l modes evaluated at two directions is P_l of the
     * direction cosine. Atmospheric Legendre coefficients therefore multiply
     * the analyzed modes directly.
     */
    class ScalarAngularBasis {
      public:
        ScalarAngularBasis(const sasktran2::math::UnitSphere& incoming,
                           const sasktran2::math::UnitSphere& outgoing,
                           int num_coefficients);

        int input_size() const { return static_cast<int>(m_analysis.cols()); }
        int output_size() const { return static_cast<int>(m_synthesis.rows()); }
        int num_coefficients() const { return m_num_coefficients; }
        int num_modes() const {
            return static_cast<int>(m_mode_degrees.size());
        }

        std::size_t storage_bytes() const;

        /** Applies point-local scattering blocks.
         *
         * Inputs and outputs are point-major matrices. Coefficients have one
         * row per point and one column per Legendre degree.
         */
        void apply(Eigen::Ref<const Eigen::MatrixXd> incoming,
                   Eigen::Ref<const Eigen::MatrixXd> coefficients,
                   Eigen::Ref<Eigen::MatrixXd> outgoing,
                   Eigen::MatrixXd& moment_workspace) const;

        void apply_active(Eigen::Ref<const Eigen::MatrixXd> incoming,
                          Eigen::Ref<const Eigen::MatrixXd> coefficients,
                          int active_coefficients,
                          Eigen::Ref<Eigen::MatrixXd> outgoing,
                          Eigen::MatrixXd& moment_workspace) const;

        /** Applies the transpose with respect to angular radiances. */
        void apply_transpose(Eigen::Ref<const Eigen::MatrixXd> outgoing,
                             Eigen::Ref<const Eigen::MatrixXd> coefficients,
                             Eigen::Ref<Eigen::MatrixXd> incoming,
                             Eigen::MatrixXd& moment_workspace) const;

        void
        apply_transpose_active(Eigen::Ref<const Eigen::MatrixXd> outgoing,
                               Eigen::Ref<const Eigen::MatrixXd> coefficients,
                               int active_coefficients,
                               Eigen::Ref<Eigen::MatrixXd> incoming,
                               Eigen::MatrixXd& moment_workspace) const;

        /** Forward product of the scattering operation. */
        void apply_jvp(Eigen::Ref<const Eigen::MatrixXd> incoming,
                       Eigen::Ref<const Eigen::MatrixXd> incoming_tangent,
                       Eigen::Ref<const Eigen::MatrixXd> coefficients,
                       Eigen::Ref<const Eigen::MatrixXd> coefficient_tangent,
                       Eigen::Ref<Eigen::MatrixXd> outgoing_tangent,
                       Eigen::MatrixXd& moment_workspace,
                       Eigen::MatrixXd& tangent_workspace) const;

        void
        apply_jvp_active(Eigen::Ref<const Eigen::MatrixXd> incoming,
                         Eigen::Ref<const Eigen::MatrixXd> incoming_tangent,
                         Eigen::Ref<const Eigen::MatrixXd> coefficients,
                         Eigen::Ref<const Eigen::MatrixXd> coefficient_tangent,
                         int active_coefficients,
                         Eigen::Ref<Eigen::MatrixXd> outgoing_tangent,
                         Eigen::MatrixXd& moment_workspace,
                         Eigen::MatrixXd& tangent_workspace) const;

        /** Reverse product of the scattering operation. */
        void apply_vjp(Eigen::Ref<const Eigen::MatrixXd> incoming,
                       Eigen::Ref<const Eigen::MatrixXd> coefficients,
                       Eigen::Ref<const Eigen::MatrixXd> outgoing_cotangent,
                       Eigen::Ref<Eigen::MatrixXd> incoming_cotangent,
                       Eigen::Ref<Eigen::MatrixXd> coefficient_gradient,
                       Eigen::MatrixXd& analyzed_workspace,
                       Eigen::MatrixXd& moment_cotangent_workspace) const;

      private:
        void validate_blocks(Eigen::Ref<const Eigen::MatrixXd> incoming,
                             Eigen::Ref<const Eigen::MatrixXd> coefficients,
                             Eigen::Ref<const Eigen::MatrixXd> outgoing) const;
        void multiply_coefficients(
            Eigen::MatrixXd& moments,
            Eigen::Ref<const Eigen::MatrixXd> coefficients) const;

        int m_num_coefficients;
        Eigen::MatrixXd m_analysis;
        Eigen::MatrixXd m_synthesis;
        std::vector<int> m_mode_degrees;
    };

    /** Reusable storage for polarized spin-harmonic transforms. */
    struct VectorAngularWorkspace {
        std::array<Eigen::MatrixXd, 3> stokes_input;
        std::array<Eigen::MatrixXd, 3> stokes_output;
        std::array<Eigen::MatrixXd, 6> moments;
        Eigen::MatrixXd stacked_input;
        Eigen::MatrixXd stacked_moments;
        Eigen::MatrixXd real_products;
        Eigen::MatrixXd imaginary_products;

        std::size_t storage_bytes() const;
    };

    /** Geometry-only spin-0/spin-2 transform for polarized I/Q/U scattering.
     *
     * Each (l,m) mode contains intensity, Q+iU, and Q-iU channels. The four
     * Greek coefficient families mix these channels locally, allowing all
     * atmospheric points to be scattered with matrix-matrix transforms rather
     * than retaining one dense angular matrix per point.
     */
    class VectorAngularBasis {
      public:
        VectorAngularBasis(const sasktran2::math::UnitSphere& incoming,
                           const sasktran2::math::UnitSphere& outgoing,
                           int num_coefficients);

        int input_directions() const {
            return static_cast<int>(m_analysis_real[0].cols());
        }
        int output_directions() const {
            return static_cast<int>(m_synthesis_real[0].rows());
        }
        int input_size() const { return 3 * input_directions(); }
        int output_size() const { return 3 * output_directions(); }
        int num_coefficients() const { return m_num_coefficients; }
        int num_modes() const {
            return static_cast<int>(m_mode_degrees.size());
        }

        std::size_t storage_bytes() const;

        void apply_active(Eigen::Ref<const Eigen::MatrixXd> incoming,
                          Eigen::Ref<const Eigen::MatrixXd> coefficients,
                          int active_coefficients,
                          Eigen::Ref<Eigen::MatrixXd> outgoing,
                          VectorAngularWorkspace& workspace) const;

        void
        apply_transpose_active(Eigen::Ref<const Eigen::MatrixXd> outgoing,
                               Eigen::Ref<const Eigen::MatrixXd> coefficients,
                               int active_coefficients,
                               Eigen::Ref<Eigen::MatrixXd> incoming,
                               VectorAngularWorkspace& workspace) const;

        void accumulate_coefficient_vjp(
            Eigen::Ref<const Eigen::MatrixXd> incoming,
            Eigen::Ref<const Eigen::MatrixXd> outgoing_cotangent,
            Eigen::Ref<Eigen::MatrixXd> coefficient_gradient,
            VectorAngularWorkspace& workspace) const;

      private:
        struct FrameCorrection {
            int input_index = 0;
            int output_index = 0;
            // [degree, Greek family, output Stokes, input Stokes]
            std::vector<double> values;
        };

        void validate_blocks(Eigen::Ref<const Eigen::MatrixXd> incoming,
                             Eigen::Ref<const Eigen::MatrixXd> coefficients,
                             Eigen::Ref<const Eigen::MatrixXd> outgoing) const;
        void analyze(Eigen::Ref<const Eigen::MatrixXd> incoming,
                     int active_modes, VectorAngularWorkspace& workspace) const;
        void
        analyze_synthesis_cotangent(Eigen::Ref<const Eigen::MatrixXd> outgoing,
                                    int active_modes,
                                    VectorAngularWorkspace& workspace) const;
        void
        apply_frame_corrections(Eigen::Ref<const Eigen::MatrixXd> incoming,
                                Eigen::Ref<const Eigen::MatrixXd> coefficients,
                                int active_coefficients,
                                Eigen::Ref<Eigen::MatrixXd> outgoing) const;
        void apply_frame_corrections_transpose(
            Eigen::Ref<const Eigen::MatrixXd> outgoing,
            Eigen::Ref<const Eigen::MatrixXd> coefficients,
            int active_coefficients,
            Eigen::Ref<Eigen::MatrixXd> incoming) const;
        void accumulate_frame_correction_vjp(
            Eigen::Ref<const Eigen::MatrixXd> incoming,
            Eigen::Ref<const Eigen::MatrixXd> outgoing_cotangent,
            Eigen::Ref<Eigen::MatrixXd> coefficient_gradient) const;

        int m_num_coefficients;
        std::array<Eigen::MatrixXd, 3> m_analysis_real;
        std::array<Eigen::MatrixXd, 3> m_analysis_imaginary;
        std::array<Eigen::MatrixXd, 3> m_synthesis_real;
        std::array<Eigen::MatrixXd, 3> m_synthesis_imaginary;
        std::vector<int> m_mode_degrees;
        std::vector<FrameCorrection> m_frame_corrections;
    };

} // namespace sasktran2::successive_orders
