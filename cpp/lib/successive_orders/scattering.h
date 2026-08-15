#pragma once

#include "angular_basis.h"

#include <Eigen/Core>

#include <cstddef>
#include <memory>
#include <vector>

namespace sasktran2::successive_orders {

    /** Point-major block layout shared by scattering and ray transport.
     *
     * Each point owns a contiguous block of directions and each direction owns
     * NSTOKES contiguous values. Atmospheric blocks precede boundary blocks.
     * Explicit offsets allow a ground hemisphere to use fewer directions than
     * an atmospheric sphere.
     */
    class ScatteringBlockLayout {
      public:
        ScatteringBlockLayout(int atmospheric_blocks, int stokes_components,
                              std::vector<int> input_offsets,
                              std::vector<int> output_offsets);

        ScatteringBlockLayout(int atmospheric_blocks, int ground_blocks,
                              int atmospheric_input_directions,
                              int atmospheric_output_directions,
                              int ground_input_directions,
                              int ground_output_directions,
                              int stokes_components);

        int blocks() const {
            return static_cast<int>(m_input_offsets.size()) - 1;
        }
        int atmospheric_blocks() const { return m_atmospheric_blocks; }
        int ground_blocks() const { return blocks() - m_atmospheric_blocks; }
        int stokes_components() const { return m_stokes_components; }
        int input_size() const { return m_input_offsets.back(); }
        int output_size() const { return m_output_offsets.back(); }

        int input_block_size(int block) const;
        int output_block_size(int block) const;
        int input_directions(int block) const {
            return input_block_size(block) / m_stokes_components;
        }
        int output_directions(int block) const {
            return output_block_size(block) / m_stokes_components;
        }

        const std::vector<int>& input_offsets() const {
            return m_input_offsets;
        }
        const std::vector<int>& output_offsets() const {
            return m_output_offsets;
        }
        std::size_t storage_bytes() const;

      private:
        int m_atmospheric_blocks;
        int m_stokes_components;
        std::vector<int> m_input_offsets;
        std::vector<int> m_output_offsets;
    };

    struct ScatteringMemoryUsage {
        std::size_t layout_bytes = 0;
        std::size_t angular_basis_bytes = 0;
        std::size_t atmospheric_value_bytes = 0;
        std::size_t boundary_value_bytes = 0;
        std::size_t workspace_bytes = 0;

        std::size_t operator_bytes() const {
            return layout_bytes + angular_basis_bytes +
                   atmospheric_value_bytes + boundary_value_bytes;
        }
        std::size_t total_bytes() const {
            return operator_bytes() + workspace_bytes;
        }
    };

    template <int NSTOKES> class ScatteringOperator;
    template <int NSTOKES> class ScatteringWorkspace;

    template <> class ScatteringWorkspace<1> {
      public:
        ScatteringWorkspace() = default;

        std::size_t storage_bytes() const;

      private:
        friend class ScatteringOperator<1>;
        void prepare(int atmospheric_blocks, int input_directions,
                     int output_directions, int modes);

        Eigen::MatrixXd m_atmospheric_input;
        Eigen::MatrixXd m_atmospheric_output;
        Eigen::MatrixXd m_auxiliary_input;
        Eigen::MatrixXd m_moments;
        Eigen::MatrixXd m_auxiliary_moments;
    };

    /** Scalar scattering with coefficient-space atmospheric blocks.
     *
     * The expensive angular basis is stored once. Each atmospheric point only
     * stores its Legendre coefficients. Ground blocks remain dense because a
     * general BRDF is not diagonal in the atmospheric angular basis.
     */
    template <> class ScatteringOperator<1> {
      public:
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;

        ScatteringOperator(ScatteringBlockLayout layout,
                           ScalarAngularBasis angular_basis);
        ScatteringOperator(
            ScatteringBlockLayout layout,
            std::shared_ptr<const ScalarAngularBasis> angular_basis);

        const ScatteringBlockLayout& layout() const { return m_layout; }
        const std::vector<int>& input_offsets() const {
            return m_layout.input_offsets();
        }
        const std::vector<int>& output_offsets() const {
            return m_layout.output_offsets();
        }
        int input_size() const { return m_layout.input_size(); }
        int output_size() const { return m_layout.output_size(); }

        const ScalarAngularBasis& angular_basis() const { return *m_basis; }
        Eigen::MatrixXd& atmospheric_coefficients() {
            // Mutable access can activate any trailing mode. Callers that know
            // the resulting exact active count may narrow it after filling.
            m_active_coefficients = m_basis->num_coefficients();
            return m_atmospheric_coefficients;
        }
        const Eigen::MatrixXd& atmospheric_coefficients() const {
            return m_atmospheric_coefficients;
        }
        Eigen::VectorXd& ground_values() { return m_ground_values; }
        const Eigen::VectorXd& ground_values() const { return m_ground_values; }

        int coefficient_value_size() const {
            return static_cast<int>(m_atmospheric_coefficients.size());
        }
        int ground_value_size() const {
            return static_cast<int>(m_ground_values.size());
        }

        void set_active_coefficients(int active_coefficients);
        int active_coefficients() const { return m_active_coefficients; }

        void set_atmospheric_coefficients(
            Eigen::Ref<const Eigen::MatrixXd> coefficients);
        void set_ground_values(Eigen::Ref<const Eigen::VectorXd> values);
        void set_ground_block(int ground_block,
                              Eigen::Ref<const Eigen::MatrixXd> values);

        Eigen::Map<RowMajorMatrix> ground_block(int ground_block);
        Eigen::Map<const RowMajorMatrix> ground_block(int ground_block) const;

        ScatteringWorkspace<1> make_workspace() const;
        ScatteringMemoryUsage
        memory_usage(const ScatteringWorkspace<1>* workspace = nullptr) const;

        void apply(Eigen::Ref<const Eigen::VectorXd> incoming,
                   Eigen::Ref<Eigen::VectorXd> outgoing,
                   ScatteringWorkspace<1>& workspace) const;
        void apply_transpose(Eigen::Ref<const Eigen::VectorXd> outgoing,
                             Eigen::Ref<Eigen::VectorXd> incoming,
                             ScatteringWorkspace<1>& workspace) const;

        void apply_jvp(Eigen::Ref<const Eigen::VectorXd> incoming,
                       Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
                       Eigen::Ref<const Eigen::MatrixXd> coefficient_tangent,
                       Eigen::Ref<const Eigen::VectorXd> ground_value_tangent,
                       Eigen::Ref<Eigen::VectorXd> outgoing_tangent,
                       ScatteringWorkspace<1>& workspace) const;

        void apply_vjp(Eigen::Ref<const Eigen::VectorXd> incoming,
                       Eigen::Ref<const Eigen::VectorXd> outgoing_cotangent,
                       Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
                       Eigen::Ref<Eigen::MatrixXd> coefficient_gradient,
                       Eigen::Ref<Eigen::VectorXd> ground_value_gradient,
                       ScatteringWorkspace<1>& workspace) const;

      private:
        void prepare_workspace(ScatteringWorkspace<1>& workspace) const;
        void validate_input_output(Eigen::Index incoming_size,
                                   Eigen::Index outgoing_size) const;
        int ground_value_offset(int ground_block) const;

        ScatteringBlockLayout m_layout;
        std::shared_ptr<const ScalarAngularBasis> m_basis;
        Eigen::MatrixXd m_atmospheric_coefficients;
        std::vector<int> m_ground_value_offsets;
        Eigen::VectorXd m_ground_values;
        int m_active_coefficients = 1;
    };

    template <> class ScatteringWorkspace<3> {
      public:
        std::size_t storage_bytes() const;

      private:
        friend class ScatteringOperator<3>;
        void prepare(int atmospheric_blocks, int input_directions,
                     int output_directions);

        Eigen::MatrixXd m_atmospheric_input;
        Eigen::MatrixXd m_atmospheric_output;
        Eigen::MatrixXd m_auxiliary_input;
        Eigen::MatrixXd m_auxiliary_output;
        VectorAngularWorkspace m_angular;
    };

    /** Coefficient-space I/Q/U atmospheric scattering with dense boundaries. */
    template <> class ScatteringOperator<3> {
      public:
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;

        ScatteringOperator(
            ScatteringBlockLayout layout,
            std::shared_ptr<const VectorAngularBasis> angular_basis);

        const ScatteringBlockLayout& layout() const { return m_layout; }
        const std::vector<int>& input_offsets() const {
            return m_layout.input_offsets();
        }
        const std::vector<int>& output_offsets() const {
            return m_layout.output_offsets();
        }
        int input_size() const { return m_layout.input_size(); }
        int output_size() const { return m_layout.output_size(); }

        Eigen::MatrixXd& atmospheric_coefficients() {
            m_active_coefficients = m_basis->num_coefficients();
            return m_atmospheric_coefficients;
        }
        const Eigen::MatrixXd& atmospheric_coefficients() const {
            return m_atmospheric_coefficients;
        }
        Eigen::VectorXd& ground_values() { return m_ground_values; }
        const Eigen::VectorXd& ground_values() const { return m_ground_values; }

        int coefficient_value_size() const {
            return static_cast<int>(m_atmospheric_coefficients.size());
        }
        int ground_value_size() const {
            return static_cast<int>(m_ground_values.size());
        }

        void set_active_coefficients(int active_coefficients);
        int active_coefficients() const { return m_active_coefficients; }
        void set_atmospheric_coefficients(
            Eigen::Ref<const Eigen::MatrixXd> coefficients);
        void set_ground_values(Eigen::Ref<const Eigen::VectorXd> values);
        void set_ground_block(int ground_block,
                              Eigen::Ref<const Eigen::MatrixXd> values);

        Eigen::Map<RowMajorMatrix> ground_block(int ground_block);
        Eigen::Map<const RowMajorMatrix> ground_block(int ground_block) const;

        ScatteringWorkspace<3> make_workspace() const;
        ScatteringMemoryUsage
        memory_usage(const ScatteringWorkspace<3>* workspace = nullptr) const;

        void apply(Eigen::Ref<const Eigen::VectorXd> incoming,
                   Eigen::Ref<Eigen::VectorXd> outgoing,
                   ScatteringWorkspace<3>& workspace) const;
        void apply_transpose(Eigen::Ref<const Eigen::VectorXd> outgoing,
                             Eigen::Ref<Eigen::VectorXd> incoming,
                             ScatteringWorkspace<3>& workspace) const;
        void apply_jvp(Eigen::Ref<const Eigen::VectorXd> incoming,
                       Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
                       Eigen::Ref<const Eigen::MatrixXd> coefficient_tangent,
                       Eigen::Ref<const Eigen::VectorXd> ground_value_tangent,
                       Eigen::Ref<Eigen::VectorXd> outgoing_tangent,
                       ScatteringWorkspace<3>& workspace) const;
        void apply_vjp(Eigen::Ref<const Eigen::VectorXd> incoming,
                       Eigen::Ref<const Eigen::VectorXd> outgoing_cotangent,
                       Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
                       Eigen::Ref<Eigen::MatrixXd> coefficient_gradient,
                       Eigen::Ref<Eigen::VectorXd> ground_value_gradient,
                       ScatteringWorkspace<3>& workspace) const;

      private:
        void validate_input_output(Eigen::Index incoming_size,
                                   Eigen::Index outgoing_size) const;
        void prepare_workspace(ScatteringWorkspace<3>& workspace) const;
        int ground_value_offset(int ground_block) const;

        ScatteringBlockLayout m_layout;
        std::shared_ptr<const VectorAngularBasis> m_basis;
        std::vector<int> m_ground_value_offsets;
        Eigen::MatrixXd m_atmospheric_coefficients;
        Eigen::VectorXd m_ground_values;
        int m_active_coefficients = 1;
    };

} // namespace sasktran2::successive_orders
