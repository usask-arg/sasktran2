#include "scattering.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace sasktran2::successive_orders {
    namespace {
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;

        int checked_product(int left, int right, const char* message) {
            const auto product = static_cast<long long>(left) * right;
            if (left < 0 || right < 0 ||
                product > std::numeric_limits<int>::max()) {
                throw std::invalid_argument(message);
            }
            return static_cast<int>(product);
        }

        std::vector<int> make_uniform_offsets(int atmospheric_blocks,
                                              int ground_blocks,
                                              int atmospheric_directions,
                                              int ground_directions,
                                              int stokes_components) {
            if (atmospheric_blocks < 0 || ground_blocks < 0 ||
                atmospheric_blocks + ground_blocks < 1 ||
                atmospheric_directions < 1 || ground_directions < 1 ||
                stokes_components < 1) {
                throw std::invalid_argument(
                    "invalid successive-orders scattering block counts");
            }
            const int atmospheric_size = checked_product(
                atmospheric_directions, stokes_components,
                "successive-orders atmospheric block size overflow");
            const int ground_size =
                checked_product(ground_directions, stokes_components,
                                "successive-orders ground block size overflow");
            std::vector<int> offsets(
                static_cast<std::size_t>(atmospheric_blocks + ground_blocks) +
                    1,
                0);
            for (int block = 0; block < atmospheric_blocks; ++block) {
                if (offsets[block] >
                    std::numeric_limits<int>::max() - atmospheric_size) {
                    throw std::invalid_argument(
                        "successive-orders scattering layout overflow");
                }
                offsets[block + 1] = offsets[block] + atmospheric_size;
            }
            for (int block = atmospheric_blocks;
                 block < atmospheric_blocks + ground_blocks; ++block) {
                if (offsets[block] >
                    std::numeric_limits<int>::max() - ground_size) {
                    throw std::invalid_argument(
                        "successive-orders scattering layout overflow");
                }
                offsets[block + 1] = offsets[block] + ground_size;
            }
            return offsets;
        }

        std::vector<int>
        make_dense_value_offsets(const ScatteringBlockLayout& layout,
                                 int first_block, int block_count) {
            if (first_block < 0 || block_count < 0 ||
                first_block + block_count > layout.blocks()) {
                throw std::invalid_argument(
                    "invalid successive-orders dense scattering block range");
            }
            std::vector<int> offsets(static_cast<std::size_t>(block_count) + 1,
                                     0);
            for (int local_block = 0; local_block < block_count;
                 ++local_block) {
                const int block = first_block + local_block;
                const int values = checked_product(
                    layout.input_block_size(block),
                    layout.output_block_size(block),
                    "successive-orders dense scattering size overflow");
                if (offsets[local_block] >
                    std::numeric_limits<int>::max() - values) {
                    throw std::invalid_argument(
                        "successive-orders dense scattering size overflow");
                }
                offsets[local_block + 1] = offsets[local_block] + values;
            }
            return offsets;
        }

        void validate_dense_tangent(const Eigen::VectorXd& values,
                                    Eigen::Ref<const Eigen::VectorXd> tangent) {
            if (tangent.size() != values.size()) {
                throw std::invalid_argument(
                    "invalid successive-orders dense scattering tangent size");
            }
        }

        void apply_dense_blocks(const ScatteringBlockLayout& layout,
                                int first_block,
                                const std::vector<int>& value_offsets,
                                const Eigen::VectorXd& values,
                                Eigen::Ref<const Eigen::VectorXd> incoming,
                                Eigen::Ref<Eigen::VectorXd> outgoing) {
            const int block_count = static_cast<int>(value_offsets.size()) - 1;
            for (int local_block = 0; local_block < block_count;
                 ++local_block) {
                const int block = first_block + local_block;
                const int rows = layout.output_block_size(block);
                const int columns = layout.input_block_size(block);
                const Eigen::Map<const RowMajorMatrix> matrix(
                    values.data() + value_offsets[local_block], rows, columns);
                outgoing.segment(layout.output_offsets()[block], rows)
                    .noalias() =
                    matrix *
                    incoming.segment(layout.input_offsets()[block], columns);
            }
        }

        void
        apply_dense_blocks_transpose(const ScatteringBlockLayout& layout,
                                     int first_block,
                                     const std::vector<int>& value_offsets,
                                     const Eigen::VectorXd& values,
                                     Eigen::Ref<const Eigen::VectorXd> outgoing,
                                     Eigen::Ref<Eigen::VectorXd> incoming) {
            const int block_count = static_cast<int>(value_offsets.size()) - 1;
            for (int local_block = 0; local_block < block_count;
                 ++local_block) {
                const int block = first_block + local_block;
                const int rows = layout.output_block_size(block);
                const int columns = layout.input_block_size(block);
                const Eigen::Map<const RowMajorMatrix> matrix(
                    values.data() + value_offsets[local_block], rows, columns);
                incoming.segment(layout.input_offsets()[block], columns)
                    .noalias() =
                    matrix.transpose() *
                    outgoing.segment(layout.output_offsets()[block], rows);
            }
        }

        void apply_dense_blocks_jvp(
            const ScatteringBlockLayout& layout, int first_block,
            const std::vector<int>& value_offsets,
            const Eigen::VectorXd& values,
            Eigen::Ref<const Eigen::VectorXd> value_tangent,
            Eigen::Ref<const Eigen::VectorXd> incoming,
            Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
            Eigen::Ref<Eigen::VectorXd> outgoing_tangent) {
            const int block_count = static_cast<int>(value_offsets.size()) - 1;
            for (int local_block = 0; local_block < block_count;
                 ++local_block) {
                const int block = first_block + local_block;
                const int rows = layout.output_block_size(block);
                const int columns = layout.input_block_size(block);
                const int value_offset = value_offsets[local_block];
                const Eigen::Map<const RowMajorMatrix> matrix(
                    values.data() + value_offset, rows, columns);
                const Eigen::Map<const RowMajorMatrix> matrix_tangent(
                    value_tangent.data() + value_offset, rows, columns);
                auto result = outgoing_tangent.segment(
                    layout.output_offsets()[block], rows);
                result.noalias() =
                    matrix * incoming_tangent.segment(
                                 layout.input_offsets()[block], columns);
                result.noalias() +=
                    matrix_tangent *
                    incoming.segment(layout.input_offsets()[block], columns);
            }
        }

        void apply_dense_blocks_vjp(
            const ScatteringBlockLayout& layout, int first_block,
            const std::vector<int>& value_offsets,
            const Eigen::VectorXd& values,
            Eigen::Ref<const Eigen::VectorXd> incoming,
            Eigen::Ref<const Eigen::VectorXd> outgoing_cotangent,
            Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
            Eigen::Ref<Eigen::VectorXd> value_gradient) {
            const int block_count = static_cast<int>(value_offsets.size()) - 1;
            for (int local_block = 0; local_block < block_count;
                 ++local_block) {
                const int block = first_block + local_block;
                const int rows = layout.output_block_size(block);
                const int columns = layout.input_block_size(block);
                const int value_offset = value_offsets[local_block];
                const Eigen::Map<const RowMajorMatrix> matrix(
                    values.data() + value_offset, rows, columns);
                const auto block_cotangent = outgoing_cotangent.segment(
                    layout.output_offsets()[block], rows);
                const auto block_input =
                    incoming.segment(layout.input_offsets()[block], columns);
                incoming_cotangent
                    .segment(layout.input_offsets()[block], columns)
                    .noalias() = matrix.transpose() * block_cotangent;
                Eigen::Map<RowMajorMatrix> gradient(
                    value_gradient.data() + value_offset, rows, columns);
                gradient.noalias() = block_cotangent * block_input.transpose();
            }
        }

        void pack_atmospheric_blocks(const ScatteringBlockLayout& layout,
                                     Eigen::Ref<const Eigen::VectorXd> values,
                                     Eigen::MatrixXd& packed) {
            for (int point = 0; point < layout.atmospheric_blocks(); ++point) {
                packed.row(point) = values
                                        .segment(layout.input_offsets()[point],
                                                 layout.input_block_size(point))
                                        .transpose();
            }
        }

        void
        pack_atmospheric_output_blocks(const ScatteringBlockLayout& layout,
                                       Eigen::Ref<const Eigen::VectorXd> values,
                                       Eigen::MatrixXd& packed) {
            for (int point = 0; point < layout.atmospheric_blocks(); ++point) {
                packed.row(point) =
                    values
                        .segment(layout.output_offsets()[point],
                                 layout.output_block_size(point))
                        .transpose();
            }
        }

        void unpack_atmospheric_blocks(const ScatteringBlockLayout& layout,
                                       const Eigen::MatrixXd& packed,
                                       Eigen::Ref<Eigen::VectorXd> values) {
            for (int point = 0; point < layout.atmospheric_blocks(); ++point) {
                values.segment(layout.output_offsets()[point],
                               layout.output_block_size(point)) =
                    packed.row(point).transpose();
            }
        }

        void
        unpack_atmospheric_input_blocks(const ScatteringBlockLayout& layout,
                                        const Eigen::MatrixXd& packed,
                                        Eigen::Ref<Eigen::VectorXd> values) {
            for (int point = 0; point < layout.atmospheric_blocks(); ++point) {
                values.segment(layout.input_offsets()[point],
                               layout.input_block_size(point)) =
                    packed.row(point).transpose();
            }
        }
    } // namespace

    ScatteringBlockLayout::ScatteringBlockLayout(
        int atmospheric_blocks, int stokes_components,
        std::vector<int> input_offsets, std::vector<int> output_offsets)
        : m_atmospheric_blocks(atmospheric_blocks),
          m_stokes_components(stokes_components),
          m_input_offsets(std::move(input_offsets)),
          m_output_offsets(std::move(output_offsets)) {
        if (m_stokes_components < 1 || m_input_offsets.size() < 2 ||
            m_input_offsets.size() != m_output_offsets.size() ||
            m_atmospheric_blocks < 0 || m_atmospheric_blocks > blocks() ||
            m_input_offsets.front() != 0 || m_output_offsets.front() != 0 ||
            !std::is_sorted(m_input_offsets.begin(), m_input_offsets.end()) ||
            !std::is_sorted(m_output_offsets.begin(), m_output_offsets.end())) {
            throw std::invalid_argument(
                "invalid successive-orders scattering block offsets");
        }
        for (int block = 0; block < blocks(); ++block) {
            if (input_block_size(block) < 1 || output_block_size(block) < 1 ||
                input_block_size(block) % m_stokes_components != 0 ||
                output_block_size(block) % m_stokes_components != 0) {
                throw std::invalid_argument(
                    "invalid successive-orders scattering block dimensions");
            }
        }
    }

    ScatteringBlockLayout::ScatteringBlockLayout(
        int atmospheric_blocks, int ground_blocks,
        int atmospheric_input_directions, int atmospheric_output_directions,
        int ground_input_directions, int ground_output_directions,
        int stokes_components)
        : ScatteringBlockLayout(
              atmospheric_blocks, stokes_components,
              make_uniform_offsets(atmospheric_blocks, ground_blocks,
                                   atmospheric_input_directions,
                                   ground_input_directions, stokes_components),
              make_uniform_offsets(atmospheric_blocks, ground_blocks,
                                   atmospheric_output_directions,
                                   ground_output_directions,
                                   stokes_components)) {}

    int ScatteringBlockLayout::input_block_size(int block) const {
        if (block < 0 || block >= blocks()) {
            throw std::out_of_range(
                "successive-orders scattering input block out of range");
        }
        return m_input_offsets[block + 1] - m_input_offsets[block];
    }

    int ScatteringBlockLayout::output_block_size(int block) const {
        if (block < 0 || block >= blocks()) {
            throw std::out_of_range(
                "successive-orders scattering output block out of range");
        }
        return m_output_offsets[block + 1] - m_output_offsets[block];
    }

    std::size_t ScatteringBlockLayout::storage_bytes() const {
        return (m_input_offsets.capacity() + m_output_offsets.capacity()) *
               sizeof(int);
    }

    void ScatteringWorkspace<1>::prepare(int atmospheric_blocks,
                                         int input_directions,
                                         int output_directions, int modes) {
        m_atmospheric_input.resize(atmospheric_blocks, input_directions);
        m_atmospheric_output.resize(atmospheric_blocks, output_directions);
        m_auxiliary_input.resize(atmospheric_blocks, input_directions);
        m_moments.resize(atmospheric_blocks, modes);
        m_auxiliary_moments.resize(atmospheric_blocks, modes);
    }

    std::size_t ScatteringWorkspace<1>::storage_bytes() const {
        return static_cast<std::size_t>(
                   m_atmospheric_input.size() + m_atmospheric_output.size() +
                   m_auxiliary_input.size() + m_moments.size() +
                   m_auxiliary_moments.size()) *
               sizeof(double);
    }

    ScatteringOperator<1>::ScatteringOperator(ScatteringBlockLayout layout,
                                              ScalarAngularBasis angular_basis)
        : ScatteringOperator(
              std::move(layout),
              std::make_shared<ScalarAngularBasis>(std::move(angular_basis))) {}

    ScatteringOperator<1>::ScatteringOperator(
        ScatteringBlockLayout layout,
        std::shared_ptr<const ScalarAngularBasis> angular_basis)
        : m_layout(std::move(layout)), m_basis(std::move(angular_basis)),
          m_atmospheric_coefficients(
              m_layout.atmospheric_blocks(),
              m_basis == nullptr ? 0 : m_basis->num_coefficients()),
          m_ground_value_offsets(
              make_dense_value_offsets(m_layout, m_layout.atmospheric_blocks(),
                                       m_layout.ground_blocks())),
          m_ground_values(m_ground_value_offsets.back()) {
        if (m_basis == nullptr) {
            throw std::invalid_argument(
                "scalar successive-orders scattering requires an angular "
                "basis");
        }
        if (m_layout.stokes_components() != 1) {
            throw std::invalid_argument(
                "scalar successive-orders scattering requires NSTOKES=1");
        }
        for (int block = 0; block < m_layout.atmospheric_blocks(); ++block) {
            if (m_layout.input_block_size(block) != m_basis->input_size() ||
                m_layout.output_block_size(block) != m_basis->output_size()) {
                throw std::invalid_argument(
                    "scalar successive-orders atmospheric block does not "
                    "match its angular basis");
            }
        }
        m_atmospheric_coefficients.setZero();
        m_ground_values.setZero();
    }

    void
    ScatteringOperator<1>::set_active_coefficients(int active_coefficients) {
        if (active_coefficients < 1 ||
            active_coefficients > m_basis->num_coefficients()) {
            throw std::invalid_argument(
                "invalid active scalar scattering coefficient count");
        }
        m_active_coefficients = active_coefficients;
    }

    void ScatteringOperator<1>::set_atmospheric_coefficients(
        Eigen::Ref<const Eigen::MatrixXd> coefficients) {
        if (coefficients.rows() != m_atmospheric_coefficients.rows() ||
            coefficients.cols() != m_atmospheric_coefficients.cols() ||
            !coefficients.allFinite()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders atmospheric coefficients");
        }
        m_atmospheric_coefficients = coefficients;
        int active_coefficients = 1;
        for (int degree = coefficients.cols() - 1; degree >= 1; --degree) {
            if (!coefficients.col(degree).isZero(0.0)) {
                active_coefficients = degree + 1;
                break;
            }
        }
        m_active_coefficients = active_coefficients;
    }

    void ScatteringOperator<1>::set_ground_values(
        Eigen::Ref<const Eigen::VectorXd> values) {
        if (values.size() != m_ground_values.size() || !values.allFinite()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders ground scattering values");
        }
        m_ground_values = values;
    }

    void ScatteringOperator<1>::set_ground_block(
        int ground_block, Eigen::Ref<const Eigen::MatrixXd> values) {
        auto target = this->ground_block(ground_block);
        if (values.rows() != target.rows() || values.cols() != target.cols() ||
            !values.allFinite()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders ground scattering block");
        }
        target = values;
    }

    int ScatteringOperator<1>::ground_value_offset(int ground_block) const {
        if (ground_block < 0 || ground_block >= m_layout.ground_blocks()) {
            throw std::out_of_range(
                "successive-orders ground scattering block out of range");
        }
        return m_ground_value_offsets[ground_block];
    }

    Eigen::Map<ScatteringOperator<1>::RowMajorMatrix>
    ScatteringOperator<1>::ground_block(int ground_block) {
        const int block = m_layout.atmospheric_blocks() + ground_block;
        return {m_ground_values.data() + ground_value_offset(ground_block),
                m_layout.output_block_size(block),
                m_layout.input_block_size(block)};
    }

    Eigen::Map<const ScatteringOperator<1>::RowMajorMatrix>
    ScatteringOperator<1>::ground_block(int ground_block) const {
        const int block = m_layout.atmospheric_blocks() + ground_block;
        return {m_ground_values.data() + ground_value_offset(ground_block),
                m_layout.output_block_size(block),
                m_layout.input_block_size(block)};
    }

    void ScatteringOperator<1>::prepare_workspace(
        ScatteringWorkspace<1>& workspace) const {
        workspace.prepare(m_layout.atmospheric_blocks(), m_basis->input_size(),
                          m_basis->output_size(),
                          m_active_coefficients * m_active_coefficients);
    }

    ScatteringWorkspace<1> ScatteringOperator<1>::make_workspace() const {
        ScatteringWorkspace<1> workspace;
        prepare_workspace(workspace);
        return workspace;
    }

    ScatteringMemoryUsage ScatteringOperator<1>::memory_usage(
        const ScatteringWorkspace<1>* workspace) const {
        ScatteringMemoryUsage result;
        result.layout_bytes = m_layout.storage_bytes() +
                              m_ground_value_offsets.capacity() * sizeof(int);
        result.angular_basis_bytes = m_basis->storage_bytes();
        result.atmospheric_value_bytes =
            static_cast<std::size_t>(m_atmospheric_coefficients.size()) *
            sizeof(double);
        result.boundary_value_bytes =
            static_cast<std::size_t>(m_ground_values.size()) * sizeof(double);
        result.workspace_bytes =
            workspace == nullptr ? 0 : workspace->storage_bytes();
        return result;
    }

    void ScatteringOperator<1>::validate_input_output(
        Eigen::Index incoming_size, Eigen::Index outgoing_size) const {
        if (incoming_size != input_size() || outgoing_size != output_size()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders scattering vector sizes");
        }
    }

    void
    ScatteringOperator<1>::apply(Eigen::Ref<const Eigen::VectorXd> incoming,
                                 Eigen::Ref<Eigen::VectorXd> outgoing,
                                 ScatteringWorkspace<1>& workspace) const {
        validate_input_output(incoming.size(), outgoing.size());
        prepare_workspace(workspace);
        if (m_layout.atmospheric_blocks() != 0) {
            pack_atmospheric_blocks(m_layout, incoming,
                                    workspace.m_atmospheric_input);
            m_basis->apply_active(
                workspace.m_atmospheric_input, m_atmospheric_coefficients,
                m_active_coefficients, workspace.m_atmospheric_output,
                workspace.m_moments);
            unpack_atmospheric_blocks(m_layout, workspace.m_atmospheric_output,
                                      outgoing);
        }
        apply_dense_blocks(m_layout, m_layout.atmospheric_blocks(),
                           m_ground_value_offsets, m_ground_values, incoming,
                           outgoing);
    }

    void ScatteringOperator<1>::apply_transpose(
        Eigen::Ref<const Eigen::VectorXd> outgoing,
        Eigen::Ref<Eigen::VectorXd> incoming,
        ScatteringWorkspace<1>& workspace) const {
        validate_input_output(incoming.size(), outgoing.size());
        prepare_workspace(workspace);
        if (m_layout.atmospheric_blocks() != 0) {
            pack_atmospheric_output_blocks(m_layout, outgoing,
                                           workspace.m_atmospheric_output);
            m_basis->apply_transpose_active(
                workspace.m_atmospheric_output, m_atmospheric_coefficients,
                m_active_coefficients, workspace.m_atmospheric_input,
                workspace.m_moments);
            unpack_atmospheric_input_blocks(
                m_layout, workspace.m_atmospheric_input, incoming);
        }
        apply_dense_blocks_transpose(m_layout, m_layout.atmospheric_blocks(),
                                     m_ground_value_offsets, m_ground_values,
                                     outgoing, incoming);
    }

    void ScatteringOperator<1>::apply_jvp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
        Eigen::Ref<const Eigen::MatrixXd> coefficient_tangent,
        Eigen::Ref<const Eigen::VectorXd> ground_value_tangent,
        Eigen::Ref<Eigen::VectorXd> outgoing_tangent,
        ScatteringWorkspace<1>& workspace) const {
        validate_input_output(incoming.size(), outgoing_tangent.size());
        if (incoming_tangent.size() != incoming.size() ||
            coefficient_tangent.rows() != m_atmospheric_coefficients.rows() ||
            coefficient_tangent.cols() != m_atmospheric_coefficients.cols()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders scattering JVP sizes");
        }
        validate_dense_tangent(m_ground_values, ground_value_tangent);
        if (coefficient_tangent.isZero(0.0) &&
            ground_value_tangent.isZero(0.0)) {
            apply(incoming_tangent, outgoing_tangent, workspace);
            return;
        }
        prepare_workspace(workspace);
        if (m_layout.atmospheric_blocks() != 0) {
            pack_atmospheric_blocks(m_layout, incoming,
                                    workspace.m_atmospheric_input);
            pack_atmospheric_blocks(m_layout, incoming_tangent,
                                    workspace.m_auxiliary_input);
            int active_coefficients = m_active_coefficients;
            for (int degree = m_basis->num_coefficients() - 1;
                 degree >= m_active_coefficients; --degree) {
                if (!coefficient_tangent.col(degree).isZero(0.0)) {
                    active_coefficients = degree + 1;
                    break;
                }
            }
            m_basis->apply_jvp_active(
                workspace.m_atmospheric_input, workspace.m_auxiliary_input,
                m_atmospheric_coefficients, coefficient_tangent,
                active_coefficients, workspace.m_atmospheric_output,
                workspace.m_moments, workspace.m_auxiliary_moments);
            unpack_atmospheric_blocks(m_layout, workspace.m_atmospheric_output,
                                      outgoing_tangent);
        }
        apply_dense_blocks_jvp(m_layout, m_layout.atmospheric_blocks(),
                               m_ground_value_offsets, m_ground_values,
                               ground_value_tangent, incoming, incoming_tangent,
                               outgoing_tangent);
    }

    void ScatteringOperator<1>::apply_vjp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> outgoing_cotangent,
        Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
        Eigen::Ref<Eigen::MatrixXd> coefficient_gradient,
        Eigen::Ref<Eigen::VectorXd> ground_value_gradient,
        ScatteringWorkspace<1>& workspace) const {
        validate_input_output(incoming.size(), outgoing_cotangent.size());
        if (incoming_cotangent.size() != incoming.size() ||
            coefficient_gradient.rows() != m_atmospheric_coefficients.rows() ||
            coefficient_gradient.cols() != m_atmospheric_coefficients.cols() ||
            ground_value_gradient.size() != m_ground_values.size()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders scattering VJP sizes");
        }
        prepare_workspace(workspace);
        incoming_cotangent.setZero();
        coefficient_gradient.setZero();
        ground_value_gradient.setZero();
        if (m_layout.atmospheric_blocks() != 0) {
            pack_atmospheric_blocks(m_layout, incoming,
                                    workspace.m_atmospheric_input);
            pack_atmospheric_output_blocks(m_layout, outgoing_cotangent,
                                           workspace.m_atmospheric_output);
            m_basis->apply_vjp(
                workspace.m_atmospheric_input, m_atmospheric_coefficients,
                workspace.m_atmospheric_output, workspace.m_auxiliary_input,
                coefficient_gradient, workspace.m_moments,
                workspace.m_auxiliary_moments);
            unpack_atmospheric_input_blocks(
                m_layout, workspace.m_auxiliary_input, incoming_cotangent);
        }
        apply_dense_blocks_vjp(m_layout, m_layout.atmospheric_blocks(),
                               m_ground_value_offsets, m_ground_values,
                               incoming, outgoing_cotangent, incoming_cotangent,
                               ground_value_gradient);
    }

    void ScatteringWorkspace<3>::prepare(int atmospheric_blocks,
                                         int input_directions,
                                         int output_directions) {
        m_atmospheric_input.resize(atmospheric_blocks, 3 * input_directions);
        m_atmospheric_output.resize(atmospheric_blocks, 3 * output_directions);
        m_auxiliary_input.resize(atmospheric_blocks, 3 * input_directions);
        m_auxiliary_output.resize(atmospheric_blocks, 3 * output_directions);
    }

    std::size_t ScatteringWorkspace<3>::storage_bytes() const {
        return static_cast<std::size_t>(
                   m_atmospheric_input.size() + m_atmospheric_output.size() +
                   m_auxiliary_input.size() + m_auxiliary_output.size()) *
                   sizeof(double) +
               m_angular.storage_bytes();
    }

    ScatteringOperator<3>::ScatteringOperator(
        ScatteringBlockLayout layout,
        std::shared_ptr<const VectorAngularBasis> angular_basis)
        : m_layout(std::move(layout)), m_basis(std::move(angular_basis)),
          m_ground_value_offsets(
              make_dense_value_offsets(m_layout, m_layout.atmospheric_blocks(),
                                       m_layout.ground_blocks())),
          m_atmospheric_coefficients(
              m_layout.atmospheric_blocks(),
              m_basis == nullptr ? 0 : 4 * m_basis->num_coefficients()),
          m_ground_values(m_ground_value_offsets.back()) {
        if (m_basis == nullptr) {
            throw std::invalid_argument(
                "vector successive-orders scattering requires an angular "
                "basis");
        }
        if (m_layout.stokes_components() != 3) {
            throw std::invalid_argument(
                "vector successive-orders scattering requires NSTOKES=3");
        }
        for (int block = 0; block < m_layout.atmospheric_blocks(); ++block) {
            if (m_layout.input_block_size(block) != m_basis->input_size() ||
                m_layout.output_block_size(block) != m_basis->output_size()) {
                throw std::invalid_argument(
                    "vector successive-orders atmospheric block does not "
                    "match its angular basis");
            }
        }
        m_atmospheric_coefficients.setZero();
        m_ground_values.setZero();
    }

    void
    ScatteringOperator<3>::set_active_coefficients(int active_coefficients) {
        if (active_coefficients < 1 ||
            active_coefficients > m_basis->num_coefficients()) {
            throw std::invalid_argument(
                "invalid active vector scattering coefficient count");
        }
        m_active_coefficients = active_coefficients;
    }

    void ScatteringOperator<3>::set_atmospheric_coefficients(
        Eigen::Ref<const Eigen::MatrixXd> coefficients) {
        if (coefficients.rows() != m_atmospheric_coefficients.rows() ||
            coefficients.cols() != m_atmospheric_coefficients.cols() ||
            !coefficients.allFinite()) {
            throw std::invalid_argument(
                "invalid vector successive-orders atmospheric coefficients");
        }
        m_atmospheric_coefficients = coefficients;
        int active_coefficients = 1;
        for (int degree = m_basis->num_coefficients() - 1; degree >= 1;
             --degree) {
            if (!coefficients.middleCols(4 * degree, 4).isZero(0.0)) {
                active_coefficients = degree + 1;
                break;
            }
        }
        m_active_coefficients = active_coefficients;
    }

    void ScatteringOperator<3>::set_ground_values(
        Eigen::Ref<const Eigen::VectorXd> values) {
        if (values.size() != m_ground_values.size() || !values.allFinite()) {
            throw std::invalid_argument(
                "invalid vector successive-orders ground values");
        }
        m_ground_values = values;
    }

    int ScatteringOperator<3>::ground_value_offset(int ground_block) const {
        if (ground_block < 0 || ground_block >= m_layout.ground_blocks()) {
            throw std::out_of_range(
                "successive-orders ground scattering block out of range");
        }
        return m_ground_value_offsets[ground_block];
    }

    Eigen::Map<ScatteringOperator<3>::RowMajorMatrix>
    ScatteringOperator<3>::ground_block(int ground_block) {
        const int block = m_layout.atmospheric_blocks() + ground_block;
        return {m_ground_values.data() + ground_value_offset(ground_block),
                m_layout.output_block_size(block),
                m_layout.input_block_size(block)};
    }

    Eigen::Map<const ScatteringOperator<3>::RowMajorMatrix>
    ScatteringOperator<3>::ground_block(int ground_block) const {
        const int block = m_layout.atmospheric_blocks() + ground_block;
        return {m_ground_values.data() + ground_value_offset(ground_block),
                m_layout.output_block_size(block),
                m_layout.input_block_size(block)};
    }

    void ScatteringOperator<3>::set_ground_block(
        int ground_block, Eigen::Ref<const Eigen::MatrixXd> values) {
        auto target = this->ground_block(ground_block);
        if (values.rows() != target.rows() || values.cols() != target.cols() ||
            !values.allFinite()) {
            throw std::invalid_argument(
                "invalid vector successive-orders ground block");
        }
        target = values;
    }

    ScatteringMemoryUsage ScatteringOperator<3>::memory_usage(
        const ScatteringWorkspace<3>* workspace) const {
        ScatteringMemoryUsage result;
        result.layout_bytes = m_layout.storage_bytes() +
                              m_ground_value_offsets.capacity() * sizeof(int);
        result.angular_basis_bytes = m_basis->storage_bytes();
        result.atmospheric_value_bytes =
            static_cast<std::size_t>(m_atmospheric_coefficients.size()) *
            sizeof(double);
        result.boundary_value_bytes =
            static_cast<std::size_t>(m_ground_values.size()) * sizeof(double);
        result.workspace_bytes =
            workspace == nullptr ? 0 : workspace->storage_bytes();
        return result;
    }

    void ScatteringOperator<3>::validate_input_output(
        Eigen::Index incoming_size, Eigen::Index outgoing_size) const {
        if (incoming_size != input_size() || outgoing_size != output_size()) {
            throw std::invalid_argument(
                "invalid vector successive-orders scattering vector sizes");
        }
    }

    void ScatteringOperator<3>::prepare_workspace(
        ScatteringWorkspace<3>& workspace) const {
        workspace.prepare(m_layout.atmospheric_blocks(),
                          m_basis->input_directions(),
                          m_basis->output_directions());
    }

    ScatteringWorkspace<3> ScatteringOperator<3>::make_workspace() const {
        ScatteringWorkspace<3> workspace;
        prepare_workspace(workspace);
        return workspace;
    }

    void
    ScatteringOperator<3>::apply(Eigen::Ref<const Eigen::VectorXd> incoming,
                                 Eigen::Ref<Eigen::VectorXd> outgoing,
                                 ScatteringWorkspace<3>& workspace) const {
        validate_input_output(incoming.size(), outgoing.size());
        prepare_workspace(workspace);
        pack_atmospheric_blocks(m_layout, incoming,
                                workspace.m_atmospheric_input);
        m_basis->apply_active(workspace.m_atmospheric_input,
                              m_atmospheric_coefficients, m_active_coefficients,
                              workspace.m_atmospheric_output,
                              workspace.m_angular);
        unpack_atmospheric_blocks(m_layout, workspace.m_atmospheric_output,
                                  outgoing);
        apply_dense_blocks(m_layout, m_layout.atmospheric_blocks(),
                           m_ground_value_offsets, m_ground_values, incoming,
                           outgoing);
    }

    void ScatteringOperator<3>::apply_transpose(
        Eigen::Ref<const Eigen::VectorXd> outgoing,
        Eigen::Ref<Eigen::VectorXd> incoming,
        ScatteringWorkspace<3>& workspace) const {
        validate_input_output(incoming.size(), outgoing.size());
        prepare_workspace(workspace);
        pack_atmospheric_output_blocks(m_layout, outgoing,
                                       workspace.m_atmospheric_output);
        m_basis->apply_transpose_active(
            workspace.m_atmospheric_output, m_atmospheric_coefficients,
            m_active_coefficients, workspace.m_atmospheric_input,
            workspace.m_angular);
        unpack_atmospheric_input_blocks(m_layout, workspace.m_atmospheric_input,
                                        incoming);
        apply_dense_blocks_transpose(m_layout, m_layout.atmospheric_blocks(),
                                     m_ground_value_offsets, m_ground_values,
                                     outgoing, incoming);
    }

    void ScatteringOperator<3>::apply_jvp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
        Eigen::Ref<const Eigen::MatrixXd> coefficient_tangent,
        Eigen::Ref<const Eigen::VectorXd> ground_value_tangent,
        Eigen::Ref<Eigen::VectorXd> outgoing_tangent,
        ScatteringWorkspace<3>& workspace) const {
        validate_input_output(incoming.size(), outgoing_tangent.size());
        if (incoming_tangent.size() != incoming.size() ||
            coefficient_tangent.rows() != m_atmospheric_coefficients.rows() ||
            coefficient_tangent.cols() != m_atmospheric_coefficients.cols()) {
            throw std::invalid_argument(
                "invalid vector successive-orders scattering JVP sizes");
        }
        validate_dense_tangent(m_ground_values, ground_value_tangent);
        prepare_workspace(workspace);
        pack_atmospheric_blocks(m_layout, incoming_tangent,
                                workspace.m_auxiliary_input);
        m_basis->apply_active(workspace.m_auxiliary_input,
                              m_atmospheric_coefficients, m_active_coefficients,
                              workspace.m_atmospheric_output,
                              workspace.m_angular);
        if (!coefficient_tangent.isZero(0.0)) {
            int active_coefficients = m_active_coefficients;
            for (int degree = m_basis->num_coefficients() - 1;
                 degree >= m_active_coefficients; --degree) {
                if (!coefficient_tangent.middleCols(4 * degree, 4)
                         .isZero(0.0)) {
                    active_coefficients = degree + 1;
                    break;
                }
            }
            pack_atmospheric_blocks(m_layout, incoming,
                                    workspace.m_atmospheric_input);
            m_basis->apply_active(workspace.m_atmospheric_input,
                                  coefficient_tangent, active_coefficients,
                                  workspace.m_auxiliary_output,
                                  workspace.m_angular);
            workspace.m_atmospheric_output += workspace.m_auxiliary_output;
        }
        unpack_atmospheric_blocks(m_layout, workspace.m_atmospheric_output,
                                  outgoing_tangent);
        apply_dense_blocks_jvp(m_layout, m_layout.atmospheric_blocks(),
                               m_ground_value_offsets, m_ground_values,
                               ground_value_tangent, incoming, incoming_tangent,
                               outgoing_tangent);
    }

    void ScatteringOperator<3>::apply_vjp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> outgoing_cotangent,
        Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
        Eigen::Ref<Eigen::MatrixXd> coefficient_gradient,
        Eigen::Ref<Eigen::VectorXd> ground_value_gradient,
        ScatteringWorkspace<3>& workspace) const {
        validate_input_output(incoming.size(), outgoing_cotangent.size());
        if (incoming_cotangent.size() != incoming.size() ||
            coefficient_gradient.rows() != m_atmospheric_coefficients.rows() ||
            coefficient_gradient.cols() != m_atmospheric_coefficients.cols() ||
            ground_value_gradient.size() != m_ground_values.size()) {
            throw std::invalid_argument(
                "invalid vector successive-orders scattering VJP sizes");
        }
        incoming_cotangent.setZero();
        coefficient_gradient.setZero();
        ground_value_gradient.setZero();
        prepare_workspace(workspace);
        pack_atmospheric_blocks(m_layout, incoming,
                                workspace.m_atmospheric_input);
        pack_atmospheric_output_blocks(m_layout, outgoing_cotangent,
                                       workspace.m_atmospheric_output);
        m_basis->apply_transpose_active(
            workspace.m_atmospheric_output, m_atmospheric_coefficients,
            m_active_coefficients, workspace.m_auxiliary_input,
            workspace.m_angular);
        unpack_atmospheric_input_blocks(m_layout, workspace.m_auxiliary_input,
                                        incoming_cotangent);
        m_basis->accumulate_coefficient_vjp(
            workspace.m_atmospheric_input, workspace.m_atmospheric_output,
            coefficient_gradient, workspace.m_angular);
        apply_dense_blocks_vjp(m_layout, m_layout.atmospheric_blocks(),
                               m_ground_value_offsets, m_ground_values,
                               incoming, outgoing_cotangent, incoming_cotangent,
                               ground_value_gradient);
    }

} // namespace sasktran2::successive_orders
