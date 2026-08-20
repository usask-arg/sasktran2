#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace sasktran2::successive_orders {

    /** Fixed CSR topology for ray transport.
     *
     * Geometry construction owns the topology; only the values are rebuilt
     * for a changing atmosphere. Rows are incoming angular radiances and
     * columns are outgoing source samples.
     */
    class TransportSparsity {
      public:
        TransportSparsity() = default;
        TransportSparsity(int columns, std::vector<int> row_offsets,
                          std::vector<int> column_indices)
            : m_columns(columns), m_row_offsets(std::move(row_offsets)),
              m_column_indices(std::move(column_indices)) {
            validate();
        }

        int rows() const {
            return m_row_offsets.empty()
                       ? 0
                       : static_cast<int>(m_row_offsets.size()) - 1;
        }
        int columns() const { return m_columns; }
        int nonzeros() const {
            return static_cast<int>(m_column_indices.size());
        }
        const std::vector<int>& row_offsets() const { return m_row_offsets; }
        const std::vector<int>& column_indices() const {
            return m_column_indices;
        }
        std::size_t storage_bytes() const {
            return m_row_offsets.capacity() * sizeof(int) +
                   m_column_indices.capacity() * sizeof(int);
        }

      private:
        void validate() const {
            if (m_columns < 0 || m_row_offsets.empty() ||
                m_row_offsets.front() != 0 ||
                m_row_offsets.back() !=
                    static_cast<int>(m_column_indices.size()) ||
                !std::is_sorted(m_row_offsets.begin(), m_row_offsets.end())) {
                throw std::invalid_argument(
                    "invalid successive-orders transport CSR offsets");
            }
            for (int row = 0; row < rows(); ++row) {
                const auto begin =
                    m_column_indices.begin() + m_row_offsets[row];
                const auto end =
                    m_column_indices.begin() + m_row_offsets[row + 1];
                if (!std::is_sorted(begin, end) ||
                    std::adjacent_find(begin, end) != end ||
                    std::any_of(begin, end, [&](int column) {
                        return column < 0 || column >= m_columns;
                    })) {
                    throw std::invalid_argument(
                        "invalid successive-orders transport CSR columns");
                }
            }
        }

        int m_columns = 0;
        std::vector<int> m_row_offsets{0};
        std::vector<int> m_column_indices;
    };

    class TransportOperator {
      public:
        explicit TransportOperator(const TransportSparsity& sparsity)
            : m_sparsity(&sparsity),
              m_values(Eigen::VectorXd::Zero(sparsity.nonzeros())) {}

        const TransportSparsity& sparsity() const { return *m_sparsity; }
        Eigen::VectorXd& values() { return m_values; }
        const Eigen::VectorXd& values() const { return m_values; }

        std::size_t storage_bytes() const {
            return static_cast<std::size_t>(m_values.size()) * sizeof(double);
        }

        void apply(Eigen::Ref<const Eigen::VectorXd> state,
                   Eigen::Ref<Eigen::VectorXd> incoming) const {
            validate_vectors(state, incoming);
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                double value = 0.0;
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    value += m_values(index) * state(columns[index]);
                }
                incoming(row) = value;
            }
        }

        /** Applies one geometry operator to interleaved Stokes channels
         * without duplicating CSR values or column indices. */
        template <int NSTOKES>
        void apply_stokes(Eigen::Ref<const Eigen::VectorXd> state,
                          Eigen::Ref<Eigen::VectorXd> incoming) const {
            if constexpr (NSTOKES == 1) {
                apply(state, incoming);
                return;
            }
            validate_stokes_vectors<NSTOKES>(state, incoming);
            incoming.setZero();
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    const int column = columns[index];
                    for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                        incoming(row * NSTOKES + stokes) +=
                            m_values(index) * state(column * NSTOKES + stokes);
                    }
                }
            }
        }

        void apply_affine(Eigen::Ref<const Eigen::VectorXd> state,
                          Eigen::Ref<const Eigen::VectorXd> forcing,
                          Eigen::Ref<Eigen::VectorXd> incoming) const {
            if (forcing.size() != m_sparsity->rows()) {
                throw std::invalid_argument(
                    "invalid successive-orders transport forcing size");
            }
            apply(state, incoming);
            incoming += forcing;
        }

        void apply_transpose(Eigen::Ref<const Eigen::VectorXd> incoming,
                             Eigen::Ref<Eigen::VectorXd> state) const {
            if (incoming.size() != m_sparsity->rows() ||
                state.size() != m_sparsity->columns()) {
                throw std::invalid_argument(
                    "invalid successive-orders transpose transport sizes");
            }
            state.setZero();
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                const double row_value = incoming(row);
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    state(columns[index]) += m_values(index) * row_value;
                }
            }
        }

        template <int NSTOKES>
        void apply_transpose_stokes(Eigen::Ref<const Eigen::VectorXd> incoming,
                                    Eigen::Ref<Eigen::VectorXd> state) const {
            if constexpr (NSTOKES == 1) {
                apply_transpose(incoming, state);
                return;
            }
            validate_stokes_vectors<NSTOKES>(state, incoming);
            state.setZero();
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    const int column = columns[index];
                    for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                        state(column * NSTOKES + stokes) +=
                            m_values(index) * incoming(row * NSTOKES + stokes);
                    }
                }
            }
        }

        void apply_jvp(Eigen::Ref<const Eigen::VectorXd> state,
                       Eigen::Ref<const Eigen::VectorXd> state_tangent,
                       Eigen::Ref<const Eigen::VectorXd> value_tangent,
                       Eigen::Ref<Eigen::VectorXd> incoming_tangent) const {
            if (state_tangent.size() != state.size() ||
                value_tangent.size() != m_values.size()) {
                throw std::invalid_argument(
                    "invalid successive-orders transport JVP sizes");
            }
            validate_vectors(state, incoming_tangent);
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                double value = 0.0;
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    const int column = columns[index];
                    value += m_values(index) * state_tangent(column) +
                             value_tangent(index) * state(column);
                }
                incoming_tangent(row) = value;
            }
        }

        /** Applies only the direct value derivative, dT * state. */
        void
        apply_value_jvp(Eigen::Ref<const Eigen::VectorXd> state,
                        Eigen::Ref<const Eigen::VectorXd> value_tangent,
                        Eigen::Ref<Eigen::VectorXd> incoming_tangent) const {
            if (value_tangent.size() != m_values.size()) {
                throw std::invalid_argument(
                    "invalid successive-orders transport value JVP size");
            }
            validate_vectors(state, incoming_tangent);
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                double value = 0.0;
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    value += value_tangent(index) * state(columns[index]);
                }
                incoming_tangent(row) = value;
            }
        }

        template <int NSTOKES>
        void apply_value_jvp_stokes(
            Eigen::Ref<const Eigen::VectorXd> state,
            Eigen::Ref<const Eigen::VectorXd> value_tangent,
            Eigen::Ref<Eigen::VectorXd> incoming_tangent) const {
            if constexpr (NSTOKES == 1) {
                apply_value_jvp(state, value_tangent, incoming_tangent);
                return;
            }
            validate_stokes_vectors<NSTOKES>(state, incoming_tangent);
            if (value_tangent.size() != m_values.size()) {
                throw std::invalid_argument(
                    "invalid successive-orders Stokes transport value JVP "
                    "size");
            }
            incoming_tangent.setZero();
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    const int column = columns[index];
                    for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                        incoming_tangent(row * NSTOKES + stokes) +=
                            value_tangent(index) *
                            state(column * NSTOKES + stokes);
                    }
                }
            }
        }

        template <int NSTOKES>
        void
        apply_jvp_stokes(Eigen::Ref<const Eigen::VectorXd> state,
                         Eigen::Ref<const Eigen::VectorXd> state_tangent,
                         Eigen::Ref<const Eigen::VectorXd> value_tangent,
                         Eigen::Ref<Eigen::VectorXd> incoming_tangent) const {
            if constexpr (NSTOKES == 1) {
                apply_jvp(state, state_tangent, value_tangent,
                          incoming_tangent);
                return;
            }
            validate_stokes_vectors<NSTOKES>(state, incoming_tangent);
            if (state_tangent.size() != state.size() ||
                value_tangent.size() != m_values.size()) {
                throw std::invalid_argument(
                    "invalid successive-orders Stokes transport JVP sizes");
            }
            incoming_tangent.setZero();
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    const int column = columns[index];
                    for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                        incoming_tangent(row * NSTOKES + stokes) +=
                            m_values(index) *
                                state_tangent(column * NSTOKES + stokes) +
                            value_tangent(index) *
                                state(column * NSTOKES + stokes);
                    }
                }
            }
        }

        void apply_vjp(Eigen::Ref<const Eigen::VectorXd> state,
                       Eigen::Ref<const Eigen::VectorXd> incoming_cotangent,
                       Eigen::Ref<Eigen::VectorXd> state_cotangent,
                       Eigen::Ref<Eigen::VectorXd> value_gradient) const {
            if (state.size() != m_sparsity->columns() ||
                incoming_cotangent.size() != m_sparsity->rows() ||
                state_cotangent.size() != m_sparsity->columns() ||
                value_gradient.size() != m_values.size()) {
                throw std::invalid_argument(
                    "invalid successive-orders transport VJP sizes");
            }
            state_cotangent.setZero();
            value_gradient.setZero();
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                const double row_cotangent = incoming_cotangent(row);
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    const int column = columns[index];
                    state_cotangent(column) += m_values(index) * row_cotangent;
                    value_gradient(index) = state(column) * row_cotangent;
                }
            }
        }

        template <int NSTOKES>
        void
        apply_vjp_stokes(Eigen::Ref<const Eigen::VectorXd> state,
                         Eigen::Ref<const Eigen::VectorXd> incoming_cotangent,
                         Eigen::Ref<Eigen::VectorXd> state_cotangent,
                         Eigen::Ref<Eigen::VectorXd> value_gradient) const {
            if constexpr (NSTOKES == 1) {
                apply_vjp(state, incoming_cotangent, state_cotangent,
                          value_gradient);
                return;
            }
            validate_stokes_vectors<NSTOKES>(state, incoming_cotangent);
            if (state_cotangent.size() != state.size() ||
                value_gradient.size() != m_values.size()) {
                throw std::invalid_argument(
                    "invalid successive-orders Stokes transport VJP sizes");
            }
            state_cotangent.setZero();
            value_gradient.setZero();
            const auto& offsets = m_sparsity->row_offsets();
            const auto& columns = m_sparsity->column_indices();
            for (int row = 0; row < m_sparsity->rows(); ++row) {
                for (int index = offsets[row]; index < offsets[row + 1];
                     ++index) {
                    const int column = columns[index];
                    for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                        const double cotangent =
                            incoming_cotangent(row * NSTOKES + stokes);
                        state_cotangent(column * NSTOKES + stokes) +=
                            m_values(index) * cotangent;
                        value_gradient(index) +=
                            state(column * NSTOKES + stokes) * cotangent;
                    }
                }
            }
        }

      private:
        void validate_vectors(Eigen::Ref<const Eigen::VectorXd> state,
                              Eigen::Ref<Eigen::VectorXd> incoming) const {
            if (state.size() != m_sparsity->columns() ||
                incoming.size() != m_sparsity->rows()) {
                throw std::invalid_argument(
                    "invalid successive-orders transport vector sizes");
            }
        }

        template <int NSTOKES>
        void validate_stokes_vectors(
            Eigen::Ref<const Eigen::VectorXd> state,
            Eigen::Ref<const Eigen::VectorXd> incoming) const {
            static_assert(NSTOKES > 0);
            if (state.size() != m_sparsity->columns() * NSTOKES ||
                incoming.size() != m_sparsity->rows() * NSTOKES) {
                throw std::invalid_argument(
                    "invalid successive-orders Stokes transport vector sizes");
            }
        }

        const TransportSparsity* m_sparsity;
        Eigen::VectorXd m_values;
    };

} // namespace sasktran2::successive_orders
