#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace sasktran2::successive_orders {

    template <typename Index> class TransportTypedColumnView {
      public:
        TransportTypedColumnView(const Index* values, std::size_t size)
            : m_values(values), m_size(size) {}
        const Index* data() const { return m_values; }
        std::size_t size() const { return m_size; }
        const Index* begin() const { return m_values; }
        const Index* end() const {
            return m_size == 0 ? m_values : m_values + m_size;
        }
        int operator[](std::size_t index) const { return m_values[index]; }

      private:
        const Index* m_values;
        std::size_t m_size;
    };

    /** Non-owning columns from one immutable CSR generation. Hot consumers
     * visit once around their loop to select the actual index width. */
    class TransportColumnView {
      public:
        TransportColumnView() = default;
        TransportColumnView(const int* values, std::size_t size)
            : m_values32(values), m_size(size) {}
        TransportColumnView(const std::uint16_t* values, std::size_t size)
            : m_values16(values), m_size(size), m_compact(true) {}

        std::size_t size() const { return m_size; }
        bool empty() const { return m_size == 0; }
        bool is_compact() const { return m_compact; }
        std::size_t element_bytes() const {
            return m_compact ? sizeof(std::uint16_t) : sizeof(int);
        }
        const void* data() const {
            return m_compact ? static_cast<const void*>(m_values16)
                             : static_cast<const void*>(m_values32);
        }
        int operator[](std::size_t index) const {
            return m_compact ? m_values16[index] : m_values32[index];
        }
        TransportColumnView subview(std::size_t offset,
                                    std::size_t size) const {
            if (offset > m_size || size > m_size - offset) {
                throw std::out_of_range("invalid transport column subview");
            }
            if (m_compact) {
                return {offset == 0 ? m_values16 : m_values16 + offset, size};
            }
            return {offset == 0 ? m_values32 : m_values32 + offset, size};
        }
        template <typename Visitor>
        decltype(auto) visit(Visitor&& visitor) const {
            if (m_compact) {
                return std::forward<Visitor>(visitor)(
                    TransportTypedColumnView<std::uint16_t>(m_values16,
                                                            m_size));
            }
            return std::forward<Visitor>(visitor)(
                TransportTypedColumnView<int>(m_values32, m_size));
        }
        bool operator==(const std::vector<int>& values) const {
            if (m_size == 0) {
                return values.empty();
            }
            return values.size() == m_size && visit([&](const auto& columns) {
                       return std::equal(columns.begin(), columns.end(),
                                         values.begin());
                   });
        }
        bool operator!=(const std::vector<int>& values) const {
            return !(*this == values);
        }
        /** Compatibility copy for owning constructors during construction. */
        std::vector<int> to_vector() const {
            if (empty()) {
                return {};
            }
            return visit([](const auto& columns) {
                return std::vector<int>(columns.begin(), columns.end());
            });
        }

      private:
        const int* m_values32 = nullptr;
        const std::uint16_t* m_values16 = nullptr;
        std::size_t m_size = 0;
        bool m_compact = false;
    };

    /** Fixed CSR topology for ray transport.
     *
     * Geometry construction creates one immutable topology generation, shared
     * by geometry and its transport maps. Only the values are rebuilt for a
     * changing atmosphere. Rows are incoming angular radiances and columns are
     * outgoing source samples. Copies keep the generation alive without
     * copying either CSR array. Column indices use 16 bits when every possible
     * column is representable; wider grids retain the original int storage.
     */
    class TransportSparsity {
      public:
        TransportSparsity() : TransportSparsity(0, {0}, {}) {}
        TransportSparsity(int columns, std::vector<int> row_offsets,
                          std::vector<int> column_indices)
            : m_data(std::make_shared<const Data>(
                  columns, std::move(row_offsets), std::move(column_indices))) {
        }

        int rows() const {
            return static_cast<int>(m_data->row_offsets.size() - 1);
        }
        int columns() const { return m_data->columns; }
        int nonzeros() const {
            return static_cast<int>(column_indices().size());
        }
        const std::vector<int>& row_offsets() const {
            return m_data->row_offsets;
        }
        TransportColumnView column_indices() const {
            if (m_data->compact) {
                return {m_data->column_indices16.data(),
                        m_data->column_indices16.size()};
            }
            return {m_data->column_indices32.data(),
                    m_data->column_indices32.size()};
        }
        bool compact_column_indices() const { return m_data->compact; }
        std::size_t column_index_bytes() const {
            return m_data->column_indices16.capacity() * sizeof(std::uint16_t) +
                   m_data->column_indices32.capacity() * sizeof(int);
        }
        /** Payload size of this generation; shared handles do not add it. */
        std::size_t storage_bytes() const {
            return m_data->row_offsets.capacity() * sizeof(int) +
                   column_index_bytes();
        }

      private:
        struct Data {
            Data(int num_columns, std::vector<int> offsets,
                 std::vector<int> indices)
                : columns(num_columns), row_offsets(std::move(offsets)) {
                if (columns < 0 || row_offsets.empty() ||
                    row_offsets.size() - 1 >
                        static_cast<std::size_t>(
                            std::numeric_limits<int>::max()) ||
                    indices.size() > static_cast<std::size_t>(
                                         std::numeric_limits<int>::max()) ||
                    row_offsets.front() != 0 ||
                    row_offsets.back() != static_cast<int>(indices.size()) ||
                    !std::is_sorted(row_offsets.begin(), row_offsets.end())) {
                    throw std::invalid_argument(
                        "invalid successive-orders transport CSR offsets");
                }
                for (std::size_t row = 0; row + 1 < row_offsets.size(); ++row) {
                    const auto begin = indices.begin() + row_offsets[row];
                    const auto end = indices.begin() + row_offsets[row + 1];
                    if (!std::is_sorted(begin, end) ||
                        std::adjacent_find(begin, end) != end ||
                        std::any_of(begin, end, [&](int column) {
                            return column < 0 || column >= columns;
                        })) {
                        throw std::invalid_argument(
                            "invalid successive-orders transport CSR columns");
                    }
                }
                compact =
                    columns <= static_cast<int>(
                                   std::numeric_limits<std::uint16_t>::max()) +
                                   1;
                if (compact) {
                    column_indices16.reserve(indices.size());
                    for (const int column : indices) {
                        column_indices16.push_back(
                            static_cast<std::uint16_t>(column));
                    }
                } else {
                    column_indices32 = std::move(indices);
                }
            }

            int columns;
            std::vector<int> row_offsets;
            std::vector<std::uint16_t> column_indices16;
            std::vector<int> column_indices32;
            bool compact = false;
        };

        std::shared_ptr<const Data> m_data;
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
                for (int row = 0; row < m_sparsity->rows(); ++row) {
                    double value = 0.0;
                    for (int index = offsets[row]; index < offsets[row + 1];
                         ++index) {
                        value += m_values(index) * state(columns[index]);
                    }
                    incoming(row) = value;
                }
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
                for (int row = 0; row < m_sparsity->rows(); ++row) {
                    for (int index = offsets[row]; index < offsets[row + 1];
                         ++index) {
                        const int column = columns[index];
                        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                            incoming(row * NSTOKES + stokes) +=
                                m_values(index) *
                                state(column * NSTOKES + stokes);
                        }
                    }
                }
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
                for (int row = 0; row < m_sparsity->rows(); ++row) {
                    const double row_value = incoming(row);
                    for (int index = offsets[row]; index < offsets[row + 1];
                         ++index) {
                        state(columns[index]) += m_values(index) * row_value;
                    }
                }
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
                for (int row = 0; row < m_sparsity->rows(); ++row) {
                    for (int index = offsets[row]; index < offsets[row + 1];
                         ++index) {
                        const int column = columns[index];
                        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                            state(column * NSTOKES + stokes) +=
                                m_values(index) *
                                incoming(row * NSTOKES + stokes);
                        }
                    }
                }
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
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
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
                for (int row = 0; row < m_sparsity->rows(); ++row) {
                    double value = 0.0;
                    for (int index = offsets[row]; index < offsets[row + 1];
                         ++index) {
                        value += value_tangent(index) * state(columns[index]);
                    }
                    incoming_tangent(row) = value;
                }
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
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
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
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
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
                for (int row = 0; row < m_sparsity->rows(); ++row) {
                    const double row_cotangent = incoming_cotangent(row);
                    for (int index = offsets[row]; index < offsets[row + 1];
                         ++index) {
                        const int column = columns[index];
                        state_cotangent(column) +=
                            m_values(index) * row_cotangent;
                        value_gradient(index) = state(column) * row_cotangent;
                    }
                }
            });
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
            m_sparsity->column_indices().visit([&](const auto& columns) {
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
            });
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
