#pragma once

#include "interpolation.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

namespace sasktran2::successive_orders {

    /** Four original endpoint weights on one verified structured cell. */
    class StructuredEndpointStencilView {
      public:
        StructuredEndpointStencilView(int base, int altitude_stride,
                                      const std::array<double, 4>& weights)
            : m_base(base), m_altitude_stride(altitude_stride),
              m_weights(&weights) {}

        std::size_t size() const { return 4; }
        std::pair<int, double> operator[](std::size_t index) const {
            if (index >= 4) {
                throw std::out_of_range(
                    "Endpoint stencil index is out of range");
            }
            return {m_base + static_cast<int>(index & 1) +
                        static_cast<int>(index >> 1) * m_altitude_stride,
                    (*m_weights)[index]};
        }

      private:
        int m_base;
        int m_altitude_stride;
        const std::array<double, 4>* m_weights;
    };

    /** Dispatch selects the endpoint representation before its weight loop. */
    class EndpointStencilView {
      public:
        using GenericView = InterpolationView<InterpolationWeight>;

        explicit EndpointStencilView(GenericView values) : m_values(values) {}
        explicit EndpointStencilView(StructuredEndpointStencilView values)
            : m_values(values) {}

        template <typename Callback>
        decltype(auto) visit(Callback&& callback) const {
            return std::visit(std::forward<Callback>(callback), m_values);
        }

      private:
        std::variant<GenericView, StructuredEndpointStencilView> m_values;
    };

    /** Immutable endpoint stencils with exact structured-grid compaction.
     *
     * Uniform four-corner stencils keep aligned, original double weights and
     * one cell base per stencil. Corner indices and offsets are implicit.
     * Other geometries retain their original indices, offsets, and weights.
     */
    class EndpointStencilStorage {
      public:
        struct GenericStorage {
            std::vector<int> offsets;
            std::vector<InterpolationWeight> weights;

            std::size_t size() const {
                return offsets.empty() ? 0 : offsets.size() - 1;
            }
            InterpolationView<InterpolationWeight>
            operator[](std::size_t slot) const {
                return {weights, static_cast<std::size_t>(offsets[slot]),
                        static_cast<std::size_t>(offsets[slot + 1] -
                                                 offsets[slot])};
            }
            std::size_t storage_bytes() const {
                return offsets.capacity() * sizeof(int) +
                       weights.capacity() * sizeof(InterpolationWeight);
            }
        };

        template <typename Base> struct StructuredStorage {
            int altitude_stride = 0;
            std::vector<Base> bases;
            std::vector<std::array<double, 4>> weights;

            std::size_t size() const { return bases.size(); }
            StructuredEndpointStencilView operator[](std::size_t slot) const {
                return {static_cast<int>(bases[slot]), altitude_stride,
                        weights[slot]};
            }
            std::size_t storage_bytes() const {
                return bases.capacity() * sizeof(Base) +
                       weights.capacity() * sizeof(std::array<double, 4>);
            }
        };

        template <typename Callback>
        decltype(auto) visit(Callback&& callback) const {
            return std::visit(std::forward<Callback>(callback), m_values);
        }

        void clear() { m_values = GenericStorage{}; }
        std::size_t size() const {
            return visit([](const auto& values) { return values.size(); });
        }
        bool empty() const { return size() == 0; }
        bool is_structured() const { return m_values.index() != 0; }
        int base_bits() const {
            return m_values.index() == 1 ? 16 : m_values.index() == 2 ? 32 : 0;
        }
        std::size_t storage_bytes() const {
            return visit(
                [](const auto& values) { return values.storage_bytes(); });
        }
        std::size_t weight_count() const {
            if (const auto* generic = std::get_if<GenericStorage>(&m_values)) {
                return generic->weights.size();
            }
            return 4 * size();
        }

        EndpointStencilView view(std::size_t slot) const {
            return visit([slot](const auto& values) {
                if (slot >= values.size()) {
                    throw std::out_of_range(
                        "Endpoint stencil slot is out of range");
                }
                return EndpointStencilView(values[slot]);
            });
        }

        void assign(std::vector<int>&& offsets,
                    std::vector<InterpolationWeight>&& weights,
                    int altitude_stride, int geometry_size) {
            if (offsets.empty()) {
                if (!weights.empty()) {
                    throw std::invalid_argument(
                        "Endpoint stencils require their original offsets");
                }
                clear();
                return;
            }
            if (offsets.front() != 0 || offsets.back() < 0 ||
                static_cast<std::size_t>(offsets.back()) != weights.size()) {
                throw std::invalid_argument("Invalid endpoint stencil offsets");
            }
            bool structured = altitude_stride > 0 &&
                              geometry_size > altitude_stride &&
                              geometry_size % altitude_stride == 0;
            bool narrow = true;
            for (std::size_t slot = 0; slot + 1 < offsets.size(); ++slot) {
                const int begin = offsets[slot];
                const int end = offsets[slot + 1];
                if (begin < 0 || end < begin ||
                    static_cast<std::size_t>(end) > weights.size()) {
                    throw std::invalid_argument(
                        "Invalid endpoint stencil offsets");
                }
                if (!structured || end - begin != 4) {
                    structured = false;
                    continue;
                }
                const int base = weights[begin].index;
                const bool corner =
                    base >= 0 && base < geometry_size &&
                    base / altitude_stride + 1 <
                        geometry_size / altitude_stride &&
                    base % altitude_stride + 1 < altitude_stride &&
                    weights[begin + 1].index == base + 1 &&
                    weights[begin + 2].index == base + altitude_stride &&
                    weights[begin + 3].index == base + altitude_stride + 1;
                structured = structured && corner;
                narrow =
                    narrow && base <= std::numeric_limits<std::uint16_t>::max();
            }
            if (!structured) {
                m_values =
                    GenericStorage{std::move(offsets), std::move(weights)};
                return;
            }
            if (narrow) {
                assign_structured<std::uint16_t>(offsets, weights,
                                                 altitude_stride);
            } else {
                assign_structured<std::uint32_t>(offsets, weights,
                                                 altitude_stride);
            }
            // The original construction arrays are no longer needed.
            std::vector<int>().swap(offsets);
            std::vector<InterpolationWeight>().swap(weights);
        }

      private:
        template <typename Base>
        void assign_structured(const std::vector<int>& offsets,
                               const std::vector<InterpolationWeight>& weights,
                               int altitude_stride) {
            StructuredStorage<Base> values;
            values.altitude_stride = altitude_stride;
            const std::size_t stencils = offsets.size() - 1;
            values.bases.resize(stencils);
            values.weights.resize(stencils);
            for (std::size_t slot = 0; slot < stencils; ++slot) {
                const int begin = offsets[slot];
                values.bases[slot] = static_cast<Base>(weights[begin].index);
                for (int corner = 0; corner < 4; ++corner) {
                    values.weights[slot][corner] =
                        weights[begin + corner].weight();
                }
            }
            m_values = std::move(values);
        }

        std::variant<GenericStorage, StructuredStorage<std::uint16_t>,
                     StructuredStorage<std::uint32_t>>
            m_values;
    };

} // namespace sasktran2::successive_orders
