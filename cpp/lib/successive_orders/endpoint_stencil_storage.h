#pragma once

#include "interpolation.h"

#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

namespace sasktran2::successive_orders {

    /** Preserve the ray tracer's six separate floating-point operations. */
    inline std::array<double, 4>
    expand_endpoint_factors(const std::array<double, 2>& factors) {
        const double altitude_lower = 1.0 - factors[0];
        const double horizontal_lower = 1.0 - factors[1];
        return {horizontal_lower * altitude_lower,
                horizontal_lower * factors[0], factors[1] * altitude_lower,
                factors[1] * factors[0]};
    }

    /** Cold recovery executes in the selected original floating environment.
     * Keep this strict leaf non-template: Clang does not consistently retain
     * lexical FENV_ACCESS in an instantiated template's floating operations.
     */
    inline std::array<double, 4> expand_endpoint_factors_in_environment(
        const std::array<double, 2>& factors) {
#if defined(__clang__)
#pragma STDC FENV_ACCESS ON
        const double altitude_lower = 1.0 - factors[0];
        const double horizontal_lower = 1.0 - factors[1];
        return {horizontal_lower * altitude_lower,
                horizontal_lower * factors[0], factors[1] * altitude_lower,
                factors[1] * factors[0]};
#else
        volatile double altitude = factors[0];
        volatile double horizontal = factors[1];
        volatile double altitude_lower = 1.0 - altitude;
        volatile double horizontal_lower = 1.0 - horizontal;
        volatile double weight0 = horizontal_lower * altitude_lower;
        volatile double weight1 = horizontal_lower * altitude;
        volatile double weight2 = horizontal * altitude_lower;
        volatile double weight3 = horizontal * altitude;
        return {weight0, weight1, weight2, weight3};
#endif
    }

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

    /** Expanded coefficients belong to this view, including after a move. */
    class FactoredEndpointStencilView {
      public:
        FactoredEndpointStencilView(int base, int altitude_stride,
                                    const std::array<double, 2>& factors)
            : m_base(base), m_altitude_stride(altitude_stride),
              m_weights(expand_endpoint_factors(factors)) {}

        std::size_t size() const { return 4; }
        std::pair<int, double> operator[](std::size_t index) const {
            if (index >= 4) {
                throw std::out_of_range(
                    "Endpoint stencil index is out of range");
            }
            return {m_base + static_cast<int>(index & 1) +
                        static_cast<int>(index >> 1) * m_altitude_stride,
                    m_weights[index]};
        }

      private:
        int m_base;
        int m_altitude_stride;
        std::array<double, 4> m_weights;
    };

    /** Dispatch selects the endpoint representation before its weight loop. */
    class EndpointStencilView {
      public:
        using GenericView = InterpolationView<InterpolationWeight>;

        explicit EndpointStencilView(GenericView values) : m_values(values) {}
        explicit EndpointStencilView(StructuredEndpointStencilView values)
            : m_values(values) {}
        explicit EndpointStencilView(FactoredEndpointStencilView values)
            : m_values(values) {}

        template <typename Callback>
        decltype(auto) visit(Callback&& callback) const {
            // Keep hot endpoint dispatch direct as representations are added.
            if (const auto* generic = std::get_if<GenericView>(&m_values)) {
                return callback(*generic);
            }
            if (const auto* structured =
                    std::get_if<StructuredEndpointStencilView>(&m_values)) {
                return callback(*structured);
            }
            return callback(std::get<FactoredEndpointStencilView>(m_values));
        }

      private:
        std::variant<GenericView, StructuredEndpointStencilView,
                     FactoredEndpointStencilView>
            m_values;
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

        template <typename Base> struct FactoredStorage {
            int altitude_stride = 0;
            int rounding_mode = FE_TONEAREST;
            std::vector<Base> bases;
            std::vector<std::array<double, 2>> factors;

            std::size_t size() const { return bases.size(); }
            FactoredEndpointStencilView operator[](std::size_t slot) const {
                return {static_cast<int>(bases[slot]), altitude_stride,
                        factors[slot]};
            }
            std::size_t storage_bytes() const {
                return bases.capacity() * sizeof(Base) +
                       factors.capacity() * sizeof(std::array<double, 2>);
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
        bool is_factored() const { return m_values.index() >= 3; }
        int base_bits() const {
            return m_values.index() == 1 || m_values.index() == 3   ? 16
                   : m_values.index() == 2 || m_values.index() == 4 ? 32
                                                                    : 0;
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

        template <typename Callback>
        decltype(auto) visit_view(std::size_t slot, Callback&& callback) const {
            // Pass the concrete view directly into its consumer. In particular,
            // a factored view owns its expanded weights throughout the call.
            if (const auto* values =
                    std::get_if<FactoredStorage<std::uint16_t>>(&m_values)) {
                return visit_checked_view(*values, slot,
                                          std::forward<Callback>(callback));
            }
            if (const auto* values =
                    std::get_if<FactoredStorage<std::uint32_t>>(&m_values)) {
                return visit_checked_view(*values, slot,
                                          std::forward<Callback>(callback));
            }
            if (const auto* values =
                    std::get_if<StructuredStorage<std::uint16_t>>(&m_values)) {
                return visit_checked_view(*values, slot,
                                          std::forward<Callback>(callback));
            }
            if (const auto* values =
                    std::get_if<StructuredStorage<std::uint32_t>>(&m_values)) {
                return visit_checked_view(*values, slot,
                                          std::forward<Callback>(callback));
            }
            return visit_checked_view(std::get<GenericStorage>(m_values), slot,
                                      std::forward<Callback>(callback));
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

        void
        assign(std::vector<int>&& offsets,
               std::vector<InterpolationWeight>&& weights, int altitude_stride,
               int geometry_size,
               std::vector<std::array<double, 2>>&& original_factors = {}) {
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
            const int rounding_mode = std::fegetround();
            bool factored = original_factors.size() == offsets.size() - 1 &&
                            !original_factors.empty() && rounding_mode != -1;
            // Adoption is whole-provider and checks original coefficient bits;
            // factors are never inferred from already rounded coefficients.
            for (std::size_t slot = 0;
                 factored && slot < original_factors.size(); ++slot) {
                const auto& factors = original_factors[slot];
                if (!std::isfinite(factors[0]) || !std::isfinite(factors[1]) ||
                    factors[0] < 0.0 || factors[0] > 1.0 || factors[1] < 0.0 ||
                    factors[1] > 1.0) {
                    factored = false;
                    break;
                }
                const auto expanded = expand_endpoint_factors(factors);
                for (int corner = 0; corner < 4; ++corner) {
                    const double original =
                        weights[offsets[slot] + corner].weight();
                    if (std::memcmp(&expanded[corner], &original,
                                    sizeof(double)) != 0) {
                        factored = false;
                        break;
                    }
                }
            }
            if (narrow) {
                if (factored) {
                    assign_factored<std::uint16_t>(
                        offsets, weights, altitude_stride, rounding_mode,
                        std::move(original_factors));
                } else {
                    assign_structured<std::uint16_t>(offsets, weights,
                                                     altitude_stride);
                }
            } else {
                if (factored) {
                    assign_factored<std::uint32_t>(
                        offsets, weights, altitude_stride, rounding_mode,
                        std::move(original_factors));
                } else {
                    assign_structured<std::uint32_t>(offsets, weights,
                                                     altitude_stride);
                }
            }
            // The original construction arrays are no longer needed.
            std::vector<int>().swap(offsets);
            std::vector<InterpolationWeight>().swap(weights);
            std::vector<std::array<double, 2>>().swap(original_factors);
        }

        /** Call once before a bulk operation, while no views are outstanding.
         *
         * Factored providers run on one calling thread. If that thread's
         * rounding mode changes, recover the verified original coefficients
         * once and permanently use four-double storage. Physical calculations
         * continue in the caller's rounding mode, unchanged.
         */
        void prepare_for_current_rounding() {
            if (!is_factored()) {
                return;
            }
            const int current_mode = std::fegetround();
            if (const auto* values =
                    std::get_if<FactoredStorage<std::uint16_t>>(&m_values)) {
                if (current_mode != values->rounding_mode) {
                    materialize(*values, current_mode);
                }
            } else {
                const auto& values32 =
                    std::get<FactoredStorage<std::uint32_t>>(m_values);
                if (current_mode != values32.rounding_mode) {
                    materialize(values32, current_mode);
                }
            }
        }

      private:
        template <typename Storage, typename Callback>
        static decltype(auto) visit_checked_view(const Storage& values,
                                                 std::size_t slot,
                                                 Callback&& callback) {
            if (slot >= values.size()) {
                throw std::out_of_range(
                    "Endpoint stencil slot is out of range");
            }
            return std::forward<Callback>(callback)(values[slot]);
        }

        template <typename Base>
        void
        assign_factored(const std::vector<int>& offsets,
                        const std::vector<InterpolationWeight>& weights,
                        int altitude_stride, int rounding_mode,
                        std::vector<std::array<double, 2>>&& original_factors) {
            FactoredStorage<Base> values;
            values.altitude_stride = altitude_stride;
            values.rounding_mode = rounding_mode;
            values.bases.resize(offsets.size() - 1);
            for (std::size_t slot = 0; slot < values.bases.size(); ++slot) {
                values.bases[slot] =
                    static_cast<Base>(weights[offsets[slot]].index);
            }
            values.factors = std::move(original_factors);
            m_values = std::move(values);
        }

        template <typename Base>
        void materialize(const FactoredStorage<Base>& original,
                         int current_mode) {
            StructuredStorage<Base> values;
            values.altitude_stride = original.altitude_stride;
            values.bases = original.bases;
            values.weights.resize(original.size());
            ScopedRoundingMode rounding(original.rounding_mode, current_mode);
            // Each strict leaf finishes before restoring the caller's mode.
            // The usual hot reconstruction uses the existing compiler rules.
            for (std::size_t slot = 0; slot < original.size(); ++slot) {
                values.weights[slot] = expand_endpoint_factors_in_environment(
                    original.factors[slot]);
            }
            rounding.restore();
            m_values = std::move(values);
        }

        class ScopedRoundingMode {
          public:
            ScopedRoundingMode(int target, int original)
                : m_original(original) {
                if (original == -1 || std::fesetround(target) != 0) {
                    throw std::runtime_error(
                        "Cannot recover endpoint interpolation rounding mode");
                }
            }
            ~ScopedRoundingMode() {
                if (m_original != -1) {
                    std::fesetround(m_original);
                }
            }
            void restore() {
                if (std::fesetround(m_original) != 0) {
                    throw std::runtime_error(
                        "Cannot restore endpoint interpolation rounding mode");
                }
                m_original = -1;
            }

          private:
            int m_original;
        };

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
                     StructuredStorage<std::uint32_t>,
                     FactoredStorage<std::uint16_t>,
                     FactoredStorage<std::uint32_t>>
            m_values;
    };

} // namespace sasktran2::successive_orders
