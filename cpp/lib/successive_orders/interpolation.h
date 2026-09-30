#pragma once

#include "source_weight_storage.h"

#include <sasktran2/raytracing.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace sasktran2 {
    class Geometry;
}

namespace sasktran2::grids {
    class SourceLocationInterpolator;
}

namespace sasktran2::successive_orders {
    class SourcePoint;

    /** Unit tangent direction toward the local solar azimuth.
     *
     * The coordinate basis supplies a deterministic tangent fallback at exact
     * solar zenith and nadir. The construction is equivariant under rigid
     * rotations of the complete geometry.
     */
    Eigen::Vector3d
    solar_horizontal_reference(const Eigen::Vector3d& local_up,
                               const sasktran2::Geometry& geometry);

    /** Small C++17-compatible view over immutable packed interpolation data. */
    template <typename T> class InterpolationView {
      public:
        InterpolationView() = default;
        explicit InterpolationView(const std::vector<T>& values)
            : InterpolationView(values, 0, values.size()) {}
        InterpolationView(const std::vector<T>& values, std::size_t offset,
                          std::size_t size) {
            if (offset > values.size() || size > values.size() - offset) {
                throw std::out_of_range(
                    "Successive-orders interpolation view is out of range");
            }
            m_data = size == 0 ? nullptr : values.data() + offset;
            m_size = size;
        }

        const T* data() const { return m_data; }
        const T* begin() const { return m_data; }
        const T* end() const {
            return m_data == nullptr ? nullptr : m_data + m_size;
        }
        std::size_t size() const { return m_size; }
        bool empty() const { return m_size == 0; }
        const T& operator[](std::size_t index) const {
            if (index >= m_size) {
                throw std::out_of_range(
                    "Successive-orders interpolation index is out of range");
            }
            return m_data[index];
        }

      private:
        const T* m_data = nullptr;
        std::size_t m_size = 0;
    };

    /** One geometry-only interpolation coefficient.
     *
     * Store the double's complete representation in byte storage so the
     * record has no alignment padding. memcpy gives aligned double values to
     * callers and preserves their bits on platforms requiring aligned loads.
     */
    struct InterpolationWeight {
        int index = 0;

        InterpolationWeight() { set_weight(0.0); }
        InterpolationWeight(int grid_index, double value) : index(grid_index) {
            set_weight(value);
        }

        double weight() const {
            double value;
            std::memcpy(&value, m_weight.data(), sizeof(value));
            return value;
        }
        void set_weight(double value) {
            std::memcpy(m_weight.data(), &value, sizeof(value));
        }
        void add_weight(double value) { set_weight(weight() + value); }

      private:
        std::array<std::byte, sizeof(double)> m_weight;
    };

    static_assert(sizeof(double) == 8,
                  "Successive-orders interpolation requires 64-bit doubles");
    static_assert(sizeof(InterpolationWeight) == 12,
                  "Successive-orders interpolation weights must stay compact");

    /** Interpolation from the global outgoing-source vector.
     *
     * The index names a global outgoing source during construction. Once the
     * owning ray's transport row is compiled, it names that source's local
     * CSR slot instead. Global indices are then read from the shared columns.
     * No second index or floating-point alignment padding is retained.
     */
    struct SourceInterpolationWeight {
        SourceInterpolationWeight() = default;
        SourceInterpolationWeight(int index, double interpolation_weight)
            : m_entry(index, interpolation_weight) {}

        /** Valid only before the owning ray's transport row is compiled. */
        int source_index() const { return m_entry.index; }
        /** Valid only after the owning ray's transport row is compiled. */
        std::uint32_t row_inner_index() const {
            return static_cast<std::uint32_t>(m_entry.index);
        }
        void set_row_inner_index(std::uint32_t slot) {
            if (slot >
                static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
                throw std::length_error(
                    "Successive-orders source slot exceeds its index range");
            }
            m_entry.index = static_cast<int>(slot);
        }
        double weight() const { return m_entry.weight(); }
        void add_weight(double value) { m_entry.add_weight(value); }

      private:
        InterpolationWeight m_entry;
    };

    static_assert(sizeof(SourceInterpolationWeight) == 12,
                  "Successive-orders source weights must stay compact");

    using SourceInterpolationStorage =
        SourceWeightStorage<SourceInterpolationWeight>;
    using SourceInterpolationView = SourceWeightView<SourceInterpolationWeight>;

    /** Geometry metadata required to integrate one traced layer. */
    struct LayerInterpolation {
        std::uint32_t atmosphere_offset = 0;
        std::uint32_t atmosphere_count = 0;
        std::uint32_t source_offset = 0;
        std::uint32_t source_count = 0;
        std::uint32_t optical_depth_offset = 0;
        std::uint32_t optical_depth_count = 0;
    };

    /** A verified structured 2D cell, with double coefficients stored
     * separately.
     *
     * OD corners are base, base+1, base+altitude_stride,
     * base+altitude_stride+1. The midpoint mask selects an ordered subset of
     * those corners. Source counts are derived from the following layer's
     * offset or the array end.
     */
    struct StructuredLayerInterpolation {
        std::uint16_t atmosphere_offset = 0;
        std::uint16_t source_offset = 0;
        std::uint16_t cell_base = 0;
        std::uint8_t atmosphere_mask = 0;
        std::uint8_t reserved = 0;
    };

    static_assert(sizeof(StructuredLayerInterpolation) == 8,
                  "Structured layer descriptors must stay compact");

    /** A structured cell whose midpoint offset fits twelve bits.
     *
     * The source offset and cell base keep their full sixteen-bit ranges. The
     * four corner-mask bits share a naturally aligned uint16_t with the
     * midpoint offset; no unaligned loads or floating-point decoding are used.
     */
    class CompactStructuredLayerInterpolation {
      public:
        static constexpr std::uint16_t maximum_atmosphere_offset = 4095;

        explicit CompactStructuredLayerInterpolation(
            const StructuredLayerInterpolation& layer)
            : m_source_offset(layer.source_offset),
              m_cell_base(layer.cell_base) {
            if (layer.atmosphere_offset > maximum_atmosphere_offset ||
                layer.atmosphere_mask > 15 || layer.reserved != 0) {
                throw std::out_of_range("Compact structured layer overflow");
            }
            m_atmosphere_offset_mask = static_cast<std::uint16_t>(
                layer.atmosphere_offset |
                (static_cast<std::uint16_t>(layer.atmosphere_mask) << 12));
        }
        StructuredLayerInterpolation expanded() const {
            return {static_cast<std::uint16_t>(m_atmosphere_offset_mask & 4095),
                    m_source_offset, m_cell_base,
                    static_cast<std::uint8_t>(m_atmosphere_offset_mask >> 12),
                    0};
        }

      private:
        std::uint16_t m_atmosphere_offset_mask = 0;
        std::uint16_t m_source_offset = 0;
        std::uint16_t m_cell_base = 0;
    };

    static_assert(sizeof(CompactStructuredLayerInterpolation) == 6,
                  "Compact structured layer descriptors must use six bytes");
    static_assert(alignof(CompactStructuredLayerInterpolation) == 2,
                  "Compact structured layer descriptors use aligned fields");

    /** Construction uses wide descriptors; immutable reads decode either form.
     */
    class LayerInterpolationStorage {
      public:
        using value_type = LayerInterpolation;

        LayerInterpolationStorage&
        operator=(std::initializer_list<LayerInterpolation> values) {
            m_values = std::vector<LayerInterpolation>(values);
            return *this;
        }
        std::size_t size() const {
            return std::visit([](const auto& values) { return values.size(); },
                              m_values);
        }
        bool empty() const { return size() == 0; }
        bool is_structured() const { return m_values.index() != 0; }
        bool is_compact_structured() const { return m_values.index() == 2; }
        std::size_t element_bytes() const {
            return std::visit(
                [](const auto& values) {
                    using Value =
                        typename std::decay_t<decltype(values)>::value_type;
                    return sizeof(Value);
                },
                m_values);
        }
        std::size_t capacity_bytes() const {
            return std::visit(
                [](const auto& values) {
                    using Value =
                        typename std::decay_t<decltype(values)>::value_type;
                    return values.capacity() * sizeof(Value);
                },
                m_values);
        }
        void resize(std::size_t size) { wide_values().resize(size); }
        void shrink_to_fit() {
            std::visit([](auto& values) { values.shrink_to_fit(); }, m_values);
        }
        std::vector<LayerInterpolation>& wide_values() {
            auto* values =
                std::get_if<std::vector<LayerInterpolation>>(&m_values);
            if (values == nullptr) {
                throw std::logic_error(
                    "Structured layer descriptors are immutable");
            }
            return *values;
        }
        StructuredLayerInterpolation structured_layer(std::size_t index) const {
            if (const auto* values = std::get_if<
                    std::vector<CompactStructuredLayerInterpolation>>(
                    &m_values)) {
                return (*values)[index].expanded();
            }
            return std::get<std::vector<StructuredLayerInterpolation>>(
                m_values)[index];
        }
        LayerInterpolation operator[](std::size_t index) const {
            if (const auto* values = std::get_if<
                    std::vector<CompactStructuredLayerInterpolation>>(
                    &m_values)) {
                return structured_descriptor(*values, index);
            }
            if (const auto* values =
                    std::get_if<std::vector<StructuredLayerInterpolation>>(
                        &m_values)) {
                return structured_descriptor(*values, index);
            }
            return std::get<std::vector<LayerInterpolation>>(m_values)[index];
        }
        void
        assign_structured(std::vector<StructuredLayerInterpolation>&& values,
                          std::uint32_t source_weight_count) {
            const bool compact = std::all_of(
                values.begin(), values.end(), [](const auto& layer) {
                    return layer.atmosphere_offset <=
                               CompactStructuredLayerInterpolation::
                                   maximum_atmosphere_offset &&
                           layer.atmosphere_mask <= 15 && layer.reserved == 0;
                });
            if (compact) {
                std::vector<CompactStructuredLayerInterpolation> compact_values;
                compact_values.reserve(values.size());
                for (const auto& layer : values) {
                    compact_values.emplace_back(layer);
                }
                m_values = std::move(compact_values);
            } else {
                m_values = std::move(values);
            }
            m_source_weight_count = source_weight_count;
        }

      private:
        template <typename Descriptor>
        LayerInterpolation
        structured_descriptor(const std::vector<Descriptor>& values,
                              std::size_t index) const {
            const auto current = values.cbegin() + index;
            const auto next = current + 1;
            const auto layer = expanded_descriptor(*current);
            const std::uint32_t source_end =
                next == values.cend()
                    ? m_source_weight_count
                    : expanded_descriptor(*next).source_offset;
            constexpr std::array<std::uint8_t, 16> mask_counts = {
                0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
            return {layer.atmosphere_offset,
                    mask_counts[layer.atmosphere_mask],
                    layer.source_offset,
                    source_end - layer.source_offset,
                    static_cast<std::uint32_t>(index * 4),
                    4};
        }
        static StructuredLayerInterpolation
        expanded_descriptor(const StructuredLayerInterpolation& layer) {
            return layer;
        }
        static StructuredLayerInterpolation
        expanded_descriptor(const CompactStructuredLayerInterpolation& layer) {
            return layer.expanded();
        }
        std::variant<std::vector<LayerInterpolation>,
                     std::vector<StructuredLayerInterpolation>,
                     std::vector<CompactStructuredLayerInterpolation>>
            m_values;
        std::uint32_t m_source_weight_count = 0;
    };

    /** Midpoint view preserving the generic coefficient and index interface. */
    class AtmosphereInterpolationView {
      public:
        AtmosphereInterpolationView() = default;
        AtmosphereInterpolationView(
            const std::vector<InterpolationWeight>& values, std::size_t offset,
            std::size_t size)
            : m_generic(values, offset, size), m_size(size) {}
        AtmosphereInterpolationView(const double* values, int cell_base,
                                    int altitude_stride, std::uint8_t mask,
                                    std::size_t size)
            : m_values(values), m_base(cell_base), m_stride(altitude_stride),
              m_mask(mask), m_size(size) {}

        std::size_t size() const { return m_size; }
        bool empty() const { return m_size == 0; }
        InterpolationWeight operator[](std::size_t index) const {
            if (index >= m_size) {
                throw std::out_of_range(
                    "Successive-orders midpoint index is out of range");
            }
            if (m_values == nullptr) {
                return m_generic[index];
            }
            constexpr std::array<std::array<std::uint8_t, 4>, 16> corners = {
                {{0, 0, 0, 0},
                 {0, 0, 0, 0},
                 {1, 0, 0, 0},
                 {0, 1, 0, 0},
                 {2, 0, 0, 0},
                 {0, 2, 0, 0},
                 {1, 2, 0, 0},
                 {0, 1, 2, 0},
                 {3, 0, 0, 0},
                 {0, 3, 0, 0},
                 {1, 3, 0, 0},
                 {0, 1, 3, 0},
                 {2, 3, 0, 0},
                 {0, 2, 3, 0},
                 {1, 2, 3, 0},
                 {0, 1, 2, 3}}};
            const int corner = corners[m_mask][index];
            return {m_base + (corner & 1) + (corner >> 1) * m_stride,
                    m_values[index]};
        }

        class Iterator {
          public:
            using iterator_category = std::input_iterator_tag;
            using value_type = InterpolationWeight;
            using difference_type = std::ptrdiff_t;
            using reference = InterpolationWeight;
            using pointer = void;
            Iterator() = default;
            Iterator(const AtmosphereInterpolationView* view, std::size_t index)
                : m_data(view->m_values == nullptr
                             ? static_cast<const void*>(view->m_generic.data())
                             : static_cast<const void*>(view->m_values)),
                  m_size(view->m_size), m_index(index), m_base(view->m_base),
                  m_stride(view->m_stride), m_mask(view->m_mask),
                  m_structured(view->m_values != nullptr) {}
            InterpolationWeight operator*() const {
                if (m_index >= m_size) {
                    throw std::out_of_range(
                        "Successive-orders midpoint iterator is out of range");
                }
                if (!m_structured) {
                    return static_cast<const InterpolationWeight*>(
                        m_data)[m_index];
                }
                return AtmosphereInterpolationView(
                    static_cast<const double*>(m_data), m_base, m_stride,
                    m_mask, m_size)[m_index];
            }
            Iterator& operator++() {
                ++m_index;
                return *this;
            }
            Iterator operator++(int) {
                auto copy = *this;
                ++*this;
                return copy;
            }
            bool operator==(const Iterator& other) const {
                return m_data == other.m_data && m_index == other.m_index &&
                       m_base == other.m_base && m_stride == other.m_stride &&
                       m_mask == other.m_mask &&
                       m_structured == other.m_structured;
            }
            bool operator!=(const Iterator& other) const {
                return !(*this == other);
            }

          private:
            const void* m_data = nullptr;
            std::size_t m_size = 0;
            std::size_t m_index = 0;
            int m_base = 0;
            int m_stride = 0;
            std::uint8_t m_mask = 0;
            bool m_structured = false;
        };
        Iterator begin() const { return {this, 0}; }
        Iterator end() const { return {this, m_size}; }

      private:
        InterpolationView<InterpolationWeight> m_generic;
        const double* m_values = nullptr;
        int m_base = 0;
        int m_stride = 0;
        std::uint8_t m_mask = 0;
        std::size_t m_size = 0;
    };

    /** OD view over existing indices or four implicitly indexed cell corners.
     */
    class OpticalDepthInterpolationView {
      public:
        OpticalDepthInterpolationView(
            sasktran2::raytracing::GridWeightStencilView generic = {})
            : m_generic(generic) {}
        OpticalDepthInterpolationView(const double* values, int cell_base,
                                      int altitude_stride)
            : m_values(values), m_base(cell_base), m_stride(altitude_stride) {}
        std::size_t size() const {
            return m_values == nullptr ? m_generic.size() : 4;
        }
        bool empty() const { return size() == 0; }
        std::pair<int, double> operator[](std::size_t index) const {
            if (m_values == nullptr) {
                return m_generic[index];
            }
            assert(index < 4);
            return {m_base + static_cast<int>(index & 1) +
                        static_cast<int>(index >> 1) * m_stride,
                    m_values[index]};
        }

      private:
        sasktran2::raytracing::GridWeightStencilView m_generic;
        const double* m_values = nullptr;
        int m_base = 0;
        int m_stride = 0;
    };

    /** Compiled source interpolation for one traced ray.
     *
     * Optical-depth stencils are retained directly so transport calculations
     * do not need the much larger traced-layer geometry after setup.
     */
    struct RayInterpolation {
        const sasktran2::raytracing::TracedRay* traced_ray = nullptr;
        LayerInterpolationStorage layers;
        std::vector<InterpolationWeight> atmosphere_weights;
        std::vector<double> structured_atmosphere_weights;
        int structured_altitude_stride = 0;
        SourceInterpolationStorage source_weights;
        std::vector<int> optical_depth_indices;
        std::vector<double> optical_depth_weights;
        SourceInterpolationStorage ground_weights;
        std::vector<std::pair<int, double>> ground_horizontal_weights;
        bool ground_hit = false;
        bool transport_compiled = false;

        /** Offset of this row in SourceGeometry1D::transport_column_indices. */
        std::size_t transport_value_offset = 0;
        std::uint32_t transport_row_nnz = 0;

        AtmosphereInterpolationView
        atmosphere_for_layer(std::size_t layer_index) const {
            const auto layer = layers[layer_index];
            if (layers.is_structured()) {
                const auto& structured = layers.structured_layer(layer_index);
                return {layer.atmosphere_count == 0
                            ? nullptr
                            : structured_atmosphere_weights.data() +
                                  layer.atmosphere_offset,
                        structured.cell_base, structured_altitude_stride,
                        structured.atmosphere_mask, layer.atmosphere_count};
            }
            return {atmosphere_weights, layer.atmosphere_offset,
                    layer.atmosphere_count};
        }
        SourceInterpolationView
        source_for_layer(std::size_t layer_index) const {
            const auto layer = layers[layer_index];
            return source_weights.view(layer.source_offset, layer.source_count);
        }
        OpticalDepthInterpolationView
        optical_depth_for_layer(std::size_t layer_index) const {
            const auto layer = layers[layer_index];
            if (layers.is_structured()) {
                return {optical_depth_weights.data() +
                            layer.optical_depth_offset,
                        layers.structured_layer(layer_index).cell_base,
                        structured_altitude_stride};
            }
            if (!optical_depth_indices.empty()) {
                return {sasktran2::raytracing::GridWeightStencilView{
                    optical_depth_indices.data() + layer.optical_depth_offset,
                    optical_depth_weights.data() + layer.optical_depth_offset,
                    layer.optical_depth_count}};
            }
            if (traced_ray != nullptr) {
                return traced_ray->optical_depth_weights(layer_index);
            }
            return {};
        }
        SourceInterpolationView ground() const {
            return ground_weights.view(0, ground_weights.size());
        }

        const std::vector<std::pair<int, double>>& ground_horizontal() const {
            return ground_horizontal_weights;
        }

        bool ground_is_hit() const { return ground_hit; }
    };

    /** Reusable temporary storage for compiling ray interpolation. */
    struct InterpolationScratch {
        std::vector<std::pair<int, double>> location;
        std::vector<std::pair<int, double>> direction;
        std::vector<std::pair<int, double>> atmosphere;
        std::vector<InterpolationWeight> compiled_atmosphere;
        std::vector<SourceInterpolationWeight> compiled_source;
    };

    /** Compiles deterministic interpolation metadata for one traced ray. */
    void compile_ray_interpolation(
        const sasktran2::raytracing::TracedRay& ray,
        const sasktran2::Geometry& geometry,
        sasktran2::grids::SourceLocationInterpolator& location_interpolator,
        const std::vector<SourcePoint>& source_points, RayInterpolation& result,
        InterpolationScratch& scratch);

    /** Moves the OD stencil buffers from a construction-time traced ray into
     * its compact runtime interpolation record. */
    void adopt_optical_depth_storage(sasktran2::raytracing::TracedRay& ray,
                                     RayInterpolation& interpolation);

    /** Remove construction capacity from immutable interpolation buffers.
     *
     * Call only before views of these buffers escape geometry construction.
     * Returns the number of bytes released without changing any coefficients,
     * stencil order, or layer offsets.
     */
    std::size_t compact_ray_interpolation(RayInterpolation& interpolation);

    /** Collect global columns without changing construction-time indices. */
    void collect_transport_row(const RayInterpolation& interpolation,
                               std::vector<int>& columns);

    /** Builds one sorted unique CSR row and replaces global indices with
     * local slots. This finalization is performed once per ray. */
    void compile_transport_row(RayInterpolation& interpolation,
                               std::vector<int>& columns);

} // namespace sasktran2::successive_orders
