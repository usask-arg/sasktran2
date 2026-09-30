#pragma once

#include <sasktran2/raytracing.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
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

    /** Geometry metadata required to integrate one traced layer. */
    struct LayerInterpolation {
        std::uint32_t atmosphere_offset = 0;
        std::uint32_t atmosphere_count = 0;
        std::uint32_t source_offset = 0;
        std::uint32_t source_count = 0;
        std::uint32_t optical_depth_offset = 0;
        std::uint32_t optical_depth_count = 0;
    };

    /** Compiled source interpolation for one traced ray.
     *
     * Optical-depth stencils are retained directly so transport calculations
     * do not need the much larger traced-layer geometry after setup.
     */
    struct RayInterpolation {
        const sasktran2::raytracing::TracedRay* traced_ray = nullptr;
        std::vector<LayerInterpolation> layers;
        std::vector<InterpolationWeight> atmosphere_weights;
        std::vector<SourceInterpolationWeight> source_weights;
        std::vector<int> optical_depth_indices;
        std::vector<double> optical_depth_weights;
        std::vector<SourceInterpolationWeight> ground_weights;
        std::vector<std::pair<int, double>> ground_horizontal_weights;
        bool ground_hit = false;
        bool transport_compiled = false;

        /** Offset of this row in SourceGeometry1D::transport_column_indices. */
        std::size_t transport_value_offset = 0;
        std::uint32_t transport_row_nnz = 0;

        InterpolationView<InterpolationWeight>
        atmosphere_for_layer(std::size_t layer_index) const {
            const auto& layer = layers[layer_index];
            return {atmosphere_weights, layer.atmosphere_offset,
                    layer.atmosphere_count};
        }
        InterpolationView<SourceInterpolationWeight>
        source_for_layer(std::size_t layer_index) const {
            const auto& layer = layers[layer_index];
            return {source_weights, layer.source_offset, layer.source_count};
        }
        sasktran2::raytracing::GridWeightStencilView
        optical_depth_for_layer(std::size_t layer_index) const {
            const auto& layer = layers[layer_index];
            if (!optical_depth_indices.empty()) {
                return {
                    optical_depth_indices.data() + layer.optical_depth_offset,
                    optical_depth_weights.data() + layer.optical_depth_offset,
                    layer.optical_depth_count};
            }
            if (traced_ray != nullptr) {
                return traced_ray->optical_depth_weights(layer_index);
            }
            return {};
        }
        InterpolationView<SourceInterpolationWeight> ground() const {
            return InterpolationView<SourceInterpolationWeight>(ground_weights);
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
