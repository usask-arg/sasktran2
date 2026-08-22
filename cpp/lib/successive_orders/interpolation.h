#pragma once

#include <sasktran2/raytracing.h>

#include <cstddef>
#include <cstdint>
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

    /** One geometry-only interpolation coefficient. */
    struct InterpolationWeight {
        int index = 0;
        double weight = 0.0;
    };

    /** Interpolation from the global outgoing-source vector.
     *
     * `row_inner_index` is the position of `source_index` in the owning ray's
     * sorted `transport_columns` array. It lets the atmosphere-dependent
     * transport assembler accumulate directly into an existing CSR row.
     */
    struct SourceInterpolationWeight {
        double weight = 0.0;
        int source_index = 0;
        std::uint32_t row_inner_index = 0;

        SourceInterpolationWeight() = default;
        SourceInterpolationWeight(int index, double interpolation_weight,
                                  std::uint32_t inner_index = 0)
            : weight(interpolation_weight), source_index(index),
              row_inner_index(inner_index) {}
    };

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

    /** Builds one sorted unique CSR row and assigns local slots to its source
     * interpolation entries. */
    void compile_transport_row(RayInterpolation& interpolation,
                               std::vector<int>& columns);

} // namespace sasktran2::successive_orders
