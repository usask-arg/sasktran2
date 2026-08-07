#pragma once

#include <sasktran2/geometry.h>
#include <sasktran2/raytracing.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace sasktran2::successive_orders {
    /** Packed source interpolation for one traced diffuse ray.
     *
     * This representation is deliberately private to the Rust
     * successive-orders adapter. The legacy C++ HR source keeps its original
     * interpolation representation and ABI.
     */
    template <int NSTOKES> struct RaySourceInterpolationWeights {
        struct Layer {
            std::uint32_t atmosphere_offset = 0;
            std::uint32_t source_offset = 0;
            std::uint32_t atmosphere_count = 0;
            std::uint32_t source_count = 0;
        };

        std::vector<Layer> interior_weights;
        std::vector<double> atmosphere_weights;
        std::vector<int> source_indices;
        std::vector<double> source_weights;
        std::vector<std::uint16_t> source_accumulation_inner_indices;

        std::vector<int> ground_source_indices;
        std::vector<double> ground_source_weights;
        std::vector<std::uint16_t> ground_accumulation_inner_indices;
        std::uint32_t accumulation_row_offset = 0;
        std::uint32_t accumulation_row_nnz = 0;
        bool ground_is_hit = false;
    };

    /** Geometry-only transport stencils transferred once to Rust. */
    struct RayTransportGeometry {
        std::vector<std::uint32_t> ray_layer_offsets;
        std::vector<std::uint32_t> layer_atmosphere_offsets;
        std::vector<std::uint32_t> atmosphere_indices;
        std::vector<double> optical_depth_weights;
        std::vector<double> albedo_weights;
        std::vector<double> entrance_weights;
        std::vector<double> exit_weights;
        std::vector<double> layer_distance;
        std::vector<double> layer_start_fraction;
        std::vector<double> layer_end_fraction;
        std::vector<double> ray_scattering_cosine;
        std::vector<double> ray_phase_q_factor;
        std::vector<double> ray_phase_u_factor;
        std::vector<std::uint32_t> ray_transport_value_offsets;
        std::vector<std::uint32_t> ray_transport_row_nnz;
        std::vector<std::uint32_t> layer_source_offsets;
        std::vector<std::uint16_t> source_value_inner_indices;
        std::vector<double> source_weights;
        std::vector<std::uint32_t> ray_ground_offsets;
        std::vector<std::uint16_t> ground_value_inner_indices;
        std::vector<double> ground_weights;
    };

    /** Temporary C++ owner for traced geometry until Rust takes ownership.
     *
     * In 2D, each ray is compacted immediately after tracing so peak setup
     * memory does not scale with the full C++ TracedRay representation.
     */
    template <int NSTOKES> class TransportGeometry {
      private:
        struct CompactLayer {
            double layer_distance;
            double od_quad_start;
            double od_quad_end;
            double od_quad_start_fraction;
            double od_quad_end_fraction;
            std::uint32_t grid_weight_offset;
            std::uint8_t grid_weight_count;
        };

        struct CompactGridWeights {
            std::vector<int> indices;
            std::vector<double> entrance;
            std::vector<double> exit;
            std::vector<double> optical_depth;
        };

        struct CompactRay {
            std::vector<CompactLayer> layers;
            CompactGridWeights weights;
        };

        const std::vector<sasktran2::raytracing::TracedRay>* m_traced_rays =
            nullptr;
        bool m_use_compact_geometry = false;
        bool m_released = false;
        std::vector<CompactRay> m_compact_rays;

        int num_layers(int rayidx) const;
        const sasktran2::raytracing::TracedLayer&
        layer(int rayidx, int layeridx,
              sasktran2::raytracing::TracedLayer& scratch) const;
        sasktran2::raytracing::GridWeightStencilView
        entrance_weights(int rayidx, int layeridx) const;
        sasktran2::raytracing::GridWeightStencilView
        exit_weights(int rayidx, int layeridx) const;
        sasktran2::raytracing::GridWeightStencilView
        optical_depth_weights(int rayidx, int layeridx) const;

      public:
        using SInterpolator =
            std::vector<RaySourceInterpolationWeights<NSTOKES>>;

        void initialize(
            const std::vector<sasktran2::raytracing::TracedRay>& traced_rays);
        void begin_compact_2d(
            std::vector<sasktran2::raytracing::TracedRay>& traced_rays,
            const sasktran2::Geometry& geometry);
        void compact_2d_ray(int rayidx,
                            sasktran2::raytracing::TracedRay& traced_ray);
        void finalize_compact_2d();

        RayTransportGeometry
        pack(const SInterpolator& source_interpolator,
             const sasktran2::Geometry& model_geometry) const;

        void release();
    };
} // namespace sasktran2::successive_orders
