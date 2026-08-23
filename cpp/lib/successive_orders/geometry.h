#pragma once

#include "interpolation.h"

#include <sasktran2/geometry.h>
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/viewinggeometry_internal.h>

#include <memory>
#include <vector>

namespace sasktran2::grids {
    class SourceLocationInterpolator;
}

namespace sasktran2::successive_orders {
    /** Geometry-only configuration for the successive-orders source grid. */
    struct SourceGeometrySettings {
        int num_incoming = 110;
        int num_outgoing = 110;
        int num_sza = 1;
        int num_threads = 1;
        bool include_refraction = false;

        /** Empty selects one source altitude at each atmosphere-layer
         * midpoint. */
        std::vector<double> altitude_grid_m;

        /** Empty selects an evenly spaced Geometry2D horizontal grid with
         * num_sza points. Values are local Geometry2D angles in radians. */
        std::vector<double> horizontal_angle_grid_radians;

        void validate() const;
    };

    /** One location on the successive-orders source grid. */
    class SourcePoint {
      public:
        const sasktran2::Location& location() const { return m_location; }
        bool is_ground() const { return m_is_ground; }

        const sasktran2::math::UnitSphere& incoming_sphere() const {
            return *m_incoming_sphere;
        }
        const sasktran2::math::UnitSphere& outgoing_sphere() const {
            return *m_outgoing_sphere;
        }

        int num_incoming() const { return m_incoming_sphere->num_points(); }
        int num_outgoing() const { return m_outgoing_sphere->num_points(); }
        int incoming_offset() const { return m_incoming_offset; }
        int outgoing_offset() const { return m_outgoing_offset; }

        const std::vector<InterpolationWeight>& atmosphere_weights() const {
            return m_atmosphere_weights;
        }

      private:
        friend class SourceGeometry1D;

        sasktran2::Location m_location;
        const sasktran2::math::UnitSphere* m_incoming_sphere = nullptr;
        const sasktran2::math::UnitSphere* m_outgoing_sphere = nullptr;
        std::vector<InterpolationWeight> m_atmosphere_weights;
        int m_incoming_offset = 0;
        int m_outgoing_offset = 0;
        bool m_is_ground = false;
    };

    /** Immutable source geometry and interpolation topology.
     *
     * Atmospheric optical properties are intentionally absent. Incoming rays
     * are traced once and both incoming-ray and LOS source interpolation are
     * compiled into deterministic metadata suitable for a changing
     * atmosphere and either scalar or vector radiance storage.
     */
    class SourceGeometry1D {
      public:
        SourceGeometry1D(const sasktran2::raytracing::RayTracerBase& raytracer,
                         const sasktran2::Geometry1D& geometry);
#ifdef SKTRAN_RUST_SUPPORT
        SourceGeometry1D(
            const sasktran2::raytracing::RustRayTracer2D& raytracer,
            const sasktran2::Geometry2D& geometry);
#endif
        ~SourceGeometry1D();

        SourceGeometry1D(const SourceGeometry1D&) = delete;
        SourceGeometry1D& operator=(const SourceGeometry1D&) = delete;
        SourceGeometry1D(SourceGeometry1D&&) = delete;
        SourceGeometry1D& operator=(SourceGeometry1D&&) = delete;

        void
        initialize(const sasktran2::viewinggeometry::InternalViewingGeometry&
                       internal_viewing,
                   const SourceGeometrySettings& settings);

        /** Recompile only observer-LOS interpolation and transport topology.
         * Source points and diffuse incoming rays are unchanged. */
        void
        refresh_los(const sasktran2::viewinggeometry::InternalViewingGeometry&
                        internal_viewing);

        const SourceGeometrySettings& settings() const { return m_settings; }
        const std::vector<double>& source_altitudes_m() const {
            return m_source_altitudes_m;
        }
        const std::vector<double>& source_cos_sza() const {
            return m_source_cos_sza;
        }
        const std::vector<double>& source_horizontal_angles_rad() const {
            return m_source_horizontal_angles_rad;
        }

        int num_interior_points() const { return m_num_interior_points; }
        int num_ground_points() const { return m_num_ground_points; }
        int num_points() const {
            return static_cast<int>(m_source_points.size());
        }
        int total_num_incoming() const {
            return m_incoming_point_offsets.empty()
                       ? 0
                       : m_incoming_point_offsets.back();
        }
        int total_num_outgoing() const {
            return m_outgoing_point_offsets.empty()
                       ? 0
                       : m_outgoing_point_offsets.back();
        }

        const std::vector<SourcePoint>& source_points() const {
            return m_source_points;
        }
        const SourcePoint& source_point(int index) const {
            return m_source_points.at(index);
        }
        const std::vector<std::pair<int, double>>&
        ground_horizontal_weights(int ground_index) const {
            return m_ground_horizontal_weights.at(ground_index);
        }
        const std::vector<int>& incoming_point_offsets() const {
            return m_incoming_point_offsets;
        }
        const std::vector<int>& outgoing_point_offsets() const {
            return m_outgoing_point_offsets;
        }

        const sasktran2::viewinggeometry::InternalViewingGeometry&
        incoming_viewing_geometry() const {
            return m_incoming_viewing;
        }
        const std::vector<sasktran2::raytracing::TracedRay>&
        incoming_rays() const {
            return m_incoming_viewing.traced_rays;
        }
        const std::vector<RayInterpolation>& incoming_interpolation() const {
            return m_incoming_interpolation;
        }
        const std::vector<RayInterpolation>& los_interpolation() const {
            return m_los_interpolation;
        }

        /** Releases the construction-only diffuse traced layers after all
         * scalar consumers have compiled their compact geometry. */
        void release_incoming_traced_rays();

        /** CSR topology for outgoing-source transport to incoming rays. */
        const std::vector<int>& transport_row_offsets() const {
            return m_transport_row_offsets;
        }
        const std::vector<int>& transport_column_indices() const {
            return m_transport_column_indices;
        }
        InterpolationView<int>
        transport_columns_for_ray(std::size_t ray_index) const {
            if (ray_index >= m_incoming_interpolation.size()) {
                throw std::out_of_range(
                    "Successive-orders transport ray is out of range");
            }
            const auto& ray = m_incoming_interpolation[ray_index];
            return {m_transport_column_indices, ray.transport_value_offset,
                    ray.transport_row_nnz};
        }
        const std::vector<int>& los_transport_row_offsets() const {
            return m_los_transport_row_offsets;
        }
        const std::vector<int>& los_transport_column_indices() const {
            return m_los_transport_column_indices;
        }
        InterpolationView<int>
        los_transport_columns_for_ray(std::size_t ray_index) const {
            if (ray_index >= m_los_interpolation.size()) {
                throw std::out_of_range(
                    "Successive-orders LOS ray is out of range");
            }
            const auto& ray = m_los_interpolation[ray_index];
            return {m_los_transport_column_indices, ray.transport_value_offset,
                    ray.transport_row_nnz};
        }

      private:
        struct AngularGridPair {
            std::unique_ptr<const sasktran2::math::UnitSphere> incoming;
            std::unique_ptr<const sasktran2::math::UnitSphere> outgoing;
        };

        sasktran2::grids::AltitudeGrid make_altitude_grid();
        sasktran2::grids::Grid make_cos_sza_grid(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing);
        void construct_source_points();
        void trace_and_compile_incoming();
        void compile_los_interpolation(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing);
        void
        compile_transport_topology(std::vector<RayInterpolation>& interpolation,
                                   std::vector<int>& row_offsets,
                                   std::vector<int>& column_indices);
        void
        trace_ray(const sasktran2::viewinggeometry::ViewingRay& viewing_ray,
                  sasktran2::raytracing::TracedRay& traced_ray) const;

        const sasktran2::Geometry& m_geometry;
        const sasktran2::Geometry1D* m_geometry_1d = nullptr;
        const sasktran2::Geometry2D* m_geometry_2d = nullptr;
        const sasktran2::raytracing::RayTracerBase* m_raytracer_1d = nullptr;
#ifdef SKTRAN_RUST_SUPPORT
        const sasktran2::raytracing::RustRayTracer2D* m_raytracer_2d = nullptr;
#endif
        SourceGeometrySettings m_settings;

        std::unique_ptr<sasktran2::grids::SourceLocationInterpolator>
            m_location_interpolator;
        std::vector<std::unique_ptr<AngularGridPair>> m_angular_grids;
        std::vector<SourcePoint> m_source_points;
        std::vector<int> m_incoming_point_offsets;
        std::vector<int> m_outgoing_point_offsets;
        std::vector<double> m_source_altitudes_m;
        std::vector<double> m_source_cos_sza;
        std::vector<double> m_source_horizontal_angles_rad;
        std::vector<std::vector<std::pair<int, double>>>
            m_ground_horizontal_weights;
        int m_num_interior_points = 0;
        int m_num_ground_points = 0;

        sasktran2::viewinggeometry::InternalViewingGeometry m_incoming_viewing;
        std::vector<RayInterpolation> m_incoming_interpolation;
        std::vector<RayInterpolation> m_los_interpolation;
        std::vector<int> m_transport_row_offsets{0};
        std::vector<int> m_transport_column_indices;
        std::vector<int> m_los_transport_row_offsets{0};
        std::vector<int> m_los_transport_column_indices;
    };

} // namespace sasktran2::successive_orders
