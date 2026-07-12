#include <sasktran2/raytracing.h>

#ifdef SKTRAN_RUST_SUPPORT

#include "sasktran2-core/src/raytracer/cxx.rs.h"

namespace {
    std::vector<double> eigen_vector_to_std(const Eigen::VectorXd& values) {
        std::vector<double> result;
        result.reserve(values.size());
        for (int i = 0; i < values.size(); ++i) {
            result.push_back(values(i));
        }
        return result;
    }

    Eigen::Vector3d make_vector(double x, double y, double z) {
        return Eigen::Vector3d(x, y, z);
    }

    void set_location(sasktran2::Location& location,
                      const Eigen::Vector3d& position, bool on_exact_altitude,
                      int lower_alt_index) {
        location.position = position;
        location.on_exact_altitude = on_exact_altitude;
        location.lower_alt_index = lower_alt_index;
    }

    void set_interpolation_weights(
        sasktran2::Location& location, double altitude,
        const sasktran2::Geometry1D& geometry,
        sasktran2::raytracing::InterpolationStencil1D& weights) {
        const int lower = location.lower_alt_index;

        if (location.on_exact_altitude && lower >= 0) {
            weights.cell_index = std::min(lower, geometry.size() - 2);
            weights.upper_weight = lower == weights.cell_index ? 0.0 : 1.0;
            return;
        }

        const auto& altitude_grid = geometry.altitude_grid();
        if (lower < 0 || lower + 1 >= altitude_grid.grid().size()) {
            std::vector<std::pair<int, double>> workspace;
            geometry.assign_interpolation_weights(location, workspace);
            weights.assign(workspace, geometry.size());
            return;
        }

        weights.cell_index = lower;
        if (altitude_grid.interpolation_method() ==
            sasktran2::grids::interpolation::lower) {
            weights.upper_weight = 0.0;
            return;
        }

        if (altitude_grid.interpolation_method() ==
            sasktran2::grids::interpolation::shell) {
            weights.upper_weight = 0.5;
            return;
        }

        const double lower_altitude = altitude_grid.grid()(lower);
        const double upper_altitude = altitude_grid.grid()(lower + 1);
        weights.upper_weight =
            (altitude - lower_altitude) / (upper_altitude - lower_altitude);
    }
} // namespace

namespace sasktran2::rust::raytracer {

    class CppTraceResult {
      public:
        CppTraceResult(sasktran2::raytracing::TracedRay& result,
                       const sasktran2::Geometry1D& geometry,
                       const sasktran2::viewinggeometry::ViewingRay& ray)
            : m_result(result), m_geometry(geometry), m_ray(ray) {}

        void prepare(const RustTraceSummary& summary) {
            m_result.observer_and_look = m_ray;
            m_result.ground_is_hit = summary.ground_is_hit;
            m_result.is_straight = summary.is_straight;
            m_result.tangent_radius = summary.tangent_radius;
            m_result.layers.resize(summary.num_layers);
        }

        void set_layer(std::size_t index, const RustTraceLayer& rust_layer) {
            auto& layer = m_result.layers[index];

            layer.type = static_cast<sasktran2::raytracing::LayerType>(
                rust_layer.layer_type);
            set_location(layer.entrance,
                         make_vector(rust_layer.entrance_x,
                                     rust_layer.entrance_y,
                                     rust_layer.entrance_z),
                         rust_layer.entrance_on_exact_altitude,
                         rust_layer.entrance_lower_alt_index);
            set_location(layer.exit,
                         make_vector(rust_layer.exit_x, rust_layer.exit_y,
                                     rust_layer.exit_z),
                         rust_layer.exit_on_exact_altitude,
                         rust_layer.exit_lower_alt_index);

            layer.r_entrance = m_geometry.coordinates().earth_radius() +
                               rust_layer.entrance_altitude;
            layer.r_exit = m_geometry.coordinates().earth_radius() +
                           rust_layer.exit_altitude;

            layer.layer_distance = rust_layer.layer_distance;
            if (m_result.is_straight || layer.layer_distance == 0.0) {
                layer.average_look_away = m_ray.look_away;
            } else {
                layer.average_look_away =
                    (layer.exit.position - layer.entrance.position) /
                    layer.layer_distance;
            }
            layer.curvature_factor = rust_layer.curvature_factor;
            layer.od_quad_start = rust_layer.od_quad_start;
            layer.od_quad_end = rust_layer.od_quad_end;
            layer.od_quad_start_fraction = rust_layer.od_quad_start_fraction;
            layer.od_quad_end_fraction = rust_layer.od_quad_end_fraction;
            layer.cos_sza_entrance = rust_layer.cos_sza_entrance;
            layer.cos_sza_exit = rust_layer.cos_sza_exit;
            layer.saz_entrance = rust_layer.saz_entrance;
            layer.saz_exit = rust_layer.saz_exit;

            set_interpolation_weights(layer.entrance,
                                      rust_layer.entrance_altitude, m_geometry,
                                      layer.entrance_interpolation_weights);
            set_interpolation_weights(layer.exit, rust_layer.exit_altitude,
                                      m_geometry,
                                      layer.exit_interpolation_weights);
            sasktran2::raytracing::add_integrated_od_weights(layer, m_geometry);
        }

        void set_layers(::rust::Slice<const RustTraceLayer> layers) {
            for (std::size_t i = 0; i < layers.size(); ++i) {
                set_layer(i, layers[i]);
            }
        }

      private:
        sasktran2::raytracing::TracedRay& m_result;
        const sasktran2::Geometry1D& m_geometry;
        const sasktran2::viewinggeometry::ViewingRay& m_ray;
    };

    class CppTraceResult2D {
      public:
        CppTraceResult2D(sasktran2::raytracing::TracedRay2D& result,
                         const sasktran2::Geometry2D& geometry,
                         const sasktran2::viewinggeometry::ViewingRay& ray)
            : m_result(result), m_geometry(geometry), m_ray(ray) {}

        void prepare(const RustTraceSummary& summary) {
            m_result.observer_and_look = m_ray;
            m_result.ground_is_hit = summary.ground_is_hit;
            m_result.is_straight = summary.is_straight;
            m_result.tangent_radius = summary.tangent_radius;
            m_result.layers.resize(summary.num_layers);
        }

        void set_layer(std::size_t index, const RustTraceLayer& rust_layer) {
            auto& layer = m_result.layers[index];

            layer.type = static_cast<sasktran2::raytracing::LayerType>(
                rust_layer.layer_type);
            set_location(layer.entrance,
                         make_vector(rust_layer.entrance_x,
                                     rust_layer.entrance_y,
                                     rust_layer.entrance_z),
                         rust_layer.entrance_on_exact_altitude, -1);
            set_location(layer.exit,
                         make_vector(rust_layer.exit_x, rust_layer.exit_y,
                                     rust_layer.exit_z),
                         rust_layer.exit_on_exact_altitude, -1);

            layer.r_entrance = m_geometry.coordinates().earth_radius() +
                               rust_layer.entrance_altitude;
            layer.r_exit = m_geometry.coordinates().earth_radius() +
                           rust_layer.exit_altitude;
            layer.layer_distance = rust_layer.layer_distance;
            if (layer.layer_distance == 0.0) {
                layer.average_look_away = m_ray.look_away;
            } else {
                layer.average_look_away =
                    (layer.exit.position - layer.entrance.position) /
                    layer.layer_distance;
            }
            layer.curvature_factor = rust_layer.curvature_factor;
            layer.od_quad_start = rust_layer.od_quad_start;
            layer.od_quad_end = rust_layer.od_quad_end;
            layer.od_quad_start_fraction = rust_layer.od_quad_start_fraction;
            layer.od_quad_end_fraction = rust_layer.od_quad_end_fraction;
            layer.cos_sza_entrance = rust_layer.cos_sza_entrance;
            layer.cos_sza_exit = rust_layer.cos_sza_exit;
            layer.saz_entrance = rust_layer.saz_entrance;
            layer.saz_exit = rust_layer.saz_exit;
            layer.altitude_cell = rust_layer.cell_altitude_index;
            layer.horizontal_cell = rust_layer.cell_horizontal_index;
            const auto entrance_coordinates =
                m_geometry.cell_interpolation_coordinates(
                    layer.entrance, layer.altitude_cell, layer.horizontal_cell);
            const auto exit_coordinates =
                m_geometry.cell_interpolation_coordinates(
                    layer.exit, layer.altitude_cell, layer.horizontal_cell);
            layer.entrance_interpolation.altitude_upper_weight =
                entrance_coordinates.first;
            layer.entrance_interpolation.horizontal_upper_weight =
                entrance_coordinates.second;
            layer.exit_interpolation.altitude_upper_weight =
                exit_coordinates.first;
            layer.exit_interpolation.horizontal_upper_weight =
                exit_coordinates.second;
            layer.integrated_od.weights = {rust_layer.integrated_od_weight_0,
                                           rust_layer.integrated_od_weight_1,
                                           rust_layer.integrated_od_weight_2,
                                           rust_layer.integrated_od_weight_3};
            layer.entrance.lower_alt_index = layer.altitude_cell;
            layer.exit.lower_alt_index = layer.altitude_cell;
        }

      private:
        sasktran2::raytracing::TracedRay2D& m_result;
        const sasktran2::Geometry2D& m_geometry;
        const sasktran2::viewinggeometry::ViewingRay& m_ray;
    };

    void prepare_trace_result(CppTraceResult& result,
                              const RustTraceSummary& summary) {
        result.prepare(summary);
    }

    void set_trace_layer(CppTraceResult& result, std::size_t index,
                         const RustTraceLayer& layer) {
        result.set_layer(index, layer);
    }

    void prepare_trace_result_2d(CppTraceResult2D& result,
                                 const RustTraceSummary& summary) {
        result.prepare(summary);
    }

    void set_trace_layer_2d(CppTraceResult2D& result, std::size_t index,
                            const RustTraceLayer& layer) {
        result.set_layer(index, layer);
    }

    void set_trace_layers(CppTraceResult& result,
                          ::rust::Slice<const RustTraceLayer> layers) {
        result.set_layers(layers);
    }

} // namespace sasktran2::rust::raytracer

namespace sasktran2::raytracing {

    class RustRayTracerImpl {
      public:
        explicit RustRayTracerImpl(
            ::rust::Box<sasktran2::rust::raytracer::RustVerticalTracer>
                rust_tracer)
            : rust_tracer(std::move(rust_tracer)) {}

        ::rust::Box<sasktran2::rust::raytracer::RustVerticalTracer> rust_tracer;
    };

    class RustRayTracer2DImpl {
      public:
        explicit RustRayTracer2DImpl(
            ::rust::Box<sasktran2::rust::raytracer::RustStructuredTracer2D>
                rust_tracer)
            : rust_tracer(std::move(rust_tracer)) {}

        ::rust::Box<sasktran2::rust::raytracer::RustStructuredTracer2D>
            rust_tracer;
    };

    RustRayTracer::RustRayTracer(const sasktran2::Geometry1D& geometry)
        : m_geometry(geometry) {
        auto altitudes = eigen_vector_to_std(geometry.altitude_grid().grid());
        auto refractive_index =
            eigen_vector_to_std(geometry.refractive_index());
        const auto& sun = geometry.coordinates().sun_unit();
        auto rust_tracer = sasktran2::rust::raytracer::new_vertical_tracer(
            geometry.coordinates().earth_radius(), altitudes, refractive_index,
            static_cast<int>(geometry.coordinates().geometry_type()),
            static_cast<int>(geometry.altitude_grid().interpolation_method()),
            sun.x(), sun.y(), sun.z());

        m_impl = std::make_unique<RustRayTracerImpl>(std::move(rust_tracer));
    }

    RustRayTracer::~RustRayTracer() = default;

    RustRayTracer2D::RustRayTracer2D(const sasktran2::Geometry2D& geometry)
        : m_geometry(geometry) {
        auto altitudes = eigen_vector_to_std(geometry.altitude_grid().grid());
        auto horizontal_angles =
            eigen_vector_to_std(geometry.horizontal_angle_grid());
        const auto& coordinates = geometry.coordinates();
        const auto& reference_x = coordinates.reference_x();
        const auto& reference_z = coordinates.reference_z();
        const auto& sun = coordinates.sun_unit();
        auto rust_tracer = sasktran2::rust::raytracer::new_structured_tracer_2d(
            coordinates.earth_radius(), altitudes, horizontal_angles,
            static_cast<int>(geometry.altitude_grid().interpolation_method()),
            reference_x.x(), reference_x.y(), reference_x.z(), reference_z.x(),
            reference_z.y(), reference_z.z(), sun.x(), sun.y(), sun.z());

        m_impl = std::make_unique<RustRayTracer2DImpl>(std::move(rust_tracer));
    }

    RustRayTracer2D::~RustRayTracer2D() = default;

    void
    RustRayTracer::trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray,
                             TracedRay& result, bool include_refraction) const {
        sasktran2::rust::raytracer::CppTraceResult cpp_result(result,
                                                              m_geometry, ray);
        sasktran2::rust::raytracer::trace_vertical_ray_into_cpp_result(
            *m_impl->rust_tracer, ray.observer.position.x(),
            ray.observer.position.y(), ray.observer.position.z(),
            ray.look_away.x(), ray.look_away.y(), ray.look_away.z(),
            include_refraction, cpp_result);
    }

    void RustRayTracer2D::trace_ray(
        const sasktran2::viewinggeometry::ViewingRay& ray,
        TracedRay2D& result) const {
        trace_ray_impl(ray, nullptr, result);
    }

    void RustRayTracer2D::trace_ray(
        const sasktran2::viewinggeometry::ViewingRay& ray,
        const Eigen::VectorXd& refractive_index, TracedRay2D& result) const {
        trace_ray_impl(ray, &refractive_index, result);
    }

    void RustRayTracer2D::trace_ray_impl(
        const sasktran2::viewinggeometry::ViewingRay& ray,
        const Eigen::VectorXd* refractive_index, TracedRay2D& result) const {
        if (refractive_index != nullptr &&
            refractive_index->size() != m_geometry.num_altitudes()) {
            throw std::invalid_argument(
                "The per-ray refractive-index profile must have one value "
                "per Geometry2D altitude");
        }

        const ::rust::Slice<const double> profile(
            refractive_index == nullptr ? nullptr : refractive_index->data(),
            refractive_index == nullptr
                ? 0
                : static_cast<std::size_t>(refractive_index->size()));
        sasktran2::rust::raytracer::CppTraceResult2D cpp_result(
            result, m_geometry, ray);
        sasktran2::rust::raytracer::trace_structured_ray_2d_into_cpp_result(
            *m_impl->rust_tracer, ray.observer.position.x(),
            ray.observer.position.y(), ray.observer.position.z(),
            ray.look_away.x(), ray.look_away.y(), ray.look_away.z(), profile,
            cpp_result);
    }

} // namespace sasktran2::raytracing

#endif
