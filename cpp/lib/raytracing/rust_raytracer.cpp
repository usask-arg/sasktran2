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
                      const Eigen::Vector3d& position,
                      bool on_exact_altitude, int lower_alt_index) {
        location.position = position;
        location.on_exact_altitude = on_exact_altitude;
        location.lower_alt_index = lower_alt_index;
        location.interpolation_weights.clear();
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
            m_result.reset();
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

            layer.r_entrance = layer.entrance.radius();
            layer.r_exit = layer.exit.radius();

            const Eigen::Vector3d delta =
                layer.exit.position - layer.entrance.position;
            const double delta_norm = delta.norm();
            layer.average_look_away =
                delta_norm == 0.0 ? m_ray.look_away : delta / delta_norm;

            layer.layer_distance = rust_layer.layer_distance;
            layer.curvature_factor = rust_layer.curvature_factor;
            layer.od_quad_start = rust_layer.od_quad_start;
            layer.od_quad_end = rust_layer.od_quad_end;
            layer.od_quad_start_fraction = rust_layer.od_quad_start_fraction;
            layer.od_quad_end_fraction = rust_layer.od_quad_end_fraction;
            layer.cos_sza_entrance = rust_layer.cos_sza_entrance;
            layer.cos_sza_exit = rust_layer.cos_sza_exit;
            layer.saz_entrance = rust_layer.saz_entrance;
            layer.saz_exit = rust_layer.saz_exit;

            m_geometry.assign_interpolation_weights(
                layer.entrance, layer.entrance.interpolation_weights);
            m_geometry.assign_interpolation_weights(
                layer.exit, layer.exit.interpolation_weights);
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

    void prepare_trace_result(CppTraceResult& result,
                              const RustTraceSummary& summary) {
        result.prepare(summary);
    }

    void set_trace_layer(CppTraceResult& result, std::size_t index,
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

    RustRayTracer::RustRayTracer(const sasktran2::Geometry1D& geometry)
        : m_geometry(geometry) {
        auto altitudes = eigen_vector_to_std(geometry.altitude_grid().grid());
        auto refractive_index = eigen_vector_to_std(geometry.refractive_index());
        const auto& sun = geometry.coordinates().sun_unit();
        auto rust_tracer = sasktran2::rust::raytracer::new_vertical_tracer(
            geometry.coordinates().earth_radius(), altitudes, refractive_index,
            static_cast<int>(geometry.coordinates().geometry_type()),
            static_cast<int>(geometry.altitude_grid().interpolation_method()),
            sun.x(), sun.y(), sun.z());

        m_impl = std::make_unique<RustRayTracerImpl>(std::move(rust_tracer));
    }

    RustRayTracer::~RustRayTracer() = default;

    void RustRayTracer::trace_ray(
        const sasktran2::viewinggeometry::ViewingRay& ray, TracedRay& result,
        bool include_refraction) const {
        sasktran2::rust::raytracer::CppTraceResult cpp_result(result,
                                                             m_geometry, ray);
        sasktran2::rust::raytracer::trace_vertical_ray_into_cpp_result(
            *m_impl->rust_tracer, ray.observer.position.x(),
            ray.observer.position.y(), ray.observer.position.z(),
            ray.look_away.x(), ray.look_away.y(), ray.look_away.z(),
            include_refraction, cpp_result);
    }

} // namespace sasktran2::raytracing

#endif
