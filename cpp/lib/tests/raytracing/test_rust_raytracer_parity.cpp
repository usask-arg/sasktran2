#include <sasktran2/test_helper.h>

#include <sasktran2.h>
#include <chrono>
#include <iostream>

#ifdef SKTRAN_RUST_SUPPORT

#include "sasktran2-core/src/raytracer/cxx.rs.h"

namespace {
    using RustLayer = sasktran2::rust::raytracer::RustTraceLayer;
    using RustSummary = sasktran2::rust::raytracer::RustTraceSummary;

    struct ParityGeometry {
        sasktran2::Geometry1D geometry;
        std::vector<double> altitudes;
        std::vector<double> refractive_index;

        ParityGeometry(sasktran2::Geometry1D&& geometry,
                       std::vector<double>&& altitudes,
                       std::vector<double>&& refractive_index)
            : geometry(std::move(geometry)), altitudes(std::move(altitudes)),
              refractive_index(std::move(refractive_index)) {}
    };

    std::vector<double> atmosphere_like_refractive_index() {
        return {1.00027691, 1.00025137, 1.00022759, 1.00020559, 1.00018526,
                1.00016652, 1.00014926, 1.00013342, 1.00011888, 1.0001056,
                1.0000935,  1.00008041, 1.00006916, 1.00005949, 1.00005117,
                1.00004402, 1.00003763, 1.00003217, 1.0000275,  1.00002351,
                1.0000201,  1.00001713, 1.00001461, 1.00001246, 1.00001062,
                1.00000906, 1.00000775, 1.00000664, 1.00000568, 1.00000486,
                1.00000416, 1.00000357, 1.00000306, 1.00000263, 1.00000226,
                1.00000194, 1.00000166, 1.00000143, 1.00000122, 1.00000105,
                1.0000009,  1.00000079, 1.00000069, 1.0000006,  1.00000052,
                1.00000046, 1.0000004,  1.00000035, 1.0000003,  1.00000027,
                1.00000023, 1.00000021, 1.00000018, 1.00000016, 1.00000014,
                1.00000013, 1.00000011, 1.0000001,  1.00000009, 1.00000008,
                1.00000007, 1.00000006, 1.00000005, 1.00000005, 1.00000004,
                1.00000004};
    }

    ParityGeometry
    make_geometry_with_grid(std::vector<double>&& altitudes,
                            sasktran2::grids::gridspacing grid_spacing,
                            sasktran2::grids::interpolation interpolation,
                            std::vector<double>&& refractive_index) {
        constexpr double earth_radius = 6372000.0;
        Eigen::VectorXd grid_values(altitudes.size());
        for (int i = 0; i < static_cast<int>(altitudes.size()); ++i) {
            grid_values(i) = altitudes[i];
        }

        sasktran2::grids::AltitudeGrid grid(
            std::move(grid_values), grid_spacing,
            sasktran2::grids::outofbounds::extend, interpolation);
        sasktran2::Coordinates coordinates(0.6, 0.3, earth_radius);
        sasktran2::Geometry1D geometry(std::move(coordinates), std::move(grid));

        for (int i = 0; i < static_cast<int>(refractive_index.size()); ++i) {
            geometry.refractive_index()(i) = refractive_index[i];
        }

        return ParityGeometry(std::move(geometry), std::move(altitudes),
                              std::move(refractive_index));
    }

    ParityGeometry make_geometry() {
        constexpr int num_altitudes = 66;
        constexpr double spacing = 1000.0;

        std::vector<double> altitudes(num_altitudes);
        for (int i = 0; i < num_altitudes; ++i) {
            altitudes[i] = i * spacing;
        }

        return make_geometry_with_grid(std::move(altitudes),
                                       sasktran2::grids::gridspacing::constant,
                                       sasktran2::grids::interpolation::linear,
                                       atmosphere_like_refractive_index());
    }

    double tolerance(double expected, double absolute_floor = 1e-6) {
        return std::max(absolute_floor, 1e-8 * std::abs(expected));
    }

    void require_close(const char* label, double actual, double expected,
                       double absolute_floor = 1e-6) {
        INFO(label);
        CAPTURE(actual);
        CAPTURE(expected);
        REQUIRE(std::abs(actual - expected) <=
                tolerance(expected, absolute_floor));
    }

    void require_position_close(const char* label,
                                const Eigen::Vector3d& actual, double x,
                                double y, double z) {
        INFO(label);
        require_close("x", x, actual.x());
        require_close("y", y, actual.y());
        require_close("z", z, actual.z());
    }

    void require_angle_close(const char* label, double actual, double expected,
                             double absolute_floor = 1e-8) {
        INFO(label);
        CAPTURE(actual);
        CAPTURE(expected);
        REQUIRE(std::abs(std::remainder(actual - expected, 2 * EIGEN_PI)) <=
                tolerance(expected, absolute_floor));
    }

    RustSummary trace_rust(const ParityGeometry& parity_geometry,
                           const sasktran2::viewinggeometry::ViewingRay& ray,
                           bool include_refraction,
                           std::vector<RustLayer>& rust_layers) {
        rust_layers.clear();
        const auto& sun = parity_geometry.geometry.coordinates().sun_unit();
        return sasktran2::rust::raytracer::trace_vertical_ray(
            parity_geometry.geometry.coordinates().earth_radius(),
            parity_geometry.altitudes, parity_geometry.refractive_index,
            static_cast<int>(
                parity_geometry.geometry.coordinates().geometry_type()),
            static_cast<int>(parity_geometry.geometry.altitude_grid()
                                 .interpolation_method()),
            sun.x(), sun.y(), sun.z(), ray.observer.position.x(),
            ray.observer.position.y(), ray.observer.position.z(),
            ray.look_away.x(), ray.look_away.y(), ray.look_away.z(),
            include_refraction, rust_layers);
    }

    double checksum_cpp_ray(const sasktran2::raytracing::TracedRay& ray) {
        double checksum = ray.tangent_radius * 1e-9;
        checksum += ray.ground_is_hit ? 7.0 : 11.0;
        checksum += ray.is_straight ? 13.0 : 17.0;
        checksum += static_cast<double>(ray.layers.size());
        for (const auto& layer : ray.layers) {
            checksum += layer.layer_distance * 1e-6;
            checksum += layer.curvature_factor;
            checksum += layer.od_quad_start * 1e-7;
            checksum += layer.od_quad_end * 1e-7;
            checksum += layer.cos_sza_entrance * 1e-3;
            checksum += layer.cos_sza_exit * 1e-3;
        }
        return checksum;
    }

    struct TimingResult {
        double seconds;
        double checksum;
        size_t num_traces;
        size_t total_layers;
    };

    TimingResult time_cpp_tracer(
        const sasktran2::raytracing::RayTracerBase& tracer,
        const std::vector<sasktran2::viewinggeometry::ViewingRay>& rays,
        bool include_refraction, size_t iterations) {
        sasktran2::raytracing::TracedRay traced;
        double checksum = 0.0;
        size_t total_layers = 0;
        auto start = std::chrono::steady_clock::now();
        for (size_t iteration = 0; iteration < iterations; ++iteration) {
            for (const auto& ray : rays) {
                tracer.trace_ray(ray, traced, include_refraction);
                checksum += checksum_cpp_ray(traced);
                total_layers += traced.layers.size();
            }
        }
        auto end = std::chrono::steady_clock::now();
        return {
            std::chrono::duration<double>(end - start).count(),
            checksum,
            iterations * rays.size(),
            total_layers,
        };
    }

    TimingResult time_rust_tracer(
        const ParityGeometry& parity_geometry,
        const std::vector<sasktran2::viewinggeometry::ViewingRay>& rays,
        bool include_refraction, size_t iterations) {
        std::vector<double> origins;
        std::vector<double> looks;
        origins.reserve(rays.size() * 3);
        looks.reserve(rays.size() * 3);
        for (const auto& ray : rays) {
            origins.push_back(ray.observer.position.x());
            origins.push_back(ray.observer.position.y());
            origins.push_back(ray.observer.position.z());
            looks.push_back(ray.look_away.x());
            looks.push_back(ray.look_away.y());
            looks.push_back(ray.look_away.z());
        }

        const auto& sun = parity_geometry.geometry.coordinates().sun_unit();
        auto start = std::chrono::steady_clock::now();
        auto summary =
            sasktran2::rust::raytracer::trace_vertical_ray_batch_checksum(
                parity_geometry.geometry.coordinates().earth_radius(),
                parity_geometry.altitudes, parity_geometry.refractive_index,
                static_cast<int>(
                    parity_geometry.geometry.coordinates().geometry_type()),
                static_cast<int>(parity_geometry.geometry.altitude_grid()
                                     .interpolation_method()),
                sun.x(), sun.y(), sun.z(), origins, looks, include_refraction,
                iterations);
        auto end = std::chrono::steady_clock::now();

        return {
            std::chrono::duration<double>(end - start).count(),
            summary.checksum,
            summary.num_traces,
            summary.total_layers,
        };
    }

    void print_timing(const char* label, const TimingResult& cpp,
                      const TimingResult& rust) {
        const double cpp_us_per_ray = cpp.seconds / cpp.num_traces * 1e6;
        const double rust_us_per_ray = rust.seconds / rust.num_traces * 1e6;
        std::cout << "\nRaytracer timing - " << label << "\n"
                  << "  traces: " << cpp.num_traces << "\n"
                  << "  C++:  " << cpp_us_per_ray << " us/ray"
                  << " (" << cpp.total_layers << " total layers)\n"
                  << "  Rust: " << rust_us_per_ray << " us/ray"
                  << " (" << rust.total_layers << " total layers)\n"
                  << "  Rust/C++: " << rust_us_per_ray / cpp_us_per_ray << "\n"
                  << "  checksums: C++=" << cpp.checksum
                  << ", Rust=" << rust.checksum << "\n";
    }

    std::vector<sasktran2::viewinggeometry::ViewingRay>
    make_timing_rays(const sasktran2::Coordinates& coordinates) {
        std::vector<sasktran2::viewinggeometry::ViewingRay> rays;
        rays.reserve(48);

        for (int i = 0; i < 16; ++i) {
            const double tangent_altitude = 5000.0 + 1500.0 * i;
            sasktran2::viewinggeometry::TangentAltitude policy(
                tangent_altitude, 0.2 + 0.01 * i, 600000.0);
            rays.push_back(policy.construct_ray(coordinates));
        }

        for (int i = 0; i < 16; ++i) {
            const double cos_viewing = 0.992 + 0.00045 * i;
            sasktran2::viewinggeometry::GroundViewingSolar policy(
                0.3, 0.2, cos_viewing, 200000.0);
            rays.push_back(policy.construct_ray(coordinates));
        }

        return rays;
    }

    void require_trace_parity(const sasktran2::raytracing::TracedRay& cpp_ray,
                              const RustSummary& rust_summary,
                              const std::vector<RustLayer>& rust_layers,
                              double earth_radius) {
        REQUIRE(rust_summary.ground_is_hit == cpp_ray.ground_is_hit);
        REQUIRE(rust_summary.is_straight == cpp_ray.is_straight);
        require_close("tangent radius", rust_summary.tangent_radius,
                      cpp_ray.tangent_radius);
        REQUIRE(rust_summary.num_layers == cpp_ray.layers.size());
        REQUIRE(rust_layers.size() == cpp_ray.layers.size());

        for (size_t i = 0; i < cpp_ray.layers.size(); ++i) {
            CAPTURE(i);
            const auto& cpp_layer = cpp_ray.layers[i];
            const auto& rust_layer = rust_layers[i];

            REQUIRE(rust_layer.layer_type == static_cast<int>(cpp_layer.type));
            require_position_close("entrance position",
                                   cpp_layer.entrance.position,
                                   rust_layer.entrance_x, rust_layer.entrance_y,
                                   rust_layer.entrance_z);
            require_position_close("exit position", cpp_layer.exit.position,
                                   rust_layer.exit_x, rust_layer.exit_y,
                                   rust_layer.exit_z);
            require_close("entrance altitude", rust_layer.entrance_altitude,
                          cpp_layer.entrance.radius() - earth_radius);
            require_close("exit altitude", rust_layer.exit_altitude,
                          cpp_layer.exit.radius() - earth_radius);
            require_close("layer distance", rust_layer.layer_distance,
                          cpp_layer.layer_distance);
            require_close("curvature factor", rust_layer.curvature_factor,
                          cpp_layer.curvature_factor, 1e-8);
            // Ground paths retain a few micrometers of accumulated C++/Rust
            // endpoint roundoff even after the geometry fields match.
            constexpr double od_abs_floor = 2e-5;
            require_close("od quad start", rust_layer.od_quad_start,
                          cpp_layer.od_quad_start, od_abs_floor);
            require_close("od quad end", rust_layer.od_quad_end,
                          cpp_layer.od_quad_end, od_abs_floor);
            require_close("od quad start fraction",
                          rust_layer.od_quad_start_fraction,
                          cpp_layer.od_quad_start_fraction, 5e-8);
            require_close("od quad end fraction",
                          rust_layer.od_quad_end_fraction,
                          cpp_layer.od_quad_end_fraction, 5e-8);
            require_close("cos sza entrance", rust_layer.cos_sza_entrance,
                          cpp_layer.cos_sza_entrance, 1e-10);
            require_close("cos sza exit", rust_layer.cos_sza_exit,
                          cpp_layer.cos_sza_exit, 1e-10);
            const double saz_abs_floor =
                std::abs(cpp_ray.observer_and_look.cos_viewing()) > 0.999999
                    ? 1e-4
                    : 1e-8;
            require_angle_close("saz entrance", rust_layer.saz_entrance,
                                cpp_layer.saz_entrance, saz_abs_floor);
            require_angle_close("saz exit", rust_layer.saz_exit,
                                cpp_layer.saz_exit, saz_abs_floor);
        }
    }

    void require_adapter_conversion(
        const sasktran2::raytracing::TracedRay& adapted_ray,
        const RustSummary& rust_summary,
        const std::vector<RustLayer>& rust_layers,
        const sasktran2::Geometry1D& geometry) {
        REQUIRE(adapted_ray.ground_is_hit == rust_summary.ground_is_hit);
        REQUIRE(adapted_ray.is_straight == rust_summary.is_straight);
        require_close("adapter tangent radius", adapted_ray.tangent_radius,
                      rust_summary.tangent_radius);
        REQUIRE(adapted_ray.layers.size() == rust_layers.size());

        std::vector<std::pair<int, double>> expected_weights;
        const auto require_weights = [](const auto& actual,
                                        const auto& expected) {
            std::vector<std::pair<int, double>> actual_nonzero;
            for (std::size_t index = 0; index < actual.size(); ++index) {
                if (actual[index].second != 0.0) {
                    actual_nonzero.push_back(actual[index]);
                }
            }
            REQUIRE(actual_nonzero == expected);
        };
        for (size_t i = 0; i < rust_layers.size(); ++i) {
            CAPTURE(i);
            const auto& layer = adapted_ray.layers[i];
            const auto& rust_layer = rust_layers[i];
            REQUIRE(static_cast<int>(layer.type) == rust_layer.layer_type);
            require_position_close("adapter entrance position",
                                   layer.entrance.position,
                                   rust_layer.entrance_x, rust_layer.entrance_y,
                                   rust_layer.entrance_z);
            require_position_close("adapter exit position", layer.exit.position,
                                   rust_layer.exit_x, rust_layer.exit_y,
                                   rust_layer.exit_z);

            const double earth_radius = geometry.coordinates().earth_radius();
            require_close("adapter entrance radius", layer.r_entrance,
                          earth_radius + rust_layer.entrance_altitude);
            require_close("adapter exit radius", layer.r_exit,
                          earth_radius + rust_layer.exit_altitude);
            require_close("adapter layer distance", layer.layer_distance,
                          rust_layer.layer_distance);
            require_close("adapter curvature factor", layer.curvature_factor,
                          rust_layer.curvature_factor, 1e-10);
            require_close("adapter od quad start", layer.od_quad_start,
                          rust_layer.od_quad_start);
            require_close("adapter od quad end", layer.od_quad_end,
                          rust_layer.od_quad_end);

            geometry.assign_interpolation_weights(layer.entrance,
                                                  expected_weights);
            require_weights(adapted_ray.entrance_weights(i), expected_weights);
            geometry.assign_interpolation_weights(layer.exit, expected_weights);
            require_weights(adapted_ray.exit_weights(i), expected_weights);
        }
    }

    void
    require_spherical_parity(const ParityGeometry& parity_geometry,
                             const sasktran2::viewinggeometry::ViewingRay& ray,
                             bool include_refraction) {
        sasktran2::raytracing::SphericalShellRayTracer cpp_tracer(
            parity_geometry.geometry);
        sasktran2::raytracing::TracedRay cpp_ray;
        cpp_tracer.trace_ray(ray, cpp_ray, include_refraction);

        std::vector<RustLayer> rust_layers;
        auto rust_summary =
            trace_rust(parity_geometry, ray, include_refraction, rust_layers);
        require_trace_parity(
            cpp_ray, rust_summary, rust_layers,
            parity_geometry.geometry.coordinates().earth_radius());

        sasktran2::raytracing::RustRayTracer rust_adapter(
            parity_geometry.geometry);
        sasktran2::raytracing::TracedRay adapted_ray;
        rust_adapter.trace_ray(ray, adapted_ray, include_refraction);
        require_adapter_conversion(adapted_ray, rust_summary, rust_layers,
                                   parity_geometry.geometry);
    }
} // namespace

TEST_CASE("Rust raytracer parity - spherical limb",
          "[sasktran2][raytracing][rust]") {
    auto parity_geometry = make_geometry();

    double tangent_altitude = GENERATE(5000.0, 10000.0, 30405.0);
    bool include_refraction = GENERATE(false, true);
    CAPTURE(tangent_altitude);
    CAPTURE(include_refraction);
    sasktran2::viewinggeometry::TangentAltitude ray_policy(tangent_altitude,
                                                           0.2, 600000.0);
    auto ray = ray_policy.construct_ray(parity_geometry.geometry.coordinates());
    require_spherical_parity(parity_geometry, ray, include_refraction);
}

TEST_CASE("Rust raytracer parity - spherical ground",
          "[sasktran2][raytracing][rust]") {
    auto parity_geometry = make_geometry();

    double cos_viewing_zenith = GENERATE(0.9999995, 0.999, 0.995);
    bool include_refraction = GENERATE(false, true);
    CAPTURE(cos_viewing_zenith);
    CAPTURE(include_refraction);
    sasktran2::viewinggeometry::GroundViewingSolar ray_policy(
        0.3, 0.2, cos_viewing_zenith, 200000.0);
    auto ray = ray_policy.construct_ray(parity_geometry.geometry.coordinates());
    require_spherical_parity(parity_geometry, ray, include_refraction);
}

TEST_CASE("Rust raytracer parity - spherical inside limb",
          "[sasktran2][raytracing][rust]") {
    auto parity_geometry = make_geometry();

    double tangent_altitude = GENERATE(5000.0, 10000.0);
    bool include_refraction = GENERATE(false, true);
    CAPTURE(tangent_altitude);
    CAPTURE(include_refraction);
    sasktran2::viewinggeometry::TangentAltitude ray_policy(tangent_altitude,
                                                           0.2, 30500.0);
    auto ray = ray_policy.construct_ray(parity_geometry.geometry.coordinates());
    require_spherical_parity(parity_geometry, ray, include_refraction);
}

TEST_CASE("Rust raytracer parity - spherical inside ground",
          "[sasktran2][raytracing][rust]") {
    auto parity_geometry = make_geometry();

    double cos_viewing_zenith = GENERATE(0.999, 0.995);
    bool include_refraction = GENERATE(false, true);
    CAPTURE(cos_viewing_zenith);
    CAPTURE(include_refraction);
    sasktran2::viewinggeometry::GroundViewingSolar ray_policy(
        0.3, 0.2, cos_viewing_zenith, 30500.0);
    auto ray = ray_policy.construct_ray(parity_geometry.geometry.coordinates());
    require_spherical_parity(parity_geometry, ray, include_refraction);
}

TEST_CASE("Rust raytracer parity - variable grids and interpolation",
          "[sasktran2][raytracing][rust][edge]") {
    auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                  sasktran2::grids::interpolation::shell,
                                  sasktran2::grids::interpolation::lower);
    bool include_refraction = GENERATE(false, true);
    CAPTURE(interpolation);
    CAPTURE(include_refraction);

    auto parity_geometry = make_geometry_with_grid(
        {0.0, 1000.0, 1700.0, 6300.0, 14000.0, 27000.0, 50000.0, 80000.0,
         120000.0},
        sasktran2::grids::gridspacing::variable, interpolation,
        {1.00030, 1.00024, 1.00018, 1.00013, 1.00009, 1.00006, 1.00004, 1.00002,
         1.00000});

    std::vector<sasktran2::viewinggeometry::ViewingRay> rays;
    rays.push_back(
        sasktran2::viewinggeometry::TangentAltitude(0.1, 0.2, 600000.0)
            .construct_ray(parity_geometry.geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::TangentAltitude(3500.0, 0.4, 600000.0)
            .construct_ray(parity_geometry.geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::GroundViewingSolar(0.3, 0.2, 0.99, 60000.0)
            .construct_ray(parity_geometry.geometry.coordinates()));

    for (size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        CAPTURE(ray_index);
        const auto& ray = rays[ray_index];
        require_spherical_parity(parity_geometry, ray, include_refraction);
    }
}

TEST_CASE("Rust raytracer parity - spherical boundary and empty rays",
          "[sasktran2][raytracing][rust][edge]") {
    auto parity_geometry = make_geometry();
    constexpr bool include_refraction = false;

    std::vector<sasktran2::viewinggeometry::ViewingRay> rays;
    for (double tangent_altitude :
         {-1.0, 0.0001, 50000.0, 65000.000001, 70000.0}) {
        rays.push_back(
            sasktran2::viewinggeometry::TangentAltitude(tangent_altitude, 0.2,
                                                        600000.0)
                .construct_ray(parity_geometry.geometry.coordinates()));
    }
    rays.push_back(
        sasktran2::viewinggeometry::ViewingUpSolar(0.3, 0.2, 0.9, 600000.0)
            .construct_ray(parity_geometry.geometry.coordinates()));

    for (size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        CAPTURE(ray_index);
        const auto& ray = rays[ray_index];
        require_spherical_parity(parity_geometry, ray, include_refraction);
    }
}

TEST_CASE("Rust raytracer timing - spherical batch",
          "[sasktran2][rust][.benchmark]") {
    auto parity_geometry = make_geometry();
    sasktran2::raytracing::SphericalShellRayTracer cpp_tracer(
        parity_geometry.geometry);
    sasktran2::raytracing::RustRayTracer rust_adapter(parity_geometry.geometry);
    auto rays = make_timing_rays(parity_geometry.geometry.coordinates());

    constexpr size_t straight_iterations = 3000;
    constexpr size_t refracted_iterations = 300;

    auto cpp_straight =
        time_cpp_tracer(cpp_tracer, rays, false, straight_iterations);
    auto rust_straight =
        time_rust_tracer(parity_geometry, rays, false, straight_iterations);
    REQUIRE(rust_straight.num_traces == cpp_straight.num_traces);
    REQUIRE(rust_straight.total_layers == cpp_straight.total_layers);
    REQUIRE(std::isfinite(rust_straight.checksum));
    print_timing("straight spherical", cpp_straight, rust_straight);
    auto rust_adapter_straight =
        time_cpp_tracer(rust_adapter, rays, false, straight_iterations);
    REQUIRE(rust_adapter_straight.total_layers == cpp_straight.total_layers);
    REQUIRE(std::isfinite(rust_adapter_straight.checksum));
    print_timing("straight spherical adapter", cpp_straight,
                 rust_adapter_straight);

    auto cpp_refracted =
        time_cpp_tracer(cpp_tracer, rays, true, refracted_iterations);
    auto rust_refracted =
        time_rust_tracer(parity_geometry, rays, true, refracted_iterations);
    REQUIRE(rust_refracted.num_traces == cpp_refracted.num_traces);
    REQUIRE(rust_refracted.total_layers == cpp_refracted.total_layers);
    REQUIRE(std::isfinite(rust_refracted.checksum));
    print_timing("refracted spherical", cpp_refracted, rust_refracted);
    auto rust_adapter_refracted =
        time_cpp_tracer(rust_adapter, rays, true, refracted_iterations);
    REQUIRE(rust_adapter_refracted.total_layers == cpp_refracted.total_layers);
    REQUIRE(std::isfinite(rust_adapter_refracted.checksum));
    print_timing("refracted spherical adapter", cpp_refracted,
                 rust_adapter_refracted);
}

#endif
