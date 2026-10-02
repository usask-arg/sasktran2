#include <sasktran2/test_helper.h>

#include <sasktran2.h>
#include <sasktran2/refraction.h>

namespace {
    constexpr double earth_radius_m = 6372000.0;
    constexpr double top_altitude_m = 65000.0;
    constexpr double altitude_spacing_m = 1000.0;

    // The logarithm of the refractive index is linear in altitude over the
    // whole grid, so the log-linear interpolation used by the geometry is
    // exact and the reference integration below never crosses a
    // discontinuity in the refractive index gradient. The refractive index is
    // exactly one at the top of the atmosphere.
    constexpr double surface_log_refractive_index = 2.8e-4;
    constexpr double log_refractive_index_gradient =
        -surface_log_refractive_index / top_altitude_m;

    double log_refractive_index(double altitude_m) {
        return surface_log_refractive_index *
               (1.0 - altitude_m / top_altitude_m);
    }

    double extinction_at_level(double altitude_m) {
        return 1e-5 * std::exp(-altitude_m / 7000.0);
    }

    // Extinction linearly interpolated between the grid levels
    double extinction(double altitude_m) {
        const double clamped =
            std::clamp(altitude_m, 0.0, top_altitude_m - 1e-9);
        const double lower =
            std::floor(clamped / altitude_spacing_m) * altitude_spacing_m;
        const double weight = (clamped - lower) / altitude_spacing_m;
        return (1.0 - weight) * extinction_at_level(lower) +
               weight * extinction_at_level(lower + altitude_spacing_m);
    }

    sasktran2::Geometry1D refracting_geometry() {
        const int num_levels =
            static_cast<int>(top_altitude_m / altitude_spacing_m) + 1;
        Eigen::VectorXd altitudes =
            Eigen::VectorXd::LinSpaced(num_levels, 0.0, top_altitude_m);

        sasktran2::grids::AltitudeGrid grid(
            std::move(altitudes), sasktran2::grids::gridspacing::constant,
            sasktran2::grids::outofbounds::extend,
            sasktran2::grids::interpolation::linear);
        sasktran2::Coordinates coordinates(0.5, 0, earth_radius_m);
        sasktran2::Geometry1D geometry(std::move(coordinates),
                                       std::move(grid));

        const auto& levels = geometry.altitude_grid().grid();
        for (Eigen::Index i = 0; i < levels.size(); ++i) {
            geometry.refractive_index()(i) =
                std::exp(log_refractive_index(levels(i)));
        }
        return geometry;
    }

    struct ReferenceRay {
        bool hits_ground = false;
        double minimum_radius = 0.0;
        // Angle between the start position and the ray direction after
        // leaving the atmosphere
        double asymptotic_angle = 0.0;
        double optical_depth = 0.0;
    };

    /**
     * Integrates the ray equation in the plane of the ray with RK4, using the
     * radius, the swept central angle and the local zenith angle as the
     * state. Independent of the path integrals and layering used by the ray
     * tracer.
     */
    ReferenceRay integrate_reference_ray(double radius, double zenith,
                                         double step_m = 5.0) {
        const double ground_radius = earth_radius_m;
        const double top_radius = earth_radius_m + top_altitude_m;

        // radius, swept angle, zenith, optical depth
        using State = std::array<double, 4>;
        const auto derivative = [](const State& y) {
            const double sin_zenith = std::sin(y[2]);
            return State{
                std::cos(y[2]), sin_zenith / y[0],
                -sin_zenith * (1.0 / y[0] + log_refractive_index_gradient),
                extinction(y[0] - earth_radius_m)};
        };
        const auto rk4_step = [&](const State& y, double h) {
            const auto add = [](const State& a, const State& b, double s) {
                return State{a[0] + s * b[0], a[1] + s * b[1],
                             a[2] + s * b[2], a[3] + s * b[3]};
            };
            const State k1 = derivative(y);
            const State k2 = derivative(add(y, k1, h / 2.0));
            const State k3 = derivative(add(y, k2, h / 2.0));
            const State k4 = derivative(add(y, k3, h));
            State result;
            for (int i = 0; i < 4; ++i) {
                result[i] =
                    y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            return result;
        };

        ReferenceRay result;
        State y = {radius, 0.0, zenith, 0.0};
        result.minimum_radius = radius;
        while (true) {
            State next = rk4_step(y, step_m);
            if (next[0] < ground_radius) {
                result.hits_ground = true;
                return result;
            }
            if (next[0] >= top_radius) {
                // Land on the top of the atmosphere, the remaining distance
                // is short enough that a straight line estimate converges
                for (int i = 0; i < 3; ++i) {
                    const double remaining =
                        (top_radius - y[0]) / std::cos(y[2]);
                    next = rk4_step(y, remaining);
                    if (std::abs(next[0] - top_radius) < 1e-9) {
                        break;
                    }
                }
                y = next;
                break;
            }
            y = next;
            result.minimum_radius = std::min(result.minimum_radius, y[0]);
        }
        // The refractive index is one at the top of the atmosphere, so the
        // direction is unchanged when leaving it
        result.asymptotic_angle = y[1] + y[2];
        result.optical_depth = y[3];
        return result;
    }

    double zenith_angle(const Eigen::Vector3d& position,
                        const Eigen::Vector3d& direction) {
        const Eigen::Vector3d up = position.normalized();
        return std::atan2(up.cross(direction).norm(), up.dot(direction));
    }
} // namespace

TEST_CASE("Refraction - Apparent sun direction reaches the sun",
          "[sasktran2][raytracing][refraction]") {
    const double altitude = GENERATE(0.0, 500.0, 5000.0, 20000.0, 40000.0);
    const double geometric_sza_deg =
        GENERATE(60.0, 85.0, 89.0, 90.0, 90.5, 91.0, 92.0, 93.0, 95.0);

    auto geometry = refracting_geometry();
    const Eigen::Vector3d& sun = geometry.coordinates().sun_unit();
    const double geometric_sza = geometric_sza_deg * EIGEN_PI / 180.0;
    const Eigen::Vector3d position =
        geometry.coordinates().solar_coordinate_vector(std::cos(geometric_sza),
                                                       0.0, altitude);
    const double radius = position.norm();

    Eigen::Vector3d direction_to_sun;
    std::vector<std::pair<int, double>> index_weights;
    const bool illuminated =
        sasktran2::raytracing::refraction::refracted_direction_to_sun(
            geometry, position, direction_to_sun, index_weights);

    // Largest zenith angle that clears the surface follows from the ray
    // invariant n r sin(zenith)
    const double surface_invariant =
        std::exp(log_refractive_index(0.0)) * earth_radius_m;
    const double start_invariant =
        std::exp(log_refractive_index(altitude)) * radius;
    const double grazing_zenith =
        start_invariant > surface_invariant
            ? EIGEN_PI - std::asin(surface_invariant / start_invariant)
            : EIGEN_PI / 2.0;
    const auto grazing =
        integrate_reference_ray(radius, grazing_zenith - 1e-7);

    CAPTURE(altitude, geometric_sza_deg, illuminated,
            grazing.asymptotic_angle - geometric_sza, grazing.hits_ground);
    if (!illuminated) {
        // Even the grazing ray does not bend far enough to reach the sun
        REQUIRE(grazing.asymptotic_angle < geometric_sza);
        return;
    }
    if (!grazing.hits_ground) {
        REQUIRE(grazing.asymptotic_angle >= geometric_sza - 1e-8);
    }

    // The apparent direction is in the vertical plane containing the sun,
    // and the sun appears higher than its geometric position
    REQUIRE(std::abs(direction_to_sun.norm() - 1.0) < 1e-12);
    REQUIRE(std::abs(direction_to_sun.dot(position.cross(sun).normalized())) <
            1e-12);
    const double apparent_zenith = zenith_angle(position, direction_to_sun);
    REQUIRE(apparent_zenith <= geometric_sza + 1e-12);

    const auto reference = integrate_reference_ray(radius, apparent_zenith);
    CAPTURE(altitude, geometric_sza_deg, apparent_zenith,
            reference.asymptotic_angle - geometric_sza);
    REQUIRE(!reference.hits_ground);
    REQUIRE(std::abs(reference.asymptotic_angle - geometric_sza) < 1e-8);
}

TEST_CASE("Refraction - Refracted solar ray optical depth",
          "[sasktran2][raytracing][refraction]") {
    const double altitude = GENERATE(0.0, 5000.0, 20000.0, 40000.0);
    const double geometric_sza_deg = GENERATE(70.0, 88.0, 90.0, 91.5, 93.0);

    auto geometry = refracting_geometry();
    const double geometric_sza = geometric_sza_deg * EIGEN_PI / 180.0;
    const Eigen::Vector3d position =
        geometry.coordinates().solar_coordinate_vector(std::cos(geometric_sza),
                                                       0.0, altitude);

    sasktran2::viewinggeometry::ViewingRay ray_to_sun;
    ray_to_sun.observer.position = position;
    std::vector<std::pair<int, double>> index_weights;
    if (!sasktran2::raytracing::refraction::refracted_direction_to_sun(
            geometry, position, ray_to_sun.look_away, index_weights)) {
        return;
    }

    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    sasktran2::raytracing::TracedRay traced;
    raytracer.trace_ray(ray_to_sun, traced, true);
    REQUIRE(!traced.ground_is_hit);

    const auto& levels = geometry.altitude_grid().grid();
    double optical_depth = 0.0;
    for (std::size_t layer = 0; layer < traced.layers.size(); ++layer) {
        const auto weights = traced.optical_depth_weights(layer);
        for (std::size_t i = 0; i < weights.size(); ++i) {
            optical_depth +=
                weights[i].second * extinction_at_level(levels(weights[i].first));
        }
    }

    const auto reference = integrate_reference_ray(
        position.norm(), zenith_angle(position, ray_to_sun.look_away));
    CAPTURE(altitude, geometric_sza_deg, optical_depth,
            reference.optical_depth,
            optical_depth / reference.optical_depth - 1.0);
    REQUIRE(!reference.hits_ground);
    // The ray tracer evaluates the extinction of a refracted layer along the
    // straight chord between its end points, which lies slightly below the
    // curved path. That is a ~2e-4 effect here; errors in the ray invariant
    // or the layering move the path by kilometres and change the optical
    // depth by percent.
    REQUIRE(std::abs(optical_depth - reference.optical_depth) <
            5e-4 * reference.optical_depth);

    // The traced ray leaves the atmosphere where the reference ray does
    const Eigen::Vector3d& exit = traced.layers.front().exit.position;
    REQUIRE(std::abs(exit.norm() - (earth_radius_m + top_altitude_m)) < 1e-3);
}

TEST_CASE("Refraction - Observer inside the atmosphere conserves the ray "
          "invariant",
          "[sasktran2][raytracing][refraction]") {
    const double altitude = GENERATE(500.0, 5000.0, 12345.0, 30000.0);
    const double zenith_deg = GENERATE(91.0, 93.0, 95.0);

    auto geometry = refracting_geometry();
    const double zenith = zenith_deg * EIGEN_PI / 180.0;

    sasktran2::viewinggeometry::ViewingRay ray;
    ray.observer.position =
        geometry.coordinates().solar_coordinate_vector(0.5, 0.0, altitude);
    const Eigen::Vector3d up = ray.observer.position.normalized();
    const Eigen::Vector3d horizontal =
        (geometry.coordinates().sun_unit() -
         geometry.coordinates().sun_unit().dot(up) * up)
            .normalized();
    ray.look_away = std::cos(zenith) * up + std::sin(zenith) * horizontal;

    const auto reference =
        integrate_reference_ray(ray.observer.position.norm(), zenith, 1.0);
    if (reference.hits_ground) {
        return;
    }

    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    sasktran2::raytracing::TracedRay traced;
    raytracer.trace_ray(ray, traced, true);
    REQUIRE(!traced.ground_is_hit);

    CAPTURE(altitude, zenith_deg, traced.tangent_radius,
            reference.minimum_radius);
    // The refracted tangent point matches the integrated ray
    REQUIRE(std::abs(traced.tangent_radius - reference.minimum_radius) < 0.05);

    // and the traced layers reach down to it
    double minimum_radius = std::numeric_limits<double>::max();
    for (const auto& layer : traced.layers) {
        minimum_radius = std::min({minimum_radius, layer.entrance.radius(),
                                   layer.exit.radius()});
    }
    REQUIRE(std::abs(minimum_radius - traced.tangent_radius) < 1e-3);
}

