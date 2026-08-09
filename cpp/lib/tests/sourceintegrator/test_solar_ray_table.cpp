#include <sasktran2/test_helper.h>

#include <sasktran2.h>
#include <sasktran2/solartransmission.h>

namespace {
    Eigen::VectorXd altitude_grid(std::initializer_list<double> values) {
        return Eigen::Map<const Eigen::VectorXd>(values.begin(), values.size());
    }

    Eigen::VectorXd
    direct_path_weights(const sasktran2::raytracing::TracedRay& ray,
                        const sasktran2::Geometry1D& geometry) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(geometry.size());
        for (std::size_t layer_index = 0; layer_index < ray.layers.size();
             ++layer_index) {
            const auto weights = ray.optical_depth_weights(layer_index);
            for (std::size_t index = 0; index < weights.size(); ++index) {
                const auto weight = weights[index];
                result[weight.first] += weight.second;
            }
        }
        return result;
    }

    Eigen::VectorXd table_path_weights(
        const sasktran2::solartransmission::SolarTransmissionRayTable& table,
        const sasktran2::solartransmission::SolarTransmissionRayTable::Query&
            query,
        const sasktran2::Geometry1D& geometry) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(geometry.size());
        for (std::size_t sample_index = 0; sample_index < query.count;
             ++sample_index) {
            const auto& sample = query.samples[sample_index];
            const auto& ray = table.rays()[sample.ray_index];
            for (int layer_index = sample.boundary_index;
                 layer_index < ray.layers.size(); ++layer_index) {
                const auto weights = ray.optical_depth_weights(layer_index);
                for (std::size_t index = 0; index < weights.size(); ++index) {
                    const auto weight = weights[index];
                    result[weight.first] +=
                        sample.interpolation_weight * weight.second;
                }
            }
            for (std::size_t index = 0; index < sample.partial.count; ++index) {
                result[sample.partial.index[index]] +=
                    sample.interpolation_weight * sample.partial.weight[index];
            }
        }
        return result;
    }

    sasktran2::Location
    solar_ray_location(const sasktran2::Geometry1D& geometry,
                       double impact_parameter, double altitude,
                       bool sun_facing = true) {
        const auto sun = geometry.coordinates().sun_unit().normalized();
        const auto perpendicular = sun.unitOrthogonal().normalized();
        const double radius = geometry.coordinates().earth_radius() + altitude;
        const double coordinate =
            (sun_facing ? 1.0 : -1.0) *
            std::sqrt(std::max(0.0, radius * radius -
                                        impact_parameter * impact_parameter));
        sasktran2::Location result;
        result.position = impact_parameter * perpendicular + coordinate * sun;
        return result;
    }
} // namespace

TEST_CASE("Solar ray table analytically evaluates partial shells",
          "[sourceintegrator][singlescatter][solarraytable]") {
    constexpr double earth_radius = 6'372'000.0;
    const sasktran2::Geometry1D geometry(
        0.2, 0.0, earth_radius,
        altitude_grid({0.0, 800.0, 2'500.0, 7'000.0, 15'000.0, 30'000.0,
                       50'000.0, 80'000.0}),
        sasktran2::grids::interpolation::linear,
        sasktran2::geometrytype::spherical);
    const sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    sasktran2::solartransmission::SolarTransmissionRayTable table(geometry,
                                                                  raytracer);
    table.initialize();

    // This impact parameter is an exact table ray, while the location lies
    // strictly inside another shell.  The comparison therefore isolates the
    // analytic partial-shell evaluation from impact-parameter interpolation.
    const double impact_parameter = earth_radius + 7'000.0;
    const auto location =
        solar_ray_location(geometry, impact_parameter, 21'000.0);
    const auto query = table.query(location);
    REQUIRE_FALSE(query.shadowed);
    REQUIRE(query.count == 1);

    sasktran2::viewinggeometry::ViewingRay direct_ray;
    direct_ray.observer = location;
    direct_ray.look_away = geometry.coordinates().sun_unit();
    direct_ray.relative_azimuth = 0.0;
    sasktran2::raytracing::TracedRay traced;
    raytracer.trace_ray(direct_ray, traced, false);
    REQUIRE_FALSE(traced.ground_is_hit);

    const auto expected = direct_path_weights(traced, geometry);
    const auto actual = table_path_weights(table, query, geometry);
    REQUIRE(actual.isApprox(expected, 2e-10));
}

TEST_CASE("Solar ray table interpolates impact parameter and preserves shadow",
          "[sourceintegrator][singlescatter][solarraytable]") {
    constexpr double earth_radius = 6'372'000.0;
    const sasktran2::Geometry1D geometry(
        0.2, 0.0, earth_radius,
        altitude_grid({0.0, 1'000.0, 2'000.0, 4'000.0, 8'000.0, 16'000.0,
                       32'000.0, 64'000.0, 80'000.0}),
        sasktran2::grids::interpolation::linear,
        sasktran2::geometrytype::spherical);
    const sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    sasktran2::solartransmission::SolarTransmissionRayTable table(geometry,
                                                                  raytracer);
    table.initialize();

    const double impact_parameter = earth_radius + 11'300.0;
    const auto location =
        solar_ray_location(geometry, impact_parameter, 24'000.0);
    const auto query = table.query(location);
    REQUIRE_FALSE(query.shadowed);
    REQUIRE(query.count == 2);

    sasktran2::viewinggeometry::ViewingRay direct_ray;
    direct_ray.observer = location;
    direct_ray.look_away = geometry.coordinates().sun_unit();
    direct_ray.relative_azimuth = 0.0;
    sasktran2::raytracing::TracedRay traced;
    raytracer.trace_ray(direct_ray, traced, false);
    const auto expected = direct_path_weights(traced, geometry);
    const auto actual = table_path_weights(table, query, geometry);
    const Eigen::VectorXd extinction =
        (Eigen::VectorXd(geometry.size()) << 1.0, 0.95, 0.8, 0.7, 0.5, 0.25,
         0.1, 0.03, 0.0)
            .finished();
    REQUIRE(actual.dot(extinction) ==
            Catch::Approx(expected.dot(extinction)).epsilon(3e-4));

    const auto shadowed = table.query(
        solar_ray_location(geometry, earth_radius - 20'000.0, 5'000.0, false));
    REQUIRE(shadowed.shadowed);
}

TEST_CASE("Solar ray table tracks direct spherical slant optical depth",
          "[sourceintegrator][singlescatter][solarraytable]") {
    constexpr double earth_radius = 6'372'000.0;
    Eigen::VectorXd altitudes(201);
    Eigen::VectorXd extinction(201);
    for (int index = 0; index < altitudes.size(); ++index) {
        altitudes[index] = 500.0 * index;
        extinction[index] = 1.0e-4 * std::exp(-altitudes[index] / 7'000.0);
    }
    const sasktran2::Geometry1D geometry(
        0.2, 0.0, earth_radius, Eigen::VectorXd(altitudes),
        sasktran2::grids::interpolation::linear,
        sasktran2::geometrytype::spherical);
    const sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    sasktran2::solartransmission::SolarTransmissionRayTable table(geometry,
                                                                  raytracer);
    table.initialize();

    const auto sun = geometry.coordinates().sun_unit().normalized();
    const auto perpendicular = sun.unitOrthogonal().normalized();
    double maximum_relative_error = 0.0;
    double maximum_absolute_od_error = 0.0;
    double maximum_transmission_error = 0.0;
    double worst_mu = 0.0;
    double worst_altitude = 0.0;
    for (const double altitude :
         {250.0, 2'250.0, 10'250.0, 30'250.0, 80'250.0}) {
        for (const double mu : {-0.1, -0.05, 0.05, 0.1, 0.2, 0.5, 0.9}) {
            const double radius = earth_radius + altitude;
            sasktran2::Location location;
            location.position =
                radius * (std::sqrt(1.0 - mu * mu) * perpendicular + mu * sun);
            const auto query = table.query(location);
            const double impact_parameter = radius * std::sqrt(1.0 - mu * mu);
            if (mu < 0.0 && impact_parameter < earth_radius) {
                REQUIRE(query.shadowed);
                continue;
            }
            REQUIRE_FALSE(query.shadowed);

            sasktran2::viewinggeometry::ViewingRay direct_ray;
            direct_ray.observer = location;
            direct_ray.look_away = sun;
            sasktran2::raytracing::TracedRay traced;
            raytracer.trace_ray(direct_ray, traced, false);
            const double expected =
                direct_path_weights(traced, geometry).dot(extinction);
            const double actual =
                table_path_weights(table, query, geometry).dot(extinction);
            const double relative_error =
                std::abs(actual - expected) / std::max(expected, 1e-12);
            maximum_absolute_od_error = std::max(maximum_absolute_od_error,
                                                 std::abs(actual - expected));
            maximum_transmission_error =
                std::max(maximum_transmission_error,
                         std::abs(std::exp(-actual) - std::exp(-expected)));
            if (relative_error > maximum_relative_error) {
                maximum_relative_error = relative_error;
                worst_mu = mu;
                worst_altitude = altitude;
            }
        }
    }
    INFO("maximum relative OD error: " << maximum_relative_error);
    INFO("maximum absolute OD error: " << maximum_absolute_od_error);
    INFO("maximum transmission error: " << maximum_transmission_error);
    INFO("worst mu/altitude: " << worst_mu << " / " << worst_altitude);
    REQUIRE(maximum_transmission_error < 2e-4);
}
