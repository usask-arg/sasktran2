#include <sasktran2/test_helper.h>

#include <sasktran2.h>
#include <sasktran2/solartransmission.h>

namespace {
    constexpr double earth_radius = 10.0;

    Eigen::VectorXd vector(std::initializer_list<double> values) {
        Eigen::VectorXd result(values.size());
        int index = 0;
        for (const double value : values) {
            result[index++] = value;
        }
        return result;
    }

    sasktran2::Geometry2D geometry2d() {
        return {0.6, 0.3, earth_radius, vector({0.0, 10.0, 20.0}),
                vector({-0.5, 0.0, 0.5})};
    }

    Eigen::VectorXd
    direct_path_weights(const sasktran2::raytracing::TracedRay& ray,
                        const sasktran2::Geometry2D& geometry) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(geometry.size());
        for (std::size_t layer_index = 0; layer_index < ray.layers.size();
             ++layer_index) {
            const auto weights = ray.optical_depth_weights(layer_index);
            for (std::size_t index = 0; index < weights.size(); ++index) {
                result[weights[index].first] += weights[index].second;
            }
        }
        return result;
    }

    sasktran2::viewinggeometry::ViewingRay
    viewing_ray(const Eigen::Vector3d& observer, const Eigen::Vector3d& look) {
        sasktran2::viewinggeometry::ViewingRay result;
        result.observer.position = observer;
        result.look_away = look.normalized();
        result.relative_azimuth = 0.0;
        return result;
    }
} // namespace

TEST_CASE("Solar geometry storage uses compact indices with a wide fallback",
          "[sourceintegrator][singlescatter]") {
    for (const Eigen::Index columns : {Eigen::Index{10}, Eigen::Index{65537}}) {
        sasktran2::solartransmission::SolarGeometryMatrix matrix;
        matrix.initialize_exact(2, columns, 4);
        REQUIRE(matrix.is_compact() == (columns == 10));

        matrix.start_row(0);
        matrix.insert_back(0, 0, 2.0);
        matrix.insert_back(0, columns - 1, 3.0);
        matrix.start_row(1);
        matrix.insert_back(1, 1, -4.0);
        matrix.finalize();

        Eigen::VectorXd input = Eigen::VectorXd::Zero(columns);
        input[0] = 5.0;
        input[1] = 7.0;
        input[columns - 1] = 11.0;
        REQUIRE(matrix.row_dot(0, input) == Catch::Approx(43.0));

        Eigen::VectorXd product(2);
        matrix.multiply(input, product);
        REQUIRE(product[0] == Catch::Approx(43.0));
        REQUIRE(product[1] == Catch::Approx(-28.0));

        Eigen::VectorXd accumulated = Eigen::VectorXd::Zero(columns);
        matrix.accumulate_row(0, 0.5, accumulated);
        REQUIRE(accumulated[0] == Catch::Approx(1.0));
        REQUIRE(accumulated[columns - 1] == Catch::Approx(1.5));

        int entries = 0;
        for (sasktran2::solartransmission::SolarGeometryMatrix::InnerIterator
                 entry(matrix, 0);
             entry; ++entry) {
            ++entries;
        }
        REQUIRE(entries == 2);
    }
}

#ifdef SKTRAN_RUST_SUPPORT
TEST_CASE("Exact solar transmission matrix matches direct Geometry2D rays",
          "[sourceintegrator][singlescatter][geometry2d]") {
    const auto geometry = geometry2d();
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    sasktran2::raytracing::TracedRay line_of_sight;
    raytracer.trace_ray(
        viewing_ray({0.0, 0.0, earth_radius + 25.0}, {0.0, 0.0, -1.0}),
        line_of_sight);
    REQUIRE(!line_of_sight.layers.empty());

    sasktran2::solartransmission::SolarTransmissionExact transmission(
        geometry, raytracer);
    sasktran2::solartransmission::SolarGeometryMatrix matrix;
    std::vector<bool> ground_hit;
    transmission.generate_geometry_matrix({line_of_sight}, matrix, ground_hit);

    REQUIRE(matrix.rows() == line_of_sight.layers.size() + 1);
    REQUIRE(matrix.cols() == geometry.size());
    REQUIRE(ground_hit.size() == line_of_sight.layers.size() + 1);

    sasktran2::viewinggeometry::ViewingRay solar_ray;
    solar_ray.look_away = geometry.coordinates().sun_unit();
    sasktran2::raytracing::TracedRay traced_solar_ray;

    int row = 0;
    for (int layer_index = 0; layer_index < line_of_sight.layers.size();
         ++layer_index) {
        const auto check_endpoint = [&](const sasktran2::Location& endpoint) {
            solar_ray.observer = endpoint;
            raytracer.trace_ray(solar_ray, traced_solar_ray);
            REQUIRE(ground_hit[row] == traced_solar_ray.ground_is_hit);

            Eigen::VectorXd actual = Eigen::VectorXd::Zero(geometry.size());
            for (sasktran2::solartransmission::SolarGeometryMatrix::
                     InnerIterator entry(matrix, row);
                 entry; ++entry) {
                actual[entry.index()] = entry.value();
            }
            if (traced_solar_ray.ground_is_hit) {
                REQUIRE(actual.isZero(0.0));
            } else {
                const Eigen::VectorXd expected =
                    direct_path_weights(traced_solar_ray, geometry);
                REQUIRE(actual.isApprox(expected, 1e-13));
            }
            ++row;
        };

        if (layer_index == 0) {
            check_endpoint(line_of_sight.layers[layer_index].exit);
        }
        check_endpoint(line_of_sight.layers[layer_index].entrance);
    }
    REQUIRE(row == matrix.rows());
}

TEST_CASE("Geometry2D characteristic solar table follows exact endpoint OD",
          "[sourceintegrator][singlescatter][geometry2d][solartable]") {
    const auto geometry = geometry2d();
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    sasktran2::raytracing::TracedRay line_of_sight;
    line_of_sight.layers.resize(1);
    const double query_radius = earth_radius + 20.0;
    const double impact_parameter = earth_radius + 10.0;
    const double query_cos_sza =
        -std::sqrt(1.0 - impact_parameter * impact_parameter /
                             (query_radius * query_radius));
    line_of_sight.layers[0].entrance.position =
        geometry.coordinates().solar_coordinate_vector(query_cos_sza, 0.0,
                                                       20.0);
    line_of_sight.layers[0].exit = line_of_sight.layers[0].entrance;

    sasktran2::Config config;
    sasktran2::solartransmission::SolarTransmissionTable2D table(geometry,
                                                                 raytracer);
    table.initialize_config(config);
    table.initialize_geometry({line_of_sight});
    sasktran2::solartransmission::SolarTableInterpolation interpolation;
    std::vector<bool> table_ground_hit;
    table.generate_interpolation({line_of_sight}, interpolation,
                                 table_ground_hit);

    Eigen::VectorXd extinction(geometry.size());
    for (int index = 0; index < geometry.size(); ++index) {
        extinction[index] = 0.01 * (1.0 + 0.07 * index);
    }
    Eigen::VectorXd table_nodes(table.table_size());
    table.apply(extinction, table_nodes);
    Eigen::VectorXd table_od(interpolation.rows());
    interpolation.apply(table_nodes, table_od);

    sasktran2::solartransmission::SolarTransmissionExact exact(geometry,
                                                               raytracer);
    sasktran2::solartransmission::SolarGeometryMatrix exact_matrix;
    std::vector<bool> exact_ground_hit;
    exact.generate_geometry_matrix({line_of_sight}, exact_matrix,
                                   exact_ground_hit);
    Eigen::VectorXd exact_od(exact_matrix.rows());
    exact_matrix.multiply(extinction, exact_od);

    REQUIRE(table_ground_hit == exact_ground_hit);
    for (int row = 0; row < exact_od.size(); ++row) {
        if (!exact_ground_hit[row]) {
            INFO("row " << row << ", exact " << exact_od[row] << ", table "
                        << table_od[row]);
            REQUIRE(table_od[row] ==
                    Catch::Approx(exact_od[row]).epsilon(1e-11));
        }
    }

    const Eigen::VectorXd endpoint_cotangent =
        Eigen::VectorXd::LinSpaced(interpolation.rows(), -0.7, 1.1);
    Eigen::VectorXd table_cotangent(table.table_size());
    interpolation.apply_transpose(endpoint_cotangent, table_cotangent);
    REQUIRE(endpoint_cotangent.dot(table_od) ==
            Catch::Approx(table_cotangent.dot(table_nodes)).epsilon(1e-12));

    Eigen::VectorXd extinction_cotangent =
        Eigen::VectorXd::Zero(geometry.size());
    table.accumulate_transpose(table_cotangent, extinction_cotangent, 1.0);
    REQUIRE(table_cotangent.dot(table_nodes) ==
            Catch::Approx(extinction_cotangent.dot(extinction)).epsilon(1e-12));
}

TEST_CASE("Geometry2D solar table permits azimuth-local finite-window topology",
          "[sourceintegrator][singlescatter][geometry2d][solartable]") {
    constexpr int num_altitudes = 81;
    constexpr int num_horizontal = 68;
    Eigen::VectorXd altitudes =
        Eigen::VectorXd::LinSpaced(num_altitudes, 0.0, 80000.0);
    Eigen::VectorXd horizontal =
        Eigen::VectorXd::LinSpaced(num_horizontal, -0.9, 0.25);
    sasktran2::Geometry2D geometry(0.35, 0.8, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    std::vector<sasktran2::raytracing::TracedRay> queries;
    constexpr int num_azimuths = (num_horizontal + 1) / 2;
    for (int azimuth_index = 0; azimuth_index < num_azimuths; ++azimuth_index) {
        sasktran2::raytracing::TracedRay query;
        query.layers.resize(1);
        query.layers[0].entrance.position =
            geometry.coordinates().solar_coordinate_vector(
                -0.07908409944216922,
                2.0 * EIGEN_PI * azimuth_index / num_azimuths, 20000.0);
        query.layers[0].exit = query.layers[0].entrance;
        queries.push_back(std::move(query));
    }

    sasktran2::Config config;
    sasktran2::solartransmission::SolarTransmissionTable2D table(geometry,
                                                                 raytracer);
    table.initialize_config(config);
    table.initialize_geometry(queries);
    sasktran2::solartransmission::SolarTableInterpolation interpolation;
    std::vector<bool> ground_hit;
    REQUIRE_NOTHROW(
        table.generate_interpolation(queries, interpolation, ground_hit));
    REQUIRE(interpolation.rows() == 2 * num_azimuths);
}

TEST_CASE("Geometry2D solar table interpolates an azimuth-local shell-crossing "
          "gap",
          "[sourceintegrator][singlescatter][geometry2d][solartable]") {
    // Reduced from the high-SZA final time group of an OMPS orbital-plane
    // calculation. One off-plane characteristic misses this far-side shell
    // crossing even though the neighboring characteristics in that azimuth
    // plane remain valid interpolation nodes.
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(81, 500.0, 80500.0);
    Eigen::VectorXd horizontal(18);
    horizontal << -0.15707249909324214, -0.13962148634892896,
        -0.12216995304054604, -0.10471839510657402, -0.08726555521991633,
        -0.0698127714529716, -0.05235964029640497, -0.034906379204992954,
        -0.017453112724842415, 0.0, 0.01745299241956929, 0.03490598281398394,
        0.05235897258070354, 0.06981196208857511, 0.08726495170676239,
        0.10471794180429507, 0.12217093274961874, 0.1396239249101459;
    sasktran2::Geometry2D geometry(0.06722708902977514, 2.806877591558785,
                                   6357595.44403317, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    constexpr int num_azimuths = 9;
    sasktran2::raytracing::TracedRay query;
    query.layers.resize(1);
    query.layers[0].entrance.position =
        geometry.coordinates().solar_coordinate_vector(
            -0.0177337365238761, 2.0 * EIGEN_PI / num_azimuths, 1500.0);
    query.layers[0].exit = query.layers[0].entrance;

    sasktran2::Config config;
    sasktran2::solartransmission::SolarTransmissionTable2D table(geometry,
                                                                 raytracer);
    table.initialize_config(config);
    table.initialize_geometry({query});
    sasktran2::solartransmission::SolarTableInterpolation interpolation;
    std::vector<bool> ground_hit;
    REQUIRE_NOTHROW(
        table.generate_interpolation({query}, interpolation, ground_hit));
    REQUIRE(interpolation.rows() == 2);
    REQUIRE(std::none_of(ground_hit.begin(), ground_hit.end(),
                         [](bool value) { return value; }));

    Eigen::VectorXd extinction(geometry.size());
    for (int horizontal_index = 0;
         horizontal_index < geometry.num_horizontal_locations();
         ++horizontal_index) {
        for (int altitude_index = 0; altitude_index < geometry.num_altitudes();
             ++altitude_index) {
            const double altitude =
                geometry.altitude_grid().grid()[altitude_index];
            const double angle =
                geometry.horizontal_angle_grid()[horizontal_index];
            extinction[geometry.location_index(altitude_index,
                                               horizontal_index)] =
                1.0e-5 * std::exp(-(altitude - 500.0) / 7000.0) *
                (1.0 + 0.1 * std::sin(4.0 * angle));
        }
    }
    Eigen::VectorXd table_nodes(table.table_size());
    Eigen::VectorXd table_od(interpolation.rows());
    table.apply(extinction, table_nodes);
    interpolation.apply(table_nodes, table_od);

    sasktran2::solartransmission::SolarTransmissionExact exact(geometry,
                                                               raytracer);
    sasktran2::solartransmission::SolarGeometryMatrix exact_matrix;
    std::vector<bool> exact_ground_hit;
    exact.generate_geometry_matrix({query}, exact_matrix, exact_ground_hit);
    Eigen::VectorXd exact_od(exact_matrix.rows());
    exact_matrix.multiply(extinction, exact_od);
    REQUIRE(ground_hit == exact_ground_hit);
    const double relative_error = (table_od - exact_od).cwiseAbs().maxCoeff() /
                                  exact_od.cwiseAbs().maxCoeff();
    CAPTURE(table_od.transpose(), exact_od.transpose(), relative_error);
    REQUIRE(relative_error < 1.0e-6);

    const Eigen::Vector2d endpoint_cotangent(0.7, -0.3);
    Eigen::VectorXd table_cotangent(table.table_size());
    interpolation.apply_transpose(endpoint_cotangent, table_cotangent);
    Eigen::VectorXd table_gradient = Eigen::VectorXd::Zero(geometry.size());
    table.accumulate_transpose(table_cotangent, table_gradient, 1.0);
    Eigen::VectorXd exact_gradient = Eigen::VectorXd::Zero(geometry.size());
    for (Eigen::Index row = 0; row < exact_matrix.rows(); ++row) {
        for (sasktran2::solartransmission::SolarGeometryMatrix::InnerIterator
                 value(exact_matrix, row);
             value; ++value) {
            exact_gradient[value.index()] +=
                endpoint_cotangent[row] * value.value();
        }
    }
    const double relative_gradient_error =
        (table_gradient - exact_gradient).norm() / exact_gradient.norm();
    CAPTURE(relative_gradient_error);
    REQUIRE(relative_gradient_error < 1.0e-4);
    REQUIRE(endpoint_cotangent.dot(table_od) ==
            Catch::Approx(table_gradient.dot(extinction)).epsilon(1.0e-12));
}

TEST_CASE("Geometry2D solar table completes an azimuth-local surface tangent",
          "[sourceintegrator][singlescatter][geometry2d][solartable]") {
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(81, 500.0, 80500.0);
    Eigen::VectorXd horizontal(26);
    horizontal << -0.22687388786483126, -0.20942296752292075,
        -0.19197206009364517, -0.17452116267172288, -0.15707027233969004,
        -0.13961938617145178, -0.12216850123584493, -0.104716918569682,
        -0.08726460712167375, -0.06981241610955545, -0.05235947810698597,
        -0.034906301152543374, -0.017453217776235524, 0.0, 0.017453236429728323,
        0.034906401529258245, 0.052359269565022806, 0.06981192937150443,
        0.08726465540793593, 0.10471607914063315, 0.1221674704049722,
        0.1396181704573594, 0.15706761621702, 0.17451769934386843,
        0.19196668818341656, 0.2094140853642529;
    sasktran2::Geometry2D geometry(0.16217357907765928, 0.31908455164488897,
                                   6360006.564655005, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    constexpr int num_azimuths = 13;
    sasktran2::raytracing::TracedRay query;
    query.layers.resize(1);
    query.layers[0].entrance.position =
        geometry.coordinates().solar_coordinate_vector(
            1.0e-12, 12.0 * 2.0 * EIGEN_PI / num_azimuths, 500.0);
    query.layers[0].exit = query.layers[0].entrance;

    sasktran2::Config config;
    sasktran2::solartransmission::SolarTransmissionTable2D table(geometry,
                                                                 raytracer);
    table.initialize_config(config);
    table.initialize_geometry({query});
    sasktran2::solartransmission::SolarTableInterpolation interpolation;
    std::vector<bool> ground_hit;
    REQUIRE_NOTHROW(
        table.generate_interpolation({query}, interpolation, ground_hit));
    REQUIRE(std::none_of(ground_hit.begin(), ground_hit.end(),
                         [](bool value) { return value; }));

    const Eigen::VectorXd extinction =
        Eigen::VectorXd::Constant(geometry.size(), 1.0e-5);
    Eigen::VectorXd table_nodes(table.table_size());
    Eigen::VectorXd table_od(interpolation.rows());
    table.apply(extinction, table_nodes);
    interpolation.apply(table_nodes, table_od);
    sasktran2::solartransmission::SolarTransmissionExact exact(geometry,
                                                               raytracer);
    sasktran2::solartransmission::SolarGeometryMatrix exact_matrix;
    std::vector<bool> exact_ground_hit;
    exact.generate_geometry_matrix({query}, exact_matrix, exact_ground_hit);
    Eigen::VectorXd exact_od(exact_matrix.rows());
    exact_matrix.multiply(extinction, exact_od);
    REQUIRE(ground_hit == exact_ground_hit);
    REQUIRE(table_od.isApprox(exact_od, 1.0e-6));
}

TEST_CASE("Geometry2D characteristic solar table follows refracted endpoint "
          "OD",
          "[sourceintegrator][singlescatter][geometry2d][solartable]") {
    auto geometry = geometry2d();
    geometry.refractive_index() = vector({1.0003, 1.0001, 1.0});
    geometry.validate();
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    const double query_radius = earth_radius + 20.0;
    const double impact_parameter =
        (earth_radius + 10.0) * geometry.refractive_index()[1];
    const auto& sun = geometry.coordinates().sun_unit();
    Eigen::Vector3d impact_direction =
        geometry.coordinates().reference_z() -
        geometry.coordinates().reference_z().dot(sun) * sun;
    impact_direction.normalize();
    const double observer_radius = query_radius + 1.0;
    sasktran2::viewinggeometry::ViewingRay incoming;
    incoming.observer.position =
        impact_parameter * impact_direction +
        std::sqrt(observer_radius * observer_radius -
                  impact_parameter * impact_parameter) *
            sun;
    incoming.look_away = -sun;
    sasktran2::raytracing::TracedRay characteristic;
    raytracer.trace_ray(incoming, geometry.refractive_index(), characteristic);
    REQUIRE(!characteristic.ground_is_hit);

    sasktran2::Location far_top;
    double far_top_cos_sza = 1.0;
    for (const auto& layer : characteristic.layers) {
        const double radius_error =
            std::abs(layer.exit.radius() - query_radius);
        const double cos_sza = layer.exit.position.normalized().dot(sun);
        if (radius_error < 1.0e-8 && cos_sza < far_top_cos_sza) {
            far_top = layer.exit;
            far_top_cos_sza = cos_sza;
        }
    }
    REQUIRE(far_top_cos_sza < 0.0);

    sasktran2::raytracing::TracedRay line_of_sight;
    line_of_sight.layers.resize(1);
    line_of_sight.layers[0].entrance = far_top;
    line_of_sight.layers[0].exit = line_of_sight.layers[0].entrance;

    sasktran2::Config config;
    config.set_solar_refraction(true);
    sasktran2::solartransmission::SolarTransmissionTable2D table(geometry,
                                                                 raytracer);
    table.initialize_config(config);
    table.initialize_geometry({line_of_sight});
    sasktran2::solartransmission::SolarTableInterpolation interpolation;
    std::vector<bool> table_ground_hit;
    std::vector<Eigen::Vector3d> solar_propagation_directions;
    table.generate_interpolation({line_of_sight}, interpolation,
                                 table_ground_hit,
                                 &solar_propagation_directions);

    Eigen::VectorXd extinction(geometry.size());
    for (int index = 0; index < geometry.size(); ++index) {
        extinction[index] = 0.01 * (1.0 + 0.07 * index);
    }
    Eigen::VectorXd table_nodes(table.table_size());
    table.apply(extinction, table_nodes);
    Eigen::VectorXd table_od(interpolation.rows());
    interpolation.apply(table_nodes, table_od);

    const double exact_od =
        direct_path_weights(characteristic, geometry).dot(extinction);

    CAPTURE(far_top_cos_sza, table_ground_hit);
    REQUIRE(std::none_of(table_ground_hit.begin(), table_ground_hit.end(),
                         [](bool value) { return value; }));
    for (Eigen::Index row = 0; row < table_od.size(); ++row) {
        REQUIRE(table_od[row] == Catch::Approx(exact_od).epsilon(1e-11));
    }

    REQUIRE(solar_propagation_directions.size() == 2);
    for (const auto& direction : solar_propagation_directions) {
        REQUIRE(direction.norm() == Catch::Approx(1.0).epsilon(1e-13));
        REQUIRE(direction.dot(-sun) < 0.999999999);
    }

    sasktran2::solartransmission::SolarTransmissionExact exact(geometry,
                                                               raytracer);
    sasktran2::solartransmission::SolarGeometryMatrix exact_matrix;
    auto exact_ground_hit = table_ground_hit;
    exact.generate_refracted_geometry_matrix({line_of_sight},
                                             solar_propagation_directions,
                                             exact_matrix, exact_ground_hit);
    Eigen::VectorXd retraced_od(exact_matrix.rows());
    exact_matrix.multiply(extinction, retraced_od);
    REQUIRE(exact_ground_hit == table_ground_hit);
    for (Eigen::Index row = 0; row < retraced_od.size(); ++row) {
        REQUIRE(retraced_od[row] == Catch::Approx(exact_od).epsilon(2e-10));
    }
}

TEST_CASE("Exact Geometry2D solar shadow is stable at a grazing surface ray",
          "[sourceintegrator][singlescatter][geometry2d]") {
    constexpr double target_radius = 20.0;
    const double grazing_cos_sza = -std::sqrt(
        1.0 - earth_radius * earth_radius / (target_radius * target_radius));
    const sasktran2::Geometry2D geometry(grazing_cos_sza, 0.0, earth_radius,
                                         vector({0.0, 10.0, 20.0}),
                                         vector({-0.5, 0.5}));
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    sasktran2::raytracing::TracedRay line_of_sight;
    line_of_sight.layers.resize(1);
    line_of_sight.layers[0].entrance.position =
        Eigen::Vector3d(0.0, 0.0, target_radius);
    line_of_sight.layers[0].exit = line_of_sight.layers[0].entrance;

    sasktran2::solartransmission::SolarTransmissionExact transmission(
        geometry, raytracer);
    sasktran2::solartransmission::SolarGeometryMatrix matrix;
    std::vector<bool> ground_hit;
    transmission.generate_geometry_matrix({line_of_sight}, matrix, ground_hit);

    REQUIRE(ground_hit.size() == 2);
    REQUIRE(ground_hit[0]);
    REQUIRE(ground_hit[1]);
    REQUIRE(matrix.non_zeros() == 0);
}

TEST_CASE("Geometry2D exact single scatter is invariant across thread models",
          "[engine][singlescatter][geometry2d][threading]") {
    const auto geometry = geometry2d();
    sasktran2::viewinggeometry::ViewingGeometryContainer viewing;
    viewing.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(
            5.0, 0.2, 30.0, 0.6));
    viewing.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            0.6, -0.1, 0.7, 30.0));

    constexpr int num_wavelengths = 12;
    sasktran2::Config single_thread_config;
    single_thread_config.set_num_threads(1);
    sasktran2::atmosphere::Atmosphere<1> atmosphere(num_wavelengths, geometry,
                                                    single_thread_config, true);
    atmosphere.storage().resize_derivatives(1);
    atmosphere.storage().leg_coeff.setZero();
    for (int location = 0; location < geometry.size(); ++location) {
        for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
            atmosphere.storage().total_extinction(location, wavelength) =
                0.002 * (1.0 + 0.05 * location + 0.03 * wavelength);
            atmosphere.storage().ssa(location, wavelength) =
                0.45 + 0.01 * location;
            atmosphere.storage().solar_irradiance(wavelength) =
                1.0 + 0.1 * wavelength;
            atmosphere.storage().leg_coeff(0, location, wavelength) = 1.0;
            atmosphere.storage().leg_coeff(2, location, wavelength) = 0.25;
            atmosphere.storage().d_leg_coeff(0, location, wavelength, 0) = 0.1;
            atmosphere.storage().d_leg_coeff(2, location, wavelength, 0) =
                0.03 * (1.0 + 0.02 * location);
        }
    }

    Sasktran2<1> single_thread_engine(single_thread_config, &geometry, viewing);
    sasktran2::atmosphere::Atmosphere<1> derivative_free_atmosphere(
        num_wavelengths, geometry, single_thread_config, false);
    derivative_free_atmosphere.storage().total_extinction =
        atmosphere.storage().total_extinction;
    derivative_free_atmosphere.storage().ssa = atmosphere.storage().ssa;
    derivative_free_atmosphere.storage().solar_irradiance =
        atmosphere.storage().solar_irradiance;
    derivative_free_atmosphere.storage().leg_coeff =
        atmosphere.storage().leg_coeff;
    sasktran2::OutputIdealDense<1> derivative_free_output;
    single_thread_engine.calculate_radiance(derivative_free_atmosphere,
                                            derivative_free_output);

    sasktran2::OutputIdealDense<1> single_thread_output;
    single_thread_engine.calculate_radiance(atmosphere, single_thread_output);
    REQUIRE(single_thread_output.radiance().deriv.cols() ==
            atmosphere.num_deriv());

    for (const auto threading_model :
         {sasktran2::Config::ThreadingModel::wavelength,
          sasktran2::Config::ThreadingModel::source}) {
        sasktran2::Config threaded_config;
        threaded_config.set_num_threads(4);
        threaded_config.set_threading_model(threading_model);
        Sasktran2<1> threaded_engine(threaded_config, &geometry, viewing);
        sasktran2::OutputIdealDense<1> threaded_output;
        threaded_engine.calculate_radiance(atmosphere, threaded_output);

        REQUIRE(threaded_output.radiance().value.isApprox(
            single_thread_output.radiance().value, 1e-13));
        REQUIRE(threaded_output.radiance().deriv.isApprox(
            single_thread_output.radiance().deriv, 1e-13));
    }

    // Exercise the active-column 2D source-derivative path against finite
    // differences. Choose the strongest native extinction and SSA entries so
    // this remains independent of the exact ray-cell traversal.
    constexpr int wavelength = 5;
    constexpr int los = 0;
    const int output_index = wavelength * viewing.observer_rays().size() + los;
    Eigen::Index extinction_location;
    const double max_extinction_derivative =
        single_thread_output.radiance()
            .deriv.row(output_index)
            .head(geometry.size())
            .cwiseAbs()
            .maxCoeff(&extinction_location);
    REQUIRE(max_extinction_derivative > 0.0);

    constexpr double extinction_delta = 1e-7;
    atmosphere.storage().total_extinction(extinction_location, wavelength) +=
        extinction_delta;
    sasktran2::OutputIdealDense<1> extinction_above;
    single_thread_engine.calculate_radiance(atmosphere, extinction_above);
    atmosphere.storage().total_extinction(extinction_location, wavelength) -=
        2 * extinction_delta;
    sasktran2::OutputIdealDense<1> extinction_below;
    single_thread_engine.calculate_radiance(atmosphere, extinction_below);
    atmosphere.storage().total_extinction(extinction_location, wavelength) +=
        extinction_delta;

    const double numeric_extinction_derivative =
        (extinction_above.radiance().value(output_index) -
         extinction_below.radiance().value(output_index)) /
        (2 * extinction_delta);
    REQUIRE(numeric_extinction_derivative ==
            Catch::Approx(single_thread_output.radiance().deriv(
                              output_index, extinction_location))
                .epsilon(2e-7));

    Eigen::Index ssa_location;
    const double max_ssa_derivative =
        single_thread_output.radiance()
            .deriv.row(output_index)
            .segment(geometry.size(), geometry.size())
            .cwiseAbs()
            .maxCoeff(&ssa_location);
    REQUIRE(max_ssa_derivative > 0.0);

    constexpr double ssa_delta = 1e-6;
    atmosphere.storage().ssa(ssa_location, wavelength) += ssa_delta;
    sasktran2::OutputIdealDense<1> ssa_above;
    single_thread_engine.calculate_radiance(atmosphere, ssa_above);
    atmosphere.storage().ssa(ssa_location, wavelength) -= 2 * ssa_delta;
    sasktran2::OutputIdealDense<1> ssa_below;
    single_thread_engine.calculate_radiance(atmosphere, ssa_below);
    atmosphere.storage().ssa(ssa_location, wavelength) += ssa_delta;

    const double numeric_ssa_derivative =
        (ssa_above.radiance().value(output_index) -
         ssa_below.radiance().value(output_index)) /
        (2 * ssa_delta);
    REQUIRE(numeric_ssa_derivative ==
            Catch::Approx(single_thread_output.radiance().deriv(
                              output_index, geometry.size() + ssa_location))
                .epsilon(2e-8));

    Eigen::Index scattering_location;
    const double max_scattering_derivative =
        single_thread_output.radiance()
            .deriv.row(output_index)
            .segment(2 * geometry.size(), geometry.size())
            .cwiseAbs()
            .maxCoeff(&scattering_location);
    REQUIRE(max_scattering_derivative > 0.0);

    constexpr double scattering_delta = 1e-6;
    for (int coefficient = 0;
         coefficient < atmosphere.storage().leg_coeff.dimension(0);
         ++coefficient) {
        atmosphere.storage().leg_coeff(coefficient, scattering_location,
                                       wavelength) +=
            scattering_delta *
            atmosphere.storage().d_leg_coeff(coefficient, scattering_location,
                                             wavelength, 0);
    }
    sasktran2::OutputIdealDense<1> scattering_above;
    single_thread_engine.calculate_radiance(atmosphere, scattering_above);
    for (int coefficient = 0;
         coefficient < atmosphere.storage().leg_coeff.dimension(0);
         ++coefficient) {
        atmosphere.storage().leg_coeff(coefficient, scattering_location,
                                       wavelength) -=
            2 * scattering_delta *
            atmosphere.storage().d_leg_coeff(coefficient, scattering_location,
                                             wavelength, 0);
    }
    sasktran2::OutputIdealDense<1> scattering_below;
    single_thread_engine.calculate_radiance(atmosphere, scattering_below);
    for (int coefficient = 0;
         coefficient < atmosphere.storage().leg_coeff.dimension(0);
         ++coefficient) {
        atmosphere.storage().leg_coeff(coefficient, scattering_location,
                                       wavelength) +=
            scattering_delta *
            atmosphere.storage().d_leg_coeff(coefficient, scattering_location,
                                             wavelength, 0);
    }

    const double numeric_scattering_derivative =
        (scattering_above.radiance().value(output_index) -
         scattering_below.radiance().value(output_index)) /
        (2 * scattering_delta);
    REQUIRE(numeric_scattering_derivative ==
            Catch::Approx(
                single_thread_output.radiance().deriv(
                    output_index, 2 * geometry.size() + scattering_location))
                .epsilon(2e-8));
}
#endif
