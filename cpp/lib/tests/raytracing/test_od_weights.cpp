#include <sasktran2/test_helper.h>

#include <sasktran2.h>

namespace {
    Eigen::VectorXd grid(std::initializer_list<double> values) {
        return Eigen::Map<const Eigen::VectorXd>(values.begin(), values.size());
    }

    sasktran2::Geometry1D
    make_geometry(sasktran2::geometrytype geometry_type,
                  sasktran2::grids::interpolation interpolation) {
        return sasktran2::Geometry1D(0.4, 0.2, 6372000.0,
                                     grid({0.0, 800.0, 2500.0, 7000.0, 15000.0,
                                           30000.0, 50000.0, 80000.0}),
                                     interpolation, geometry_type);
    }

    Eigen::VectorXd
    legacy_layer_weights(const sasktran2::raytracing::SphericalLayer& layer,
                         const sasktran2::Geometry& geometry) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(geometry.size());
        std::vector<std::pair<int, double>> interpolation_weights;

        geometry.assign_interpolation_weights(layer.entrance,
                                              interpolation_weights);
        for (const auto& iw : interpolation_weights) {
            result[iw.first] += iw.second * layer.od_quad_start;
        }

        geometry.assign_interpolation_weights(layer.exit,
                                              interpolation_weights);
        for (const auto& iw : interpolation_weights) {
            result[iw.first] += iw.second * layer.od_quad_end;
        }
        return result;
    }

    Eigen::MatrixXd
    legacy_od_matrix(const sasktran2::raytracing::TracedRay& ray,
                     const sasktran2::Geometry& geometry) {
        Eigen::MatrixXd result =
            Eigen::MatrixXd::Zero(ray.layers.size(), geometry.size());
        for (int i = 0; i < ray.layers.size(); ++i) {
            result.row(i) = legacy_layer_weights(ray.layers[i], geometry);
        }
        return result;
    }

    void require_stencil_matches(
        const sasktran2::raytracing::GridWeightStencilView& actual,
        const sasktran2::Location& location,
        const sasktran2::Geometry1D& geometry) {
        std::vector<std::pair<int, double>> expected;
        geometry.assign_interpolation_weights(location, expected);
        expected.erase(std::remove_if(expected.begin(), expected.end(),
                                      [](const auto& weight) {
                                          return weight.second == 0.0;
                                      }),
                       expected.end());

        std::vector<std::pair<int, double>> actual_nonzero;
        for (std::size_t index = 0; index < actual.size(); ++index) {
            if (actual[index].second != 0.0) {
                actual_nonzero.push_back(actual[index]);
            }
        }
        REQUIRE(actual_nonzero.size() == expected.size());
        for (std::size_t index = 0; index < expected.size(); ++index) {
            const auto actual_weight = actual_nonzero[index];
            REQUIRE(actual_weight.first == expected[index].first);
            REQUIRE(actual_weight.second == expected[index].second);
        }
    }

    void require_od_weights(const sasktran2::raytracing::TracedRay& ray,
                            const sasktran2::Geometry1D& geometry) {
        const auto legacy = legacy_od_matrix(ray, geometry);

        Eigen::SparseMatrix<double, Eigen::RowMajor> sparse;
        sasktran2::raytracing::construct_od_matrix(ray, geometry, sparse);
        const Eigen::MatrixXd actual(sparse);
        REQUIRE(actual.isApprox(legacy, 2e-13));

        for (int layer_index = 0; layer_index < ray.layers.size();
             ++layer_index) {
            const auto& layer = ray.layers[layer_index];
            const auto entrance_weights = ray.entrance_weights(layer_index);
            const auto exit_weights = ray.exit_weights(layer_index);
            const auto od_weights = ray.optical_depth_weights(layer_index);
            CAPTURE(layer_index);
            REQUIRE_FALSE(od_weights.empty());
            REQUIRE(od_weights.size() <= 4);
            require_stencil_matches(entrance_weights, layer.entrance, geometry);
            require_stencil_matches(exit_weights, layer.exit, geometry);

            Eigen::VectorXd compact = Eigen::VectorXd::Zero(geometry.size());
            for (std::size_t index = 0; index < od_weights.size(); ++index) {
                const auto weight = od_weights[index];
                compact[weight.first] = weight.second;
            }
            REQUIRE(
                compact.isApprox(legacy.row(layer_index).transpose(), 2e-13));

            const double effective_distance =
                layer.layer_distance * layer.curvature_factor;
            double integrated_weight_sum = 0.0;
            for (std::size_t index = 0; index < od_weights.size(); ++index) {
                integrated_weight_sum += od_weights[index].second;
            }
            REQUIRE(
                integrated_weight_sum ==
                Catch::Approx(effective_distance).epsilon(2e-12).margin(1e-9));

            int nonzero_count = 0;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(
                     sparse, layer_index);
                 it; ++it) {
                ++nonzero_count;
                bool matches_compiled_weight = false;
                for (std::size_t index = 0; index < od_weights.size();
                     ++index) {
                    matches_compiled_weight |=
                        it.col() == od_weights[index].first;
                }
                REQUIRE(matches_compiled_weight);
            }
            REQUIRE(nonzero_count <= od_weights.size());

            Eigen::VectorXd extinction(geometry.size());
            for (int i = 0; i < geometry.size(); ++i) {
                extinction[i] = 1e-6 * (1.0 + 0.17 * i + 0.03 * i * i);
            }
            const double baseline = legacy.row(layer_index).dot(extinction);
            constexpr double perturbation = 1e-7;
            for (std::size_t local_index = 0; local_index < od_weights.size();
                 ++local_index) {
                auto perturbed = extinction;
                const int grid_index = od_weights[local_index].first;
                perturbed[grid_index] += perturbation;
                const double derivative =
                    (legacy.row(layer_index).dot(perturbed) - baseline) /
                    perturbation;
                REQUIRE(derivative ==
                        Catch::Approx(od_weights[local_index].second)
                            .epsilon(2e-11)
                            .margin(1e-8));
            }
        }

        Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(1, geometry.size());
        sasktran2::solartransmission::assign_dense_matrix_column(0, ray, dense);
        REQUIRE(dense.row(0).isApprox(legacy.colwise().sum(), 2e-13));

        Eigen::MatrixXd spectra(geometry.size(), 4);
        for (int i = 0; i < spectra.rows(); ++i) {
            for (int j = 0; j < spectra.cols(); ++j) {
                spectra(i, j) =
                    1e-7 * (1.0 + i * 0.11 + j * 0.23 + i * j * 0.013);
            }
        }
        REQUIRE((sparse * spectra).isApprox(legacy * spectra, 2e-13));
    }

    void set_refractive_index(sasktran2::Geometry1D& geometry) {
        geometry.refractive_index() << 1.00030, 1.00024, 1.00018, 1.00013,
            1.00009, 1.000055, 1.000025, 1.0;
    }
} // namespace

TEST_CASE("Dimension-independent grid-weight stencil view",
          "[sasktran2][raytracing][od_weights]") {
    STATIC_REQUIRE(sizeof(sasktran2::Location) <= 40);
    STATIC_REQUIRE(sizeof(sasktran2::raytracing::TracedLayer) <=
                   sizeof(sasktran2::raytracing::LayerGeometry) + 16);
    STATIC_REQUIRE(std::is_trivially_copyable_v<
                   sasktran2::raytracing::GridWeightStencilView>);

    const std::array<int, 2> indices = {2, 3};
    const std::array<double, 2> weights = {0.25, 0.75};
    const sasktran2::raytracing::GridWeightStencilView stencil(
        indices.data(), weights.data(), weights.size());
    REQUIRE(stencil.size() == 2);
    REQUIRE(stencil[0] == std::make_pair(2, 0.25));
    REQUIRE(stencil[1] == std::make_pair(3, 0.75));

    sasktran2::raytracing::TracedRay ray;
    ray.layers.resize(1);
    const std::array<int, 8> cube_indices = {0, 1, 2, 3, 4, 5, 6, 7};
    const std::array<double, 8> entrance = {0.125, 0.125, 0.125, 0.125,
                                            0.125, 0.125, 0.125, 0.125};
    const std::array<double, 8> exit = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const std::array<double, 8> integrated = {1.0, 2.0, 3.0, 4.0,
                                              5.0, 6.0, 7.0, 8.0};
    ray.set_layer_weights(0, cube_indices, entrance, exit, integrated);
    REQUIRE(ray.num_grid_weights() == 8);
    REQUIRE(ray.entrance_weights(0).size() == 8);
    REQUIRE(ray.exit_weights(0)[0] == std::make_pair(0, 1.0));
    REQUIRE(ray.optical_depth_weights(0)[7] == std::make_pair(7, 8.0));

    REQUIRE_THROWS_AS(ray.set_layer_weights(0, cube_indices, entrance,
                                            std::array<double, 1>{1.0},
                                            integrated),
                      std::invalid_argument);

    ray.reset();
    REQUIRE(ray.layers.empty());
    REQUIRE(ray.num_grid_weights() == 0);
}

TEST_CASE("Compiled 1D weights preserve rare three-node boundary stencils",
          "[sasktran2][raytracing][od_weights]") {
    constexpr double earth_radius = 6372000.0;
    auto geometry = sasktran2::Geometry1D(
        0.4, 0.2, earth_radius, grid({0.0, 10000.0, 30000.0}),
        sasktran2::grids::interpolation::linear,
        sasktran2::geometrytype::spherical);

    sasktran2::raytracing::TracedRay ray;
    ray.layers.resize(1);
    auto& layer = ray.layers.front();
    layer.entrance.position =
        (earth_radius + 9998.0) * Eigen::Vector3d::UnitZ();
    layer.entrance.on_exact_altitude = true;
    layer.entrance.lower_alt_index = 1;
    layer.exit.position = (earth_radius + 30000.0) * Eigen::Vector3d::UnitZ();
    layer.exit.on_exact_altitude = true;
    layer.exit.lower_alt_index = 2;
    layer.od_quad_start = 2.0;
    layer.od_quad_end = 3.0;

    std::vector<std::pair<int, double>> workspace;
    sasktran2::raytracing::add_interpolation_weights(ray, 0, geometry,
                                                     workspace);

    const auto entrance = ray.entrance_weights(0);
    const auto exit = ray.exit_weights(0);
    const auto integrated = ray.optical_depth_weights(0);
    REQUIRE(entrance.size() == 3);
    REQUIRE(exit.size() == 3);
    REQUIRE(integrated.size() == 3);
    for (std::size_t index = 0; index < integrated.size(); ++index) {
        REQUIRE(integrated[index].first == static_cast<int>(index));
    }
    REQUIRE(entrance[0].second + entrance[1].second + entrance[2].second ==
            Catch::Approx(1.0));
    REQUIRE(exit[0].second + exit[1].second + exit[2].second ==
            Catch::Approx(1.0));
    REQUIRE(integrated[0].second + integrated[1].second +
                integrated[2].second ==
            Catch::Approx(5.0));
}

TEST_CASE("Integrated 1D OD weights - spherical paths",
          "[sasktran2][raytracing][od_weights]") {
    const auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                        sasktran2::grids::interpolation::shell,
                                        sasktran2::grids::interpolation::lower);
    const bool include_refraction = GENERATE(false, true);
    CAPTURE(interpolation);
    CAPTURE(include_refraction);

    auto geometry =
        make_geometry(sasktran2::geometrytype::spherical, interpolation);
    set_refractive_index(geometry);
    sasktran2::raytracing::SphericalShellRayTracer tracer(geometry);

    std::vector<sasktran2::viewinggeometry::ViewingRay> rays;
    rays.push_back(
        sasktran2::viewinggeometry::TangentAltitude(2500.0, 0.2, 120000.0)
            .construct_ray(geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::TangentAltitude(3900.0, 0.2, 120000.0)
            .construct_ray(geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::TangentAltitude(7000.0, 0.2, 10000.0)
            .construct_ray(geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::GroundViewingSolar(0.3, 0.2, 0.92, 120000.0)
            .construct_ray(geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::GroundViewingSolar(0.3, 0.2, 0.98, 23300.0)
            .construct_ray(geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::ViewingUpSolar(0.3, 0.2, 0.75, 12300.0)
            .construct_ray(geometry.coordinates()));

    for (int ray_index = 0; ray_index < rays.size(); ++ray_index) {
        CAPTURE(ray_index);
        sasktran2::raytracing::TracedRay result;
        tracer.trace_ray(rays[ray_index], result, include_refraction);
        REQUIRE_FALSE(result.layers.empty());
        require_od_weights(result, geometry);
    }
}

TEST_CASE("Integrated 1D OD weights - plane parallel paths",
          "[sasktran2][raytracing][od_weights]") {
    const auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                        sasktran2::grids::interpolation::shell,
                                        sasktran2::grids::interpolation::lower);
    CAPTURE(interpolation);

    auto geometry =
        make_geometry(sasktran2::geometrytype::planeparallel, interpolation);
    sasktran2::raytracing::PlaneParallelRayTracer tracer(geometry);

    std::vector<sasktran2::viewinggeometry::ViewingRay> rays;
    rays.push_back(
        sasktran2::viewinggeometry::GroundViewingSolar(0.3, 0.2, 0.45, 120000.0)
            .construct_ray(geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::GroundViewingSolar(0.3, 0.2, 0.83, 23300.0)
            .construct_ray(geometry.coordinates()));
    rays.push_back(
        sasktran2::viewinggeometry::ViewingUpSolar(0.3, 0.2, 0.37, 12300.0)
            .construct_ray(geometry.coordinates()));

    for (int ray_index = 0; ray_index < rays.size(); ++ray_index) {
        CAPTURE(ray_index);
        sasktran2::raytracing::TracedRay result;
        tracer.trace_ray(rays[ray_index], result);
        REQUIRE_FALSE(result.layers.empty());
        if (ray_index == 0) {
            REQUIRE(result.layers.size() == geometry.size() - 1);
        }
        require_od_weights(result, geometry);
    }
}

TEST_CASE("Integrated 1D OD weights - traced result reuse",
          "[sasktran2][raytracing][od_weights]") {
    auto geometry = make_geometry(sasktran2::geometrytype::spherical,
                                  sasktran2::grids::interpolation::linear);
    set_refractive_index(geometry);
    sasktran2::raytracing::SphericalShellRayTracer tracer(geometry);
    sasktran2::raytracing::TracedRay result;

    auto limb =
        sasktran2::viewinggeometry::TangentAltitude(3900.0, 0.2, 120000.0)
            .construct_ray(geometry.coordinates());
    tracer.trace_ray(limb, result, true);
    const auto limb_layer_count = result.layers.size();
    REQUIRE(limb_layer_count > 1);
    require_od_weights(result, geometry);

    auto empty =
        sasktran2::viewinggeometry::ViewingUpSolar(0.3, 0.2, 0.8, 120000.0)
            .construct_ray(geometry.coordinates());
    tracer.trace_ray(empty, result, false);
    REQUIRE(result.layers.empty());
    require_od_weights(result, geometry);

    auto ground =
        sasktran2::viewinggeometry::GroundViewingSolar(0.3, 0.2, 0.97, 23300.0)
            .construct_ray(geometry.coordinates());
    tracer.trace_ray(ground, result, false);
    REQUIRE_FALSE(result.layers.empty());
    REQUIRE(result.layers.size() != limb_layer_count);
    require_od_weights(result, geometry);
}

TEST_CASE("Integrated 1D OD weights remain authoritative after copying",
          "[sasktran2][raytracing][od_weights]") {
    auto geometry = make_geometry(sasktran2::geometrytype::spherical,
                                  sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::SphericalShellRayTracer tracer(geometry);
    const auto viewing_ray =
        sasktran2::viewinggeometry::TangentAltitude(3900.0, 0.2, 120000.0)
            .construct_ray(geometry.coordinates());

    sasktran2::raytracing::TracedRay original;
    tracer.trace_ray(viewing_ray, original);
    auto copied = original;

    Eigen::SparseMatrix<double, Eigen::RowMajor> original_matrix;
    Eigen::SparseMatrix<double, Eigen::RowMajor> copied_matrix;
    sasktran2::raytracing::construct_od_matrix(original, geometry,
                                               original_matrix);
    sasktran2::raytracing::construct_od_matrix(copied, geometry, copied_matrix);
    REQUIRE(Eigen::MatrixXd(original_matrix)
                .isApprox(Eigen::MatrixXd(copied_matrix), 2e-13));

    Eigen::MatrixXd original_dense = Eigen::MatrixXd::Zero(1, geometry.size());
    Eigen::MatrixXd copied_dense = Eigen::MatrixXd::Zero(1, geometry.size());
    sasktran2::solartransmission::assign_dense_matrix_column(0, original,
                                                             original_dense);
    sasktran2::solartransmission::assign_dense_matrix_column(0, copied,
                                                             copied_dense);
    REQUIRE(original_dense.isApprox(copied_dense, 2e-13));
}

#ifdef SKTRAN_RUST_SUPPORT
TEST_CASE("Integrated 1D OD weights - Rust adapter",
          "[sasktran2][raytracing][od_weights][rust]") {
    const auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                        sasktran2::grids::interpolation::shell,
                                        sasktran2::grids::interpolation::lower);
    const bool include_refraction = GENERATE(false, true);
    CAPTURE(interpolation);
    CAPTURE(include_refraction);

    auto geometry =
        make_geometry(sasktran2::geometrytype::spherical, interpolation);
    set_refractive_index(geometry);
    sasktran2::raytracing::RustRayTracer tracer(geometry);
    sasktran2::raytracing::TracedRay result;

    for (double tangent_altitude : {2500.0, 3900.0}) {
        CAPTURE(tangent_altitude);
        const auto ray = sasktran2::viewinggeometry::TangentAltitude(
                             tangent_altitude, 0.2, 120000.0)
                             .construct_ray(geometry.coordinates());
        tracer.trace_ray(ray, result, include_refraction);
        REQUIRE_FALSE(result.layers.empty());
        require_od_weights(result, geometry);
    }

    const auto ground =
        sasktran2::viewinggeometry::GroundViewingSolar(0.3, 0.2, 0.97, 23300.0)
            .construct_ray(geometry.coordinates());
    tracer.trace_ray(ground, result, include_refraction);
    REQUIRE_FALSE(result.layers.empty());
    require_od_weights(result, geometry);
}
#endif
