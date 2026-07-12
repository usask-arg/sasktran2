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
        const sasktran2::raytracing::InterpolationStencil1D& actual,
        const sasktran2::Location& location,
        const sasktran2::Geometry1D& geometry) {
        std::vector<std::pair<int, double>> expected;
        geometry.assign_interpolation_weights(location, expected);
        expected.erase(std::remove_if(expected.begin(), expected.end(),
                                      [](const auto& weight) {
                                          return weight.second == 0.0;
                                      }),
                       expected.end());

        REQUIRE(actual.size() == expected.size());
        for (std::size_t index = 0; index < expected.size(); ++index) {
            const auto actual_weight = actual[index];
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
            CAPTURE(layer_index);
            REQUIRE(layer.integrated_od.valid());
            REQUIRE(layer.integrated_od.cell_index >= 0);
            REQUIRE(layer.integrated_od.cell_index + 1 < geometry.size());
            require_stencil_matches(layer.entrance_interpolation_weights,
                                    layer.entrance, geometry);
            require_stencil_matches(layer.exit_interpolation_weights,
                                    layer.exit, geometry);

            Eigen::VectorXd compact = Eigen::VectorXd::Zero(geometry.size());
            compact[layer.integrated_od.cell_index] =
                layer.integrated_od.weights[0];
            compact[layer.integrated_od.cell_index + 1] =
                layer.integrated_od.weights[1];
            REQUIRE(
                compact.isApprox(legacy.row(layer_index).transpose(), 2e-13));

            const double effective_distance =
                layer.layer_distance * layer.curvature_factor;
            REQUIRE(
                layer.integrated_od.weights[0] +
                    layer.integrated_od.weights[1] ==
                Catch::Approx(effective_distance).epsilon(2e-12).margin(1e-9));

            int nonzero_count = 0;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(
                     sparse, layer_index);
                 it; ++it) {
                ++nonzero_count;
                REQUIRE((it.col() == layer.integrated_od.cell_index ||
                         it.col() == layer.integrated_od.cell_index + 1));
            }
            REQUIRE(nonzero_count <= 2);

            Eigen::VectorXd extinction(geometry.size());
            for (int i = 0; i < geometry.size(); ++i) {
                extinction[i] = 1e-6 * (1.0 + 0.17 * i + 0.03 * i * i);
            }
            const double baseline = legacy.row(layer_index).dot(extinction);
            constexpr double perturbation = 1e-7;
            for (int local_index = 0; local_index < 2; ++local_index) {
                auto perturbed = extinction;
                const int grid_index =
                    layer.integrated_od.cell_index + local_index;
                perturbed[grid_index] += perturbation;
                const double derivative =
                    (legacy.row(layer_index).dot(perturbed) - baseline) /
                    perturbation;
                REQUIRE(derivative ==
                        Catch::Approx(layer.integrated_od.weights[local_index])
                            .epsilon(2e-11)
                            .margin(1e-8));
            }
        }

        Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(1, geometry.size());
        std::vector<std::pair<int, double>> workspace;
        sasktran2::solartransmission::assign_dense_matrix_column(
            0, ray, geometry, dense, workspace);
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

TEST_CASE("Compact 1D interpolation stencil storage and topology",
          "[sasktran2][raytracing][od_weights]") {
    STATIC_REQUIRE(sizeof(sasktran2::raytracing::InterpolationStencil1D) == 16);
    STATIC_REQUIRE(sizeof(sasktran2::Location) <= 40);
    STATIC_REQUIRE(std::is_trivially_copyable_v<
                   sasktran2::raytracing::InterpolationStencil1D>);

    sasktran2::raytracing::InterpolationStencil1D stencil;
    REQUIRE(stencil.empty());

    stencil.assign({{2, 0.25}, {3, 0.75}}, 5);
    REQUIRE(stencil.valid());
    REQUIRE(stencil.size() == 2);
    REQUIRE(stencil[0] == std::make_pair(2, 0.25));
    REQUIRE(stencil[1] == std::make_pair(3, 0.75));

    stencil.assign({{4, 1.0}}, 5);
    REQUIRE(stencil.size() == 1);
    REQUIRE(stencil[0] == std::make_pair(4, 1.0));

    stencil.assign({{1, 0.0}, {2, 1.0}}, 5);
    REQUIRE(stencil.size() == 1);
    REQUIRE(stencil[0] == std::make_pair(2, 1.0));

    stencil.assign({{1, 0.5}, {3, 0.5}}, 5);
    REQUIRE(stencil.empty());
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

TEST_CASE("Integrated 1D OD weights - legacy fallback",
          "[sasktran2][raytracing][od_weights]") {
    auto geometry = make_geometry(sasktran2::geometrytype::spherical,
                                  sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::SphericalShellRayTracer tracer(geometry);
    const auto viewing_ray =
        sasktran2::viewinggeometry::TangentAltitude(3900.0, 0.2, 120000.0)
            .construct_ray(geometry.coordinates());

    sasktran2::raytracing::TracedRay cached;
    tracer.trace_ray(viewing_ray, cached);
    auto legacy = cached;
    for (auto& layer : legacy.layers) {
        layer.integrated_od.reset();
    }

    Eigen::SparseMatrix<double, Eigen::RowMajor> cached_matrix;
    Eigen::SparseMatrix<double, Eigen::RowMajor> fallback_matrix;
    sasktran2::raytracing::construct_od_matrix(cached, geometry, cached_matrix);
    sasktran2::raytracing::construct_od_matrix(legacy, geometry,
                                               fallback_matrix);
    REQUIRE(Eigen::MatrixXd(cached_matrix)
                .isApprox(Eigen::MatrixXd(fallback_matrix), 2e-13));

    Eigen::MatrixXd cached_dense = Eigen::MatrixXd::Zero(1, geometry.size());
    Eigen::MatrixXd fallback_dense = Eigen::MatrixXd::Zero(1, geometry.size());
    std::vector<std::pair<int, double>> workspace;
    sasktran2::solartransmission::assign_dense_matrix_column(
        0, cached, geometry, cached_dense, workspace);
    sasktran2::solartransmission::assign_dense_matrix_column(
        0, legacy, geometry, fallback_dense, workspace);
    REQUIRE(cached_dense.isApprox(fallback_dense, 2e-13));
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
