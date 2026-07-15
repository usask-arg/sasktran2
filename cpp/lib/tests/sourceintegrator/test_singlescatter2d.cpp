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
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix;
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
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
                     entry(matrix, row);
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
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix;
    std::vector<bool> ground_hit;
    transmission.generate_geometry_matrix({line_of_sight}, matrix, ground_hit);

    REQUIRE(ground_hit.size() == 2);
    REQUIRE(ground_hit[0]);
    REQUIRE(ground_hit[1]);
    REQUIRE(matrix.nonZeros() == 0);
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
