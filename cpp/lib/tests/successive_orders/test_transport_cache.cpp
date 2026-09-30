#include <sasktran2.h>
#include <sasktran2/test_helper.h>

#include <cmath>
#include <cstring>
#include <memory>

namespace {
    constexpr int Levels = 5;
    constexpr int Wavelengths = 3;
    constexpr int Outputs = Wavelengths * 2;

    sasktran2::Config make_config(int cache_count, int iterations,
                                  int threads = 1) {
        sasktran2::Config config;
        config.set_num_threads(threads);
        config.set_threading_model(
            sasktran2::Config::ThreadingModel::wavelength);
        config.set_single_scatter_source(
            sasktran2::Config::SingleScatterSource::exact);
        config.set_multiple_scatter_source(
            sasktran2::Config::MultipleScatterSource::successive_orders);
        config.set_num_hr_incoming(6);
        config.set_num_hr_outgoing(6);
        config.set_num_hr_spherical_iterations(iterations);
        config.set_successive_orders_relative_tolerance(iterations > 2 ? 1.0e-12
                                                                       : 0.0);
        config.set_successive_orders_absolute_tolerance(iterations > 2 ? 1.0e-14
                                                                       : 0.0);
        config.set_num_do_streams(4);
        config.set_num_singlescatter_moments(4);
        config.set_apply_delta_scaling(false);
        config.set_successive_orders_transport_cache_wavelengths(cache_count);
        return config;
    }

    sasktran2::Geometry1D make_geometry() {
        sasktran2::Coordinates coordinates(0.55, 0.15, 6372000.0,
                                           sasktran2::geometrytype::spherical);
        Eigen::VectorXd altitudes =
            Eigen::VectorXd::LinSpaced(Levels, 0.0, 40000.0);
        sasktran2::grids::AltitudeGrid grid(
            std::move(altitudes), sasktran2::grids::gridspacing::constant,
            sasktran2::grids::outofbounds::extend,
            sasktran2::grids::interpolation::linear);
        return {std::move(coordinates), std::move(grid)};
    }

    sasktran2::viewinggeometry::ViewingGeometryContainer make_viewing() {
        sasktran2::viewinggeometry::ViewingGeometryContainer viewing;
        viewing.observer_rays().emplace_back(
            std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
                0.55, 0.35, 0.72, 100000.0));
        viewing.observer_rays().emplace_back(
            std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(
                10000.0, -0.4, 100000.0, 0.55));
        return viewing;
    }

    void initialize(sasktran2::atmosphere::Atmosphere<1>& atmosphere) {
        atmosphere.storage().resize_derivatives(0);
        atmosphere.storage().leg_coeff.setZero();
        for (int wavelength = 0; wavelength < Wavelengths; ++wavelength) {
            for (int altitude = 0; altitude < Levels; ++altitude) {
                atmosphere.storage().total_extinction(altitude, wavelength) =
                    (1.4e-5 * std::exp(-altitude / 1.7) + 2.0e-9) *
                    (0.8 + 0.3 * wavelength);
                atmosphere.storage().ssa(altitude, wavelength) =
                    0.90 - 0.03 * wavelength - 0.02 * altitude;
                atmosphere.storage().leg_coeff(0, altitude, wavelength) = 1.0;
                atmosphere.storage().leg_coeff(2, altitude, wavelength) =
                    0.3 + 0.02 * altitude + 0.04 * wavelength;
            }
            atmosphere.surface().brdf_args()(0, wavelength) =
                0.06 + 0.07 * wavelength;
        }
        auto& mapping =
            atmosphere.storage().get_derivative_mapping("retrieval_state");
        mapping.allocate_extinction_derivatives();
        mapping.allocate_ssa_derivatives();
        mapping.allocate_legendre_derivatives();
        mapping.native_mapping().d_extinction->setConstant(1.0e-5);
        mapping.native_mapping().d_ssa->setConstant(1.0e-2);
        mapping.native_mapping().d_legendre->setZero();
        mapping.native_mapping().d_legendre->chip(2, 0).setConstant(0.03);
        mapping.native_mapping().scat_factor->setOnes();
        atmosphere.storage().finalize_scattering_derivatives(0);
        REQUIRE(atmosphere.num_scattering_deriv_groups() == 1);
        atmosphere.mark_changed();
    }

    struct Products {
        Eigen::VectorXd forward;
        Eigen::VectorXd radiance;
        Eigen::VectorXd gradient;
        Eigen::VectorXd jvp;
    };

    Products calculate(Sasktran2<1>& engine,
                       sasktran2::atmosphere::Atmosphere<1>& atmosphere) {
        Products result;
        sasktran2::OutputIdealDense<1> forward;
        engine.calculate_radiance(atmosphere, forward);
        result.forward = forward.radiance().value;
        result.radiance = Eigen::VectorXd::Zero(Outputs);
        result.gradient = Eigen::VectorXd::Zero(Levels);
        const Eigen::VectorXd cotangent =
            Eigen::VectorXd::LinSpaced(Outputs, 0.4, 1.1);
        Eigen::Map<Eigen::VectorXd> radiance_map(result.radiance.data(),
                                                 Outputs);
        Eigen::Map<const Eigen::VectorXd> cotangent_map(cotangent.data(),
                                                        Outputs);
        Eigen::Map<Eigen::VectorXd> gradient_map(result.gradient.data(),
                                                 Levels);
        sasktran2::OutputVJP<1> output(radiance_map, cotangent_map);
        output.set_derivative_gradient_memory("retrieval_state", gradient_map);
        engine.calculate_vjp(atmosphere, output);
        output.finalize();
        result.jvp = Eigen::VectorXd::Zero(Outputs);
        Eigen::Map<Eigen::VectorXd> jvp_map(result.jvp.data(), Outputs);
        sasktran2::OutputJVP<1> directional(radiance_map, jvp_map);
        const Eigen::VectorXd tangent =
            Eigen::VectorXd::LinSpaced(Levels, -0.7, 1.1);
        directional.set_derivative_tangent("retrieval_state", tangent);
        engine.calculate_jvp(atmosphere, directional);
        return result;
    }

    template <typename Values>
    void require_bits(const Values& actual, const Values& expected) {
        REQUIRE(actual.size() == expected.size());
        REQUIRE(actual.allFinite());
        REQUIRE(std::memcmp(actual.data(), expected.data(),
                            actual.size() * sizeof(double)) == 0);
    }

    void require_products(const Products& actual, const Products& expected) {
        require_bits(actual.forward, expected.forward);
        require_bits(actual.radiance, expected.radiance);
        require_bits(actual.gradient, expected.gradient);
        require_bits(actual.jvp, expected.jvp);
    }
} // namespace

TEST_CASE("Pinned transport cache preserves native products and true "
          "surface-only updates",
          "[successive_orders][engine][linearization][transport_cache]") {
    for (const int cache_count : {1, Wavelengths, Wavelengths + 2}) {
        for (const int iterations : {0, 2, 60}) {
            const auto config = make_config(cache_count, iterations);
            const auto reference_config = make_config(0, iterations);
            auto geometry = make_geometry();
            const auto viewing = make_viewing();
            Sasktran2<1> engine(config, &geometry, viewing);
            Sasktran2<1> reference(reference_config, &geometry, viewing);
            sasktran2::atmosphere::Atmosphere<1> atmosphere(
                Wavelengths, geometry, config, true);
            initialize(atmosphere);
            const Eigen::MatrixXd extinction =
                atmosphere.storage().total_extinction;
            const Eigen::MatrixXd ssa = atmosphere.storage().ssa;
            const Eigen::MatrixXd surface = atmosphere.surface().brdf_args();
            for (int evaluation = 0; evaluation < 4; ++evaluation) {
                const auto volume_revision = atmosphere.volume_revision();
                if (evaluation == 1) {
                    atmosphere.surface().brdf_args().array() += 0.11;
                    atmosphere.mark_surface_changed();
                    REQUIRE(atmosphere.volume_revision() == volume_revision);
                } else if (evaluation == 2) {
                    atmosphere.storage().total_extinction *= 1.25;
                    atmosphere.storage().ssa.array() -= 0.04;
                    atmosphere.mark_changed();
                } else if (evaluation == 3) {
                    atmosphere.storage().total_extinction = extinction;
                    atmosphere.storage().ssa = ssa;
                    atmosphere.surface().brdf_args() = surface;
                    atmosphere.mark_changed();
                }
                const auto actual = calculate(engine, atmosphere);
                const auto expected = calculate(reference, atmosphere);
                require_products(actual, expected);
                require_products(calculate(engine, atmosphere),
                                 calculate(reference, atmosphere));
            }
        }
    }
}

TEST_CASE("Pinned transport cache remains per worker with multiple "
          "wavelength workers",
          "[successive_orders][engine][linearization][transport_cache]") {
    const auto config = make_config(Wavelengths, 2, 2);
    const auto reference_config = make_config(0, 2, 2);
    auto geometry = make_geometry();
    const auto viewing = make_viewing();
    Sasktran2<1> engine(config, &geometry, viewing);
    Sasktran2<1> reference(reference_config, &geometry, viewing);
    sasktran2::atmosphere::Atmosphere<1> atmosphere(Wavelengths, geometry,
                                                    config, true);
    initialize(atmosphere);
    for (int evaluation = 0; evaluation < 3; ++evaluation) {
        require_products(calculate(engine, atmosphere),
                         calculate(reference, atmosphere));
    }
}

TEST_CASE("Negative native transport cache count is rejected",
          "[successive_orders][engine][transport_cache]") {
    const auto config = make_config(-1, 2);
    auto geometry = make_geometry();
    const auto viewing = make_viewing();
    REQUIRE_THROWS_AS(Sasktran2<1>(config, &geometry, viewing),
                      std::invalid_argument);
}
