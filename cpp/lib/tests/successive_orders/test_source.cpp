#include <sasktran2.h>
#include <sasktran2/test_helper.h>

#include <cmath>
#include <memory>
#include <optional>
#include <utility>

namespace {
    constexpr int NumAltitudes = 9;
    constexpr int NumWavelengths = 2;
    constexpr int NumLos = 3;

    sasktran2::Config successive_orders_config() {
        sasktran2::Config config;
        config.set_num_threads(2);
        config.set_threading_model(
            sasktran2::Config::ThreadingModel::wavelength);
        config.set_single_scatter_source(
            sasktran2::Config::SingleScatterSource::exact);
        config.set_multiple_scatter_source(
            sasktran2::Config::MultipleScatterSource::successive_orders_cpp);
        config.set_num_hr_incoming(14);
        config.set_num_hr_outgoing(14);
        config.set_num_hr_spherical_iterations(30);
        config.set_successive_orders_relative_tolerance(1.0e-8);
        config.set_successive_orders_absolute_tolerance(1.0e-12);
        config.set_successive_orders_anderson_depth(3);
        config.set_num_do_streams(8);
        config.set_num_singlescatter_moments(8);
        config.set_apply_delta_scaling(false);
        return config;
    }

    sasktran2::Geometry1D successive_orders_geometry() {
        sasktran2::Coordinates coordinates(0.55, 0.15, 6372000.0,
                                           sasktran2::geometrytype::spherical);
        Eigen::VectorXd altitudes(NumAltitudes);
        for (int altitude = 0; altitude < NumAltitudes; ++altitude) {
            altitudes(altitude) = altitude * 5000.0;
        }
        sasktran2::grids::AltitudeGrid altitude_grid(
            std::move(altitudes), sasktran2::grids::gridspacing::constant,
            sasktran2::grids::outofbounds::extend,
            sasktran2::grids::interpolation::linear);
        return {std::move(coordinates), std::move(altitude_grid)};
    }

    sasktran2::viewinggeometry::ViewingGeometryContainer
    successive_orders_viewing() {
        sasktran2::viewinggeometry::ViewingGeometryContainer viewing;
        viewing.observer_rays().emplace_back(
            std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
                0.55, 0.35, 0.72, 100000.0));
        viewing.observer_rays().emplace_back(
            std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(
                10000.0, -0.4, 100000.0, 0.55));
        viewing.observer_rays().emplace_back(
            std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(
                22500.0, 0.6, 100000.0, 0.55));
        return viewing;
    }

    void initialize_atmosphere(sasktran2::atmosphere::Atmosphere<1>& atmosphere,
                               int atmosphere_lifetime) {
        atmosphere.storage().resize_derivatives(atmosphere_lifetime);
        for (int wavelength = 0; wavelength < NumWavelengths; ++wavelength) {
            for (int altitude = 0; altitude < NumAltitudes; ++altitude) {
                atmosphere.storage().total_extinction(altitude, wavelength) =
                    (1.4e-5 * std::exp(-altitude / 1.7) + 2.0e-9) *
                    (0.9 + 0.2 * wavelength) *
                    (1.0 + 0.15 * atmosphere_lifetime);
                atmosphere.storage().ssa(altitude, wavelength) =
                    0.90 - 0.008 * wavelength - 0.002 * altitude -
                    0.03 * atmosphere_lifetime;
            }
        }
        atmosphere.storage().leg_coeff.chip(0, 0).setConstant(1.0);
        atmosphere.storage().leg_coeff.chip(1, 0).setConstant(0.08);
        atmosphere.storage().leg_coeff.chip(2, 0).setConstant(0.5);
        atmosphere.surface().brdf_args().row(0).setConstant(
            0.12 + 0.05 * atmosphere_lifetime);

        auto& mapping =
            atmosphere.storage().get_derivative_mapping("retrieval_state");
        mapping.allocate_extinction_derivatives();
        mapping.allocate_ssa_derivatives();
        mapping.native_mapping().d_extinction->setConstant(1.0e-5);
        mapping.native_mapping().d_ssa->setConstant(1.0e-2);
    }
} // namespace

TEST_CASE("Successive-orders revision-zero atmosphere remains directly "
          "mutable",
          "[successive_orders][engine][linearization][cache]") {
    const auto config = successive_orders_config();
    auto geometry = successive_orders_geometry();
    const auto viewing = successive_orders_viewing();
    Sasktran2<1> reused_engine(config, &geometry, viewing);
    sasktran2::atmosphere::Atmosphere<1> atmosphere(NumWavelengths, geometry,
                                                    config, true);
    initialize_atmosphere(atmosphere, 0);
    REQUIRE(atmosphere.revision() == 0);

    const auto calculate = [&](Sasktran2<1>& engine) {
        Eigen::VectorXd radiance =
            Eigen::VectorXd::Zero(NumWavelengths * NumLos);
        const Eigen::VectorXd cotangent =
            Eigen::VectorXd::LinSpaced(NumWavelengths * NumLos, 0.4, 1.1);
        Eigen::VectorXd gradient = Eigen::VectorXd::Zero(NumAltitudes);
        Eigen::Map<Eigen::VectorXd> radiance_map(radiance.data(),
                                                 radiance.size());
        Eigen::Map<const Eigen::VectorXd> cotangent_map(cotangent.data(),
                                                        cotangent.size());
        Eigen::Map<Eigen::VectorXd> gradient_map(gradient.data(),
                                                 gradient.size());
        sasktran2::OutputVJP<1> output(radiance_map, cotangent_map);
        output.set_derivative_gradient_memory("retrieval_state", gradient_map);
        engine.calculate_vjp(atmosphere, output);
        output.finalize();
        return std::make_pair(std::move(radiance), std::move(gradient));
    };

    const auto initial = calculate(reused_engine);
    atmosphere.storage().total_extinction *= 1.25;
    atmosphere.storage().ssa.array() -= 0.04;
    const auto reused = calculate(reused_engine);
    Sasktran2<1> fresh_engine(config, &geometry, viewing);
    const auto fresh = calculate(fresh_engine);

    REQUIRE_FALSE(reused.first.isApprox(initial.first, 1.0e-7));
    REQUIRE(reused.first.isApprox(fresh.first, 1.0e-7));
    REQUIRE(reused.second.isApprox(fresh.second, 1.0e-7));
}

TEST_CASE("Successive-orders scalar VJP survives repeated atmosphere and "
          "engine lifetimes",
          "[successive_orders][engine][linearization][lifetime]") {
    const auto config = successive_orders_config();
    auto geometry = successive_orders_geometry();
    const auto viewing = successive_orders_viewing();

    constexpr int evaluations_per_atmosphere = 32;
    constexpr int atmosphere_lifetimes = 2;
    constexpr int engine_lifetimes = 2;
    constexpr int output_size = NumWavelengths * NumLos;

    for (int engine_lifetime = 0; engine_lifetime < engine_lifetimes;
         ++engine_lifetime) {
        Sasktran2<1> engine(config, &geometry, viewing);
        std::optional<sasktran2::atmosphere::Atmosphere<1>> atmosphere_slot;

        for (int atmosphere_lifetime = 0;
             atmosphere_lifetime < atmosphere_lifetimes;
             ++atmosphere_lifetime) {
            atmosphere_slot.emplace(NumWavelengths, geometry, config, true);
            auto& atmosphere = *atmosphere_slot;
            initialize_atmosphere(atmosphere, atmosphere_lifetime);
            const Eigen::MatrixXd base_extinction =
                atmosphere.storage().total_extinction;
            const Eigen::MatrixXd base_ssa = atmosphere.storage().ssa;

            // Reproduce allocator-address reuse with the same local revision
            // as the preceding atmosphere. A cache key made only from object
            // address and revision cannot distinguish these two lifetimes.
            if (atmosphere_lifetime != 0) {
                for (int revision = 1; revision < evaluations_per_atmosphere;
                     ++revision) {
                    atmosphere.mark_changed();
                }
            }

            Eigen::VectorXd radiance = Eigen::VectorXd::Zero(output_size);
            Eigen::VectorXd cotangent =
                Eigen::VectorXd::LinSpaced(output_size, 0.4, 1.1);
            Eigen::VectorXd gradient = Eigen::VectorXd::Zero(NumAltitudes);

            for (int evaluation = 0; evaluation < evaluations_per_atmosphere;
                 ++evaluation) {
                const double extinction_scale =
                    1.0 + 2.0e-3 * std::sin(0.37 * evaluation);
                const double ssa_shift = 5.0e-4 * std::cos(0.23 * evaluation);
                atmosphere.storage().total_extinction =
                    extinction_scale * base_extinction;
                atmosphere.storage().ssa =
                    (base_ssa.array() + ssa_shift).matrix();
                atmosphere.mark_changed();

                gradient.setZero();
                Eigen::Map<Eigen::VectorXd> radiance_map(radiance.data(),
                                                         radiance.size());
                Eigen::Map<const Eigen::VectorXd> cotangent_map(
                    cotangent.data(), cotangent.size());
                Eigen::Map<Eigen::VectorXd> gradient_map(gradient.data(),
                                                         gradient.size());
                sasktran2::OutputVJP<1> output(radiance_map, cotangent_map);
                output.set_derivative_gradient_memory("retrieval_state",
                                                      gradient_map);
                engine.calculate_vjp(atmosphere, output);
                output.finalize();

                REQUIRE(radiance.allFinite());
                REQUIRE(gradient.allFinite());

                if (evaluation == 0) {
                    Sasktran2<1> reference_engine(config, &geometry, viewing);
                    Eigen::VectorXd reference_radiance =
                        Eigen::VectorXd::Zero(output_size);
                    Eigen::VectorXd reference_gradient =
                        Eigen::VectorXd::Zero(NumAltitudes);
                    Eigen::Map<Eigen::VectorXd> reference_radiance_map(
                        reference_radiance.data(), reference_radiance.size());
                    Eigen::Map<Eigen::VectorXd> reference_gradient_map(
                        reference_gradient.data(), reference_gradient.size());
                    sasktran2::OutputVJP<1> reference_output(
                        reference_radiance_map, cotangent_map);
                    reference_output.set_derivative_gradient_memory(
                        "retrieval_state", reference_gradient_map);
                    reference_engine.calculate_vjp(atmosphere,
                                                   reference_output);
                    reference_output.finalize();

                    REQUIRE(radiance.isApprox(reference_radiance, 1.0e-7));
                    REQUIRE(gradient.isApprox(reference_gradient, 1.0e-7));
                }
            }
        }
    }
}
