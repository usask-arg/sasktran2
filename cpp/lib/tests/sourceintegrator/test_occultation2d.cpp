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

    sasktran2::raytracing::StructuredLayer2D
    layer(int altitude_cell, int horizontal_cell,
          std::array<double, 4> weights) {
        sasktran2::raytracing::StructuredLayer2D result;
        result.altitude_cell = altitude_cell;
        result.horizontal_cell = horizontal_cell;
        result.integrated_od.weights = weights;
        return result;
    }

    double direct_optical_depth(const sasktran2::raytracing::TracedRay2D& ray,
                                const sasktran2::Geometry2D& geometry,
                                Eigen::Ref<const Eigen::VectorXd> extinction) {
        double result = 0.0;
        for (const auto& traced_layer : ray.layers) {
            const std::array<int, 4> indices = {
                geometry.location_index(traced_layer.altitude_cell,
                                        traced_layer.horizontal_cell),
                geometry.location_index(traced_layer.altitude_cell + 1,
                                        traced_layer.horizontal_cell),
                geometry.location_index(traced_layer.altitude_cell,
                                        traced_layer.horizontal_cell + 1),
                geometry.location_index(traced_layer.altitude_cell + 1,
                                        traced_layer.horizontal_cell + 1)};
            for (int local_index = 0; local_index < 4; ++local_index) {
                result += traced_layer.integrated_od.weights[local_index] *
                          extinction[indices[local_index]];
            }
        }
        return result;
    }

    Eigen::VectorXd
    direct_path_weights(const sasktran2::raytracing::TracedRay2D& ray,
                        const sasktran2::Geometry2D& geometry) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(geometry.size());
        for (const auto& traced_layer : ray.layers) {
            const std::array<int, 4> indices = {
                geometry.location_index(traced_layer.altitude_cell,
                                        traced_layer.horizontal_cell),
                geometry.location_index(traced_layer.altitude_cell + 1,
                                        traced_layer.horizontal_cell),
                geometry.location_index(traced_layer.altitude_cell,
                                        traced_layer.horizontal_cell + 1),
                geometry.location_index(traced_layer.altitude_cell + 1,
                                        traced_layer.horizontal_cell + 1)};
            for (int local_index = 0; local_index < 4; ++local_index) {
                result[indices[local_index]] +=
                    traced_layer.integrated_od.weights[local_index];
            }
        }
        return result;
    }

    Eigen::Vector3d radial_direction(double angle) {
        return {std::sin(angle), 0.0, std::cos(angle)};
    }

    sasktran2::viewinggeometry::ViewingRay
    viewing_ray(const Eigen::Vector3d& observer, const Eigen::Vector3d& look) {
        sasktran2::viewinggeometry::ViewingRay result;
        result.observer.position = observer;
        result.look_away = look.normalized();
        result.relative_azimuth = 0.0;
        return result;
    }

    class InteriorTestSource : public SourceTermInterface<1> {
      public:
        void integrated_source(
            int, int, int, int, int,
            const sasktran2::raytracing::SphericalLayer&,
            const sasktran2::SparseODDualView&,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&,
            IntegrationDirection) const override {}

        void end_of_ray_source(
            int, int, int, int,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&)
            const override {}

        void start_of_ray_source(
            int, int, int, int,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&)
            const override {}
    };
} // namespace

TEST_CASE("Atmosphere allocates flattened storage from Geometry2D",
          "[sourceintegrator][occultation][geometry2d]") {
    const auto geometry = geometry2d();
    sasktran2::Config config;
    constexpr int num_wavelengths = 3;

    sasktran2::atmosphere::Atmosphere<1> atmosphere(num_wavelengths, geometry,
                                                    config, true);

    REQUIRE(atmosphere.storage().total_extinction.rows() == geometry.size());
    REQUIRE(atmosphere.storage().total_extinction.cols() == num_wavelengths);
    REQUIRE(atmosphere.storage().ssa.rows() == geometry.size());
    REQUIRE(atmosphere.storage().emission_source.rows() == geometry.size());
    REQUIRE(atmosphere.num_deriv() >= 2 * geometry.size());
}

TEST_CASE("SourceIntegrator applies 2D occultation transmission and native "
          "derivatives",
          "[sourceintegrator][occultation][geometry2d]") {
    const auto geometry = geometry2d();
    sasktran2::Config config;
    constexpr int num_wavelengths = 3;
    sasktran2::atmosphere::Atmosphere<1> atmosphere(num_wavelengths, geometry,
                                                    config, true);

    for (int location = 0; location < geometry.size(); ++location) {
        for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
            atmosphere.storage().total_extinction(location, wavelength) =
                0.002 * (1.0 + 0.17 * location + 0.31 * wavelength);
        }
    }

    std::vector<sasktran2::raytracing::TracedRay2D> rays(3);
    rays[1].layers = {
        layer(0, 0, {1.0, 2.0, 3.0, 4.0}),
        layer(1, 1, {0.5, 1.5, 2.5, 3.5}),
    };
    rays[2].ground_is_hit = true;
    rays[2].layers = rays[1].layers;

    sasktran2::SourceIntegrator<1> integrator(true);
    integrator.initialize_geometry(rays, geometry);
    integrator.initialize_atmosphere(atmosphere);

    sasktran2::solartransmission::OccultationSource<1> source;
    source.initialize_geometry(rays);
    source.initialize_atmosphere(atmosphere);
    std::vector<SourceTermInterface<1>*> sources = {&source};

    Eigen::MatrixXd optical_depth(num_wavelengths, rays.size());
    integrator.integrate_optical_depth(optical_depth);

    for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
        for (int ray_index = 0; ray_index < rays.size(); ++ray_index) {
            sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>
                transmission(1, atmosphere.num_deriv(), true);
            integrator.integrate(transmission, sources, wavelength, ray_index,
                                 0, 0);

            const double od = direct_optical_depth(
                rays[ray_index], geometry,
                atmosphere.storage().total_extinction.col(wavelength));
            REQUIRE(optical_depth(wavelength, ray_index) ==
                    Catch::Approx(od).margin(1e-14));

            if (rays[ray_index].ground_is_hit) {
                REQUIRE(transmission.value[0] == 0.0);
                REQUIRE(transmission.deriv.isZero(0.0));
                continue;
            }

            const double expected = std::exp(-od);
            REQUIRE(transmission.value[0] ==
                    Catch::Approx(expected).epsilon(1e-13));
            const Eigen::VectorXd path =
                direct_path_weights(rays[ray_index], geometry);
            for (int location = 0; location < geometry.size(); ++location) {
                REQUIRE(
                    transmission.deriv(0, location) ==
                    Catch::Approx(-expected * path[location]).margin(1e-13));
            }
            REQUIRE(transmission.deriv
                        .rightCols(transmission.deriv.cols() - geometry.size())
                        .isZero(0.0));
        }
    }
}

TEST_CASE("OccultationSource blocks ground-terminated 1D and 2D rays",
          "[sourceintegrator][occultation]") {
    sasktran2::solartransmission::OccultationSource<1> source;

    SECTION("1D") {
        sasktran2::viewinggeometry::InternalViewingGeometry viewing;
        viewing.traced_rays.resize(2);
        viewing.traced_rays[1].ground_is_hit = true;
        source.initialize_geometry(viewing);
    }

    SECTION("2D") {
        std::vector<sasktran2::raytracing::TracedRay2D> rays(2);
        rays[1].ground_is_hit = true;
        source.initialize_geometry(rays);
    }

    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> space(1, 0, true);
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> ground(1, 0,
                                                                     true);
    source.end_of_ray_source(0, 0, 0, 0, space);
    source.end_of_ray_source(0, 1, 0, 0, ground);
    REQUIRE(space.value[0] == 1.0);
    REQUIRE(ground.value[0] == 0.0);
}

TEST_CASE("2D occultation transmission handles extreme OD and validates "
          "atmosphere size",
          "[sourceintegrator][occultation][geometry2d]") {
    const auto geometry = geometry2d();
    std::vector<sasktran2::raytracing::TracedRay2D> rays(1);
    rays[0].layers = {layer(0, 0, {1.0, 2.0, 3.0, 4.0})};
    sasktran2::SourceIntegrator<1> integrator(false);
    integrator.initialize_geometry(rays, geometry);

    sasktran2::Config config;
    sasktran2::atmosphere::Atmosphere<1> atmosphere(1, geometry, config, false);
    sasktran2::solartransmission::OccultationSource<1> source;
    source.initialize_geometry(rays);
    std::vector<SourceTermInterface<1>*> sources = {&source};

    SECTION("zero extinction") {
        atmosphere.storage().total_extinction.setZero();
        integrator.initialize_atmosphere(atmosphere);
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> transmission(
            1, 0, true);
        integrator.integrate(transmission, sources, 0, 0, 0, 0);
        REQUIRE(transmission.value[0] == 1.0);
    }

    SECTION("optically thick") {
        atmosphere.storage().total_extinction.setConstant(1e3);
        integrator.initialize_atmosphere(atmosphere);
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> transmission(
            1, 0, true);
        integrator.integrate(transmission, sources, 0, 0, 0, 0);
        REQUIRE(transmission.value[0] == 0.0);
        REQUIRE(std::isfinite(transmission.value[0]));
    }

    SECTION("mismatched atmosphere") {
        sasktran2::atmosphere::AtmosphereGridStorageFull<1> storage(
            1, geometry.size() - 1, config.num_singlescatter_moments());
        sasktran2::atmosphere::Surface<1> surface(1);
        sasktran2::atmosphere::Atmosphere<1> mismatched(
            std::move(storage), std::move(surface), false);
        REQUIRE_THROWS_AS(integrator.initialize_atmosphere(mismatched),
                          std::invalid_argument);
    }
}

TEST_CASE("Regular 2D integration rejects interior sources before mutation",
          "[sourceintegrator][geometry2d]") {
    const auto geometry = geometry2d();
    std::vector<sasktran2::raytracing::TracedRay2D> rays(1);
    rays[0].layers = {layer(0, 0, {1.0, 2.0, 3.0, 4.0})};
    sasktran2::Config config;
    sasktran2::atmosphere::Atmosphere<1> atmosphere(1, geometry, config, false);
    atmosphere.storage().total_extinction.setConstant(0.1);
    sasktran2::SourceIntegrator<1> integrator(false);
    integrator.initialize_geometry(rays, geometry);
    integrator.initialize_atmosphere(atmosphere);

    InteriorTestSource source;
    std::vector<SourceTermInterface<1>*> sources = {&source};
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> radiance(1, 0,
                                                                       true);
    REQUIRE_THROWS_AS(integrator.integrate(radiance, sources, 0, 0, 0, 0),
                      std::invalid_argument);
    REQUIRE(radiance.value[0] == 0.0);
}

#ifdef SKTRAN_RUST_SUPPORT

TEST_CASE("RustRayTracer2D feeds refracted multi-wavelength occultation "
          "transmission",
          "[sourceintegrator][occultation][geometry2d][rust]") {
    const auto geometry = geometry2d();
    sasktran2::raytracing::RustRayTracer2D tracer(geometry);
    std::vector<sasktran2::raytracing::TracedRay2D> rays(3);

    tracer.trace_ray(
        viewing_ray(15.0 * radial_direction(-0.25), Eigen::Vector3d::UnitX()),
        rays[0]);
    const double tangent_radius = 15.0;
    const double observer_radius = 40.0;
    tracer.trace_ray(
        viewing_ray(
            {0.0, 0.0, observer_radius},
            {tangent_radius / observer_radius, 0.0,
             -std::sqrt(1.0 - std::pow(tangent_radius / observer_radius, 2))}),
        vector({1.02, 1.01, 1.0}), rays[1]);
    tracer.trace_ray(
        viewing_ray(40.0 * radial_direction(0.75), -radial_direction(0.75)),
        rays[2]);

    REQUIRE_FALSE(rays[0].ground_is_hit);
    REQUIRE_FALSE(rays[1].ground_is_hit);
    REQUIRE(rays[2].ground_is_hit);
    REQUIRE_FALSE(rays[1].is_straight);

    sasktran2::Config config;
    constexpr int num_wavelengths = 4;
    sasktran2::atmosphere::Atmosphere<1> atmosphere(num_wavelengths, geometry,
                                                    config, true);
    for (int location = 0; location < geometry.size(); ++location) {
        for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
            atmosphere.storage().total_extinction(location, wavelength) =
                1e-3 * (1.0 + 0.11 * location + 0.23 * wavelength);
        }
    }

    sasktran2::SourceIntegrator<1> integrator(true);
    integrator.initialize_geometry(rays, geometry);
    integrator.initialize_atmosphere(atmosphere);
    sasktran2::solartransmission::OccultationSource<1> source;
    source.initialize_geometry(rays);
    std::vector<SourceTermInterface<1>*> sources = {&source};

    for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
        for (int ray_index = 0; ray_index < rays.size(); ++ray_index) {
            sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>
                transmission(1, atmosphere.num_deriv(), true);
            integrator.integrate(transmission, sources, wavelength, ray_index,
                                 0, 0);
            const double od = direct_optical_depth(
                rays[ray_index], geometry,
                atmosphere.storage().total_extinction.col(wavelength));
            const double expected =
                rays[ray_index].ground_is_hit ? 0.0 : std::exp(-od);
            REQUIRE(transmission.value[0] ==
                    Catch::Approx(expected).margin(1e-13));
        }
    }

    const int perturbed_location = geometry.location_index(1, 1);
    constexpr int wavelength = 2;
    constexpr double delta = 1e-7;
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> analytic(
        1, atmosphere.num_deriv(), true);
    integrator.integrate(analytic, sources, wavelength, 1, 0, 0);
    atmosphere.storage().total_extinction(perturbed_location, wavelength) +=
        delta;
    integrator.initialize_atmosphere(atmosphere);
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> perturbed(
        1, atmosphere.num_deriv(), true);
    integrator.integrate(perturbed, sources, wavelength, 1, 0, 0);
    REQUIRE((perturbed.value[0] - analytic.value[0]) / delta ==
            Catch::Approx(analytic.deriv(0, perturbed_location)).epsilon(2e-6));
}

TEST_CASE("Horizontally constant 2D occultation transmission matches 1D",
          "[sourceintegrator][occultation][geometry2d][rust][parity]") {
    const auto geometry_2d = geometry2d();
    sasktran2::Geometry1D geometry_1d(0.6, 0.3, earth_radius,
                                      vector({0.0, 10.0, 20.0}));
    const auto ray =
        viewing_ray(15.0 * radial_direction(-0.25), Eigen::Vector3d::UnitX());

    std::vector<sasktran2::raytracing::TracedRay> rays_1d(1);
    sasktran2::raytracing::RustRayTracer tracer_1d(geometry_1d);
    tracer_1d.trace_ray(ray, rays_1d[0]);
    std::vector<sasktran2::raytracing::TracedRay2D> rays_2d(1);
    sasktran2::raytracing::RustRayTracer2D tracer_2d(geometry_2d);
    tracer_2d.trace_ray(ray, rays_2d[0]);

    sasktran2::Config config;
    constexpr int num_wavelengths = 3;
    sasktran2::atmosphere::Atmosphere<1> atmosphere_1d(
        num_wavelengths, geometry_1d, config, false);
    sasktran2::atmosphere::Atmosphere<1> atmosphere_2d(
        num_wavelengths, geometry_2d, config, false);
    for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
        for (int altitude = 0; altitude < geometry_2d.num_altitudes();
             ++altitude) {
            const double extinction =
                1e-3 * (1.0 + 0.4 * altitude + 0.2 * wavelength);
            atmosphere_1d.storage().total_extinction(altitude, wavelength) =
                extinction;
            for (int horizontal = 0;
                 horizontal < geometry_2d.num_horizontal_locations();
                 ++horizontal) {
                atmosphere_2d.storage().total_extinction(
                    geometry_2d.location_index(altitude, horizontal),
                    wavelength) = extinction;
            }
        }
    }

    sasktran2::SourceIntegrator<1> integrator_1d(false);
    integrator_1d.initialize_geometry(rays_1d, geometry_1d);
    integrator_1d.initialize_atmosphere(atmosphere_1d);
    sasktran2::SourceIntegrator<1> integrator_2d(false);
    integrator_2d.initialize_geometry(rays_2d, geometry_2d);
    integrator_2d.initialize_atmosphere(atmosphere_2d);

    sasktran2::viewinggeometry::InternalViewingGeometry internal_1d;
    internal_1d.traced_rays = rays_1d;
    sasktran2::solartransmission::OccultationSource<1> source_1d;
    source_1d.initialize_geometry(internal_1d);
    sasktran2::solartransmission::OccultationSource<1> source_2d;
    source_2d.initialize_geometry(rays_2d);
    std::vector<SourceTermInterface<1>*> sources_2d = {&source_2d};

    for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>
            transmission_1d(1, 0, true);
        std::vector<SourceTermInterface<1>*> sources = {&source_1d};
        integrator_1d.integrate(transmission_1d, sources, wavelength, 0, 0, 0);
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>
            transmission_2d(1, 0, true);
        integrator_2d.integrate(transmission_2d, sources_2d, wavelength, 0, 0,
                                0);
        REQUIRE(transmission_2d.value[0] ==
                Catch::Approx(transmission_1d.value[0]).epsilon(5e-8));
    }
}

#endif
