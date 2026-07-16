#include <sasktran2/test_helper.h>

#include <sasktran2.h>
#include <sasktran2/solartransmission.h>
#include <numeric>

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

    void add_layer(sasktran2::raytracing::TracedRay& ray,
                   const sasktran2::Geometry2D& geometry, int altitude_cell,
                   int horizontal_cell, const std::array<double, 4>& od_weights,
                   const std::array<double, 4>& entrance_weights = {0.25, 0.25,
                                                                    0.25, 0.25},
                   const std::array<double, 4>& exit_weights = {0.25, 0.25,
                                                                0.25, 0.25}) {
        sasktran2::raytracing::TracedLayer result;
        result.layer_distance =
            std::accumulate(od_weights.begin(), od_weights.end(), 0.0);
        result.od_quad_start_fraction = 0.5;
        result.od_quad_end_fraction = 0.5;
        ray.layers.push_back(result);
        const std::array<int, 4> indices = {
            geometry.location_index(altitude_cell, horizontal_cell),
            geometry.location_index(altitude_cell + 1, horizontal_cell),
            geometry.location_index(altitude_cell, horizontal_cell + 1),
            geometry.location_index(altitude_cell + 1, horizontal_cell + 1)};
        ray.set_layer_weights(ray.layers.size() - 1, indices, entrance_weights,
                              exit_weights, od_weights);
    }

    template <typename Source>
    void initialize_source_geometry(
        Source& source,
        const std::vector<sasktran2::raytracing::TracedRay>& rays) {
        sasktran2::viewinggeometry::InternalViewingGeometry viewing;
        viewing.traced_rays = rays;
        source.initialize_geometry(viewing);
    }

    double direct_optical_depth(const sasktran2::raytracing::TracedRay& ray,
                                const sasktran2::Geometry2D& geometry,
                                Eigen::Ref<const Eigen::VectorXd> extinction) {
        double result = 0.0;
        for (std::size_t layer_index = 0; layer_index < ray.layers.size();
             ++layer_index) {
            const auto weights = ray.optical_depth_weights(layer_index);
            for (std::size_t index = 0; index < weights.size(); ++index) {
                result +=
                    weights[index].second * extinction[weights[index].first];
            }
        }
        return result;
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
        void
        integrated_source(const sasktran2::WavelengthBlock<>&, int, int, int,
                          int, const sasktran2::raytracing::TracedLayer&,
                          const sasktran2::raytracing::GridWeightStencilView&,
                          const sasktran2::raytracing::GridWeightStencilView&,
                          const sasktran2::WavelengthBlockODView&,
                          sasktran2::WavelengthBlockDual<1>&,
                          IntegrationDirection) const override {}

        bool supports_geometry_dimension(int) const override { return false; }

        void
        end_of_ray_source(const sasktran2::WavelengthBlock<>&, int, int, int,
                          sasktran2::WavelengthBlockDual<1>&) const override {}

        void
        start_of_ray_source(const sasktran2::WavelengthBlock<>&, int, int, int,
                            sasktran2::WavelengthBlockDual<1>&) const override {
        }
    };

    template <int NSTOKES>
    void
    integrate_scalar(sasktran2::SourceIntegrator<NSTOKES>& integrator,
                     sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& radiance,
                     const std::vector<SourceTermInterface<NSTOKES>*>& sources,
                     int wavelength, int ray_index, int wavelength_thread_index,
                     int thread_index) {
        const sasktran2::WavelengthBlock<> block{wavelength, 1};
        sasktran2::WavelengthBlockDual<NSTOKES> block_radiance;
        block_radiance.resize(1, radiance.deriv.cols(), true);
        integrator.integrate(block_radiance, sources, block, ray_index,
                             wavelength_thread_index, thread_index);
        radiance.value = block_radiance.value.col(0);
        for (int derivative = 0; derivative < radiance.deriv.cols();
             ++derivative) {
            radiance.deriv.col(derivative) =
                block_radiance.derivative(derivative, 1).col(0);
        }
    }
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

    std::vector<sasktran2::raytracing::TracedRay> rays(3);
    add_layer(rays[1], geometry, 0, 0, {1.0, 2.0, 3.0, 4.0});
    add_layer(rays[1], geometry, 1, 1, {0.5, 1.5, 2.5, 3.5});
    rays[2] = rays[1];
    rays[2].ground_is_hit = true;

    sasktran2::SourceIntegrator<1> integrator(true);
    integrator.initialize_geometry(rays, geometry);
    integrator.initialize_atmosphere(atmosphere);

    sasktran2::solartransmission::OccultationSource<1> source;
    initialize_source_geometry(source, rays);
    source.initialize_atmosphere(atmosphere);
    std::vector<SourceTermInterface<1>*> sources = {&source};

    Eigen::MatrixXd optical_depth(num_wavelengths, rays.size());
    integrator.integrate_optical_depth(optical_depth);

    for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
        for (int ray_index = 0; ray_index < rays.size(); ++ray_index) {
            sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>
                transmission(1, atmosphere.num_deriv(), true);
            integrate_scalar(integrator, transmission, sources, wavelength,
                             ray_index, 0, 0);

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

TEST_CASE("SourceIntegrator restores derivatives after a derivative-free "
          "atmosphere",
          "[sourceintegrator][occultation][geometry2d][derivatives]") {
    const auto geometry = geometry2d();
    sasktran2::Config config;
    std::vector<sasktran2::raytracing::TracedRay> rays(1);
    add_layer(rays[0], geometry, 0, 0, {1.0, 2.0, 3.0, 4.0});

    sasktran2::atmosphere::Atmosphere<1> derivative_free(1, geometry, config,
                                                         false);
    derivative_free.storage().total_extinction.setConstant(0.01);
    sasktran2::atmosphere::Atmosphere<1> with_derivatives(1, geometry, config,
                                                          true);
    with_derivatives.storage().total_extinction.setConstant(0.01);

    sasktran2::SourceIntegrator<1> integrator(true);
    integrator.initialize_geometry(rays, geometry);
    sasktran2::solartransmission::OccultationSource<1> source;
    initialize_source_geometry(source, rays);
    std::vector<SourceTermInterface<1>*> sources = {&source};

    integrator.initialize_atmosphere(derivative_free);
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> first(1, 0, true);
    integrate_scalar(integrator, first, sources, 0, 0, 0, 0);

    integrator.initialize_atmosphere(with_derivatives);
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> second(
        1, with_derivatives.num_deriv(), true);
    integrate_scalar(integrator, second, sources, 0, 0, 0, 0);

    REQUIRE(second.deriv.leftCols(geometry.size()).cwiseAbs().maxCoeff() > 0.0);
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
        std::vector<sasktran2::raytracing::TracedRay> rays(2);
        rays[1].ground_is_hit = true;
        initialize_source_geometry(source, rays);
    }

    sasktran2::WavelengthBlockDual<1> space;
    sasktran2::WavelengthBlockDual<1> ground;
    space.resize(1, 0, true);
    ground.resize(1, 0, true);
    const sasktran2::WavelengthBlock<> block{0, 1};
    source.end_of_ray_source(block, 0, 0, 0, space);
    source.end_of_ray_source(block, 1, 0, 0, ground);
    REQUIRE(space.value(0, 0) == 1.0);
    REQUIRE(ground.value(0, 0) == 0.0);
}

TEST_CASE("2D occultation transmission handles extreme OD and validates "
          "atmosphere size",
          "[sourceintegrator][occultation][geometry2d]") {
    const auto geometry = geometry2d();
    std::vector<sasktran2::raytracing::TracedRay> rays(1);
    add_layer(rays[0], geometry, 0, 0, {1.0, 2.0, 3.0, 4.0});
    sasktran2::SourceIntegrator<1> integrator(false);
    integrator.initialize_geometry(rays, geometry);

    sasktran2::Config config;
    sasktran2::atmosphere::Atmosphere<1> atmosphere(1, geometry, config, false);
    sasktran2::solartransmission::OccultationSource<1> source;
    initialize_source_geometry(source, rays);
    std::vector<SourceTermInterface<1>*> sources = {&source};

    SECTION("zero extinction") {
        atmosphere.storage().total_extinction.setZero();
        integrator.initialize_atmosphere(atmosphere);
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> transmission(
            1, 0, true);
        integrate_scalar(integrator, transmission, sources, 0, 0, 0, 0);
        REQUIRE(transmission.value[0] == 1.0);
    }

    SECTION("optically thick") {
        atmosphere.storage().total_extinction.setConstant(1e3);
        integrator.initialize_atmosphere(atmosphere);
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> transmission(
            1, 0, true);
        integrate_scalar(integrator, transmission, sources, 0, 0, 0, 0);
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
    std::vector<sasktran2::raytracing::TracedRay> rays(1);
    add_layer(rays[0], geometry, 0, 0, {1.0, 2.0, 3.0, 4.0});
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
    REQUIRE_THROWS_AS(
        integrate_scalar(integrator, radiance, sources, 0, 0, 0, 0),
        std::invalid_argument);
    REQUIRE(radiance.value[0] == 0.0);
}

TEST_CASE("SourceIntegrator evaluates structured 2D volume emission and "
          "derivatives",
          "[sourceintegrator][emission][geometry2d]") {
    const auto geometry = geometry2d();
    std::vector<sasktran2::raytracing::TracedRay> rays(1);
    add_layer(rays[0], geometry, 0, 0, {1.0, 2.0, 3.0, 4.0},
              {1.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 1.0});
    auto& emission_layer = rays[0].layers[0];
    emission_layer.od_quad_start_fraction = 0.25;
    emission_layer.od_quad_end_fraction = 0.75;

    sasktran2::Config config;
    config.set_emission_source(
        sasktran2::Config::EmissionSource::volume_emission_rate);
    sasktran2::atmosphere::Atmosphere<1> atmosphere(1, geometry, config, true);
    atmosphere.storage().total_extinction.setZero();
    atmosphere.storage().emission_source.setZero();
    const int entrance_index = geometry.location_index(0, 0);
    const int exit_index = geometry.location_index(1, 1);
    atmosphere.storage().emission_source(entrance_index, 0) = 2.0;
    atmosphere.storage().emission_source(exit_index, 0) = 6.0;

    sasktran2::SourceIntegrator<1> integrator(true);
    integrator.initialize_geometry(rays, geometry);
    integrator.initialize_atmosphere(atmosphere);
    sasktran2::emission::EmissionSource<
        1, sasktran2::Config::EmissionSource::volume_emission_rate>
        source;
    initialize_source_geometry(source, rays);
    source.initialize_atmosphere(atmosphere);
    std::vector<SourceTermInterface<1>*> sources = {&source};
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> radiance(
        1, atmosphere.num_deriv(), true);

    integrate_scalar(integrator, radiance, sources, 0, 0, 0, 0);

    const double distance = emission_layer.layer_distance;
    REQUIRE(radiance.value[0] ==
            Catch::Approx(distance * (0.25 * 2.0 + 0.75 * 6.0)));
    const auto d_emission = radiance.d_emission(
        geometry.size(), atmosphere.num_scattering_deriv_groups());
    REQUIRE(d_emission(0, entrance_index) == Catch::Approx(0.25 * distance));
    REQUIRE(d_emission(0, exit_index) == Catch::Approx(0.75 * distance));
}

TEST_CASE("EmissionSource rejects missing atmospheric emission derivatives",
          "[sourceintegrator][emission]") {
    const auto geometry = geometry2d();
    sasktran2::Config config;
    sasktran2::atmosphere::Atmosphere<1> atmosphere(1, geometry, config, true);
    sasktran2::emission::EmissionSource<
        1, sasktran2::Config::EmissionSource::volume_emission_rate>
        source;

    REQUIRE_THROWS_AS(source.initialize_atmosphere(atmosphere),
                      std::invalid_argument);
}

#ifdef SKTRAN_RUST_SUPPORT

TEST_CASE("RustRayTracer2D feeds refracted multi-wavelength occultation "
          "transmission",
          "[sourceintegrator][occultation][geometry2d][rust]") {
    const auto geometry = geometry2d();
    sasktran2::raytracing::RustRayTracer2D tracer(geometry);
    std::vector<sasktran2::raytracing::TracedRay> rays(3);

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
    initialize_source_geometry(source, rays);
    std::vector<SourceTermInterface<1>*> sources = {&source};

    for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
        for (int ray_index = 0; ray_index < rays.size(); ++ray_index) {
            sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>
                transmission(1, atmosphere.num_deriv(), true);
            integrate_scalar(integrator, transmission, sources, wavelength,
                             ray_index, 0, 0);
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
    integrate_scalar(integrator, analytic, sources, wavelength, 1, 0, 0);
    atmosphere.storage().total_extinction(perturbed_location, wavelength) +=
        delta;
    integrator.initialize_atmosphere(atmosphere);
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> perturbed(
        1, atmosphere.num_deriv(), true);
    integrate_scalar(integrator, perturbed, sources, wavelength, 1, 0, 0);
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
    std::vector<sasktran2::raytracing::TracedRay> rays_2d(1);
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
    initialize_source_geometry(source_2d, rays_2d);
    std::vector<SourceTermInterface<1>*> sources_2d = {&source_2d};

    for (int wavelength = 0; wavelength < num_wavelengths; ++wavelength) {
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>
            transmission_1d(1, 0, true);
        std::vector<SourceTermInterface<1>*> sources = {&source_1d};
        integrate_scalar(integrator_1d, transmission_1d, sources, wavelength, 0,
                         0, 0);
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>
            transmission_2d(1, 0, true);
        integrate_scalar(integrator_2d, transmission_2d, sources_2d, wavelength,
                         0, 0, 0);
        REQUIRE(transmission_2d.value[0] ==
                Catch::Approx(transmission_1d.value[0]).epsilon(5e-8));
    }
}

#endif
