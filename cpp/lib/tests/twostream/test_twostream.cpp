#include <sasktran2/test_helper.h>
#include <sasktran2.h>

#ifdef SKTRAN_CATCH2_VERSION3

TEST_CASE("twostreamBenchmark", "[sktran_do][lowlevel][benchmark]") {
#ifdef USE_OMP
    omp_set_num_threads(1);
#endif
    // Construct the geometry
    sasktran2::Coordinates coords(0.6, 0, 6371000,
                                  sasktran2::geometrytype::planeparallel);

    int nlyr = 20;
    int nwavel = 220000;
    int nstr = 2;

    Eigen::VectorXd grid_values(nlyr + 1);
    for (int i = 0; i < nlyr + 1; ++i) {
        grid_values(i) = i;
    }
    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::lower);

    sasktran2::Geometry1D geo(std::move(coords), std::move(grid));

    // Construct the Atmosphere
    sasktran2::atmosphere::AtmosphereGridStorageFull<1> storage(nwavel,
                                                                geo.size(), 16);
    sasktran2::atmosphere::Surface<1> surface(nwavel);
    surface.brdf_args().setConstant(0.5);

    sasktran2::atmosphere::Atmosphere<1> atmo(std::move(storage),
                                              std::move(surface), true);

    atmo.storage().total_extinction.setConstant(0.01);

    atmo.storage().ssa.setConstant(1.0);

    atmo.storage().leg_coeff.chip(0, 0).setConstant(1);
    atmo.storage().leg_coeff.chip(2, 0).setConstant(0.5);

    // Construct the Viewing rays
    sasktran2::viewinggeometry::ViewingGeometryContainer viewing_geometry;
    auto& los = viewing_geometry.observer_rays();

    los.emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            0.6, 60 * EIGEN_PI / 180.0, 0.6, 200000));

    // Construct the config
    sasktran2::Config config;
    config.set_multiple_scatter_source(
        sasktran2::Config::MultipleScatterSource::twostream);
    config.set_single_scatter_source(
        sasktran2::Config::SingleScatterSource::none);

    config.set_num_do_streams(nstr);
    config.set_do_backprop(true);

    // Make the engine
    Sasktran2<1> engine(config, &geo, viewing_geometry);

    sasktran2::OutputIdealDense<1> output;
    engine.calculate_radiance(atmo, output);
}

TEST_CASE("twostream_thermal_Benchmark", "[sktran_do][lowlevel][benchmark]") {
#ifdef USE_OMP
    omp_set_num_threads(1);
#endif
    // Construct the geometry
    sasktran2::Coordinates coords(0.6, 0, 6371000,
                                  sasktran2::geometrytype::planeparallel);

    int nlyr = 20;
    int nwavel = 220000;
    int nstr = 2;

    Eigen::VectorXd grid_values(nlyr + 1);
    for (int i = 0; i < nlyr + 1; ++i) {
        grid_values(i) = i;
    }
    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::lower);

    sasktran2::Geometry1D geo(std::move(coords), std::move(grid));

    // Construct the Atmosphere
    sasktran2::atmosphere::AtmosphereGridStorageFull<1> storage(nwavel,
                                                                geo.size(), 16);
    sasktran2::atmosphere::Surface<1> surface(nwavel);
    surface.brdf_args().setConstant(0.5);

    sasktran2::atmosphere::Atmosphere<1> atmo(std::move(storage),
                                              std::move(surface), true, true);

    atmo.storage().total_extinction.setConstant(0.01);

    atmo.storage().ssa.setConstant(0.8);

    atmo.storage().leg_coeff.chip(0, 0).setConstant(1);
    atmo.storage().leg_coeff.chip(2, 0).setConstant(0.5);

    // Construct the Viewing rays
    sasktran2::viewinggeometry::ViewingGeometryContainer viewing_geometry;
    auto& los = viewing_geometry.observer_rays();

    los.emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            0.6, 60 * EIGEN_PI / 180.0, 0.6, 200000));

    // Construct the config
    sasktran2::Config config;
    config.set_emission_source(sasktran2::Config::EmissionSource::twostream);
    config.set_single_scatter_source(
        sasktran2::Config::SingleScatterSource::none);

    config.set_num_do_streams(nstr);
    config.set_do_backprop(true);

    // Make the engine
    Sasktran2<1> engine(config, &geo, viewing_geometry);

    sasktran2::OutputIdealDense<1> output;
    engine.calculate_radiance(atmo, output);
}

TEST_CASE("twostream_engine_ground_extinction_wf_fd",
          "[twostream][wf][extinction][.]") {
#ifdef USE_OMP
    omp_set_num_threads(1);
#endif
    const double csz = 0.6;
    const int nlyr = 13;
    const int nwavel = 1;

    sasktran2::Coordinates coords(csz, 0, 6372000,
                                  sasktran2::geometrytype::planeparallel);

    Eigen::VectorXd grid_values(nlyr + 1);
    for (int i = 0; i < nlyr + 1; ++i) {
        grid_values(i) = i * 5000.0;
    }
    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::linear);

    sasktran2::Geometry1D geo(std::move(coords), std::move(grid));

    sasktran2::atmosphere::AtmosphereGridStorageFull<1> storage(nwavel,
                                                                geo.size(), 16);
    sasktran2::atmosphere::Surface<1> surface(nwavel);
    surface.brdf_args().setConstant(0.5 / EIGEN_PI);

    sasktran2::atmosphere::Atmosphere<1> atmo(std::move(storage),
                                              std::move(surface), true);

    for (int i = 0; i < geo.size(); ++i) {
        atmo.storage().total_extinction(i, 0) =
            0.02 * std::exp(-0.25 * i) + 1e-12;
    }
    atmo.storage().ssa.setConstant(0.8);
    atmo.storage().leg_coeff.chip(0, 0).setConstant(1);
    atmo.storage().leg_coeff.chip(1, 0).setConstant(0.0);
    atmo.storage().leg_coeff.chip(2, 0).setConstant(0.5);

    sasktran2::viewinggeometry::ViewingGeometryContainer viewing_geometry;
    auto& los = viewing_geometry.observer_rays();
    los.emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            csz, 0.0, 0.8, 200000));

    sasktran2::Config config;
    config.set_multiple_scatter_source(
        sasktran2::Config::MultipleScatterSource::twostream);
    config.set_single_scatter_source(
        sasktran2::Config::SingleScatterSource::none);
    config.set_num_do_streams(2);

    Sasktran2<1> engine(config, &geo, viewing_geometry);

    sasktran2::OutputIdealDense<1> output, output_p, output_m;
    engine.calculate_radiance(atmo, output);

    Eigen::VectorXd analytic(geo.size());
    Eigen::VectorXd numeric_reuse(geo.size());
    Eigen::VectorXd numeric_fresh(geo.size());

    auto calc_with_fresh_engine = [&]() {
        Sasktran2<1> engine_local(config, &geo, viewing_geometry);
        sasktran2::OutputIdealDense<1> out_local;
        engine_local.calculate_radiance(atmo, out_local);
        return out_local.radiance().value(0);
    };

    const double dfrac = 1e-4;
    for (int i = 0; i < geo.size(); ++i) {
        const double dk = dfrac * atmo.storage().total_extinction(i, 0);

        atmo.storage().total_extinction(i, 0) += dk;
        engine.calculate_radiance(atmo, output_p);

        atmo.storage().total_extinction(i, 0) -= 2 * dk;
        engine.calculate_radiance(atmo, output_m);

        atmo.storage().total_extinction(i, 0) += dk;

        analytic(i) = output.radiance().deriv(0, i);
        numeric_reuse(i) =
            (output_p.radiance().value(0) - output_m.radiance().value(0)) /
            (2 * dk);

        atmo.storage().total_extinction(i, 0) += dk;
        const double fp_fresh = calc_with_fresh_engine();
        atmo.storage().total_extinction(i, 0) -= 2 * dk;
        const double fm_fresh = calc_with_fresh_engine();
        atmo.storage().total_extinction(i, 0) += dk;

        numeric_fresh(i) = (fp_fresh - fm_fresh) / (2 * dk);
    }

    const double max_abs_analytic = analytic.cwiseAbs().maxCoeff();
    REQUIRE(max_abs_analytic > 0);

    for (int i = 0; i < geo.size(); ++i) {
        const double rel_reuse =
            std::abs(analytic(i) - numeric_reuse(i)) / max_abs_analytic;
        const double rel_fresh =
            std::abs(analytic(i) - numeric_fresh(i)) / max_abs_analytic;
        CAPTURE(i, analytic(i), numeric_reuse(i), numeric_fresh(i), rel_reuse,
                rel_fresh);
        REQUIRE(std::abs(numeric_reuse(i) - numeric_fresh(i)) <
                1e-8 * std::max(1.0, std::abs(numeric_fresh(i))));
        REQUIRE(rel_reuse < 1e-5);
    }
}

#endif
