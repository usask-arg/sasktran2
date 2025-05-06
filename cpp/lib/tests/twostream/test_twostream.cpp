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
#endif
