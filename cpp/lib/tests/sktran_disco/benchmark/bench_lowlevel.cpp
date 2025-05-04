#include <sasktran2/test_helper.h>
#include <sasktran2.h>

#ifdef SKTRAN_CATCH2_VERSION3

/*
TEST_CASE("2StreamBenchmark", "[sktran_do][lowlevel][benchmark]") {
    // TODO: 2 stream single layer sometimes fails with derivatives and no full
    // compile? Check this

    int nwavel = 1000;
    int nlyr = 20;
    int nderiv = 4;

    sasktran_disco_lowlevel::CPPApi cppapi(2, nlyr, nwavel, 1, 1,
                                           nderiv * nlyr);

    for (int i = 0; i < nwavel; ++i) {
        for (int j = 0; j < nlyr; ++j) {
            cppapi.od(j, i) = 0.2;
            cppapi.ssa(j, i) = 0.8;
            cppapi.a1(0, j, i) = 1;
            // cppapi.a1(2, j, i) = 0.5;

            for (int k = 0; k < nderiv; ++k) {
                cppapi.d_layerindex(k * nlyr + j) = j;
                cppapi.d_ssa(k * nlyr + j, i) = 1;
            }
        }
    }
    // Layer Construction
    double top_alt = 0.001;
    for (int j = 0; j < nlyr; ++j) {
        cppapi.layerboundaryaltitude(j) =
            top_alt - j * (top_alt / double(nlyr));
    }

    cppapi.earth_radius() = 6372000;

    // Coordinates construction
    cppapi.cos_sza() = 0.8;
    cppapi.cos_vza(0) = 0.7;
    cppapi.saa(0) = 0;

    sasktran_disco_lowlevel::Atmosphere atmosphere;
    sasktran_disco_lowlevel::Config config;
    sasktran_disco_lowlevel::WeightingFunctions weightingfunctions;
    sasktran_disco_lowlevel::ViewingGeometry geometry;
    sasktran_disco_lowlevel::Output output;
    cppapi.initialize_c_api(&atmosphere, &config, &geometry,
                            &weightingfunctions, &output);

    config.nthreads = 1;

    BENCHMARK("Test") {
        return sasktran_disco_lowlevel::calculate(
            &atmosphere, &config, &weightingfunctions, &geometry, &output);
    };
}
*/
TEST_CASE("8StreamBenchmark", "[sktran_do][lowlevel][benchmark]") {
#ifdef USE_OMP
    omp_set_num_threads(1);
#endif
    // Construct the geometry
    sasktran2::Coordinates coords(0.6, 0, 6371000,
                                  sasktran2::geometrytype::planeparallel);

    int nlyr = 20;
    int nwavel = 40;
    int nstr = 2;

    Eigen::VectorXd grid_values(nlyr + 1);
    for (int i = 0; i < nlyr + 1; ++i) {
        grid_values(i) = i;
    }
    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::linear);

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
        sasktran2::Config::MultipleScatterSource::discrete_ordinates);
    config.set_single_scatter_source(
        sasktran2::Config::SingleScatterSource::discrete_ordinates);

    config.set_num_do_streams(nstr);
    config.set_do_backprop(true);

    // Make the engine
    Sasktran2<1> engine(config, &geo, viewing_geometry);

    sasktran2::OutputIdealDense<1> output;
    engine.calculate_radiance(atmo, output);
}
#endif
