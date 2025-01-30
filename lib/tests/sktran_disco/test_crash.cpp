#include "sasktran2/viewinggeometry.h"
#include <sasktran2.h>

#include <sasktran2/test_helper.h>

#include <sasktran2/testing/do_test_util.h>

TEST_CASE("Homogeneous Precision Issues", "[sktran_do][vector]") {
    // Test case where a crash was observed due to precision issues in the
    // EigenSolver
#ifdef USE_OMP
    omp_set_num_threads(1);
#endif
    // Construct the geometry
    sasktran2::Coordinates coords(0, sasktran2::math::PiOver2, 6372000,
                                  sasktran2::geometrytype::pseudospherical);

    Eigen::VectorXd grid_values(2);
    for (int i = 0; i < 2; ++i) {
        grid_values(i) = i * 1000;
    }
    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::lower);

    sasktran2::Geometry1D geo(std::move(coords), std::move(grid));

    // Construct the Atmosphere
    int nwavel = 1;
    sasktran2::atmosphere::AtmosphereGridStorageFull<3> storage(nwavel,
                                                                geo.size(), 16);
    sasktran2::atmosphere::Surface<3> surface(nwavel);

    sasktran2::atmosphere::Atmosphere<3> atmo(std::move(storage),
                                              std::move(surface), true);

    std::vector<double> extinction{1, 1};

    for (int i = 0; i < nwavel; ++i) {
        atmo.storage().total_extinction(Eigen::all, i) =
            Eigen::Map<Eigen::MatrixXd>(&extinction[0], 2, 1);
    }

    // Set the legendre coefficients based on
    /*
    [2025-01-23 09:11:34.685] [error] Failed to compute the eigensolution for
    layer 21, order 7, ssa 0.999999999 [2025-01-23 09:11:34.685] [error] Layer
    Legendre Coefficients: [2025-01-23 09:11:34.685] [error] 0: a1: 1, a2: 0,
    a3: 0, b1: 0 [2025-01-23 09:11:34.685] [error] 1:
    a1: 2.5655745156039553e-06, a2: 0, a3: 0, b1: 0 [2025-01-23 09:11:34.685]
    [error] 2: a1: 0.47976740901438475, a2: 2.8786036663749366,
    a3: 4.4636822946333296e-06, b1: -1.1751848992583032 [2025-01-23
    09:11:34.685] [error] 3: a1: 5.100310444691835e-07,
    a2: 1.676954749440442e-06, a3: 1.970119437952731e-07, b1:
    -9.07930802490974e-07 [2025-01-23 09:11:34.685] [error] 4:
    a1: 2.4344319757063098e-08, a2: 5.984733780462433e-08,
    a3: 1.0451459967105075e-08, b1: -3.7007234589531765e-08 [2025-01-23
    09:11:34.685] [error] 5: a1: 1.349586446891329e-09,
    a2: 2.7807818658017294e-09, a3: 6.395828400613546e-10, b1:
    -1.8526303537670522e-09 [2025-01-23 09:11:34.685] [error] 6:
    a1: 8.492757754763808e-11, a2: 1.5546314407109339e-10,
    a3: 4.2633227968269347e-11, b1: -1.0867339253613116e-10 [2025-01-23
    09:11:34.685] [error] 7: a1: 5.723724395523332e-12,
    a2: 9.645012157829066e-12, a3: 2.9412201197112983e-12, b1:
    -7.008347795204861e-12
    */
    atmo.storage().leg_coeff.setZero();

    atmo.storage().leg_coeff(0, 0, 0) = 1;
    atmo.storage().leg_coeff(4 * 1, 0, 0) = 2.5655745156039553e-06;
    atmo.storage().leg_coeff(4 * 2, 0, 0) = 0.47976740901438475;
    atmo.storage().leg_coeff(4 * 2 + 1, 0, 0) = 2.8786036663749366;
    atmo.storage().leg_coeff(4 * 2 + 2, 0, 0) = 4.4636822946333296e-06;
    atmo.storage().leg_coeff(4 * 2 + 3, 0, 0) = 1.1751848992583032;
    atmo.storage().leg_coeff(4 * 3, 0, 0) = 5.100310444691835e-07;
    atmo.storage().leg_coeff(4 * 3 + 1, 0, 0) = 1.676954749440442e-06;
    atmo.storage().leg_coeff(4 * 3 + 2, 0, 0) = 1.970119437952731e-07;
    atmo.storage().leg_coeff(4 * 3 + 3, 0, 0) = 9.07930802490974e-07;
    atmo.storage().leg_coeff(4 * 4, 0, 0) = 2.4344319757063098e-08;
    atmo.storage().leg_coeff(4 * 4 + 1, 0, 0) = 5.984733780462433e-08;
    atmo.storage().leg_coeff(4 * 4 + 2, 0, 0) = 1.0451459967105075e-08;
    atmo.storage().leg_coeff(4 * 4 + 3, 0, 0) = 3.7007234589531765e-08;
    atmo.storage().leg_coeff(4 * 5, 0, 0) = 1.349586446891329e-09;
    atmo.storage().leg_coeff(4 * 5 + 1, 0, 0) = 2.7807818658017294e-09;
    atmo.storage().leg_coeff(4 * 5 + 2, 0, 0) = 6.395828400613546e-10;
    atmo.storage().leg_coeff(4 * 5 + 3, 0, 0) = 1.8526303537670522e-09;
    atmo.storage().leg_coeff(4 * 6, 0, 0) = 8.492757754763808e-11;
    atmo.storage().leg_coeff(4 * 6 + 1, 0, 0) = 1.5546314407109339e-10;
    atmo.storage().leg_coeff(4 * 6 + 2, 0, 0) = 4.2633227968269347e-11;
    atmo.storage().leg_coeff(4 * 6 + 3, 0, 0) = 1.0867339253613116e-10;
    atmo.storage().leg_coeff(4 * 7, 0, 0) = 5.723724395523332e-12;
    atmo.storage().leg_coeff(4 * 7 + 1, 0, 0) = 9.645012157829066e-12;
    atmo.storage().leg_coeff(4 * 7 + 2, 0, 0) = 2.9412201197112983e-12;
    atmo.storage().leg_coeff(4 * 7 + 3, 0, 0) = 7.008347795204861e-12;

    atmo.storage().ssa.setConstant(1.0);

    // Construct the Viewing rays
    sasktran2::viewinggeometry::ViewingGeometryContainer viewing_geometry;
    auto& los = viewing_geometry.observer_rays();

    los.emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            0.6, 0.0, 0.6, 200000));

    // Construct the config
    sasktran2::Config config;

    config.set_num_do_streams(8);

    config.set_multiple_scatter_source(
        sasktran2::Config::MultipleScatterSource::discrete_ordinates);

    config.set_single_scatter_source(
        sasktran2::Config::SingleScatterSource::discrete_ordinates);
    config.set_num_do_forced_azimuth(8);

    // Make the engine
    Sasktran2<3> engine(config, &geo, viewing_geometry);

    sasktran2::OutputIdealDense<3> output, output_perturb_above,
        output_perturb_below;
    engine.calculate_radiance(atmo, output);
}
