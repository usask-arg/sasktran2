#include "sasktran2/atmosphere/surface.h"
#include "sasktran2/output.h"
#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("Validate Surface BRDF WF's in Single Scatter/DO",
          "[sasktran2][engine]") {
#ifdef USE_OMP
    omp_set_num_threads(1);
#endif
    // Test BRDFs
    auto brdf = GENERATE(
        std::shared_ptr<sasktran2::atmosphere::brdf::BRDF<1>>(
            std::make_shared<sasktran2::atmosphere::brdf::Lambertian<1>>()),
        std::shared_ptr<sasktran2::atmosphere::brdf::BRDF<1>>(
            std::make_shared<
                sasktran2::atmosphere::brdf::SnowKokhanovsky<1>>()),
        std::shared_ptr<sasktran2::atmosphere::brdf::BRDF<1>>(
            std::make_shared<sasktran2::atmosphere::brdf::MODIS<1>>()));

    // Construct the geometry
    sasktran2::Coordinates coords(0.6, 0, 6371000,
                                  sasktran2::geometrytype::planeparallel);

    Eigen::VectorXd grid_values(101);
    for (int i = 0; i < 101; ++i) {
        grid_values(i) = i * 1000;
    }
    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::linear);

    sasktran2::Geometry1D geo(std::move(coords), std::move(grid));

    // Construct the Atmosphere
    int nwavel = 1;
    sasktran2::atmosphere::AtmosphereGridStorageFull<1> storage(nwavel,
                                                                geo.size(), 16);
    sasktran2::atmosphere::Surface<1> surface(nwavel);

    sasktran2::atmosphere::Atmosphere<1> atmo(std::move(storage),
                                              std::move(surface), true);

    std::vector<double> extinction{
        7.07906113e-05, 6.46250950e-05, 5.86431083e-05, 5.29850715e-05,
        4.77339013e-05, 4.29288557e-05, 3.85773022e-05, 3.46642865e-05,
        3.11600517e-05, 2.80258050e-05, 2.52180748e-05, 2.26734428e-05,
        2.02816648e-05, 1.79778464e-05, 1.57467704e-05, 1.36034281e-05,
        1.15882231e-05, 9.77118267e-06, 8.18898344e-06, 6.84554061e-06,
        5.72584994e-06, 4.80319926e-06, 4.04164882e-06, 3.41027519e-06,
        2.88467502e-06, 2.44547520e-06, 2.07720643e-06, 1.76744635e-06,
        1.50616412e-06, 1.28521627e-06, 1.09795686e-06, 9.38934589e-07,
        8.03656135e-07, 6.88499846e-07, 5.90874787e-07, 5.08080278e-07,
        4.37762355e-07, 3.77949584e-07, 3.26990467e-07, 2.83500936e-07,
        2.46320382e-07, 2.14474894e-07, 1.87146548e-07, 1.63647786e-07,
        1.43400047e-07, 1.25915952e-07, 1.10782493e-07, 9.76450474e-08,
        8.62059194e-08, 7.62163724e-08, 6.74679955e-08, 5.97856954e-08,
        5.30219866e-08, 4.70523226e-08, 4.17712677e-08, 3.70893448e-08,
        3.29295331e-08, 2.92235209e-08, 2.59140539e-08, 2.29536982e-08,
        2.03028465e-08, 1.79281189e-08, 1.58010825e-08, 1.38964913e-08,
        1.21932261e-08, 1.06731597e-08, 9.31932242e-09, 8.11636527e-09,
        7.05027714e-09, 6.10817608e-09, 5.27815905e-09, 4.54919569e-09,
        3.91105464e-09, 3.35204107e-09, 2.86870538e-09, 2.45077836e-09,
        2.09082449e-09, 1.78181885e-09, 1.51726918e-09, 1.29127795e-09,
        1.09856137e-09, 9.35341764e-10, 7.95620813e-10, 6.76822327e-10,
        5.75848867e-10, 4.90034206e-10, 4.17093171e-10, 3.55073944e-10,
        3.02313970e-10, 2.57400129e-10, 2.19133412e-10, 1.86465279e-10,
        1.58426278e-10, 1.34262300e-10, 1.13404218e-10, 9.54097598e-11,
        7.99225004e-11, 6.66436434e-11, 5.53133171e-11, 4.56988383e-11,
        3.75879135e-11};

    atmo.storage().total_extinction =
        Eigen::Map<Eigen::MatrixXd>(&extinction[0], 101, 1);

    atmo.storage().ssa.setConstant(1.0);

    atmo.storage().leg_coeff.chip(0, 0).setConstant(1);
    atmo.storage().leg_coeff.chip(2, 0).setConstant(0.5);

    // Construct the Viewing rays
    sasktran2::viewinggeometry::ViewingGeometryContainer viewing_geometry;
    auto& los = viewing_geometry.observer_rays();

    los.emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            0.6, 60 * EIGEN_PI / 180.0, 0.6, 200000));

    atmo.surface().set_brdf_object(std::move(brdf));

    atmo.surface().brdf_args().setConstant(0.1);

    // Construct the config
    sasktran2::Config config;
    config.set_multiple_scatter_source(
        sasktran2::Config::MultipleScatterSource::discrete_ordinates);
    config.set_num_do_streams(4);

    // Make the engine
    Sasktran2<1> engine(config, &geo, viewing_geometry);

    sasktran2::OutputIdealDense<1> output;
    engine.calculate_radiance(atmo, output);

    sasktran2::OutputIdealDense<1> output_perturbed_above,
        output_perturbed_below;

    double dx = 1e-4;
    for (int i = 0; i < atmo.surface().d_brdf_args().size(); ++i) {
        atmo.surface().brdf_args() += atmo.surface().d_brdf_args()[i] * dx;

        engine.calculate_radiance(atmo, output_perturbed_above);
        atmo.surface().brdf_args() -=
            atmo.surface().d_brdf_args()[i] * (2.0 * dx);
        engine.calculate_radiance(atmo, output_perturbed_below);
        atmo.surface().brdf_args() += atmo.surface().d_brdf_args()[i] * dx;

        // Validate WF
        double analytic_wf = output.radiance().deriv(2 * extinction.size() + i);
        double numeric_wf = (output_perturbed_above.radiance().value(0) -
                             output_perturbed_below.radiance().value(0)) /
                            (2.0 * dx);

        double diff = std::abs(analytic_wf - numeric_wf);

        if (diff > 1e-6) {
            std::cout << "Analytic: " << analytic_wf << std::endl;
            std::cout << "Numeric: " << numeric_wf << std::endl;
        }

        REQUIRE(diff < 1e-6);
    }

    double x = 5;
}
