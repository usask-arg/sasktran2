#include "sktran_disco/twostream/meta.h"
#include <sasktran2/test_helper.h>
#include <sasktran2.h>
#include <sktran_disco/twostream/twostream.h>

#ifdef SKTRAN_CATCH2_VERSION3

void construct_atmo_geometry(
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>>& atmo,
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>>& layer_geo,
    double csz = 0.6) {
#ifdef USE_OMP
    omp_set_num_threads(1);
#endif
    // Construct the geometry
    sasktran2::Coordinates coords(csz, 0, 6371000,
                                  sasktran2::geometrytype::pseudospherical);

    int nlyr = 20;
    int nwavel = 1;
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

    atmo = std::make_unique<sasktran2::atmosphere::Atmosphere<1>>(
        std::move(storage), std::move(surface), true);

    atmo->storage().total_extinction.setConstant(0.01);
    atmo->storage()
        .total_extinction(Eigen::seq(10, 13), Eigen::placeholders::all)
        .setConstant(0.02);

    atmo->storage().ssa.setConstant(0.8);

    atmo->storage().leg_coeff.chip(0, 0).setConstant(1);
    atmo->storage().leg_coeff.chip(1, 0).setConstant(0.5);

    sasktran2::Config config;
    sasktran_disco::PersistentConfiguration<1> pconfig;
    std::vector<sasktran2::raytracing::TracedRay> traced_rays;
    sasktran_disco::SKTRAN_DO_UserSpec specs;
    pconfig.configure(specs, config, csz, nlyr, traced_rays);

    layer_geo =
        std::make_unique<sasktran_disco::GeometryLayerArray<1>>(pconfig, geo);
}

void construct_atmo_geometry_linear_lowlevel(
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>>& atmo,
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>>& layer_geo,
    double csz = 0.6) {
#ifdef USE_OMP
    omp_set_num_threads(1);
#endif
    // Match lowlevel two-stream test geometry as closely as possible.
    sasktran2::Coordinates coords(csz, 0, 6372000,
                                  sasktran2::geometrytype::planeparallel);

    int nlyr = 13;
    int nwavel = 1;

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

    atmo = std::make_unique<sasktran2::atmosphere::Atmosphere<1>>(
        std::move(storage), std::move(surface), true);

    for (int i = 0; i < geo.size(); ++i) {
        atmo->storage().total_extinction(i, 0) =
            0.02 * std::exp(-0.25 * i) + 1e-12;
    }
    atmo->storage().ssa.setConstant(0.8);
    atmo->storage().leg_coeff.chip(0, 0).setConstant(1);
    atmo->storage().leg_coeff.chip(1, 0).setConstant(0.0);
    atmo->storage().leg_coeff.chip(2, 0).setConstant(0.5);

    sasktran2::Config config;
    sasktran_disco::PersistentConfiguration<1> pconfig;
    std::vector<sasktran2::raytracing::TracedRay> traced_rays;
    sasktran_disco::SKTRAN_DO_UserSpec specs;
    pconfig.configure(specs, config, csz, nlyr, traced_rays);

    layer_geo =
        std::make_unique<sasktran_disco::GeometryLayerArray<1>>(pconfig, geo);
}

TEST_CASE("Backprop_OpticalDepth", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * od
    Eigen::RowVectorXd weights = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.od(i) += eps;

        numerical_grad(i) =
            (weights.dot(input_p.od) - weights.dot(input.od)) / eps;
    }

    sasktran2::twostream::backprop::od(input, weights, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_extinction(i) - numerical_grad(i)) < 1e-6);
    }
}

TEST_CASE("Backprop_SSA", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * ssa
    Eigen::RowVectorXd weights = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.ssa(i) += eps;

        numerical_grad(i) =
            (weights.dot(input_p.ssa) - weights.dot(input.ssa)) / eps;
    }

    sasktran2::twostream::backprop::ssa(input, weights, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_ssa(i) - numerical_grad(i)) < 1e-6);
    }
}

TEST_CASE("Backprop_Homog", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::RowVectorXd, 2> weights;
    weights[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Xp;
    weights_Xp[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Xp[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Xm;
    weights_Xm[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Xm[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_omega;
    weights_omega[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_omega[1] = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve_layers(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.ssa(i) += eps;
        sasktran2::twostream::solve_layers(input_p, soln_p);

        numerical_grad(i) = (weights[0].dot(soln_p.homog[0].k.value) -
                             weights[0].dot(soln.homog[0].k.value)) /
                            eps;
        numerical_grad(i) += (weights[1].dot(soln_p.homog[1].k.value) -
                              weights[1].dot(soln.homog[1].k.value)) /
                             eps;

        numerical_grad(i) += (weights_Xp[0].dot(soln_p.homog[0].X_plus.value) -
                              weights_Xp[0].dot(soln.homog[0].X_plus.value)) /
                             eps;
        numerical_grad(i) += (weights_Xp[1].dot(soln_p.homog[1].X_plus.value) -
                              weights_Xp[1].dot(soln.homog[1].X_plus.value)) /
                             eps;

        numerical_grad(i) += (weights_Xm[0].dot(soln_p.homog[0].X_minus.value) -
                              weights_Xm[0].dot(soln.homog[0].X_minus.value)) /
                             eps;
        numerical_grad(i) += (weights_Xm[1].dot(soln_p.homog[1].X_minus.value) -
                              weights_Xm[1].dot(soln.homog[1].X_minus.value)) /
                             eps;

        numerical_grad(i) +=
            (weights_omega[0].dot(soln_p.homog[0].omega.value) -
             weights_omega[0].dot(soln.homog[0].omega.value)) /
            eps;
        numerical_grad(i) +=
            (weights_omega[1].dot(soln_p.homog[1].omega.value) -
             weights_omega[1].dot(soln.homog[1].omega.value)) /
            eps;
    }

    sasktran2::twostream::backprop::homog_k(input, soln.homog, weights,
                                            grad_map);
    sasktran2::twostream::backprop::homog_X(input, soln.homog, weights_Xp,
                                            weights_Xm, grad_map);
    sasktran2::twostream::backprop::homog_omega(input, soln.homog,
                                                weights_omega, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_ssa(i) - numerical_grad(i)) < 1e-4);
    }
}

TEST_CASE("Backprop_Particular", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::RowVectorXd, 2> weights_Gplus_top;
    weights_Gplus_top[0] = Eigen::RowVectorXd::Random(natmo - 1) +
                           Eigen::RowVectorXd::Ones(natmo - 1);
    weights_Gplus_top[1] = Eigen::RowVectorXd::Random(natmo - 1) +
                           Eigen::RowVectorXd::Ones(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gplus_bottom;
    weights_Gplus_bottom[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gplus_bottom[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gminus_top;
    weights_Gminus_top[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gminus_top[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gminus_bottom;
    weights_Gminus_bottom[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gminus_bottom[1] = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve_layers(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.ssa(i) += eps;
        input_p.csz = 0.6;
        sasktran2::twostream::solve_layers(input_p, soln_p);

        numerical_grad(i) =
            (weights_Gplus_top[0].dot(soln_p.particular[0].G_plus_top.value) -
             weights_Gplus_top[0].dot(soln.particular[0].G_plus_top.value)) /
            (eps);
        numerical_grad(i) +=
            (weights_Gplus_top[1].dot(soln_p.particular[1].G_plus_top.value) -
             weights_Gplus_top[1].dot(soln.particular[1].G_plus_top.value)) /
            (eps);

        numerical_grad(i) += (weights_Gplus_bottom[0].dot(
                                  soln_p.particular[0].G_plus_bottom.value) -
                              weights_Gplus_bottom[0].dot(
                                  soln.particular[0].G_plus_bottom.value)) /
                             eps;
        numerical_grad(i) += (weights_Gplus_bottom[1].dot(
                                  soln_p.particular[1].G_plus_bottom.value) -
                              weights_Gplus_bottom[1].dot(
                                  soln.particular[1].G_plus_bottom.value)) /
                             eps;

        numerical_grad(i) +=
            (weights_Gminus_top[0].dot(soln_p.particular[0].G_minus_top.value) -
             weights_Gminus_top[0].dot(soln.particular[0].G_minus_top.value)) /
            eps;
        numerical_grad(i) +=
            (weights_Gminus_top[1].dot(soln_p.particular[1].G_minus_top.value) -
             weights_Gminus_top[1].dot(soln.particular[1].G_minus_top.value)) /
            eps;

        numerical_grad(i) += (weights_Gminus_bottom[0].dot(
                                  soln_p.particular[0].G_minus_bottom.value) -
                              weights_Gminus_bottom[0].dot(
                                  soln.particular[0].G_minus_bottom.value)) /
                             eps;
        numerical_grad(i) += (weights_Gminus_bottom[1].dot(
                                  soln_p.particular[1].G_minus_bottom.value) -
                              weights_Gminus_bottom[1].dot(
                                  soln.particular[1].G_minus_bottom.value)) /
                             eps;
    }

    sasktran2::twostream::backprop::particular_G_plus_top(
        input, soln.particular, weights_Gplus_top, grad_map);
    sasktran2::twostream::backprop::particular_G_plus_bottom(
        input, soln.particular, weights_Gplus_bottom, grad_map);
    sasktran2::twostream::backprop::particular_G_minus_top(
        input, soln.particular, weights_Gminus_top, grad_map);
    sasktran2::twostream::backprop::particular_G_minus_bottom(
        input, soln.particular, weights_Gminus_bottom, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_ssa(i) - numerical_grad(i)) < 1e-4);
    }
}

TEST_CASE("Backprop_BVP", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::MatrixXd, 2> weights;
    weights[0].resize(2 * (natmo - 1), 1);
    weights[0].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    weights[1].resize(2 * (natmo - 1), 1);
    weights[1].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.ssa(i) += eps;
        input_p.csz = 0.6;
        sasktran2::twostream::solve_layers(input_p, soln_p);
        sasktran2::twostream::solve_bvp(input_p, soln_p);

        numerical_grad(i) =
            (weights[0].col(0).dot(soln_p.bvp_coeffs[0].rhs.col(0)) -
             weights[0].col(0).dot(soln.bvp_coeffs[0].rhs.col(0))) /
            eps;
        numerical_grad(i) +=
            (weights[1].col(0).dot(soln_p.bvp_coeffs[1].rhs.col(0)) -
             weights[1].col(0).dot(soln.bvp_coeffs[1].rhs.col(0))) /
            eps;
    }

    sasktran2::twostream::backprop::bvp(input, soln.bvp_coeffs, soln.homog,
                                        soln.particular, weights, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_ssa(i) - numerical_grad(i)) < 1e-4);
    }
}

TEST_CASE("Backprop_Full_SSA", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    double viewing_zenith = 0.6;
    double azimuth = 0.0;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    std::array<Eigen::MatrixXd, 2> bvp_storage;
    bvp_storage[0].resize(2 * (natmo - 1), 1);
    bvp_storage[1].resize(2 * (natmo - 1), 1);

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    Eigen::RowVectorXd weights;
    weights = Eigen::RowVectorXd::Random((natmo - 1));

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    sasktran2::twostream::Sources<sasktran2::twostream::SourceType::ONLY_SOLAR>
        sources;
    soln.init(natmo - 1);
    sources.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);
    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            sources_p;
        sources_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.ssa(i) += eps;
        input_p.csz = 0.6;
        sasktran2::twostream::solve_layers(input_p, soln_p);
        sasktran2::twostream::solve_bvp(input_p, soln_p);
        sasktran2::twostream::post_process(input_p, viewing_zenith, azimuth,
                                           soln_p, sources_p);

        numerical_grad(i) = (weights.dot(sources_p.source.value) -
                             weights.dot(sources.source.value)) /
                            eps;
    }

    sasktran2::twostream::backprop::full(input, soln, sources, weights,
                                         bvp_storage, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_ssa(i) - numerical_grad(i)) < 1e-6);
    }
}

TEST_CASE("Backprop_Transmission", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * od
    Eigen::RowVectorXd weights = Eigen::RowVectorXd::Random(natmo);

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-8;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        input_p.calculate_base(0);
        input_p.od(i) += eps;
        input_p.calculate_derived(0);

        numerical_grad(i) = (weights.dot(input_p.transmission) -
                             weights.dot(input.transmission)) /
                            eps;
    }

    sasktran2::twostream::backprop::transmission(input, weights, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_extinction(i) - numerical_grad(i)) < 1e-6);
    }
}

TEST_CASE("Backprop_Secant", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo, 0.01);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * od
    Eigen::RowVectorXd weights = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-8;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        input_p.calculate_base(0);
        input_p.od(i) += eps;
        input_p.calculate_derived(0);
        input_p.od(i) -= eps;

        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_pm;
        input_pm.atmosphere = atmo.get();
        input_pm.geometry_layers = layer_geo.get();
        input_pm.init(natmo - 1);

        input_pm.calculate_base(0);
        input_pm.od(i) -= eps;
        input_pm.calculate_derived(0);
        input_pm.od(i) += eps; // reset

        numerical_grad(i) = (weights.dot(input_p.average_secant) -
                             weights.dot(input_pm.average_secant)) /
                            (2.0 * eps);
    }

    sasktran2::twostream::backprop::secant(input, weights, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_extinction(i) - numerical_grad(i)) < 1e-4);
    }
}

TEST_CASE("Backprop_Particular_ext", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo, 0.01);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::RowVectorXd, 2> weights_Gplus_top;
    weights_Gplus_top[0] = Eigen::RowVectorXd::Ones(natmo - 1);
    weights_Gplus_top[1] = Eigen::RowVectorXd::Ones(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gplus_bottom;
    weights_Gplus_bottom[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gplus_bottom[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gminus_top;
    weights_Gminus_top[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gminus_top[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gminus_bottom;
    weights_Gminus_bottom[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gminus_bottom[1] = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.01;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve_layers(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-7;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate_base(0);
        input_p.od(i) += eps;
        input_p.calculate_derived(0);
        input_p.csz = 0.01;

        sasktran2::twostream::solve_layers(input_p, soln_p);

        numerical_grad(i) =
            (weights_Gplus_top[0].dot(soln_p.particular[0].G_plus_top.value) -
             weights_Gplus_top[0].dot(soln.particular[0].G_plus_top.value)) /
            (eps);

        numerical_grad(i) +=
            (weights_Gplus_top[1].dot(soln_p.particular[1].G_plus_top.value) -
             weights_Gplus_top[1].dot(soln.particular[1].G_plus_top.value)) /
            (eps);
        numerical_grad(i) += (weights_Gplus_bottom[0].dot(
                                  soln_p.particular[0].G_plus_bottom.value) -
                              weights_Gplus_bottom[0].dot(
                                  soln.particular[0].G_plus_bottom.value)) /
                             eps;
        numerical_grad(i) += (weights_Gplus_bottom[1].dot(
                                  soln_p.particular[1].G_plus_bottom.value) -
                              weights_Gplus_bottom[1].dot(
                                  soln.particular[1].G_plus_bottom.value)) /
                             eps;

        numerical_grad(i) +=
            (weights_Gminus_top[0].dot(soln_p.particular[0].G_minus_top.value) -
             weights_Gminus_top[0].dot(soln.particular[0].G_minus_top.value)) /
            eps;
        numerical_grad(i) +=
            (weights_Gminus_top[1].dot(soln_p.particular[1].G_minus_top.value) -
             weights_Gminus_top[1].dot(soln.particular[1].G_minus_top.value)) /
            eps;

        numerical_grad(i) += (weights_Gminus_bottom[0].dot(
                                  soln_p.particular[0].G_minus_bottom.value) -
                              weights_Gminus_bottom[0].dot(
                                  soln.particular[0].G_minus_bottom.value)) /
                             eps;
        numerical_grad(i) += (weights_Gminus_bottom[1].dot(
                                  soln_p.particular[1].G_minus_bottom.value) -
                              weights_Gminus_bottom[1].dot(
                                  soln.particular[1].G_minus_bottom.value)) /
                             eps;
    }

    sasktran2::twostream::backprop::particular_G_plus_top(
        input, soln.particular, weights_Gplus_top, grad_map);
    sasktran2::twostream::backprop::particular_G_plus_bottom(
        input, soln.particular, weights_Gplus_bottom, grad_map);
    sasktran2::twostream::backprop::particular_G_minus_top(
        input, soln.particular, weights_Gminus_top, grad_map);
    sasktran2::twostream::backprop::particular_G_minus_bottom(
        input, soln.particular, weights_Gminus_bottom, grad_map);
    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_extinction(i) - numerical_grad(i)) < 1e-6);
    }
}

TEST_CASE("Backprop_BVP_ext", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());
    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::MatrixXd, 2> weights;
    weights[0].resize(2 * (natmo - 1), 1);
    weights[0].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    weights[1].resize(2 * (natmo - 1), 1);
    weights[1].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate_base(0);
        input_p.od(i) += eps;
        input_p.calculate_derived(0);
        input_p.csz = 0.6;
        sasktran2::twostream::solve_layers(input_p, soln_p);
        sasktran2::twostream::solve_bvp(input_p, soln_p);

        numerical_grad(i) =
            (weights[0].col(0).dot(soln_p.bvp_coeffs[0].rhs.col(0)) -
             weights[0].col(0).dot(soln.bvp_coeffs[0].rhs.col(0))) /
            eps;
        numerical_grad(i) +=
            (weights[1].col(0).dot(soln_p.bvp_coeffs[1].rhs.col(0)) -
             weights[1].col(0).dot(soln.bvp_coeffs[1].rhs.col(0))) /
            eps;
    }

    sasktran2::twostream::backprop::bvp(input, soln.bvp_coeffs, soln.homog,
                                        soln.particular, weights, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_extinction(i) - numerical_grad(i)) < 1e-4);
    }
}

TEST_CASE("Backprop_Full_ext", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    double viewing_zenith = 0.6;
    double azimuth = 0.0;

    construct_atmo_geometry(atmo, layer_geo, 0.01);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    std::array<Eigen::MatrixXd, 2> bvp_storage;
    bvp_storage[0].resize(2 * (natmo - 1), 1);
    bvp_storage[1].resize(2 * (natmo - 1), 1);

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    Eigen::RowVectorXd weights;
    weights = Eigen::RowVectorXd::Random((natmo - 1));

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.01;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    sasktran2::twostream::Sources<sasktran2::twostream::SourceType::ONLY_SOLAR>
        sources;
    soln.init(natmo - 1);
    sources.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);
    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            sources_p;
        sources_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate_base(0);
        input_p.od(i) += eps;
        input_p.calculate_derived(0);
        input_p.csz = 0.01;

        sasktran2::twostream::solve_layers(input_p, soln_p);
        sasktran2::twostream::solve_bvp(input_p, soln_p);
        sasktran2::twostream::post_process(input_p, viewing_zenith, azimuth,
                                           soln_p, sources_p);

        numerical_grad(i) = (weights.dot(sources_p.source.value) -
                             weights.dot(sources.source.value)) /
                            eps;
    }

    sasktran2::twostream::backprop::full(input, soln, sources, weights,
                                         bvp_storage, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_extinction(i) - numerical_grad(i)) < 1e-6);
    }
}

TEST_CASE("Backprop_Full_b1", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    double viewing_zenith = 0.6;
    double azimuth = 0.0;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(natmo * 3);
    grad.setZero();

    std::array<Eigen::MatrixXd, 2> bvp_storage;
    bvp_storage[0].resize(2 * (natmo - 1), 1);
    bvp_storage[1].resize(2 * (natmo - 1), 1);

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    Eigen::RowVectorXd weights;
    weights = Eigen::RowVectorXd::Random((natmo - 1));

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    sasktran2::twostream::Sources<sasktran2::twostream::SourceType::ONLY_SOLAR>
        sources;
    soln.init(natmo - 1);
    sources.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);
    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            sources_p;
        sources_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.b1(i) += eps;
        input_p.csz = 0.6;
        sasktran2::twostream::solve_layers(input_p, soln_p);
        sasktran2::twostream::solve_bvp(input_p, soln_p);
        sasktran2::twostream::post_process(input_p, viewing_zenith, azimuth,
                                           soln_p, sources_p);

        numerical_grad(i) = (weights.dot(sources_p.source.value) -
                             weights.dot(sources.source.value)) /
                            eps;
    }

    sasktran2::twostream::backprop::full(input, soln, sources, weights,
                                         bvp_storage, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_b1(i) - numerical_grad(i)) < 1e-6);
    }
}

TEST_CASE("Backprop_Particular_b1", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(natmo * 3);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::RowVectorXd, 2> weights_Gplus_top;
    weights_Gplus_top[0] = Eigen::RowVectorXd::Random(natmo - 1) +
                           Eigen::RowVectorXd::Ones(natmo - 1);
    weights_Gplus_top[1] = Eigen::RowVectorXd::Random(natmo - 1) +
                           Eigen::RowVectorXd::Ones(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gplus_bottom;
    weights_Gplus_bottom[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gplus_bottom[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gminus_top;
    weights_Gminus_top[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gminus_top[1] = Eigen::RowVectorXd::Random(natmo - 1);

    std::array<Eigen::RowVectorXd, 2> weights_Gminus_bottom;
    weights_Gminus_bottom[0] = Eigen::RowVectorXd::Random(natmo - 1);
    weights_Gminus_bottom[1] = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve_layers(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.b1(i) += eps;
        input_p.csz = 0.6;
        sasktran2::twostream::solve_layers(input_p, soln_p);

        numerical_grad(i) =
            (weights_Gplus_top[0].dot(soln_p.particular[0].G_plus_top.value) -
             weights_Gplus_top[0].dot(soln.particular[0].G_plus_top.value)) /
            (eps);
        numerical_grad(i) +=
            (weights_Gplus_top[1].dot(soln_p.particular[1].G_plus_top.value) -
             weights_Gplus_top[1].dot(soln.particular[1].G_plus_top.value)) /
            (eps);

        numerical_grad(i) += (weights_Gplus_bottom[0].dot(
                                  soln_p.particular[0].G_plus_bottom.value) -
                              weights_Gplus_bottom[0].dot(
                                  soln.particular[0].G_plus_bottom.value)) /
                             eps;
        numerical_grad(i) += (weights_Gplus_bottom[1].dot(
                                  soln_p.particular[1].G_plus_bottom.value) -
                              weights_Gplus_bottom[1].dot(
                                  soln.particular[1].G_plus_bottom.value)) /
                             eps;

        numerical_grad(i) +=
            (weights_Gminus_top[0].dot(soln_p.particular[0].G_minus_top.value) -
             weights_Gminus_top[0].dot(soln.particular[0].G_minus_top.value)) /
            eps;
        numerical_grad(i) +=
            (weights_Gminus_top[1].dot(soln_p.particular[1].G_minus_top.value) -
             weights_Gminus_top[1].dot(soln.particular[1].G_minus_top.value)) /
            eps;

        numerical_grad(i) += (weights_Gminus_bottom[0].dot(
                                  soln_p.particular[0].G_minus_bottom.value) -
                              weights_Gminus_bottom[0].dot(
                                  soln.particular[0].G_minus_bottom.value)) /
                             eps;
        numerical_grad(i) += (weights_Gminus_bottom[1].dot(
                                  soln_p.particular[1].G_minus_bottom.value) -
                              weights_Gminus_bottom[1].dot(
                                  soln.particular[1].G_minus_bottom.value)) /
                             eps;
    }

    sasktran2::twostream::backprop::particular_G_plus_top(
        input, soln.particular, weights_Gplus_top, grad_map);
    sasktran2::twostream::backprop::particular_G_plus_bottom(
        input, soln.particular, weights_Gplus_bottom, grad_map);
    sasktran2::twostream::backprop::particular_G_minus_top(
        input, soln.particular, weights_Gminus_top, grad_map);
    sasktran2::twostream::backprop::particular_G_minus_bottom(
        input, soln.particular, weights_Gminus_bottom, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_b1(i) - numerical_grad(i)) < 1e-4);
    }
}

TEST_CASE("Backprop_BVP_b1", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(natmo * 3);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::MatrixXd, 2> weights;
    weights[0].resize(2 * (natmo - 1), 1);
    weights[0].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    weights[1].resize(2 * (natmo - 1), 1);
    weights[1].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate(0);
        input_p.b1(i) += eps;
        input_p.csz = 0.6;
        sasktran2::twostream::solve_layers(input_p, soln_p);
        sasktran2::twostream::solve_bvp(input_p, soln_p);

        numerical_grad(i) =
            (weights[0].col(0).dot(soln_p.bvp_coeffs[0].rhs.col(0)) -
             weights[0].col(0).dot(soln.bvp_coeffs[0].rhs.col(0))) /
            eps;
        numerical_grad(i) +=
            (weights[1].col(0).dot(soln_p.bvp_coeffs[1].rhs.col(0)) -
             weights[1].col(0).dot(soln.bvp_coeffs[1].rhs.col(0))) /
            eps;
    }

    sasktran2::twostream::backprop::bvp(input, soln.bvp_coeffs, soln.homog,
                                        soln.particular, weights, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_b1(i) - numerical_grad(i)) < 1e-4);
    }
}

TEST_CASE("Backprop_Full_ext_midcsz", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    double viewing_zenith = 0.8;
    double azimuth = 0.0;

    construct_atmo_geometry(atmo, layer_geo, 0.6);
    atmo->storage().leg_coeff.chip(1, 0).setConstant(0.0);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    std::array<Eigen::MatrixXd, 2> bvp_storage;
    bvp_storage[0].resize(2 * (natmo - 1), 1);
    bvp_storage[1].resize(2 * (natmo - 1), 1);

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    Eigen::RowVectorXd weights;
    weights = Eigen::RowVectorXd::Random((natmo - 1));

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);
    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    sasktran2::twostream::Sources<sasktran2::twostream::SourceType::ONLY_SOLAR>
        sources;
    soln.init(natmo - 1);
    sources.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);
    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            sources_p;
        sources_p.init(natmo - 1);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        soln_p.init(natmo - 1);

        input_p.calculate_base(0);
        input_p.od(i) += eps;
        input_p.calculate_derived(0);
        input_p.csz = 0.6;

        sasktran2::twostream::solve_layers(input_p, soln_p);
        sasktran2::twostream::solve_bvp(input_p, soln_p);
        sasktran2::twostream::post_process(input_p, viewing_zenith, azimuth,
                                           soln_p, sources_p);

        numerical_grad(i) = (weights.dot(sources_p.source.value) -
                             weights.dot(sources.source.value)) /
                            eps;
    }

    sasktran2::twostream::backprop::full(input, soln, sources, weights,
                                         bvp_storage, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        CAPTURE(i, grad_map.d_extinction(i), numerical_grad(i));
        REQUIRE(abs(grad_map.d_extinction(i) - numerical_grad(i)) < 1e-5);
    }
}

TEST_CASE("Backprop_MapToAtmosphere_ext", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo, 0.6);

    int natmo = atmo->storage().total_extinction.rows();
    int nlyr = natmo - 1;

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.init(nlyr);
    input.calculate(0);

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();
    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    Eigen::RowVectorXd w = Eigen::RowVectorXd::Random(nlyr);
    sasktran2::twostream::backprop::od(input, w, grad_map);

    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> source(
        1, atmo->num_deriv(), true);
    sasktran2::twostream::backprop::map_to_atmosphere(input, grad_map, source);

    auto base_obj = [&]() { return w.dot(input.od); };
    const double f0 = base_obj();

    const double eps = 1e-6;
    for (int i = 0; i < natmo; ++i) {
        const double oldk = atmo->storage().total_extinction(i, 0);

        atmo->storage().total_extinction(i, 0) = oldk + eps;
        input.calculate(0);
        const double fp = base_obj();

        atmo->storage().total_extinction(i, 0) = oldk - eps;
        input.calculate(0);
        const double fm = base_obj();

        atmo->storage().total_extinction(i, 0) = oldk;
        input.calculate(0);

        const double numerical = (fp - fm) / (2 * eps);
        const double analytic = source.deriv(i);

        CAPTURE(i, analytic, numerical);
        REQUIRE(abs(analytic - numerical) < 1e-6);
    }

    input.calculate(0);
    REQUIRE(abs(base_obj() - f0) < 1e-12);
}

TEST_CASE("Backprop_StartOfRay_ext_midcsz", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    const double viewing_zenith = 0.8;
    const double azimuth = 0.0;

    construct_atmo_geometry(atmo, layer_geo, 0.6);
    atmo->storage().leg_coeff.chip(1, 0).setConstant(0.0);

    const int natmo = atmo->storage().total_extinction.rows();
    const int nlyr = natmo - 1;

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;
    input.init(nlyr);
    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    sasktran2::twostream::Sources<sasktran2::twostream::SourceType::ONLY_SOLAR>
        sources;
    soln.init(nlyr);
    sources.init(nlyr);
    sasktran2::twostream::solve(input, soln);

    Eigen::MatrixXd atten(nlyr, nlyr);
    atten.setZero();
    const double viewing_secant = 1.0 / viewing_zenith;
    for (int j = 1; j < nlyr; ++j) {
        atten(j, Eigen::seq(0, j - 1)).setConstant(viewing_secant);
    }

    sources.final_weight_factors.noalias() = (-atten * input.od);
    sources.final_weight_factors = sources.final_weight_factors.array().exp();

    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    auto objective = [&](const auto& in, const auto& so, const auto& src) {
        const double ground =
            (so.particular[0].G_plus_bottom.value(Eigen::placeholders::last) +
             so.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                 so.homog[0].X_plus.value(Eigen::placeholders::last) *
                 so.homog[0].omega.value(Eigen::placeholders::last) +
             so.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
                 so.homog[0].X_minus.value(Eigen::placeholders::last)) *
            2 * in.mu * in.albedo *
            src.final_weight_factors(Eigen::placeholders::last) *
            src.beamtrans.value(Eigen::placeholders::last);

        return ground + src.source.value.dot(src.final_weight_factors);
    };

    const double ground_weight =
        2 * input.mu * input.albedo *
        sources.final_weight_factors(Eigen::placeholders::last) *
        sources.beamtrans.value(Eigen::placeholders::last);

    const double ground_source =
        (soln.particular[0].G_plus_bottom.value(Eigen::placeholders::last) +
         soln.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
             soln.homog[0].X_plus.value(Eigen::placeholders::last) *
             soln.homog[0].omega.value(Eigen::placeholders::last) +
         soln.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
             soln.homog[0].X_minus.value(Eigen::placeholders::last)) *
        ground_weight;

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();
    std::array<Eigen::MatrixXd, 2> bvp_storage;
    bvp_storage[0].resize(2 * nlyr, 1);
    bvp_storage[1].resize(2 * nlyr, 1);

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    sasktran2::twostream::backprop::full(input, soln, sources,
                                         sources.final_weight_factors,
                                         bvp_storage, grad_map, ground_weight);

    for (int i = 0; i < nlyr; ++i) {
        sasktran2::twostream::backprop::od(
            input,
            (-sources.final_weight_factors(i) * sources.source.value(i) *
             atten.row(i)),
            grad_map);
    }

    sasktran2::twostream::backprop::od(
        input,
        (-ground_source * Eigen::RowVectorXd::Ones(nlyr) / viewing_zenith),
        grad_map);

    const double eps = 1e-6;
    for (int i = 0; i < nlyr; ++i) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.csz = 0.6;
        input_p.init(nlyr);
        input_p.calculate(0);

        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_m;
        input_m.atmosphere = atmo.get();
        input_m.geometry_layers = layer_geo.get();
        input_m.csz = 0.6;
        input_m.init(nlyr);
        input_m.calculate(0);

        input_p.od(i) += eps;
        input_p.calculate_derived(0);
        input_m.od(i) -= eps;
        input_m.calculate_derived(0);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_m;
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            sources_p;
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            sources_m;

        soln_p.init(nlyr);
        soln_m.init(nlyr);
        sources_p.init(nlyr);
        sources_m.init(nlyr);

        sasktran2::twostream::solve(input_p, soln_p);
        sasktran2::twostream::solve(input_m, soln_m);

        sources_p.final_weight_factors.noalias() = (-atten * input_p.od);
        sources_p.final_weight_factors =
            sources_p.final_weight_factors.array().exp();
        sources_m.final_weight_factors.noalias() = (-atten * input_m.od);
        sources_m.final_weight_factors =
            sources_m.final_weight_factors.array().exp();

        sasktran2::twostream::post_process(input_p, viewing_zenith, azimuth,
                                           soln_p, sources_p);
        sasktran2::twostream::post_process(input_m, viewing_zenith, azimuth,
                                           soln_m, sources_m);

        const double fp = objective(input_p, soln_p, sources_p);
        const double fm = objective(input_m, soln_m, sources_m);
        const double numerical = (fp - fm) / (2 * eps);
        const double analytic = grad_map.d_extinction(i);

        CAPTURE(i, analytic, numerical);
        REQUIRE(abs(analytic - numerical) < 1e-5);
    }
}

TEST_CASE("Backprop_MapToAtmosphere_ext_linear", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry_linear_lowlevel(atmo, layer_geo, 0.6);

    int natmo = atmo->storage().total_extinction.rows();
    int nlyr = natmo - 1;

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.init(nlyr);
    input.calculate(0);

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();
    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    Eigen::RowVectorXd w = Eigen::RowVectorXd::Random(nlyr);
    sasktran2::twostream::backprop::od(input, w, grad_map);

    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> source(
        1, atmo->num_deriv(), true);
    sasktran2::twostream::backprop::map_to_atmosphere(input, grad_map, source);

    auto base_obj = [&]() { return w.dot(input.od); };

    const double eps = 1e-6;
    for (int i = 0; i < natmo; ++i) {
        const double oldk = atmo->storage().total_extinction(i, 0);

        atmo->storage().total_extinction(i, 0) = oldk + eps;
        input.calculate(0);
        const double fp = base_obj();

        atmo->storage().total_extinction(i, 0) = oldk - eps;
        input.calculate(0);
        const double fm = base_obj();

        atmo->storage().total_extinction(i, 0) = oldk;
        input.calculate(0);

        const double numerical = (fp - fm) / (2 * eps);
        const double analytic = source.deriv(i);

        CAPTURE(i, analytic, numerical);
        REQUIRE(abs(analytic - numerical) < 1e-5);
    }
}

TEST_CASE("Backprop_MapToAtmosphere_fullk_linear", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry_linear_lowlevel(atmo, layer_geo, 0.6);

    int natmo = atmo->storage().total_extinction.rows();
    int nlyr = natmo - 1;

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.init(nlyr);
    input.calculate(0);

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();
    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    Eigen::RowVectorXd w_od = Eigen::RowVectorXd::Random(nlyr);
    Eigen::RowVectorXd w_ssa = Eigen::RowVectorXd::Random(nlyr);
    Eigen::RowVectorXd w_b1 = Eigen::RowVectorXd::Random(nlyr);

    grad_map.d_extinction = w_od;
    grad_map.d_ssa = w_ssa;
    grad_map.d_b1 = w_b1;

    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> source(
        1, atmo->num_deriv(), true);
    sasktran2::twostream::backprop::map_to_atmosphere(input, grad_map, source);

    auto base_obj = [&]() {
        return w_od.dot(input.od) + w_ssa.dot(input.ssa) + w_b1.dot(input.b1);
    };

    const double eps = 1e-6;
    for (int i = 0; i < natmo; ++i) {
        const double oldk = atmo->storage().total_extinction(i, 0);

        atmo->storage().total_extinction(i, 0) = oldk + eps;
        input.calculate(0);
        const double fp = base_obj();

        atmo->storage().total_extinction(i, 0) = oldk - eps;
        input.calculate(0);
        const double fm = base_obj();

        atmo->storage().total_extinction(i, 0) = oldk;
        input.calculate(0);

        const double numerical = (fp - fm) / (2 * eps);
        const double analytic = source.deriv(i);

        CAPTURE(i, analytic, numerical);
        REQUIRE(abs(analytic - numerical) < 1e-5);
    }
}

TEST_CASE("Backprop_StartOfRay_ext_linear_ground", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    const double viewing_zenith = 0.8;
    const double azimuth = 0.0;

    construct_atmo_geometry_linear_lowlevel(atmo, layer_geo, 0.6);
    atmo->storage().leg_coeff.chip(1, 0).setConstant(0.0);
    atmo->surface().brdf_args().setConstant(0.2);
    for (int i = 0; i < atmo->storage().total_extinction.rows(); ++i) {
        const double z = i * 5000.0;
        atmo->storage().total_extinction(i, 0) = 2.5e-8 * std::exp(-z / 7000.0);
    }
    atmo->storage().ssa.setConstant(1.0);

    const int natmo = atmo->storage().total_extinction.rows();
    const int nlyr = natmo - 1;

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;
    input.init(nlyr);
    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    sasktran2::twostream::Sources<sasktran2::twostream::SourceType::ONLY_SOLAR>
        sources;
    soln.init(nlyr);
    sources.init(nlyr);
    sasktran2::twostream::solve(input, soln);

    Eigen::MatrixXd atten(nlyr, nlyr);
    atten.setZero();
    const double viewing_secant = 1.0 / viewing_zenith;
    for (int j = 1; j < nlyr; ++j) {
        atten(j, Eigen::seq(0, j - 1)).setConstant(viewing_secant);
    }

    sources.final_weight_factors.noalias() = (-atten * input.od);
    sources.final_weight_factors = sources.final_weight_factors.array().exp();

    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    auto objective = [&](const auto& in, const auto& so, const auto& src) {
        const double ground =
            (so.particular[0].G_plus_bottom.value(Eigen::placeholders::last) +
             so.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                 so.homog[0].X_plus.value(Eigen::placeholders::last) *
                 so.homog[0].omega.value(Eigen::placeholders::last) +
             so.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
                 so.homog[0].X_minus.value(Eigen::placeholders::last)) *
            2 * in.mu * in.albedo *
            src.final_weight_factors(Eigen::placeholders::last) *
            src.beamtrans.value(Eigen::placeholders::last);

        return ground + src.source.value.dot(src.final_weight_factors);
    };

    const double ground_weight =
        2 * input.mu * input.albedo *
        sources.final_weight_factors(Eigen::placeholders::last) *
        sources.beamtrans.value(Eigen::placeholders::last);

    const double ground_source =
        (soln.particular[0].G_plus_bottom.value(Eigen::placeholders::last) +
         soln.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
             soln.homog[0].X_plus.value(Eigen::placeholders::last) *
             soln.homog[0].omega.value(Eigen::placeholders::last) +
         soln.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
             soln.homog[0].X_minus.value(Eigen::placeholders::last)) *
        ground_weight;

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();
    std::array<Eigen::MatrixXd, 2> bvp_storage;
    bvp_storage[0].resize(2 * nlyr, 1);
    bvp_storage[1].resize(2 * nlyr, 1);

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    sasktran2::twostream::backprop::full(input, soln, sources,
                                         sources.final_weight_factors,
                                         bvp_storage, grad_map, ground_weight);

    for (int i = 0; i < nlyr; ++i) {
        sasktran2::twostream::backprop::od(
            input,
            (-sources.final_weight_factors(i) * sources.source.value(i) *
             atten.row(i)),
            grad_map);
    }

    sasktran2::twostream::backprop::od(
        input,
        (-ground_source * Eigen::RowVectorXd::Ones(nlyr) / viewing_zenith),
        grad_map);

    const double eps = 1e-6;
    for (int i = 0; i < nlyr; ++i) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.csz = 0.6;
        input_p.init(nlyr);
        input_p.calculate(0);

        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            input_m;
        input_m.atmosphere = atmo.get();
        input_m.geometry_layers = layer_geo.get();
        input_m.csz = 0.6;
        input_m.init(nlyr);
        input_m.calculate(0);

        input_p.od(i) += eps;
        input_p.calculate_derived(0);
        input_m.od(i) -= eps;
        input_m.calculate_derived(0);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_p;
        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            soln_m;
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            sources_p;
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            sources_m;

        soln_p.init(nlyr);
        soln_m.init(nlyr);
        sources_p.init(nlyr);
        sources_m.init(nlyr);

        sasktran2::twostream::solve(input_p, soln_p);
        sasktran2::twostream::solve(input_m, soln_m);

        sources_p.final_weight_factors.noalias() = (-atten * input_p.od);
        sources_p.final_weight_factors =
            sources_p.final_weight_factors.array().exp();
        sources_m.final_weight_factors.noalias() = (-atten * input_m.od);
        sources_m.final_weight_factors =
            sources_m.final_weight_factors.array().exp();

        sasktran2::twostream::post_process(input_p, viewing_zenith, azimuth,
                                           soln_p, sources_p);
        sasktran2::twostream::post_process(input_m, viewing_zenith, azimuth,
                                           soln_m, sources_m);

        const double fp = objective(input_p, soln_p, sources_p);
        const double fm = objective(input_m, soln_m, sources_m);
        const double numerical = (fp - fm) / (2 * eps);
        const double analytic = grad_map.d_extinction(i);

        CAPTURE(i, analytic, numerical);
        REQUIRE(abs(analytic - numerical) < 1e-5);
    }
}

TEST_CASE("Backprop_StartOfRay_toAtmosphere_ext_linear_ground",
          "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    const double viewing_zenith = 0.8;
    const double azimuth = 0.0;

    construct_atmo_geometry_linear_lowlevel(atmo, layer_geo, 0.6);
    atmo->storage().leg_coeff.chip(1, 0).setConstant(0.0);

    const int natmo = atmo->storage().total_extinction.rows();
    const int nlyr = natmo - 1;

    auto run_objective = [&](sasktran2::atmosphere::Atmosphere<1>& atmo_ref) {
        sasktran2::twostream::Input<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            in;
        in.atmosphere = &atmo_ref;
        in.geometry_layers = layer_geo.get();
        in.csz = 0.6;
        in.init(nlyr);
        in.calculate(0);

        sasktran2::twostream::Solution<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            so;
        sasktran2::twostream::Sources<
            sasktran2::twostream::SourceType::ONLY_SOLAR>
            src;
        so.init(nlyr);
        src.init(nlyr);
        sasktran2::twostream::solve(in, so);

        Eigen::MatrixXd atten(nlyr, nlyr);
        atten.setZero();
        const double viewing_secant = 1.0 / viewing_zenith;
        for (int j = 1; j < nlyr; ++j) {
            atten(j, Eigen::seq(0, j - 1)).setConstant(viewing_secant);
        }

        src.final_weight_factors.noalias() = (-atten * in.od);
        src.final_weight_factors = src.final_weight_factors.array().exp();

        sasktran2::twostream::post_process(in, viewing_zenith, azimuth, so,
                                           src);

        const double ground =
            (so.particular[0].G_plus_bottom.value(Eigen::placeholders::last) +
             so.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                 so.homog[0].X_plus.value(Eigen::placeholders::last) *
                 so.homog[0].omega.value(Eigen::placeholders::last) +
             so.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
                 so.homog[0].X_minus.value(Eigen::placeholders::last)) *
            2 * in.mu * in.albedo *
            src.final_weight_factors(Eigen::placeholders::last) *
            src.beamtrans.value(Eigen::placeholders::last);

        return ground + src.source.value.dot(src.final_weight_factors);
    };

    sasktran2::twostream::Input<sasktran2::twostream::SourceType::ONLY_SOLAR>
        input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;
    input.init(nlyr);
    input.calculate(0);

    sasktran2::twostream::Solution<sasktran2::twostream::SourceType::ONLY_SOLAR>
        soln;
    sasktran2::twostream::Sources<sasktran2::twostream::SourceType::ONLY_SOLAR>
        sources;
    soln.init(nlyr);
    sources.init(nlyr);
    sasktran2::twostream::solve(input, soln);

    Eigen::MatrixXd atten(nlyr, nlyr);
    atten.setZero();
    const double viewing_secant = 1.0 / viewing_zenith;
    for (int j = 1; j < nlyr; ++j) {
        atten(j, Eigen::seq(0, j - 1)).setConstant(viewing_secant);
    }

    sources.final_weight_factors.noalias() = (-atten * input.od);
    sources.final_weight_factors = sources.final_weight_factors.array().exp();

    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    const double ground_weight =
        2 * input.mu * input.albedo *
        sources.final_weight_factors(Eigen::placeholders::last) *
        sources.beamtrans.value(Eigen::placeholders::last);

    const double ground_source =
        (soln.particular[0].G_plus_bottom.value(Eigen::placeholders::last) +
         soln.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
             soln.homog[0].X_plus.value(Eigen::placeholders::last) *
             soln.homog[0].omega.value(Eigen::placeholders::last) +
         soln.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
             soln.homog[0].X_minus.value(Eigen::placeholders::last)) *
        ground_weight;

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();
    std::array<Eigen::MatrixXd, 2> bvp_storage;
    bvp_storage[0].resize(2 * nlyr, 1);
    bvp_storage[1].resize(2 * nlyr, 1);

    sasktran2::twostream::backprop::GradientMap<
        sasktran2::twostream::SourceType::ONLY_SOLAR>
        grad_map(*atmo, grad.data());

    sasktran2::twostream::backprop::full(input, soln, sources,
                                         sources.final_weight_factors,
                                         bvp_storage, grad_map, ground_weight);

    for (int i = 0; i < nlyr; ++i) {
        sasktran2::twostream::backprop::od(
            input,
            (-sources.final_weight_factors(i) * sources.source.value(i) *
             atten.row(i)),
            grad_map);
    }

    sasktran2::twostream::backprop::od(
        input,
        (-ground_source * Eigen::RowVectorXd::Ones(nlyr) / viewing_zenith),
        grad_map);

    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> out(
        1, atmo->num_deriv(), true);
    sasktran2::twostream::backprop::map_to_atmosphere(input, grad_map, out);

    for (int i = 0; i < natmo; ++i) {
        const double oldk = atmo->storage().total_extinction(i, 0);
        const double eps = 1e-4 * oldk;

        atmo->storage().total_extinction(i, 0) = oldk + eps;
        const double fp = run_objective(*atmo);

        atmo->storage().total_extinction(i, 0) = oldk - eps;
        const double fm = run_objective(*atmo);

        atmo->storage().total_extinction(i, 0) = oldk;

        const double numerical = (fp - fm) / (2 * eps);
        const double analytic = out.deriv(i);

        CAPTURE(i, analytic, numerical);
        REQUIRE(abs(analytic - numerical) <
                2e-2 * std::max(1.0, abs(analytic)));
    }
}

#endif
