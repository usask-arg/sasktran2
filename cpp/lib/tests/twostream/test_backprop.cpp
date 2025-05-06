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
        .total_extinction(Eigen::seq(10, 13), Eigen::all)
        .setConstant(0.02);

    atmo->storage().ssa.setConstant(0.8);

    atmo->storage().leg_coeff.chip(0, 0).setConstant(1);
    atmo->storage().leg_coeff.chip(1, 0).setConstant(0.5);

    sasktran2::Config config;
    sasktran_disco::PersistentConfiguration<1> pconfig;
    std::vector<sasktran2::raytracing::TracedRay> traced_rays;
    sasktran_disco::SKTRAN_DO_UserSpec specs;
    pconfig.configure(specs, config, 0.6, nlyr, traced_rays);

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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * od
    Eigen::RowVectorXd weights = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * ssa
    Eigen::RowVectorXd weights = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

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

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve_layers(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

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

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve_layers(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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

    std::vector<sasktran2::twostream::backprop::GradientMap> grad_map;
    grad_map.emplace_back(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::MatrixXd, 2> weights;
    weights[0].resize(2 * (natmo - 1), 1);
    weights[0].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    weights[1].resize(2 * (natmo - 1), 1);
    weights[1].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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
        REQUIRE(abs(grad_map[0].d_ssa(i) - numerical_grad(i)) < 1e-4);
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    Eigen::RowVectorXd weights;
    weights = Eigen::RowVectorXd::Random((natmo - 1));

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    sasktran2::twostream::Sources sources;
    soln.init(natmo - 1);
    sources.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);
    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);
        sasktran2::twostream::Sources sources_p;
        sources_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * od
    Eigen::RowVectorXd weights = Eigen::RowVectorXd::Random(natmo);

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-8;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * od
    Eigen::RowVectorXd weights = Eigen::RowVectorXd::Random(natmo - 1);

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();

    input.init(natmo - 1);

    input.calculate(0);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        input_p.calculate_base(0);
        input_p.od(i) += eps;
        input_p.calculate_derived(0);

        numerical_grad(i) = (weights.dot(input_p.average_secant) -
                             weights.dot(input.average_secant)) /
                            eps;
    }

    sasktran2::twostream::backprop::secant(input, weights, grad_map);

    for (int i = 0; i < natmo - 1; i++) {
        REQUIRE(abs(grad_map.d_extinction(i) - numerical_grad(i)) < 1e-6);
    }
}

TEST_CASE("Backprop_Particular_ext", "[twostream][backprop]") {
    std::unique_ptr<sasktran2::atmosphere::Atmosphere<1>> atmo;
    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> layer_geo;

    construct_atmo_geometry(atmo, layer_geo, 0.01);

    int natmo = atmo->storage().total_extinction.rows();

    Eigen::RowVectorXd grad(3 * natmo);
    grad.setZero();

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

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

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.01;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve_layers(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-7;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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

    std::vector<sasktran2::twostream::backprop::GradientMap> grad_map;
    grad_map.emplace_back(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::MatrixXd, 2> weights;
    weights[0].resize(2 * (natmo - 1), 1);
    weights[0].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    weights[1].resize(2 * (natmo - 1), 1);
    weights[1].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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
        REQUIRE(abs(grad_map[0].d_extinction(i) - numerical_grad(i)) < 1e-4);
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    Eigen::RowVectorXd weights;
    weights = Eigen::RowVectorXd::Random((natmo - 1));

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.01;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    sasktran2::twostream::Sources sources;
    soln.init(natmo - 1);
    sources.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);
    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);
        sasktran2::twostream::Sources sources_p;
        sources_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    Eigen::RowVectorXd weights;
    weights = Eigen::RowVectorXd::Random((natmo - 1));

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    sasktran2::twostream::Sources sources;
    soln.init(natmo - 1);
    sources.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);
    sasktran2::twostream::post_process(input, viewing_zenith, azimuth, soln,
                                       sources);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);
        sasktran2::twostream::Sources sources_p;
        sources_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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

    sasktran2::twostream::backprop::GradientMap grad_map(*atmo, grad.data());

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

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve_layers(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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

    std::vector<sasktran2::twostream::backprop::GradientMap> grad_map;
    grad_map.emplace_back(*atmo, grad.data());

    // Let our output quantity be random weights * homog.k and Xplus/Xminus
    std::array<Eigen::MatrixXd, 2> weights;
    weights[0].resize(2 * (natmo - 1), 1);
    weights[0].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    weights[1].resize(2 * (natmo - 1), 1);
    weights[1].col(0) = Eigen::RowVectorXd::Random(2 * (natmo - 1));

    sasktran2::twostream::Input input;
    input.atmosphere = atmo.get();
    input.geometry_layers = layer_geo.get();
    input.csz = 0.6;

    input.init(natmo - 1);

    input.calculate(0);

    sasktran2::twostream::Solution soln;
    soln.init(natmo - 1);
    sasktran2::twostream::solve(input, soln);

    Eigen::VectorXd numerical_grad(natmo);

    double eps = 1e-6;
    for (int i = 0; i < natmo - 1; i++) {
        sasktran2::twostream::Input input_p;
        input_p.atmosphere = atmo.get();
        input_p.geometry_layers = layer_geo.get();
        input_p.init(natmo - 1);

        sasktran2::twostream::Solution soln_p;
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
        REQUIRE(abs(grad_map[0].d_b1(i) - numerical_grad(i)) < 1e-4);
    }
}

#endif
