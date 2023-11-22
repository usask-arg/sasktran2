#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("LinMie construction", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();

    auto size_param = Eigen::VectorXd::LinSpaced(10, 0.1, 1.0);
    auto refractive_index = std::complex<double>(1.5, 0.0);
    auto cos_angles = Eigen::VectorXd::LinSpaced(10, -1.0, 1.0);

    auto result = mie.calculate(size_param, refractive_index, cos_angles, true);
}

TEST_CASE("LinMie Dn test", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    auto size_param = Eigen::VectorXd::LinSpaced(1, 62, 62);
    auto refractive_index = std::complex<double>(1.28, -1.37);
    int nstop = 50;
    Eigen::MatrixXcd Dn_matrix;
    mie.Dn(Dn_matrix, refractive_index, size_param, nstop);

    double ans_real = fabs(Dn_matrix(9, 0).real() - 0.004087);
    double ans_imag = fabs(Dn_matrix(9, 0).imag() - 1.0002620);
    REQUIRE(ans_real < 0.00001);
    REQUIRE(ans_imag < 0.00001);
}

TEST_CASE("LinMie An Bn test", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    auto size_param = Eigen::VectorXd::LinSpaced(1, 50, 50);
    auto refractive_index = std::complex<double>(4.0 / 3.0, 0);
    double x = size_param.maxCoeff(); // largest size parameter, use to do the
                                      // calculations for highest N
    int N = int(x + 4.05 * pow(x, 0.33333) + 2.0) + 1;
    Eigen::MatrixXcd An_matrix;
    Eigen::MatrixXcd Bn_matrix;
    mie.An_Bn(refractive_index, size_param, N, An_matrix, Bn_matrix);

    REQUIRE(fabs(An_matrix(0, 0).real() - 0.5311058892948411929) < 0.00000001);
    REQUIRE(fabs(An_matrix(0, 0).imag() - 0.4990314856310943073) < 0.00000001);
    REQUIRE(fabs(Bn_matrix(0, 0).real() - 0.7919244759352004773) < 0.00001);
    REQUIRE(fabs(Bn_matrix(0, 0).imag() - 0.4059311522289938238) < 0.00001);

    auto size_param_2 = Eigen::VectorXd::LinSpaced(1, 2, 2);
    refractive_index = std::complex<double>(1.5, -1);
    x = size_param_2.maxCoeff(); // largest size parameter, use to do the
                                 // calculations for highest N
    N = int(x + 4.05 * pow(x, 0.33333) + 2.0) + 1;
    mie.An_Bn(refractive_index, size_param_2, N, An_matrix, Bn_matrix);

    REQUIRE(fabs(An_matrix(0, 0).real() - 0.5465202033970914511) < 0.00000001);
    REQUIRE(fabs(An_matrix(0, 0).imag() - 0.1523738572575972279) < 0.00000001);
    REQUIRE(fabs(Bn_matrix(0, 0).real() - 0.3897147278879423235) < 0.00001);
    REQUIRE(fabs(Bn_matrix(0, 0).imag() + 0.2278960752564908264) < 0.00001);

    refractive_index = std::complex<double>(1.1, -25);
    x = size_param_2.maxCoeff(); // largest size parameter, use to do the
                                 // calculations for highest N
    N = int(x + 4.05 * pow(x, 0.33333) + 2.0) + 1;
    mie.An_Bn(refractive_index, size_param_2, N, An_matrix, Bn_matrix);

    REQUIRE(fabs(An_matrix(1, 0).real() - 0.324433578437) < 0.0001);
    REQUIRE(fabs(An_matrix(1, 0).imag() - 0.465627763266) < 0.0001);
    REQUIRE(fabs(Bn_matrix(1, 0).real() - 0.060464399088) < 0.0001);
    REQUIRE(fabs(Bn_matrix(1, 0).imag() + 0.236805417045) < 0.0001);
}

TEST_CASE("LinMie Qext Qsca test non absorbing (miepython)",
          "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    double lambda0 = 0.6328;
    double radius = 0.525;

    auto size_param_lambda = Eigen::VectorXd::LinSpaced(
        1, 2 * EIGEN_PI * radius / lambda0, 2 * EIGEN_PI * radius / lambda0);
    auto refractive_index = std::complex<double>(1.55, 0);
    auto cos_angles = Eigen::VectorXd::LinSpaced(10, -1.0, 1.0);

    auto result =
        mie.calculate(size_param_lambda, refractive_index, cos_angles, true);

    REQUIRE(fabs(result.values.Qext(0) - 3.10543) < 0.00001);
    REQUIRE(fabs(result.values.Qsca(0) - 3.10543) < 0.00001);

    // MIEV0 Test Case 5
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_2 = Eigen::VectorXd::LinSpaced(1, 0.099, 0.099);
    result = mie.calculate(size_param_2, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000007) < 1e-6);

    // MIEV0 Test Case 6
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_3 = Eigen::VectorXd::LinSpaced(1, 0.101, 0.101);
    result = mie.calculate(size_param_3, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000008) < 1e-6);

    // MIEV0 Test Case 7
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_10 = Eigen::VectorXd::LinSpaced(1, 10, 10);
    result = mie.calculate(size_param_10, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.232265) < 1e-6);

    // MIEV0 Test Case 8
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_1000 = Eigen::VectorXd::LinSpaced(1, 1000, 1000);
    result = mie.calculate(size_param_1000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.997908) < 1e-6);

    // OLD MIEV0 Test Case 1
    refractive_index = std::complex<double>(1.5, 0);
    result = mie.calculate(size_param_10, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.8820) < 1e-4);

    // OLD MIEV0 Test Case 2
    auto size_param_100 = Eigen::VectorXd::LinSpaced(1, 100, 100);
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.0944) < 1e-4);

    // OLD MIEV0 Test Case 3
    result = mie.calculate(size_param_1000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.0139) < 1e-4);

    // OLD MIEV0 Test Case 4
    auto size_param_5000 = Eigen::VectorXd::LinSpaced(1, 5000, 5000);
    result = mie.calculate(size_param_5000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.0086) < 1e-4);

    // test non-dialectric
    refractive_index = std::complex<double>(1.55, -0.1);
    result =
        mie.calculate(size_param_lambda, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 2.86165188243) < 1e-7);
    REQUIRE(fabs(result.values.Qsca(0) - 1.66424911991) < 1e-7);
}

TEST_CASE("LinMie Qext Qsca test absorbing (miepython)", "[sasktran2][mie]") {
    // MIEV0 Test Case 9
    auto mie = sasktran2::mie::LinearizedMie();
    auto size_param_1 = Eigen::VectorXd::LinSpaced(1, 1, 1);
    auto refractive_index = std::complex<double>(1.33, -0.00001);
    auto cos_angles = Eigen::VectorXd::LinSpaced(10, -1.0, 1.0);
    auto result =
        mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.093923) < 1e-6);

    // MIEV0 Test Case 10
    auto size_param_100 = Eigen::VectorXd::LinSpaced(1, 100, 100);
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.096594) < 1e-6);

    // MIEV0 Test Case 11
    auto size_param_10000 = Eigen::VectorXd::LinSpaced(1, 10000, 10000);
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.723857) < 1e-6);

    // MIEV0 Test Case 12
    refractive_index = std::complex<double>(1.5, -1);
    auto size_param_55 = Eigen::VectorXd::LinSpaced(1, 0.055, 0.055);
    result = mie.calculate(size_param_55, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000011) < 1e-6);

    // MIEV0 Test Case 13
    auto size_param_56 = Eigen::VectorXd::LinSpaced(1, 0.056, 0.056);
    result = mie.calculate(size_param_56, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000012) < 1e-6);

    // MIEV0 Test Case 14
    result = mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.6634538) < 1e-6);

    // MIEV0 Test Case 15
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.283697) < 1e-3);
    REQUIRE(fabs(result.values.Qext(0) - 2.097502) < 1e-2);

    // MIEV0 Test Case 16
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.236575) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.004368) < 1e-6);

    // MIEV0 Test Case 17
    refractive_index = std::complex<double>(10, -10);
    result = mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.049405) < 1e-6);

    // MIEV0 Test Case 18
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.836785) < 1e-6);

    // MIEV0 Test Case 19
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.795393) < 1e-6);

    // test single non-magnetic
    refractive_index = std::complex<double>(1.5, -0.5);
    auto size_param_25 = Eigen::VectorXd::LinSpaced(1, 2.5, 2.5);
    result = mie.calculate(size_param_25, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.097071819088392) < 1e-7);
    REQUIRE(fabs(result.values.Qext(0) - 2.562873497454734) < 1e-7);
}

TEST_CASE("LinMie Qext Qsca test small (miepython)", "[sasktran2][mie]") {
    // MieV0 test case 5
    auto mie = sasktran2::mie::LinearizedMie();
    auto size_param_99 = Eigen::VectorXd::LinSpaced(1, 0.099, 0.099);
    auto refractive_index = std::complex<double>(0.75, 0);
    auto cos_angles = Eigen::VectorXd::LinSpaced(10, -1.0, 1.0);

    auto result =
        mie.calculate(size_param_99, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 0.000007) < 1e-6);

    // MIEV0 test case 6
    auto size_param_101 = Eigen::VectorXd::LinSpaced(1, 0.101, 0.101);
    result = mie.calculate(size_param_101, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 0.000008) < 1e-6);

    refractive_index = std::complex<double>(1.5, -1);
    auto size_param_55 = Eigen::VectorXd::LinSpaced(1, 0.055, 0.055);
    result = mie.calculate(size_param_55, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 0.101491) < 1e-6);

    auto size_param_56 = Eigen::VectorXd::LinSpaced(1, 0.056, 0.056);
    result = mie.calculate(size_param_56, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 0.103347) < 1e-6);

    refractive_index = std::complex<double>(1e-10, -1e-10);
    result = mie.calculate(size_param_99, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 6.329394366790141e-05) < 1e-6);

    result = mie.calculate(size_param_101, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext(0) - 6.853305788484561e-05) < 1e-6);
}
