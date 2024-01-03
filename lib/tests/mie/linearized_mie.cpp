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

TEST_CASE("LinMie multiple size_params", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();

    auto size_param = Eigen::VectorXd::LinSpaced(2, 0.101, 0.2);
    auto refractive_index = std::complex<double>(0.75, 0.0);
    Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(7, 0, 180.0);
    Eigen::VectorXd cos_angles = cos(angles.array() * EIGEN_PI / 180.0);

    auto result = mie.calculate(size_param, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qext.size() - 2) < 0.0001);
    REQUIRE(fabs(result.values.Qsca.size() - 2) < 0.0001);

    REQUIRE(fabs(result.values.S1.cols() - 7) < 0.0001);
    REQUIRE(fabs(result.values.S1.rows() - 2) < 0.0001);
    REQUIRE(fabs(result.values.S2.cols() - 7) < 0.0001);
    REQUIRE(fabs(result.values.S2.rows() - 2) < 0.0001);
}

TEST_CASE("LinMie Qext Qsca test non absorbing (miepython)",
          "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    double lambda0 = 0.6328;
    double radius = 0.525;

    auto size_param_lambda = Eigen::VectorXd::LinSpaced(
        1, 2 * EIGEN_PI * radius / lambda0, 2 * EIGEN_PI * radius / lambda0);
    auto refractive_index = std::complex<double>(1.55, 0);
    Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(7, 0, 180.0);
    Eigen::VectorXd cos_angles = cos(angles.array() * EIGEN_PI / 180.0);

    auto result =
        mie.calculate(size_param_lambda, refractive_index, cos_angles, true);

    REQUIRE(fabs(result.values.Qext(0) - 3.10543) < 0.00001);
    REQUIRE(fabs(result.values.Qsca(0) - 3.10543) < 0.00001);

    // MIEV0 Test Case 5
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_2 = Eigen::VectorXd::LinSpaced(1, 0.099, 0.099);
    result = mie.calculate(size_param_2, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000007) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.654225E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.654225E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.653815E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 1.574051E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.432172E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.652694E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 9.087788E-09) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() + 8.261265E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() + 1.651163E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 9.797186E-23) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() - 2.938374E-08) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() + 1.649631E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 9.087788E-09) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() - 8.250360E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() + 1.648510E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 1.574051E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() - 1.427725E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() + 1.648100E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 1.817558E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() - 1.648100E-04) < 1e-6);

    // MIEV0 Test Case 6
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_3 = Eigen::VectorXd::LinSpaced(1, 0.101, 0.101);
    result = mie.calculate(size_param_3, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000008) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 0.000008) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 2.048754E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.756419E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 2.048754E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.756419E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 2.048754E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.755965E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 1.774273E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.520629E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 2.048753E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.754726E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 1.024377E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() + 8.771198E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 2.048751E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() + 1.753033E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 1.845057E-15) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() - 3.247270E-08) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 2.048750E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() + 1.751341E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 1.024375E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() - 8.759147E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 2.048749E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() + 1.750103E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 1.774269E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() - 1.515715E-04) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 2.048749E-08) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() + 1.749650E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 2.048749E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() - 1.749650E-04) < 1e-6);

    // MIEV0 Test Case 7
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_10 = Eigen::VectorXd::LinSpaced(1, 10, 10);
    result = mie.calculate(size_param_10, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.232265) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.580662E+01) <
            1e-5); // not enough digits for 1e-5
    REQUIRE(fabs(result.values.S1(0).imag() + 9.758097E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 5.580662E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(0).imag() + 9.758097E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() + 7.672879E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 1.087317E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).real() + 1.092923E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).imag() - 9.629667E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 3.587894E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.756177E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 3.427411E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 8.082691E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() + 1.785905E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() + 5.232828E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() + 5.148748E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 7.027288E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 1.537971E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() + 8.329374E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 6.908338E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() - 2.152693E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() + 4.140427E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 1.876851E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() - 5.247557E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.923391E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() + 1.078568E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() + 3.608807E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() - 1.078568E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() - 3.608807E-02) < 1e-6);

    // MIEV0 Test Case 8
    refractive_index = std::complex<double>(0.75, 0);
    auto size_param_1000 = Eigen::VectorXd::LinSpaced(1, 1000, 1000);
    result = mie.calculate(size_param_1000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.997908) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 4.994770E+05) < 1e-1);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.336502E+04) < 1e-2);
    REQUIRE(fabs(result.values.S2(0).real() - 4.994770E+05) < 1e-1);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.336502E+04) < 1e-2);

    REQUIRE(fabs(result.values.S1(1).real() + 3.999296E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(1).imag() + 3.316361E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(1).real() + 3.946018E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.147791E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(2).real() + 5.209852E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(2).imag() + 5.776614E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(2).real() + 1.970767E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(2).imag() + 6.937470E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(3).real() + 1.600887E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(3).imag() - 1.348013E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).real() + 4.152365E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(3).imag() - 1.143000E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(4).real() - 8.431720E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(4).imag() + 1.209493E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(4).real() + 4.261732E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(4).imag() - 5.535055E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(5).real() + 7.556092E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(5).imag() + 8.134810E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).real() - 4.218303E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).imag() - 9.100831E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(6).real() - 1.705778E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(6).imag() - 4.842510E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(6).real() + 1.705778E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).imag() + 4.842510E+02) < 1e-4);

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
    Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(7, 0, 180.0);
    Eigen::VectorXd cos_angles = cos(angles.array() * EIGEN_PI / 180.0);
    auto result =
        mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.093923) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 9.395198E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 2.348800E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 2.281705E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 2.348800E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 2.281705E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 2.341722E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 2.217102E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 2.034830E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 1.938171E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 2.322408E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 2.046815E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 1.181704E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 1.075976E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 2.296081E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 1.828349E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 2.729533E-04) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() - 6.702879E-03) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 2.269820E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 1.625401E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 1.114466E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 7.646326E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 2.250635E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 1.486170E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 1.942300E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.271557E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 2.243622E-02) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 1.437106E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 2.243622E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 1.437106E-01) < 1e-6);

    // MIEV0 Test Case 10
    auto size_param_100 = Eigen::VectorXd::LinSpaced(1, 100, 100);
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.096594) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.101321E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.253302E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.243188E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(0).real() - 5.253302E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.243188E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(1).real() + 5.534573E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(1).imag() + 2.971881E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).real() + 8.467204E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.999470E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(2).real() - 1.710488E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.520096E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).real() - 3.310764E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).imag() + 2.709787E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() + 3.655758E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 8.769860E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() + 6.550512E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 4.675370E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 2.414318E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 5.380874E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() - 6.039011E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 1.169971E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(5).real() + 1.222996E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 3.283917E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).real() + 9.653812E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() - 1.474455E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(6).real() + 5.659205E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(6).imag() - 4.650974E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).real() - 5.659205E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).imag() + 4.650974E+01) < 1e-5);

    // // MIEV0 Test Case 11
    auto size_param_10000 = Eigen::VectorXd::LinSpaced(1, 10000, 10000);
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.723857) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.004089E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.010222E+07) < 1e1);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.535815E+05) < 1e-1);
    REQUIRE(fabs(result.values.S2(0).real() - 5.010222E+07) < 1e1);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.535815E+05) < 1e-1);

    REQUIRE(fabs(result.values.S1(1).real() - 3.786814E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(1).imag() + 7.654293E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).real() - 5.074755E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).imag() + 7.515986E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(2).real() + 2.731172E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(2).imag() - 1.326633E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).real() + 3.076558E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).imag() + 1.775975E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(3).real() + 1.061003E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(3).imag() + 1.930155E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).real() - 2.430920E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).imag() - 8.409836E+01) < 1e-5);

    REQUIRE(
        fabs(result.values.S1(4).real() + 1.05813972E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(4).imag() - 2.29841186E+01) <
        1e-5); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(4).real() - 5.90649314E+01) <
        1e-5); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(4).imag() + 5.37028316E+02) <
        1e-4); // different values using miepython with higher number of terms

    REQUIRE(
        fabs(result.values.S1(5).real() + 2.74885513E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(5).imag() - 2.29818128E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(5).real() + 8.03619653E+01) <
        1e-5); // different values using miepython with higher number of terms
    // the following criteria have been adjusted as confirmed by Dan Zawada - 6
    // digits matching with rounding
    REQUIRE(
        fabs(result.values.S2(5).imag() + 4.93916823E+00) <
        1e-5); // different values using miepython with higher number of terms

    REQUIRE(
        fabs(result.values.S1(6).real() + 1.82116215E+02) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(6).imag() + 9.51909674E+02) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).real() - 1.82116215E+02) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).imag() - 9.51909674E+02) <
        1e-3); // different values using miepython with higher number of terms

    // MIEV0 Test Case 12
    refractive_index = std::complex<double>(1.5, -1);
    auto size_param_55 = Eigen::VectorXd::LinSpaced(1, 0.055, 0.055);
    result = mie.calculate(size_param_55, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000011) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 1.014910E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 7.675259E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 8.343879E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 7.675259E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 8.343879E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 7.674331E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 8.343495E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 6.646948E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 7.225169E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 7.671794E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 8.342445E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 3.838246E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 4.169695E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 7.668328E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 8.341012E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 3.132066E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 2.037399E-08) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 7.664863E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 8.339578E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 3.830082E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 4.171317E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 7.662326E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 8.338529E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 6.634986E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 7.221887E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 7.661398E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 8.338145E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 7.661398E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 8.338145E-05) < 1e-6);

    // MIEV0 Test Case 13
    auto size_param_56 = Eigen::VectorXd::LinSpaced(1, 0.056, 0.056);
    result = mie.calculate(size_param_56, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.000012) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 1.033467E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 8.102381E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 8.807251E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 8.102381E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 8.807251E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 8.101364E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 8.806830E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 7.016844E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 7.626381E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 8.098587E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 8.805682E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 4.051865E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 4.401169E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 8.094795E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 8.804113E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 3.427277E-08) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 2.229631E-08) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 8.091003E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 8.802545E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 4.042932E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 4.402945E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 8.088228E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 8.801396E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 7.003755E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 7.622790E-05) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 8.087213E-05) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 8.800976E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 8.087213E-05) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 8.800976E-05) < 1e-6);

    // MIEV0 Test Case 14
    result = mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 0.6634538) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.336321E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.840802E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 1.905153E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 5.840802E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 1.905153E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 5.657020E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 1.871997E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 5.001610E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 1.456112E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 5.175251E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 1.784426E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 2.879639E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() - 4.105398E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 4.563396E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 1.671665E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 3.622847E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 6.182646E-02) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 4.002117E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 1.566427E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 1.748750E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 1.229586E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 3.621572E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 1.493910E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 3.056823E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.438460E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 3.488438E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 1.468286E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 3.488438E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 1.468286E-01) < 1e-6);

    // MIEV0 Test Case 15
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.283697) < 1e-3);
    REQUIRE(fabs(result.values.Qext(0) - 2.097502) < 1e-2);

    REQUIRE(fabs(result.values.S1(0).real() - 5.243754E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(0).imag() + 2.934167E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(0).real() - 5.243754E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(0).imag() + 2.934167E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(1).real() - 4.049055E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.898456E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).real() - 2.019198E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).imag() - 3.110731E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() + 2.646835E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(2).imag() + 1.929564E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).real() - 9.152743E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() + 7.470202E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 1.268890E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(3).imag() - 2.397474E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(3).real() + 1.232914E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(3).imag() + 7.823167E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 5.149886E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 2.290736E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(4).real() + 7.173357E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 1.655464E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(5).real() + 1.605395E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(5).imag() - 1.418642E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).real() - 1.448052E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.393594E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(6).real() + 2.029360E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(6).imag() - 4.384435E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() - 2.029360E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).imag() + 4.384435E+00) < 1e-6);

    // MIEV0 Test Case 16
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.236575) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.004368) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.010919E+07) < 1e1);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.753404E+05) < 1e-1);
    REQUIRE(fabs(result.values.S2(0).real() - 5.010919E+07) < 1e1);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.753404E+05) < 1e-1);

    REQUIRE(fabs(result.values.S1(1).real() + 3.690394E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.573897E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).real() + 9.333175E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(1).imag() + 1.839736E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(2).real() - 2.391551E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(2).imag() - 3.247786E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).real() + 1.202951E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).imag() + 1.899647E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(3).real() + 2.607463E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(3).imag() - 7.414859E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).real() - 1.013073E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(3).imag() + 1.064666E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(4).real() + 6.183154E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(4).imag() - 2.264970E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(4).real() - 1.334826E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(4).imag() + 1.800859E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(5).real() + 3.368019E+02) < 1e-4);
    REQUIRE(fabs(result.values.S1(5).imag() - 2.115750E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(5).real() - 2.293862E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(5).imag() + 1.996754E+03) < 1e-3);

    REQUIRE(
        fabs(result.values.S1(6).real() + 2.18471714E+02) <
        1e-4); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(6).imag() + 2064.61012226) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).real() - 2.18471714E+02) <
        1e-4); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).imag() - 2064.61012226) <
        1e-3); // different values using miepython with higher number of terms

    // MIEV0 Test Case 17
    refractive_index = std::complex<double>(10, -10);
    result = mie.calculate(size_param_1, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 2.049405) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.532993E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 6.332483E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(0).imag() - 4.179305E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).real() - 6.332483E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(0).imag() - 4.179305E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(1).real() - 6.162264E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(1).imag() - 4.597163E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).real() - 5.573186E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(1).imag() - 2.954338E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() - 5.736317E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(2).imag() - 5.602514E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).real() - 3.525107E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(2).imag() + 5.921611E-03) < 1e-6);

    REQUIRE(fabs(result.values.S1(3).real() - 5.238628E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 6.675352E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).real() - 7.881172E-02) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 3.435544E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(4).real() - 4.825816E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(4).imag() - 7.434033E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).real() + 1.881212E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(4).imag() + 6.028739E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(5).real() - 4.570214E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(5).imag() - 7.809867E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() + 3.793898E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).imag() + 7.473279E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() - 4.485464E-01) < 1e-6);
    REQUIRE(fabs(result.values.S1(6).imag() - 7.912365E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).real() + 4.485464E-01) < 1e-6);
    REQUIRE(fabs(result.values.S2(6).imag() + 7.912365E-01) < 1e-6);

    // MIEV0 Test Case 18
    result = mie.calculate(size_param_100, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.836785) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.071124E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.177811E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(0).imag() + 2.633811E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(0).real() - 5.177811E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(0).imag() + 2.633811E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(1).real() - 5.227436E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(1).imag() + 1.270012E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).real() + 2.380252E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(1).imag() + 3.872567E-01) < 1e-6);

    REQUIRE(fabs(result.values.S1(2).real() + 2.705712E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(2).imag() + 3.951751E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).real() - 2.585821E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(2).imag() - 3.323624E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(3).real() - 1.008860E+00) < 1e-6);
    REQUIRE(fabs(result.values.S1(3).imag() - 4.663027E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(3).real() + 3.479935E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(3).imag() + 4.364245E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(4).real() + 1.505640E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(4).imag() - 4.333057E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(4).real() - 1.360634E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(4).imag() + 4.238302E+01) < 1e-5);

    REQUIRE(fabs(result.values.S1(5).real() + 4.510770E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(5).imag() - 5.199554E+00) < 1e-6);
    REQUIRE(fabs(result.values.S2(5).real() - 4.474564E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(5).imag() + 5.452513E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(6).real() + 4.145383E+01) < 1e-5);
    REQUIRE(fabs(result.values.S1(6).imag() + 1.821808E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).real() - 4.145383E+01) < 1e-5);
    REQUIRE(fabs(result.values.S2(6).imag() - 1.821808E+01) < 1e-5);

    // MIEV0 Test Case 19
    result =
        mie.calculate(size_param_10000, refractive_index, cos_angles, true);
    REQUIRE(fabs(result.values.Qsca(0) - 1.795393) < 1e-6);
    REQUIRE(fabs(result.values.Qext(0) - 2.005914E+00) < 1e-6);

    REQUIRE(fabs(result.values.S1(0).real() - 5.014786E+07) < 1e1);
    REQUIRE(fabs(result.values.S1(0).imag() + 1.206004E+05) < 1e-1);
    REQUIRE(fabs(result.values.S2(0).real() - 5.014786E+07) < 1e1);
    REQUIRE(fabs(result.values.S2(0).imag() + 1.206004E+05) < 1e-1);

    REQUIRE(fabs(result.values.S1(1).real() + 4.080090E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(1).imag() + 2.664399E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).real() - 3.351286E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(1).imag() - 7.291906E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(2).real() + 1.224040E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(2).imag() - 4.596569E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(2).real() - 4.497446E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(2).imag() + 4.072999E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(3).real() + 4.579490E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(3).imag() + 8.590486E+02) < 1e-4);
    REQUIRE(fabs(result.values.S2(3).real() - 4.313394E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(3).imag() - 4.969719E+02) < 1e-4);

    REQUIRE(fabs(result.values.S1(4).real() + 3.356286E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(4).imag() - 3.125121E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(4).real() - 3.171910E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(4).imag() + 3.129068E+03) < 1e-3);

    REQUIRE(fabs(result.values.S1(5).real() + 3.149584E+03) < 1e-3);
    REQUIRE(fabs(result.values.S1(5).imag() - 3.270358E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(5).real() - 3.105243E+03) < 1e-3);
    REQUIRE(fabs(result.values.S2(5).imag() + 3.269355E+03) < 1e-3);

    REQUIRE(
        fabs(result.values.S1(6).real() - 2.25248071E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S1(6).imag() + 3924.46733361) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).real() + 2.25248071E+03) <
        1e-3); // different values using miepython with higher number of terms
    REQUIRE(
        fabs(result.values.S2(6).imag() - 3924.46733361) <
        1e-3); // different values using miepython with higher number of terms

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

TEST_CASE("LinMie S1 S2 (miepython)", "[sasktran2][mie]") {
    // MieV0 test case 5
    auto mie = sasktran2::mie::LinearizedMie();
    // auto size_param_1 = Eigen::VectorXd::LinSpaced(1, 1.0, 1.0);
    auto size_param_1 = Eigen::VectorXd::LinSpaced(1, 1.0, 1.0);
    auto refractive_index = std::complex<double>(1.5, -1.0);
    auto angles = Eigen::VectorXd::LinSpaced(7, 0.0, 180.0);
    auto cos_angles = cos(angles.array() / 180.0 * ((double)M_PI));

    auto result =
        mie.calculate(size_param_1, refractive_index, cos_angles, true);
    auto S1 = result.values.S1;
    auto S2 = result.values.S2;
    // auto S1 = result.values.S1.array()*sqrt(((double)
    // M_PI)*size_param_1(0)*size_param_1(0)*result.values.Qext(0)); auto S2 =
    // result.values.S2.array()*sqrt(((double)
    // M_PI)*size_param_1(0)*size_param_1(0)*result.values.Qext(0));

    REQUIRE(fabs(S1(0).real() - 0.584080) < 1e-6);
    REQUIRE(fabs(S1(0).imag() - 0.190515) < 1e-6);
    REQUIRE(fabs(S2(0).real() - 0.584080) < 1e-6);
    REQUIRE(fabs(S2(0).imag() - 0.190515) < 1e-6);

    REQUIRE(fabs(S1(1).real() - 0.565702) < 1e-6);
    REQUIRE(fabs(S1(1).imag() - 0.187200) < 1e-6);
    REQUIRE(fabs(S2(1).real() - 0.500161) < 1e-6);
    REQUIRE(fabs(S2(1).imag() - 0.145611) < 1e-6);

    REQUIRE(fabs(S1(2).real() - 0.517525) < 1e-6);
    REQUIRE(fabs(S1(2).imag() - 0.178443) < 1e-6);
    REQUIRE(fabs(S2(2).real() - 0.287964) < 1e-6);
    REQUIRE(fabs(S2(2).imag() - 0.041054) < 1e-6);

    REQUIRE(fabs(S1(3).real() - 0.456340) < 1e-6);
    REQUIRE(fabs(S1(3).imag() - 0.167167) < 1e-6);
    REQUIRE(fabs(S2(3).real() - 0.0362285) < 1e-6);
    REQUIRE(fabs(S2(3).imag() + 0.0618265) < 1e-6);

    REQUIRE(fabs(S1(4).real() - 0.400212) < 1e-6);
    REQUIRE(fabs(S1(4).imag() - 0.156643) < 1e-6);
    REQUIRE(fabs(S2(4).real() + 0.174875) < 1e-6);
    REQUIRE(fabs(S2(4).imag() + 0.122959) < 1e-6);

    REQUIRE(fabs(S1(5).real() - 0.362157) < 1e-6);
    REQUIRE(fabs(S1(5).imag() - 0.149391) < 1e-6);
    REQUIRE(fabs(S2(5).real() + 0.305682) < 1e-6);
    REQUIRE(fabs(S2(5).imag() + 0.143846) < 1e-6);

    REQUIRE(fabs(S1(6).real() - 0.348844) < 1e-6);
    REQUIRE(fabs(S1(6).imag() - 0.146829) < 1e-6);
    REQUIRE(fabs(S2(6).real() + 0.348844) < 1e-6);
    REQUIRE(fabs(S2(6).imag() + 0.146829) < 1e-6);

    auto size_param_7086 = Eigen::VectorXd::LinSpaced(1, 0.7086, 0.7086);
    refractive_index = std::complex<double>(1.507, -0.002);
    auto cos_angles_2 = Eigen::VectorXd::LinSpaced(1, -1, -1);

    result =
        mie.calculate(size_param_7086, refractive_index, cos_angles_2, true);
    auto S1_2 = result.values.S1.array() /
                sqrt(((double)M_PI) * size_param_7086(0) * size_param_7086(0) *
                     result.values.Qext(0));
    auto S2_2 = result.values.S2.array() /
                sqrt(((double)M_PI) * size_param_7086(0) * size_param_7086(0) *
                     result.values.Qext(0));

    REQUIRE(fabs(S1_2(0).real() - 0.02452300864212876) < 1e-6);
    REQUIRE(fabs(S1_2(0).imag() - 0.29539154027629805) < 1e-6);
    REQUIRE(fabs(S2_2(0).real() + 0.02452300864212876) < 1e-6);
    REQUIRE(fabs(S2_2(0).imag() + 0.29539154027629805) < 1e-6);
}

TEST_CASE("Benchmark for 1000 x cases", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();
    auto refractive_index = std::complex<double>(1.55, 0);
    Eigen::VectorXd angles = Eigen::VectorXd::LinSpaced(7, 0, 180.0);
    Eigen::VectorXd cos_angles = cos(angles.array() * EIGEN_PI / 180.0);
    // auto size_param = Eigen::VectorXd::LinSpaced(1000, 0.1, 10000);
    auto size_param = Eigen::VectorXd::LinSpaced(1000, 0.1, 10);

    BENCHMARK("Test") {
        return mie.calculate(size_param, refractive_index, cos_angles, true);
    };
}
