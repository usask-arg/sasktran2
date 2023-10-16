#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("LinMie construction", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::LinearizedMie();

    auto size_param = Eigen::VectorXd::LinSpaced(10, 0.1, 1.0);
    auto refractive_index = std::complex<double>(1.5, 0.0);
    auto cos_angles = Eigen::VectorXd::LinSpaced(10, -1.0, 1.0);

    auto result = mie.calculate(size_param, refractive_index, cos_angles, true);
}
