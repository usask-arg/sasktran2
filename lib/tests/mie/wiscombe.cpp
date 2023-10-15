#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("Check Wiscombe Mie call", "[sasktran2][mie]") {
    auto mie = sasktran2::mie::WiscombeMie();

    Eigen::VectorXd size_param(1);
    Eigen::VectorXcd refractive_index(1);
    Eigen::VectorXd angles(1);

    size_param << 5;
    refractive_index << 1.4;

    angles << 0;

    sasktran2::mie::MieOutput output;

    mie.calculate(size_param, refractive_index, angles, output);

}