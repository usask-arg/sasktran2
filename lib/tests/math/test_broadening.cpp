#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("Run Line Broadening", "[sasktran2][spectroscopy]") {
    int n_lines = 5000;
    int n_geo = 10;

    Eigen::VectorXd line_center =
        Eigen::VectorXd::LinSpaced(n_lines, 1000, 2000);
    Eigen::VectorXd line_intensity =
        Eigen::VectorXd::LinSpaced(n_lines, 0.1, 1.0);
    Eigen::VectorXd lower_energy =
        Eigen::VectorXd::LinSpaced(n_lines, 0.1, 1.0);
    Eigen::VectorXd gamma_air = Eigen::VectorXd::LinSpaced(n_lines, 0.1, 1.0);
    Eigen::VectorXd gamma_self = Eigen::VectorXd::LinSpaced(n_lines, 0.1, 1.0);
    Eigen::VectorXd delta_air = Eigen::VectorXd::LinSpaced(n_lines, 0.1, 1.0);
    Eigen::VectorXd n_air = Eigen::VectorXd::LinSpaced(n_lines, 0.1, 1.0);
    Eigen::VectorXi iso_id = Eigen::VectorXi::LinSpaced(n_lines, 1, 1);
    Eigen::VectorXd partitions = Eigen::VectorXd::LinSpaced(1, 1, 1);
    Eigen::VectorXd molecular_mass = Eigen::VectorXd::LinSpaced(1, 1, 1);
    Eigen::VectorXd pressure = Eigen::VectorXd::LinSpaced(n_geo, 1, 1);
    Eigen::VectorXd pself = Eigen::VectorXd::LinSpaced(n_geo, 0, 0);
    Eigen::VectorXd temperature = Eigen::VectorXd::LinSpaced(n_geo, 280, 280);

    Eigen::VectorXd wavenumber_grid =
        Eigen::VectorXd::LinSpaced(1000 * 100, 1000, 2000);
    Eigen::MatrixXd result =
        Eigen::MatrixXd::Zero(wavenumber_grid.size(), n_geo);

    sasktran2::math::spectroscopy::voigt_broaden(
        line_center, line_intensity, lower_energy, gamma_air, gamma_self,
        delta_air, n_air, iso_id, partitions, molecular_mass, pressure, pself,
        temperature, wavenumber_grid, result);
}
