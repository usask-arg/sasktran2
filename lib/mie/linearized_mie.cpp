#include <sasktran2/mie/mie.h>
#include <sasktran2/mie/linearized_mie.h>

namespace sasktran2::mie {

    LinearizedMie::LinearizedMie() {}

    void LinearizedMie::internal_calculate(
        const Eigen::VectorXd& size_param,
        const std::complex<double>& refractive_index,
        const Eigen::VectorXd& cos_angles, bool calculate_derivative,
        MieOutput& output) {}

} // namespace sasktran2::mie
