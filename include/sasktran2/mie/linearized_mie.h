#include <sasktran2/internal_common.h>

namespace sasktran2::mie {
    class LinearizedMie : public MieBase {
      private:
        void internal_calculate(const Eigen::VectorXd& size_param,
                                const std::complex<double>& refractive_index,
                                const Eigen::VectorXd& cos_angles,
                                bool calculate_derivative,
                                MieOutput& output) override;

      public:
        LinearizedMie();
        ~LinearizedMie(){};
    };
} // namespace sasktran2::mie
