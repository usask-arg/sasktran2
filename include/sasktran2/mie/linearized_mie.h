#include <sasktran2/internal_common.h>

namespace sasktran2::mie {
    class LinearizedMie : public MieBase {
      private:
        void internal_calculate(const Eigen::VectorXd& size_param,
                                const std::complex<double>& refractive_index,
                                const Eigen::VectorXd& cos_angles,
                                bool calculate_derivative,
                                MieOutput& output) override;

        // void Dn(Eigen::MatrixXcd& Dn_matrix,
        //  const std::complex<double>& refractive_index,
        //  const Eigen::VectorXd& size_param, const int N);

        void Dn_upwards(const Eigen::VectorXcd& z, const int N,
                        Eigen::MatrixXcd& Dn_matrix);

        void Dn_downwards(const Eigen::VectorXcd& z, const int N,
                          Eigen::MatrixXcd& Dn_matrix);

        void Dn_Lentz(const Eigen::VectorXcd& z_array, const int N,
                      Eigen::VectorXcd& Dn_array);

        void
        LinearizedMie::Regular_Q_S(const std::complex<double>& refractive_index,
                                   const Eigen::VectorXd& size_param,
                                   Eigen::VectorXd& Qext,
                                   Eigen::VectorXd& Qsca);

        void
        LinearizedMie::Small_Q_S(const std::complex<double>& refractive_index,
                                 const Eigen::VectorXd& size_param,
                                 Eigen::VectorXd& Qext, Eigen::VectorXd& Qsca);

      public:
        LinearizedMie();
        ~LinearizedMie(){};
        void Dn(Eigen::MatrixXcd& Dn_matrix,
                const std::complex<double>& refractive_index,
                const Eigen::VectorXd& size_param, const int N);
        void LinearizedMie::An_Bn(const std::complex<double>& refractive_index,
                                  const Eigen::VectorXd& size_param,
                                  const int N, Eigen::MatrixXcd& An_matrix,
                                  Eigen::MatrixXcd& Bn_matrix);
    };
} // namespace sasktran2::mie
