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

        void Regular_Q_S(const std::complex<double>& refractive_index,
                         const Eigen::VectorXd& size_param,
                         const Eigen::VectorXd& cos_angles,
                         Eigen::VectorXd& Qext, Eigen::VectorXd& Qsca,
                         Eigen::MatrixXcd& S1, Eigen::MatrixXcd& S2);

        void Small_Q_S(const std::complex<double>& refractive_index,
                       const Eigen::VectorXd& size_param,
                       const Eigen::VectorXd& cos_angles, Eigen::VectorXd& Qext,
                       Eigen::VectorXd& Qsca, Eigen::MatrixXcd& S1,
                       Eigen::MatrixXcd& S2);

      public:
        LinearizedMie();
        ~LinearizedMie(){};
        void Dn(Eigen::MatrixXcd& Dn_matrix,
                const std::complex<double>& refractive_index,
                const Eigen::VectorXd& size_param, const int N);
        void An_Bn(const std::complex<double>& refractive_index,
                   const Eigen::VectorXd& size_param, const int N,
                   Eigen::MatrixXcd& An_matrix, Eigen::MatrixXcd& Bn_matrix);
        void Tau_Pi(const Eigen::VectorXd& size_param,
                    const Eigen::VectorXd& cos_angles, const int N,
                    Eigen::Tensor<double, 3>& tau_tensor,
                    Eigen::Tensor<double, 3>& pi_tensor);
    };
} // namespace sasktran2::mie
