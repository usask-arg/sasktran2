#include "sasktran2/mie/mie.h"
#include <sasktran2/internal_common.h>

namespace sasktran2::mie {
    class LinearizedMie_OLD : public MieBase {
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
        LinearizedMie_OLD();
        ~LinearizedMie_OLD(){};
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

    class LinearizedMie : public MieBase {
      private:
        struct ThreadStorage {
            Eigen::VectorXcd An;
            Eigen::VectorXcd Bn;

            Eigen::VectorXcd Dn;

            int local_N;
        };

        std::vector<ThreadStorage> m_thread_storage;
        Eigen::MatrixXd m_tau_matrix;
        Eigen::MatrixXd m_pi_matrix;

        void allocate(double max_x, int num_angle);

        void internal_calculate(const Eigen::VectorXd& size_param,
                                const std::complex<double>& refractive_index,
                                const Eigen::VectorXd& cos_angles,
                                bool calculate_derivative,
                                MieOutput& output) override;

        void regular_Q_S(const std::complex<double>& refractive_index,
                         const double& size_param, int size_index,
                         int thread_idx, MieData& output);

        void small_Q_S(const std::complex<double>& refractive_index,
                       const double& size_param,
                       const Eigen::VectorXd& cos_angles, int size_index,
                       int thread_idx, MieData& output);

        void An_Bn(const std::complex<double>& refractive_index,
                   const double& size_param, int thread_idx);

        void Dn(const std::complex<double>& refractive_index,
                const double& size_param, int thread_idx);

        void Dn_upwards(const std::complex<double>& refractive_index,
                        const std::complex<double>& z, int thread_idx);

        void Dn_downwards(const std::complex<double>& refractive_index,
                          const std::complex<double>& z, int thread_idx);

        void Dn_Lentz(const std::complex<double>& z, int thread_idx,
                      std::complex<double>& result);

        void tau_pi(const Eigen::VectorXd& cos_angles);

      public:
        LinearizedMie();
        ~LinearizedMie(){};
    };

} // namespace sasktran2::mie
