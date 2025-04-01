#pragma once

#include "sasktran2/mie/mie.h"
#include <sasktran2/internal_common.h>

namespace sasktran2::mie {
    class LinearizedMieWorker {
      private:
        Eigen::VectorXcd m_An;
        Eigen::VectorXcd m_Bn;

        Eigen::VectorXcd m_Dn;

        const Eigen::MatrixXd& m_tau_matrix;
        const Eigen::MatrixXd& m_pi_matrix;

        int m_local_N;

      public:
        LinearizedMieWorker(int max_N, const Eigen::MatrixXd& tau, const Eigen::MatrixXd& pi) : m_local_N(max_N),
            m_tau_matrix(tau), m_pi_matrix(pi) {
            m_An.resize(max_N);
            m_Bn.resize(max_N);
            m_Dn.resize(max_N);
        }

        void set_local_N(double size_param) {
            m_local_N = int(size_param + 4.05 * pow(size_param, 0.33333) + 2.0) + 2;

            if(m_local_N > m_An.size()) {
                m_An.resize(m_local_N);
                m_Bn.resize(m_local_N);
                m_Dn.resize(m_local_N);
            }
        }

        void regular_Q_S(const std::complex<double>& refractive_index,
                         const double& size_param, int size_index,
                        MieData& output);

        void small_Q_S(const std::complex<double>& refractive_index,
                       const double& size_param,
                       const Eigen::VectorXd& cos_angles, int size_index,
                       MieData& output);

        void An_Bn(const std::complex<double>& refractive_index,
                   const double& size_param);

        void Dn(const std::complex<double>& refractive_index,
                const double& size_param);

        void Dn_upwards(const std::complex<double>& refractive_index,
                        const std::complex<double>& z);

        void Dn_downwards(const std::complex<double>& refractive_index,
                          const std::complex<double>& z);

        void Dn_Lentz(const std::complex<double>& z,
                      std::complex<double>& result);
    };

    class LinearizedMie : public MieBase {
      private:
        std::vector<LinearizedMieWorker> m_thread_storage;
        Eigen::MatrixXd m_tau_matrix;
        Eigen::MatrixXd m_pi_matrix;

        int m_num_threads;

        void allocate(double max_x, int num_angle);

        void internal_calculate(const Eigen::VectorXd& size_param,
                                const std::complex<double>& refractive_index,
                                const Eigen::VectorXd& cos_angles,
                                bool calculate_derivative,
                                MieOutput& output) override;

        void tau_pi(const Eigen::VectorXd& cos_angles);

      public:
        LinearizedMie(int num_threads = 1);
        ~LinearizedMie(){};
    };

} // namespace sasktran2::mie
