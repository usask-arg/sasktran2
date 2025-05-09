#include "sasktran2/math/trig.h"
#include <sasktran2/mie/mie.h>
#include <sasktran2/mie/linearized_mie.h>
#ifdef SKTRAN_OPENMP_SUPPORT
#include <omp.h>
#endif

namespace sasktran2::mie {

    LinearizedMie::LinearizedMie(int num_threads)
        : m_num_threads(num_threads) {}

    void LinearizedMie::internal_calculate(
        const Eigen::VectorXd& size_param,
        const std::complex<double>& refractive_index,
        const Eigen::VectorXd& cos_angles, bool calculate_derivative,
        MieOutput& output) {
        // max x value
        double max_x = size_param.maxCoeff();

        this->allocate(max_x, cos_angles.size());

        tau_pi(cos_angles);

#pragma omp parallel for num_threads(m_num_threads)
        for (int i = 0; i < size_param.size(); ++i) {
#ifdef SKTRAN_OPENMP_SUPPORT
            auto& worker = m_thread_storage[omp_get_thread_num()];
#else
            auto& worker = m_thread_storage[0];
#endif
            worker.set_local_N(size_param(i));
            if (abs(refractive_index) * size_param(i) < 0.1) {
                worker.small_Q_S(refractive_index, size_param(i), cos_angles, i,
                                 output.values);
            } else {
                worker.regular_Q_S(refractive_index, size_param(i), i,
                                   output.values);
            }
        }
    }

    void LinearizedMie::recalculate_scattering(
        Eigen::Ref<const Eigen::VectorXd> cos_angles, MieOutput& output) {
        double max_x = output.size_param.maxCoeff();
        this->allocate(max_x, cos_angles.size());

        tau_pi(cos_angles);

#pragma omp parallel for num_threads(m_num_threads)
        for (int i = 0; i < output.size_param.size(); ++i) {
#ifdef SKTRAN_OPENMP_SUPPORT
            auto& worker = m_thread_storage[omp_get_thread_num()];
#else
            auto& worker = m_thread_storage[0];
#endif
            worker.set_local_N(output.size_param(i));
            if (abs(output.refractive_index) * output.size_param(i) < 0.1) {
                worker.small_Q_S(output.refractive_index, output.size_param(i),
                                 cos_angles, i, output.values);
            } else {
                worker.regular_Q_S_scat(output.refractive_index,
                                        output.size_param(i), i, output.values);
            }
        }
    }

    void LinearizedMieWorker::regular_Q_S(
        const std::complex<double>& refractive_index, const double& size_param,
        int size_index, MieData& output) {

        An_Bn(refractive_index, size_param);
        int N = m_local_N;

        Eigen::VectorXcd& An = m_An;
        Eigen::VectorXcd& Bn = m_Bn;

        double& Qext = output.Qext(size_index);
        double& Qsca = output.Qsca(size_index);

        Qext = 0.0;
        Qsca = 0.0;
        for (int n = 0; n < N; ++n) {
            Qext += (2.0 * (n + 1) + 1.0) * (An(n).real() + Bn(n).real());
            Qsca += (2.0 * (n + 1) + 1.0) *
                    (math::sqr(abs(An(n))) + math::sqr(abs(Bn(n))));
        }
        Qext *= 2.0 / size_param / size_param;
        Qsca *= 2.0 / size_param / size_param;

        regular_Q_S_scat(refractive_index, size_param, size_index, output);
    }

    void LinearizedMieWorker::regular_Q_S_scat(
        const std::complex<double>& refractive_index, const double& size_param,
        int size_index, MieData& output) {
        Eigen::MatrixXcd& S1 = output.S1;
        Eigen::MatrixXcd& S2 = output.S2;
        const Eigen::VectorXcd& An = m_An;
        const Eigen::VectorXcd& Bn = m_Bn;

        S1(size_index, Eigen::all).setZero();
        S2(size_index, Eigen::all).setZero();
        for (int i = 0; i < m_local_N; ++i) {
            double n_factor = (2.0 * (i + 1.0) + 1.0) / (i + 1.0) / (i + 2.0);

            S1(size_index, Eigen::all).array() +=
                n_factor * (An(i) * m_pi_matrix(i, Eigen::all).array() +
                            Bn(i) * m_tau_matrix(i, Eigen::all).array());
            S2(size_index, Eigen::all).array() +=
                n_factor * (An(i) * m_tau_matrix(i, Eigen::all).array() +
                            Bn(i) * m_pi_matrix(i, Eigen::all).array());
        }
    }

    void LinearizedMieWorker::small_Q_S(
        const std::complex<double>& refractive_index, const double& size_param,
        const Eigen::VectorXd& cos_angles, int size_index, MieData& output) {
        double& Qext = output.Qext(size_index);
        double& Qsca = output.Qsca(size_index);

        Eigen::MatrixXcd& S1 = output.S1;
        Eigen::MatrixXcd& S2 = output.S2;

        std::complex<double> m_2 = refractive_index * refractive_index;
        double x_2 = size_param * size_param;
        std::complex<double> j(0.0, 1.0);

        std::complex<double> N_1 =
            (1.0 - 0.1 * x_2 + (4.0 * m_2 + 5.0) * x_2 * x_2 / 1400.0);

        std::complex<double> D_1 =
            m_2 + 2.0 + (1.0 - 0.7 * m_2) * x_2 -
            (8.0 * m_2 * m_2 - 385.0 * m_2 + 350.0) * x_2 * x_2 / 1400.0;

        D_1 = D_1 + 2.0 * j * (m_2 - 1.0) * pow(size_param, 3) *
                        (1.0 - 0.1 * x_2) / 3.0;

        std::complex<double> a_hat1 = 2.0 * j * (m_2 - 1.0) / 3.0 * N_1 / D_1;

        std::complex<double> b_hat1 = (j * x_2 * (m_2 - 1.0) / 45.0 *
                                       (1.0 + (2.0 * m_2 - 5.0) / 70.0 * x_2)) /
                                      (1.0 - x_2 * (2.0 * m_2 - 5.0) / 30.0);

        std::complex<double> a_hat2 =
            (j * x_2 * (m_2 - 1.0) / 15.0 * (1.0 - x_2 / 14.0)) /
            (2.0 * m_2 + 3.0 - x_2 * (2.0 * m_2 - 7.0) / 14.0);

        double T = math::sqr(abs(a_hat1)) + math::sqr(abs(b_hat1)) +
                   5.0 / 3.0 * math::sqr(abs(a_hat2));

        Qsca = 6.0 * x_2 * x_2 * T;

        if (refractive_index.imag() == 0) {
            Qext = Qsca;
        } else {
            Qext = 6.0 * size_param *
                   ((a_hat1 + b_hat1).real() + 5.0 / 3.0 * a_hat2.real());
        }

        S1(size_index, Eigen::all).array() =
            3.0 / 2.0 * pow(size_param, 3.0) *
            (a_hat1 + cos_angles.array() * (b_hat1 + 5.0 / 3.0 * a_hat2));

        S2(size_index, Eigen::all).array() =
            3.0 / 2.0 * pow(size_param, 3.0) *
            (b_hat1 + a_hat1 * cos_angles.array() +
             5.0 / 3.0 * a_hat2 * (2.0 * pow(cos_angles.array(), 2) - 1.0));
    }

    void LinearizedMie::tau_pi(const Eigen::VectorXd& cos_angles) {
        Eigen::MatrixXd& tau_matrix = m_tau_matrix;
        Eigen::MatrixXd& pi_matrix = m_pi_matrix;
        int N = m_tau_matrix.rows();

        pi_matrix(0, Eigen::all).setConstant(1.0);
        tau_matrix(0, Eigen::all) = cos_angles.array();

        if (N > 1) {
            pi_matrix(1, Eigen::all) = 3.0 * cos_angles.array();
            tau_matrix(1, Eigen::all) = 2.0 * cos_angles.array().transpose() *
                                            pi_matrix(1, Eigen::all).array() -
                                        3.0;
        }

        for (int n = 3; n < N + 1; ++n) {
            pi_matrix(n - 1, Eigen::all) =
                ((2.0 * n - 1) * cos_angles.array().transpose() *
                     pi_matrix(n - 2, Eigen::all).array() -
                 double(n) * pi_matrix(n - 3, Eigen::all).array()) /
                (n - 1);

            tau_matrix(n - 1, Eigen::all) =
                double(n) * cos_angles.array().transpose() *
                    pi_matrix(n - 1, Eigen::all).array() -
                (n + 1) * pi_matrix(n - 2, Eigen::all).array();
        }
    }

    void LinearizedMie::allocate(double max_x, int num_angle) {
        int N = int(max_x + 4.05 * pow(max_x, 0.33333) + 2.0) + 2;

        m_thread_storage.clear();
        m_thread_storage.reserve(m_num_threads);
        for (int i = 0; i < m_num_threads; ++i) {
            m_thread_storage.emplace_back(N, m_tau_matrix, m_pi_matrix);
        }

        m_tau_matrix.resize(N, num_angle);
        m_pi_matrix.resize(N, num_angle);
    }

    void
    LinearizedMieWorker::An_Bn(const std::complex<double>& refractive_index,
                               const double& size_param) {
        int N = m_local_N;

        // Convenience accessors
        Eigen::VectorXcd& An = m_An;
        Eigen::VectorXcd& Bn = m_Bn;
        std::complex<double> j(0.0, 1.0);

        std::complex<double> psi_n_1 = sin(size_param); // this is n=0
        std::complex<double> psi_n =
            sin(size_param) / size_param - cos(size_param); // this is n=1
        std::complex<double> xi_n_1 =
            sin(size_param) + j * cos(size_param); // this is n=0
        std::complex<double> xi_n =
            sin(size_param) / size_param - cos(size_param) +
            j * (cos(size_param) / size_param + sin(size_param)); // this is n=1
        std::complex<double> temp;

        this->Dn(refractive_index, size_param);

        // Set the initial elements of An and Bn
        An(0) =
            ((m_Dn(0) / refractive_index + 1.0 / size_param) * psi_n -
             psi_n_1) /
            ((m_Dn(0) / refractive_index + 1.0 / size_param) * xi_n - xi_n_1);
        Bn(0) =
            ((m_Dn(0) * refractive_index + 1.0 / size_param) * psi_n -
             psi_n_1) /
            ((m_Dn(0) * refractive_index + 1.0 / size_param) * xi_n - xi_n_1);

        temp = xi_n_1;
        xi_n_1 = xi_n;
        xi_n = (2.0 + 1.0) / size_param * xi_n_1 - temp;

        psi_n_1 = psi_n;
        psi_n = xi_n.real();

        for (int n = 2; n < N + 1; ++n) {
            An(n - 1) =
                ((m_Dn(n - 1) / refractive_index + double(n) / size_param) *
                     psi_n -
                 psi_n_1) /
                ((m_Dn(n - 1) / refractive_index + double(n) / size_param) *
                     xi_n -
                 xi_n_1);
            Bn(n - 1) =
                ((m_Dn(n - 1) * refractive_index + double(n) / size_param) *
                     psi_n -
                 psi_n_1) /
                ((m_Dn(n - 1) * refractive_index + double(n) / size_param) *
                     xi_n -
                 xi_n_1);

            temp = xi_n_1;
            xi_n_1 = xi_n;
            xi_n = (double(2 * n) + 1.0) / size_param * xi_n_1 - temp;

            psi_n_1 = psi_n;
            psi_n = xi_n.real();
        }
    }

    void LinearizedMieWorker::Dn(const std::complex<double>& refractive_index,
                                 const double& size_param) {
        std::complex<double> z = size_param * refractive_index;

        if (refractive_index.real() < 1 || refractive_index.real() > 10 ||
            abs(refractive_index.imag()) > 10) {
            Dn_downwards(refractive_index, z);
        } else {
            double temp = 3.9 - 10.8 * refractive_index.real() +
                          13.78 * pow(refractive_index.real(), 2);
            if (abs(refractive_index.imag()) * size_param >= temp) {
                Dn_downwards(refractive_index, z);
            } else {
                Dn_upwards(refractive_index, z);
            }
        }
    }

    void LinearizedMieWorker::Dn_upwards(
        const std::complex<double>& refractive_index,
        const std::complex<double>& z) {
        Eigen::VectorXcd& Dn = m_Dn;
        int N = m_local_N;

        std::complex<double> j(0.0, 1.0);
        std::complex<double> temp = -2.0 * j * z;
        std::complex<double> exp_2 = exp(temp);

        Dn(0) =
            -1.0 / z + (1.0 - exp_2) / ((1.0 - exp_2) / z - j * (1.0 + exp_2));

        for (int n = 2; n < N + 1; ++n) {
            Dn(n - 1) = -double(n) / z + 1.0 / (double(n) / z - Dn(n - 2));
        }
    }

    void LinearizedMieWorker::Dn_downwards(
        const std::complex<double>& refractive_index,
        const std::complex<double>& z) {
        int N = m_local_N;

        Dn_Lentz(z, m_Dn(N - 1));

        for (int n = N; n > 1; --n) {
            m_Dn[n - 2] = double(n) / z - 1.0 / (double(n) / z + m_Dn[n - 1]);
        }
    }

    void LinearizedMieWorker::Dn_Lentz(const std::complex<double>& z,
                                       std::complex<double>& result) {
        int N = m_local_N;

        std::complex<double> z_inv = 2.0 / z;

        std::complex<double> alpha_1 = (N + 0.5) * z_inv;
        std::complex<double> a_j = -(N + 1.5) * z_inv;
        std::complex<double> alpha_j1 = a_j + 1.0 / alpha_1;
        std::complex<double> alpha_j2 = a_j;

        std::complex<double> cur_ratio = alpha_j1 / alpha_j2;
        std::complex<double> overall_ratio = alpha_1 * cur_ratio;

        while (abs(abs(cur_ratio) - 1) > 1e-12) {
            a_j = z_inv - a_j;
            alpha_j1 = a_j + 1.0 / alpha_j1;
            alpha_j2 = a_j + 1.0 / alpha_j2;

            z_inv *= -1;
            cur_ratio = alpha_j1 / alpha_j2;
            overall_ratio *= cur_ratio;
        }
        result = -double(N) / z + overall_ratio;
    }

} // namespace sasktran2::mie
