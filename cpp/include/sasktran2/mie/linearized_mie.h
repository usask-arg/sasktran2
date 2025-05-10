#pragma once

#include "rust/cxx.h"
#include "sasktran2/mie/mie.h"
#include <sasktran2/internal_common.h>
#include <sasktran2_core_cxx/mie/mod.h>

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
        LinearizedMieWorker(int max_N, const Eigen::MatrixXd& tau,
                            const Eigen::MatrixXd& pi)
            : m_local_N(max_N), m_tau_matrix(tau), m_pi_matrix(pi) {
            m_An.resize(max_N);
            m_Bn.resize(max_N);
            m_Dn.resize(max_N);
        }

        void set_local_N(double size_param) {
            m_local_N =
                int(size_param + 4.05 * pow(size_param, 0.33333) + 2.0) + 2;

            if (m_local_N > m_An.size()) {
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

        void regular_Q_S_scat(const std::complex<double>& refractive_index,
                              const double& size_param, int size_index,
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

    class RustMie : public MieBase {
      private:
        int m_num_threads;

        void internal_calculate(const Eigen::VectorXd& size_param,
                                const std::complex<double>& refractive_index,
                                const Eigen::VectorXd& cos_angles,
                                bool calculate_derivative,
                                MieOutput& output) override {
            ::rust::Slice<const rust::f64> sp(size_param.data(),
                                              size_param.size());
            ::rust::Slice<const rust::f64> ca(cos_angles.data(),
                                              cos_angles.size());

            auto rust_output = ffi::mie_c(sp, refractive_index.real(),
                                          refractive_index.imag(), ca);

            auto qext = rust_output->Qext_vec();
            auto qsca = rust_output->Qsca_vec();

            auto s1_im = rust_output->S1_im_vec();
            auto s1_re = rust_output->S1_re_vec();
            auto s2_im = rust_output->S2_im_vec();
            auto s2_re = rust_output->S2_re_vec();

            int c = 0;
            for (int i = 0; i < qext.size(); ++i) {
                output.values.Qext(i) = qext[i];
                output.values.Qsca(i) = qsca[i];

                for (int j = 0; j < cos_angles.size(); ++j) {
                    output.values.S1(i, j) =
                        std::complex<double>(s1_re[c], s1_im[c]);
                    output.values.S2(i, j) =
                        std::complex<double>(s2_re[c], s2_im[c]);
                    ++c;
                }
            }
        }

      public:
        RustMie(int num_threads = 1) : m_num_threads(num_threads){};
    };
} // namespace sasktran2::mie
