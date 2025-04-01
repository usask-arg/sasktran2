#include "sasktran2/math/wigner.h"
#include "sasktran2/mie/mie.h"
#include "sasktran2/mie/linearized_mie.h"
#include "sasktran2/mie/distribution.h"
#include <limits>
#include <sasktran2/internal_common.h>

namespace sasktran2::mie {

    class MieIntegrator {
        private:
            Eigen::Ref<const Eigen::VectorXd> m_cos_angles;

            int m_num_legendre;
            int m_num_threads;

            Eigen::MatrixXd m_lpoly_00;
            Eigen::MatrixXd m_lpoly_22;
            Eigen::MatrixXd m_lpoly_2m2;
            Eigen::MatrixXd m_lpoly_02;

            LinearizedMie m_mie_calculator;

            Eigen::VectorXd m_size_buffer;



        public:
            MieIntegrator(Eigen::Ref<const Eigen::VectorXd> cos_angles, int num_legendre, int num_threads) : m_cos_angles(cos_angles), m_num_legendre(num_legendre), m_mie_calculator(num_threads), m_num_threads(num_threads) {
                generate_lpoly();
                m_size_buffer.resize(num_threads);
            };
            ~MieIntegrator() {};

            void generate_lpoly() {
                auto wigner00 = sasktran2::math::WignerDCalculator(0, 0);
                auto wigner22 = sasktran2::math::WignerDCalculator(2, 2);
                auto wigner2m2 = sasktran2::math::WignerDCalculator(2, -2);
                auto wigner02 = sasktran2::math::WignerDCalculator(0, 2);

                m_lpoly_00.resize(m_num_legendre, m_cos_angles.size());
                m_lpoly_22.resize(m_num_legendre, m_cos_angles.size());
                m_lpoly_2m2.resize(m_num_legendre, m_cos_angles.size());
                m_lpoly_02.resize(m_num_legendre, m_cos_angles.size());

                for(size_t i = 0; i < m_cos_angles.size(); ++i) {
                    double theta = acos(m_cos_angles(i));
                    wigner00.vector_d(theta, m_lpoly_00.col(i));
                    wigner22.vector_d(theta, m_lpoly_22.col(i));
                    wigner2m2.vector_d(theta, m_lpoly_2m2.col(i));
                    wigner02.vector_d(theta, m_lpoly_02.col(i));
                }
            }

            void integrate_all(
                double wavelength,
                std::complex<double> refractive_index,
                Eigen::Ref<const Eigen::VectorXd> size_param, // [size_param]
                Eigen::Ref<const Eigen::ArrayXXd> pdf, // [size_param, distribution]
                Eigen::Ref<const Eigen::VectorXd> size_weights, // [size_param]
                Eigen::Ref<const Eigen::VectorXd> angle_weights, // [angle]
                Eigen::Ref<Eigen::ArrayXd> xs_total, // [distribution]
                Eigen::Ref<Eigen::ArrayXd> xs_scattering, // [distribution]
                Eigen::Ref<Eigen::ArrayXXd> p11, // [angle, distribution]
                Eigen::Ref<Eigen::ArrayXXd> p12, // [angle, distribution]
                Eigen::Ref<Eigen::ArrayXXd> p33, // [angle, distribution]
                Eigen::Ref<Eigen::ArrayXXd> p34, // [angle, distribution]
                Eigen::Ref<Eigen::ArrayXXd> lm_a1, // [legendre, distribution]
                Eigen::Ref<Eigen::ArrayXXd> lm_a2, // [legendre, distribution]
                Eigen::Ref<Eigen::ArrayXXd> lm_a3, // [legendre, distribution]
                Eigen::Ref<Eigen::ArrayXXd> lm_a4, // [legendre, distribution]
                Eigen::Ref<Eigen::ArrayXXd> lm_b1, // [legendre, distribution]
                Eigen::Ref<Eigen::ArrayXXd> lm_b2 // [legendre, distribution]
            ) {
                std::cout << "MIE CALC START" << std::endl;
                const auto& result = m_mie_calculator.calculate(size_param, refractive_index, m_cos_angles);
                std::cout << "MIE CALC END" << std::endl;

                const auto& S1 = result.values.S1;
                const auto& S2 = result.values.S2;

                const auto& Qext = result.values.Qext.array();
                const auto& Qsca = result.values.Qsca.array();

                double k = 2 * EIGEN_PI / wavelength;
                double c = 4 * EIGEN_PI / (2 * k*k);
                

                // For each distribution
                #pragma omp parallel for num_threads(m_num_threads)
                for(int i = 0; i < xs_total.size(); ++i) {

                    xs_total(i) = (size_weights.array() * pdf(Eigen::all, i).array() * EIGEN_PI * Qext * (size_param * wavelength / (2*EIGEN_PI)).array().square()).sum();
                    xs_scattering(i) = (size_weights.array() * pdf(Eigen::all, i).array() * EIGEN_PI * Qsca * (size_param * wavelength / (2*EIGEN_PI)).array().square()).sum();

                    p11.col(i).array() = ((c / xs_scattering(i) * size_weights.array() * pdf(Eigen::all, i).array()).matrix().transpose()
                        * (S1.cwiseAbs2() + S2.cwiseAbs2())).array();

                    p12.col(i).array() = (c / xs_scattering(i) * size_weights.array() * pdf(Eigen::all, i).array()).matrix().transpose()
                        * (S1.cwiseAbs2() - S2.cwiseAbs2());

                    p33.col(i).array() = ((c / xs_scattering(i) * size_weights.array() * pdf(Eigen::all, i).array()).matrix().transpose()
                        * (S1.cwiseProduct(S2.conjugate()) + S2.cwiseProduct(S1.conjugate()))).array().real();

                    p34.col(i).array() = ((c / xs_scattering(i) * size_weights.array() * pdf(Eigen::all, i).array()).matrix().transpose()
                        * (S1.cwiseProduct(S2.conjugate()) - S2.cwiseProduct(S1.conjugate()))).array().imag();
                    
                    for(int l = 0; l < m_num_legendre; ++l) {
                        double l_weight = 1.0 / (2.0 / (2.0 * l + 1.0));
                        for(int j = 0; j < m_cos_angles.size(); ++j) {
                            double w = angle_weights(j) * l_weight;

                            lm_a1(l, i) += w * m_lpoly_00(l, j) * p11(j, i);

                            double temp1 = w * m_lpoly_22(l, j) * (p11(j, i) + p33(j, i));
                            double temp2 = w * m_lpoly_2m2(l, j) * (p11(j, i) - p33(j, i));

                            lm_a2(l, i) += (temp1 + temp2) / 2.0;
                            lm_a3(l, i) += (temp1 - temp2) / 2.0;

                            lm_b1(l, i) += w * m_lpoly_02(l, j) * p12(j, i);
                            lm_b2(l, i) += (-1.0) * w * m_lpoly_02(l, j) * p34(j, i);
                        }
                    }
                }
            }

            std::tuple<double, double> integrate(const sasktran2::mie::MieOutput& mie_output,
                           Eigen::Ref<const Eigen::VectorXd> radii, // [size]
                           Eigen::Ref<const Eigen::VectorXd> quadrature_weights, // [size]
                           Eigen::Ref<Eigen::VectorXd> p11, // [angle]
                           Eigen::Ref<Eigen::VectorXd> p12, // [angle]
                           Eigen::Ref<Eigen::VectorXd> p33, // [angle]
                           Eigen::Ref<Eigen::VectorXd> p34, // [angle]
                           double wavelength
                        ) {
                // Efficiency assume all inputs have been zerod

                // Convenvience
                const auto& S1 = mie_output.values.S1;
                const auto& S2 = mie_output.values.S2;

                const auto& Qext = mie_output.values.Qext.array();
                const auto& Qsca = mie_output.values.Qsca.array();
                const auto& w = quadrature_weights;
                const auto& x = mie_output.size_param.array();
                const auto& r = radii.array();

                double xs_total = (w.array() * EIGEN_PI * Qext * r.square()).sum();
                double xs_scattering = (w.array() * EIGEN_PI * Qsca * r.square()).sum();

                double k = 2 * EIGEN_PI / wavelength;
                double c = 4 * EIGEN_PI / (2 * k*k * xs_scattering);

                p11.array() = (c * w.transpose() * (S1.cwiseAbs2() + S2.cwiseAbs2())).array();
                p12.array() = (c * w.transpose() * (S1.cwiseAbs2() - S2.cwiseAbs2())).array();
                p33.array() = (c * w.transpose() * (S1.cwiseProduct(S2.conjugate()) + S2.cwiseProduct(S1.conjugate()))).array().real();
                p34.array() = (c * w.transpose() * (S1.cwiseProduct(S2.conjugate()) - S2.cwiseProduct(S1.conjugate()))).array().imag();

                return std::make_tuple(xs_total, xs_scattering);
            }
    };

}