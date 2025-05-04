#pragma once

#include "sasktran2/math/wigner.h"
#include "sasktran2/mie/linearized_mie.h"
#include <sasktran2/internal_common.h>

namespace sasktran2::mie {

    /**
     *  A class which performs integration of the Mie scattering parameters over
     * a size distribution, and also computes the greek coefficients of the
     * scattering phase function. This is mostly a helper class, the Python
     * interface computes quadrature weights and the integration is done in this
     * class.
     */
    class MieIntegrator {
      private:
        Eigen::Ref<const Eigen::VectorXd>
            m_cos_angles; // Cosine of angles for the scattering phase function

        int m_num_legendre; // Number of legendre moments to calculate
        int m_num_threads;  // Number of threads to use for the integration

        // Wigner functions
        Eigen::MatrixXd m_lpoly_00;
        Eigen::MatrixXd m_lpoly_22;
        Eigen::MatrixXd m_lpoly_2m2;
        Eigen::MatrixXd m_lpoly_02;

        // Mie code to use in the calculation
        LinearizedMie m_mie_calculator;

      public:
        MieIntegrator(Eigen::Ref<const Eigen::VectorXd> cos_angles,
                      int num_legendre, int num_threads)
            : m_cos_angles(cos_angles), m_num_legendre(num_legendre),
              m_mie_calculator(num_threads), m_num_threads(num_threads) {
            generate_lpoly();
        };
        ~MieIntegrator(){};

        /**
         *  Generates the Wigner functions for the given cosine angles, called
         * by the constructor
         *
         */
        void generate_lpoly() {
            auto wigner00 = sasktran2::math::WignerDCalculator(0, 0);
            auto wigner22 = sasktran2::math::WignerDCalculator(2, 2);
            auto wigner2m2 = sasktran2::math::WignerDCalculator(2, -2);
            auto wigner02 = sasktran2::math::WignerDCalculator(0, 2);

            m_lpoly_00.resize(m_num_legendre, m_cos_angles.size());
            m_lpoly_22.resize(m_num_legendre, m_cos_angles.size());
            m_lpoly_2m2.resize(m_num_legendre, m_cos_angles.size());
            m_lpoly_02.resize(m_num_legendre, m_cos_angles.size());

            for (size_t i = 0; i < m_cos_angles.size(); ++i) {
                double theta = acos(m_cos_angles(i));
                wigner00.vector_d(theta, m_lpoly_00.col(i));
                wigner22.vector_d(theta, m_lpoly_22.col(i));
                wigner2m2.vector_d(theta, m_lpoly_2m2.col(i));
                wigner02.vector_d(theta, m_lpoly_02.col(i));
            }
        }

        void determine_cos_angles(
            double wavelength, std::complex<double> refractive_index,
            Eigen::Ref<const Eigen::VectorXd> size_param, // [size_param]
            Eigen::Ref<const Eigen::ArrayXXd> pdf, // [size_param, distribution]
            Eigen::Ref<const Eigen::VectorXd> size_weights,  // [size_param]
            Eigen::Ref<const Eigen::VectorXd> angle_weights) // [angle])
        {
            std::vector<double> cos_angles;
            std::vector<double> p11_scratch;
            Eigen::Map<Eigen::VectorXd> cos_angles_map(cos_angles.data(),
                                                       cos_angles.size());

            m_mie_calculator.calculate(size_param, refractive_index,
                                       cos_angles_map);

            for (size_t i = 0; i < m_cos_angles.size(); ++i) {
                double theta = acos(m_cos_angles(i));
                cos_angles.push_back(cos(theta));
            }

            for (size_t i = 1; i < cos_angles.size() - 1; ++i) {
                // Interpolate to the point from the left and right points
                double left = cos_angles[i - 1];
                double right = cos_angles[i + 1];

                double lval = p11_scratch[i - 1];
                double rval = p11_scratch[i + 1];

                double interp_point = cos_angles[i];
                // Interp to interp_point

                double interp_p11 = lval + (rval - lval) *
                                               (interp_point - left) /
                                               (right - left);

                if (abs(interp_p11 - p11_scratch[i]) > 1e-6) {
                    // Split intervals on either side of the point
                    double left_split = (left + interp_point) / 2.0;
                    double right_split = (interp_point + right) / 2.0;

                    cos_angles.insert(cos_angles.begin() + i, left_split);
                    cos_angles.insert(cos_angles.begin() + i + 1, right_split);

                    // Calculate p11 at the new split points and insert

                    // now we want to check the interpolation error at the left
                    // split, so we need to decrement i
                    --i;
                }
            }
        }

        /**
         * @brief Integrates over the size distribution, emplacing the results
         *
         * Integration is done for a single wavelength, but multiple size
         * distributions
         *
         * @param wavelength
         * @param refractive_index
         * @param size_param
         * @param pdf
         * @param size_weights
         * @param angle_weights
         * @param xs_total
         * @param xs_scattering
         * @param p11
         * @param p12
         * @param p33
         * @param p34
         * @param lm_a1
         * @param lm_a2
         * @param lm_a3
         * @param lm_a4
         * @param lm_b1
         * @param lm_b2
         */
        void integrate(
            double wavelength, std::complex<double> refractive_index,
            Eigen::Ref<const Eigen::VectorXd> size_param, // [size_param]
            Eigen::Ref<const Eigen::ArrayXXd> pdf, // [size_param, distribution]
            Eigen::Ref<const Eigen::VectorXd> size_weights,  // [size_param]
            Eigen::Ref<const Eigen::VectorXd> angle_weights, // [angle]
            Eigen::Ref<Eigen::ArrayXd> xs_total,             // [distribution]
            Eigen::Ref<Eigen::ArrayXd> xs_scattering,        // [distribution]
            Eigen::Ref<Eigen::ArrayXXd> p11,   // [angle, distribution]
            Eigen::Ref<Eigen::ArrayXXd> p12,   // [angle, distribution]
            Eigen::Ref<Eigen::ArrayXXd> p33,   // [angle, distribution]
            Eigen::Ref<Eigen::ArrayXXd> p34,   // [angle, distribution]
            Eigen::Ref<Eigen::ArrayXXd> lm_a1, // [legendre, distribution]
            Eigen::Ref<Eigen::ArrayXXd> lm_a2, // [legendre, distribution]
            Eigen::Ref<Eigen::ArrayXXd> lm_a3, // [legendre, distribution]
            Eigen::Ref<Eigen::ArrayXXd> lm_a4, // [legendre, distribution]
            Eigen::Ref<Eigen::ArrayXXd> lm_b1, // [legendre, distribution]
            Eigen::Ref<Eigen::ArrayXXd> lm_b2  // [legendre, distribution]
        ) {
            const auto& result = m_mie_calculator.calculate(
                size_param, refractive_index, m_cos_angles);

            const auto& S1 = result.values.S1;
            const auto& S2 = result.values.S2;

            const auto& Qext = result.values.Qext.array();
            const auto& Qsca = result.values.Qsca.array();

            double k = 2 * EIGEN_PI / wavelength;
            double c = 4 * EIGEN_PI / (2 * k * k);

// For each distribution
#pragma omp parallel for num_threads(m_num_threads)
            for (int i = 0; i < xs_total.size(); ++i) {

                xs_total(i) = (size_weights.array() *
                               pdf(Eigen::all, i).array() * EIGEN_PI * Qext *
                               (size_param * wavelength / (2 * EIGEN_PI))
                                   .array()
                                   .square())
                                  .sum();
                xs_scattering(i) =
                    (size_weights.array() * pdf(Eigen::all, i).array() *
                     EIGEN_PI * Qsca *
                     (size_param * wavelength / (2 * EIGEN_PI))
                         .array()
                         .square())
                        .sum();

                p11.col(i).array() =
                    ((c / xs_scattering(i) * size_weights.array() *
                      pdf(Eigen::all, i).array())
                         .matrix()
                         .transpose() *
                     (S1.cwiseAbs2() + S2.cwiseAbs2()))
                        .array();

                p12.col(i).array() =
                    (c / xs_scattering(i) * size_weights.array() *
                     pdf(Eigen::all, i).array())
                        .matrix()
                        .transpose() *
                    (S1.cwiseAbs2() - S2.cwiseAbs2());

                p33.col(i).array() =
                    ((c / xs_scattering(i) * size_weights.array() *
                      pdf(Eigen::all, i).array())
                         .matrix()
                         .transpose() *
                     (S1.cwiseProduct(S2.conjugate()) +
                      S2.cwiseProduct(S1.conjugate())))
                        .array()
                        .real();

                p34.col(i).array() =
                    ((c / xs_scattering(i) * size_weights.array() *
                      pdf(Eigen::all, i).array())
                         .matrix()
                         .transpose() *
                     (S1.cwiseProduct(S2.conjugate()) -
                      S2.cwiseProduct(S1.conjugate())))
                        .array()
                        .imag();

                for (int l = 0; l < m_num_legendre; ++l) {
                    double l_weight = 1.0 / (2.0 / (2.0 * l + 1.0));
                    for (int j = 0; j < m_cos_angles.size(); ++j) {
                        double w = angle_weights(j) * l_weight;

                        lm_a1(l, i) += w * m_lpoly_00(l, j) * p11(j, i);

                        double temp1 =
                            w * m_lpoly_22(l, j) * (p11(j, i) + p33(j, i));
                        double temp2 =
                            w * m_lpoly_2m2(l, j) * (p11(j, i) - p33(j, i));

                        lm_a2(l, i) += (temp1 + temp2) / 2.0;
                        lm_a3(l, i) += (temp1 - temp2) / 2.0;

                        lm_b1(l, i) += w * m_lpoly_02(l, j) * p12(j, i);
                        lm_b2(l, i) +=
                            (-1.0) * w * m_lpoly_02(l, j) * p34(j, i);
                    }
                }
            }
        }
    };

} // namespace sasktran2::mie
