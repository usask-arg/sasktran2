#pragma once

#include <sasktran2/internal_common.h>

namespace sasktran2::mie {
    struct MieData {
        Eigen::VectorXd Qext; // Extinction efficiency [nsize]
        Eigen::VectorXd Qsca; // Scattering efficiency [nsize]

        Eigen::MatrixXcd S1; // Scattering amplitude1 [nsize, nangles]
        Eigen::MatrixXcd S2; // Scattering amplitude2 [nsize, nangles]

        void resize(int nsize, int nangles) {
            Qext.resize(nsize);
            Qsca.resize(nsize);
            S1.resize(nsize, nangles);
            S2.resize(nsize, nangles);
        }
    };

    struct MieOutput {
        Eigen::VectorXd size_param; /** copy of the input size parameters */
        std::complex<double>
            refractive_index;       /** copy of the input refractive index */
        Eigen::VectorXd cos_angles; /** copy of the input cos angles */

        MieData values;           /** Mie values */
        MieData d_by_size_param;  /** Derivative of Mie values by the size
                                     parameter */
        MieData d_by_refrac_real; /** Derivative of Mie values by the real part
                                     of the refractive index */
        MieData d_by_refrac_imag; /** Derivative of Mie values by the imaginary
                                     part of the refractive index */
    };

    class MieBase {
      private:
        MieOutput m_output;

        virtual void
        internal_calculate(const Eigen::VectorXd& size_param,
                           const std::complex<double>& refractive_index,
                           const Eigen::VectorXd& cos_angles,
                           bool calculate_derivative, MieOutput& output) = 0;

      public:
        virtual ~MieBase(){};

        /**
         * Performs the Mie computation for an array of size parameters, a
         * single refractive index
         *
         * @param size_param Mie size parameter
         * @param refractive_index Mie refractive index
         * @param cos_angles Cosine of angles to calculate the scattering
         * amplitide at
         * @param calculate_derivative set to true to also calculate derivatives
         * of size parameter and refractive index
         * @param output class to store output in, will be resized
         */
        const MieOutput& calculate(const Eigen::VectorXd& size_param,
                                   const std::complex<double>& refractive_index,
                                   const Eigen::VectorXd& cos_angles,
                                   bool calculate_derivative = false) {

            m_output.size_param = size_param;
            m_output.refractive_index = refractive_index;
            m_output.cos_angles = cos_angles;
            m_output.values.resize(size_param.size(), cos_angles.size());

            internal_calculate(size_param, refractive_index, cos_angles,
                               calculate_derivative, m_output);
            return m_output;
        }
    };
} // namespace sasktran2::mie
