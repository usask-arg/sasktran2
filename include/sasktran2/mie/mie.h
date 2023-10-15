#pragma once

#include <sasktran2/internal_common.h>

namespace sasktran2::mie {
    struct MieData {
        Eigen::VectorXd Qext; // Extinction efficiency [nsize]
        Eigen::VectorXd Qsca; // Scattering efficiency [nsize]

        Eigen::VectorXd GQsc; // Asymmetry parameter [nsize]
        
        Eigen::VectorXcd SFORW; 
        Eigen::VectorXcd SBACK;

        Eigen::MatrixXcd TFORW;
        Eigen::MatrixXcd TBACK;

        Eigen::MatrixXcd S1; // Scattering matrix [nsize, nangles]
        Eigen::MatrixXcd S2; // Scattering matrix [nsize, nangles]

        void resize(int nsize, int nangles) {
            Qext.resize(nsize);
            Qsca.resize(nsize);
            GQsc.resize(nsize);
            SFORW.resize(nsize);
            SBACK.resize(nsize);
            TFORW.resize(nsize, 2);
            TBACK.resize(nsize, 2);
            S1.resize(nsize, nangles);
            S2.resize(nsize, nangles);
        }
    };

    struct MieOutput {
        Eigen::VectorXd size_param;
        Eigen::VectorXcd refractive_index;
        Eigen::VectorXd angles;

        MieData values;
        MieData d_by_size_param;
        MieData d_by_refrac_real;
        MieData d_by_refrac_imag;
    };


    class MieBase {
        private:

        public:
            virtual ~MieBase() {};

            virtual void calculate(const Eigen::VectorXd& size_param,
                                   const Eigen::VectorXcd& refractive_index,
                                   const Eigen::VectorXd& angles,
                                   MieOutput& output) = 0;
    };
}