#pragma once

#include "sasktran2/mie/mie.h"
#include "sasktran2/rext.h"
#include <sasktran2/internal_common.h>

namespace sasktran2::mie {

    class RustMie : public MieBase {
      private:
        std::unique_ptr<sasktran2::rext::mie::Mie,
                        decltype(&sasktran2::rext::mie::mie_free)>
            m_mie;

        sasktran2::rext::mie::Mie* m_mie_ptr = nullptr;

        void internal_calculate(const Eigen::VectorXd& size_param,
                                const std::complex<double>& refractive_index,
                                const Eigen::VectorXd& cos_angles,
                                bool calculate_derivative,
                                MieOutput& output) override {
            for (size_t i = 0; i < output.size_param.size(); ++i) {
                double size_param = output.size_param[i];
                double refractive_real = refractive_index.real();
                double refractive_imag = refractive_index.imag();

                sasktran2::rext::mie::mie_calculate(
                    m_mie.get(), size_param, refractive_real, refractive_imag,
                    &output.values.Qext[i], &output.values.Qsca[i]);
            }
        }

      public:
        RustMie(int num_threads = 1)
            : m_mie(sasktran2::rext::mie::mie_new(),
                    sasktran2::rext::mie::mie_free) {
            if (!m_mie) {
                throw std::runtime_error("Failed to create Mie instance");
            }
        }
        ~RustMie(){};
    };

} // namespace sasktran2::mie
