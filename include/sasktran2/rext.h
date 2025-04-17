#pragma once

#include <cstddef>

namespace sasktran2::rext {
    namespace mie {
        extern "C" {

        typedef struct Mie Mie;

        Mie* mie_new(void);

        void mie_with_cos_angles(Mie* mie_ptr, const double* angles,
                                 size_t len);

        void mie_calculate(Mie* mie_ptr, double size_param,
                           double refractive_real, double refractive_imag,
                           double* Qext, double* Qsca);

        void mie_free(Mie* mie_ptr);

#ifdef __cplusplus
        }
#endif
    } // namespace mie

} // namespace sasktran2::rext
