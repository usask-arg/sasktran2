#include "sasktran2/math/wigner.h"
#include "sasktran2-core/src/math/wigner.rs.h"

namespace sasktran2::math {

#ifdef SKTRAN_RUST_SUPPORT
    WignerDCalculator::WignerDCalculator(int m, int n)
        : m_calculator(sasktran2::rust::math::new_wigner_d_calculator(m, n)) {}

    double WignerDCalculator::d(double theta, int l) {
        return m_calculator->d(theta, l);
    }

    void WignerDCalculator::vec_d_emplace(double theta, int max_l,
                                          double* out_array) {
        ::rust::Slice<double> slice(out_array, max_l);

        m_calculator->vector_d(theta, slice);
    }

#endif
} // namespace sasktran2::math
