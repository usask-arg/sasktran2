#include "sasktran2/math/wigner.h"
#include "sasktran2-core/src/math/wigner.rs.h"

namespace sasktran2::math {

#ifdef SKTRAN_RUST_SUPPORT
    WignerDCalculator::WignerDCalculator(int m, int n)
        : m_calculator(sasktran2::rust::math::new_wigner_d_calculator(m, n)) {}

    double WignerDCalculator::d(double theta, int l) {
        return m_calculator->d(theta, l);
    }

#endif
} // namespace sasktran2::math
