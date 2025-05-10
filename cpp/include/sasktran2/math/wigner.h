#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2_core_cxx/math/wigner.h>

namespace sasktran2::math {
    /** Calculates the Wigner D functions using recurrence relations found from
     * Mischenko.  The "first" Wigner D function is the standard Legendre
     * polynomial, so this class is used to calculate Legendre polynomials as
     * well. The other Wigner D functions are only used for polarized scattering
     * calculations.
     *
     */
    class WignerDCalculator {
      private:
        ::rust::Box<ffi::WignerDCalculator> m_internal;

      public:
        /** Constructs the calculator for \f$d^l_{mn}\f$ for a given m and n
         *
         * @param m
         * @param n
         */
        WignerDCalculator(int m, int n)
            : m_internal(ffi::make_wigner_calculator(m, n)) {}
        /** Calculates \f$d^l_{mn}(\theta)\f$
         *
         * @param theta Angle in radians
         * @param l
         * @return \f$d^l_{mn}(\theta)\f$
         */
        double d(double theta, int l) { return m_internal->d(theta, l); }
    };
} // namespace sasktran2::math
