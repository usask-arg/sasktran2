#pragma once
#include "sktran_disco/sktran_do.h"
#include <sasktran2/math/wigner.h>

namespace sasktran_disco {

    /**
     * Many vector equations use a matrix D = diag(1, 1, -1, -1) that accounts
     * for the symmetry of the Legendre functions.  This implements the diagonal
     * elements of this matrix, templated to not slow down the scalar code
     *
     * @tparam NSTOKES
     * @param linear_index
     * @return double
     */
    template <int NSTOKES> double stokes_negation_factor(int linear_index) {
        int stokes_index = linear_index % NSTOKES;

        if (stokes_index >= 2) {
            return -1.0;
        } else {
            return 1.0;
        }
    }

    // Type to store a single Legendre coefficient value
    template <int NSTOKES> struct LegendreCoefficient {};

    /**
     * @brief NSTOKES = 1 just has one Legendre coefficient
     *
     * @tparam
     */
    template <> struct LegendreCoefficient<1> {
        double a1;

        LegendreCoefficient() { a1 = 0.0; }
    };

    /**
     * @brief NSTOKES = 4 has 6 Legendre coefficients
     *
     * @tparam
     */
    template <> struct LegendreCoefficient<4> {
        double a1;
        double a2;
        double a3;
        double a4;
        double b1;
        double b2;

        LegendreCoefficient() {
            a1 = 0.0;
            a2 = 0.0;
            a3 = 0.0;
            a4 = 0.0;
            b1 = 0.0;
            b2 = 0.0;
        }
    };

    /**
     * @brief NSTOKES=3 we just have 3 Legendre Coefficients
     *
     * @tparam
     */
    template <> struct LegendreCoefficient<3> {
        double a1;
        double a2;
        double a3;
        double b1;

        LegendreCoefficient() {
            a1 = 0.0;
            a2 = 0.0;
            a3 = 0.0;
            b1 = 0.0;
        }
    };

    /**
     * Class to store the Legendre functions evaluated at stream angles, see eq
     * A.6 in Rozanov et. al 2013 Storage is 3 elements, [P R T] in the general
     * case, we need these three elements for NSTOKES=3 as well so we leave this
     * as the full general case
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> struct LegendrePhaseContainer {
        Eigen::Vector3d values;

        void fill(int m, int l, double coszen) {
            auto calculatorneg = sasktran2::math::WignerDCalculator(m, -2);
            auto calculatorpos = sasktran2::math::WignerDCalculator(m, 2);
            double theta = acos(coszen);

            values(1) =
                -0.5 * (calculatorpos.d(theta, l) + calculatorneg.d(theta, l));
            values(2) =
                -0.5 * (calculatorpos.d(theta, l) - calculatorneg.d(theta, l));

            auto calculator = sasktran2::math::WignerDCalculator(m, 0);

            values(0) = calculator.d(theta, l);
        }

        double& P() { return values(0); }

        double& R() { return values(1); }

        double& T() { return values(2); }

        double P() const { return values(0); }

        double R() const { return values(1); }

        double T() const { return values(2); }
    };

    /**
     * @brief Special case for NSTOKES=1
     *
     * @tparam
     */
    template <> struct LegendrePhaseContainer<1> {
        double value;

        void fill(int m, int l, double coszen) {
            auto calculator = sasktran2::math::WignerDCalculator(m, 0);
            double theta = acos(coszen);

            value = calculator.d(theta, l);
        }
        double& P() { return value; }

        double P() const { return value; }
    };

    template <int NSTOKES, int CNSTR = -1> class LayerInputDerivative;

    /**
     * Holds the quantities calculated by the LPETripleProduct class and their
     * derivatives The value is a NSTOKES x NSTOKES matrix, where each value in
     * theory has derivatives with respect to every greek parameter.  For
     * NSTOKES=4 this would result in a 4x4 matrix where we need derivatives
     * with respect to 6 quantities for every l.  In reality the derivatives are
     * sparse, alpha1 only affects the 1,1 element of the matrix for example.
     * This class stores the NSTOKESxNOSTKES matrix values efficiently and the
     * capability to propagate derivatives of the greek parameters efficiently
     * which requires specialized instantiations for every value of NSTOKES
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES> class TripleProductDerivativeHolder {};

    /**
     * @brief Special case for NSTOKES=1
     *
     * @tparam
     */
    template <> class TripleProductDerivativeHolder<1> {
      public:
        /**
         * @brief Construct a new Triple Product Derivative Holder object
         * without allocating
         *
         */
        TripleProductDerivativeHolder() {}

        /**
         * @brief Construct a new Triple Product Derivative Holder object and
         * allocate storage
         *
         * @param nstr
         */
        TripleProductDerivativeHolder(int nstr) : nstr(nstr) { resize(nstr); }

        /**
         * @brief Allocates storage for a given nstr
         *
         * @param nstr
         */
        void resize(int nstr) {
            this->nstr = nstr;
            d_by_legendre_coeff.resize(nstr);
        }

        /**
         * @brief Calculates the triple product
         *
         * @param coeffs
         * @param lp1s
         * @param lp2s
         * @param negation
         * @param m
         */
        void calculate(const std::vector<LegendreCoefficient<1>>& coeffs,
                       const std::vector<LegendrePhaseContainer<1>>& lp1s,
                       const std::vector<LegendrePhaseContainer<1>>& lp2s,
                       bool negation, int m) {
            value = 0.0;
            d_by_legendre_coeff.setZero();
            for (int l = m; l < nstr; ++l) {
                const auto& coeff = coeffs[l];
                const auto& lp1 = lp1s[l];
                const auto& lp2 = lp2s[l];

                int negation_factor = 1;
                if (negation) {
                    if ((l - m) % 2 != 0) {
                        negation_factor *= -1;
                    }
                }
                value += coeff.a1 * lp1.P() * lp2.P() * negation_factor;
                d_by_legendre_coeff[l] += lp1.P() * lp2.P() * negation_factor;
            }
        }

        double value; /** Value of the triple product */
        double ssa;   /** Single scatter albedo */
        int nstr;     /** Number of streams */
        Eigen::VectorXd
            d_by_legendre_coeff; /** Derivative of value with respect to each
                                    Legendre coefficient */

        /**
         * @brief Reduce a user input derivative to derivative of the triple
         * product
         *
         * @param layer_deriv
         * @param deriv
         */
        void reduce(const LayerInputDerivative<1>& layer_deriv,
                    double& deriv) const;
    };

    /**
     * The values are a 4x4 matrix, the derivative factors for the 6 greek
     * constants are a1 - 11 a2 - 22, 23, 32, 33 a3 - 22, 23, 32, 33 a4 - 44 b1
     * - 12, 13, 21, 31 b2 - 24, 34, 42, 43
     *
     * This class has not been extensively validated
     *
     * @tparam
     */
    template <> class TripleProductDerivativeHolder<4> {
      public:
        Eigen::Matrix<double, 4, 4> value;
        Eigen::VectorXd a1deriv;
        Eigen::MatrixX4d a2deriv;
        Eigen::MatrixX4d a3deriv;
        Eigen::VectorXd a4deriv;
        Eigen::MatrixX4d b1deriv;
        Eigen::MatrixX4d b2deriv;
        int nstr;
        double ssa;

        TripleProductDerivativeHolder() {}

        TripleProductDerivativeHolder(int nstr) : nstr(nstr) { resize(nstr); }

        void resize(int nstr) {
            this->nstr = nstr;

            a1deriv.resize(nstr);
            a2deriv.resize(nstr, 4);
            a3deriv.resize(nstr, 4);
            a4deriv.resize(nstr);
            b1deriv.resize(nstr, 4);
            b2deriv.resize(nstr, 4);

            // Have to set to 0 because sum only goes from l=m upwards, might be
            // faster to do the full sum and not set to 0?
            a1deriv.setZero();
            a2deriv.setZero();
            a3deriv.setZero();
            a4deriv.setZero();
            b1deriv.setZero();
            b2deriv.setZero();
        }

        void calculate(const std::vector<LegendreCoefficient<4>>& coeffs,
                       const std::vector<LegendrePhaseContainer<4>>& lp1s,
                       const std::vector<LegendrePhaseContainer<4>>& lp2s,
                       bool negation, int m) {
            value.setZero();

            for (int l = m; l < nstr; ++l) {
                const auto& coeff = coeffs[l];
                const auto& lp1 = lp1s[l];
                const auto& lp2 = lp2s[l];

                int negation_factor_upper = 1;
                int negation_factor_lower = 1;
                if (negation) {
                    negation_factor_lower = -1;
                    if ((l - m) % 2 != 0) {
                        negation_factor_upper *= -1;
                        negation_factor_lower *= -1;
                    }
                }
                // Calculated product by hand
                value(0, 0) +=
                    lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
                value(0, 1) +=
                    -1.0 * lp1.P() * lp2.R() * coeff.b1 * negation_factor_upper;
                value(0, 2) +=
                    1.0 * lp1.P() * lp2.T() * coeff.b1 * negation_factor_upper;
                // 0, 3 is always 0

                value(1, 0) +=
                    -1.0 * lp1.R() * lp2.P() * coeff.b1 * negation_factor_upper;
                value(1, 1) +=
                    lp1.R() * lp2.R() * coeff.a2 * negation_factor_upper +
                    lp1.T() * lp2.T() * coeff.a3 * negation_factor_lower;
                value(1, 2) +=
                    -lp1.R() * lp2.T() * coeff.a2 * negation_factor_upper -
                    lp1.T() * lp2.R() * coeff.a3 * negation_factor_lower;
                value(1, 3) +=
                    1.0 * lp1.T() * lp2.P() * coeff.b2 * negation_factor_lower;

                value(2, 0) +=
                    lp1.T() * lp2.P() * coeff.b1 * negation_factor_upper;
                value(2, 1) +=
                    -lp1.T() * lp2.R() * coeff.a2 * negation_factor_upper -
                    lp1.R() * lp2.T() * coeff.a3 * negation_factor_lower;
                value(2, 2) +=
                    lp1.T() * lp2.T() * coeff.a2 * negation_factor_upper +
                    lp1.R() * lp2.R() * coeff.a3 * negation_factor_lower;
                value(2, 3) +=
                    -lp1.R() * lp2.P() * coeff.b2 * negation_factor_lower;

                // 3, 0 is always 0
                value(3, 1) +=
                    -lp1.P() * lp2.T() * coeff.b2 * negation_factor_lower;
                value(3, 2) +=
                    lp1.P() * lp2.R() * coeff.b2 * negation_factor_lower;
                value(3, 3) +=
                    lp1.P() * lp2.P() * coeff.a4 * negation_factor_lower;

                // Assign derivatives
                a1deriv(l) = lp1.P() * lp2.P() * negation_factor_upper;

                a2deriv(l, 0) = lp1.R() * lp2.R() * negation_factor_upper;
                a2deriv(l, 1) = -lp1.R() * lp2.T() * negation_factor_upper;
                a2deriv(l, 2) = -lp1.T() * lp2.R() * negation_factor_upper;
                a2deriv(l, 3) = lp1.T() * lp2.T() * negation_factor_upper;

                a3deriv(l, 0) = lp1.T() * lp2.T() * negation_factor_lower;
                a3deriv(l, 1) = -lp1.T() * lp2.R() * negation_factor_lower;
                a3deriv(l, 2) = -lp1.R() * lp2.T() * negation_factor_lower;
                a3deriv(l, 3) = lp1.R() * lp2.R() * negation_factor_lower;

                b1deriv(l, 0) =
                    -1.0 * lp1.P() * lp2.R() * negation_factor_upper;
                b1deriv(l, 1) = lp1.P() * lp2.T() * negation_factor_upper;
                b1deriv(l, 2) =
                    -1.0 * lp1.R() * lp2.P() * negation_factor_upper;
                b1deriv(l, 3) = lp1.T() * lp2.P() * negation_factor_upper;

                b2deriv(l, 0) = lp1.T() * lp2.P() * negation_factor_lower;
                b2deriv(l, 1) = -lp1.R() * lp2.P() * negation_factor_lower;
                b2deriv(l, 2) = -lp1.P() * lp2.T() * negation_factor_lower;
                b2deriv(l, 3) = lp1.P() * lp2.R() * negation_factor_lower;

                a4deriv(l) = lp1.P() * lp2.P() * negation_factor_lower;
            }
        }

        void reduce(const LayerInputDerivative<4>& layer_deriv,
                    Eigen::Matrix<double, 4, 4>& deriv) const;
    };

    /**
     * The values are a 3x3 matrix, the derivative factors for the 4 greek
     * constants are a1 - 11 a2 - 22, 23, 32, 33 a3 - 22, 23, 32, 33 b1 - 12,
     * 13, 21, 31
     *
     * @tparam
     */
    template <> class TripleProductDerivativeHolder<3> {
      public:
        Eigen::Matrix<double, 3, 3> value;
        Eigen::VectorXd a1deriv;
        Eigen::MatrixX4d a2deriv;
        Eigen::MatrixX4d a3deriv;
        Eigen::MatrixX4d b1deriv;
        int nstr;
        double ssa;

        TripleProductDerivativeHolder() {}

        TripleProductDerivativeHolder(int nstr) : nstr(nstr) { resize(nstr); }

        void resize(int nstr) {
            this->nstr = nstr;

            a1deriv.resize(nstr);
            a2deriv.resize(nstr, 4);
            a3deriv.resize(nstr, 4);
            b1deriv.resize(nstr, 4);
        }

        void calculate(const std::vector<LegendreCoefficient<3>>& coeffs,
                       const std::vector<LegendrePhaseContainer<3>>& lp1s,
                       const std::vector<LegendrePhaseContainer<3>>& lp2s,
                       bool negation, int m) {

            value.setZero();
            for (int l = 0; l < m; ++l) {
                a1deriv(l) = 0.0;
                a2deriv(l, Eigen::all).setZero();
                a3deriv(l, Eigen::all).setZero();
                b1deriv(l, Eigen::all).setZero();
            }

            for (int l = m; l < nstr; ++l) {
                const auto& coeff = coeffs[l];
                const auto& lp1 = lp1s[l];
                const auto& lp2 = lp2s[l];

                int negation_factor_upper = 1;
                int negation_factor_lower = 1;
                if (negation) {
                    negation_factor_lower = -1;
                    if ((l - m) % 2 != 0) {
                        negation_factor_upper *= -1;
                        negation_factor_lower *= -1;
                    }
                }

                // Calculated product by hand
                value(0, 0) +=
                    lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
                value(0, 1) +=
                    -1.0 * lp1.P() * lp2.R() * coeff.b1 * negation_factor_upper;
                value(0, 2) +=
                    1.0 * lp1.P() * lp2.T() * coeff.b1 * negation_factor_upper;
                // 0, 3 is always 0

                value(1, 0) +=
                    -1.0 * lp1.R() * lp2.P() * coeff.b1 * negation_factor_upper;
                value(1, 1) +=
                    lp1.R() * lp2.R() * coeff.a2 * negation_factor_upper +
                    lp1.T() * lp2.T() * coeff.a3 * negation_factor_lower;
                value(1, 2) +=
                    -lp1.R() * lp2.T() * coeff.a2 * negation_factor_upper -
                    lp1.T() * lp2.R() * coeff.a3 * negation_factor_lower;

                value(2, 0) +=
                    lp1.T() * lp2.P() * coeff.b1 * negation_factor_upper;
                value(2, 1) +=
                    -lp1.T() * lp2.R() * coeff.a2 * negation_factor_upper -
                    lp1.R() * lp2.T() * coeff.a3 * negation_factor_lower;
                value(2, 2) +=
                    lp1.T() * lp2.T() * coeff.a2 * negation_factor_upper +
                    lp1.R() * lp2.R() * coeff.a3 * negation_factor_lower;

                // Assign derivatives
                a1deriv(l) = lp1.P() * lp2.P() * negation_factor_upper;

                a2deriv(l, 0) = lp1.R() * lp2.R() * negation_factor_upper;
                a2deriv(l, 1) = -lp1.R() * lp2.T() * negation_factor_upper;
                a2deriv(l, 2) = -lp1.T() * lp2.R() * negation_factor_upper;
                a2deriv(l, 3) = lp1.T() * lp2.T() * negation_factor_upper;

                a3deriv(l, 0) = lp1.T() * lp2.T() * negation_factor_lower;
                a3deriv(l, 1) = -lp1.T() * lp2.R() * negation_factor_lower;
                a3deriv(l, 2) = -lp1.R() * lp2.T() * negation_factor_lower;
                a3deriv(l, 3) = lp1.R() * lp2.R() * negation_factor_lower;

                b1deriv(l, 0) =
                    -1.0 * lp1.P() * lp2.R() * negation_factor_upper;
                b1deriv(l, 1) = lp1.P() * lp2.T() * negation_factor_upper;
                b1deriv(l, 2) =
                    -1.0 * lp1.R() * lp2.P() * negation_factor_upper;
                b1deriv(l, 3) = lp1.T() * lp2.P() * negation_factor_upper;
            }
        }

        void reduce(const LayerInputDerivative<3>& layer_deriv,
                    Eigen::Matrix<double, 3, 3>& deriv) const;
    };

    /**
     * @brief Holds the value and derivative of the solar source which is size
     * NSTOKES
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class InhomogeneousSourceHolder {
      public:
        Eigen::Vector<double, NSTOKES> value; /** Value */
        Eigen::VectorXd d_by_a1; /** Derivative with respect to a1 */
        Eigen::VectorXd
            d_by_b1_first; /** First term of derivative with respect to b1 */
        Eigen::VectorXd
            d_by_b1_second; /** Second term of deriative with respect to b1 */
        Eigen::Vector<double, NSTOKES>
            d_by_ssa; /** Derivative with respect to ssa */
        int nstr;

        /**
         * @brief Construct a new Inhomogeneous Source Holder object without
         * allocating
         *
         */
        InhomogeneousSourceHolder() {}

        /**
         * @brief Construct a new Inhomogeneous Source Holder object with
         * allocating
         *
         * @param nstr
         */
        InhomogeneousSourceHolder(int nstr) : nstr(nstr) { resize(nstr); }

        /**
         * @brief Allocate memory
         *
         * @param nstr
         */
        void resize(int nstr) {
            this->nstr = nstr;
            d_by_a1.resize(nstr);
            d_by_b1_first.resize(nstr);
            d_by_b1_second.resize(nstr);
        }

        /**
         * @brief Takes a user input derivative and reduces it to the derivative
         * of the source
         *
         * @param layer_deriv
         * @param deriv
         */
        void reduce(const LayerInputDerivative<NSTOKES>& layer_deriv,
                    Eigen::Vector<double, NSTOKES>& deriv) const;
    };

    /**
     * @brief For NSTOKES=1 we only have one value
     *
     * @tparam
     */
    template <> class InhomogeneousSourceHolder<1> {
      public:
        double value;
        Eigen::VectorXd d_by_legendre_coeff;
        double d_by_ssa;
        int nstr;

        InhomogeneousSourceHolder() {}

        InhomogeneousSourceHolder(int nstr) : nstr(nstr) { resize(nstr); }

        void resize(int nstr) {
            this->nstr = nstr;
            d_by_legendre_coeff.resize(nstr);
        }

        void reduce(const LayerInputDerivative<1>& layer_deriv,
                    double& deriv) const;
    };

    // Forward declaration
    template <typename T> struct Dual;

    /**
     * Class to store a stokes vector (or scalar radiance) and it's derivative
     * We specialize the scalar case to use a double instead of an eigen object
     * (or scalars)
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> struct Radiance {
        Eigen::Vector<double, NSTOKES> value;
        Eigen::Matrix<double, -1, NSTOKES> deriv;

        Radiance() {}

        Radiance(int nderiv, bool zero = true) { resize(nderiv, zero); }

        void resize(int nderiv, bool zero = true) {
            deriv.resize(nderiv, NSTOKES);

            if (zero) {
                value.setZero();
                deriv.setZero();
            }
        }

        void setzero() {
            value.setZero();
            deriv.setZero();
        }

        void apply_transmission_factor(const Dual<double>& transmission);
        void apply_azimuth_expansion(double angle, int m);
        bool converged(double I, double epsilon);

        double I() const { return value(0); }
    };

    template <> struct Radiance<1> {
        double value;
        Eigen::VectorXd deriv;

        Radiance() {}

        Radiance(int nderiv, bool zero = true) { resize(nderiv, zero); }

        void resize(int nderiv, bool zero = true) {
            deriv.resize(nderiv);

            if (zero) {
                value = 0.0;
                deriv.setZero();
            }
        }

        void setzero() {
            value = 0.0;
            deriv.setZero();
        }

        void apply_transmission_factor(const Dual<double>& transmission);
        void apply_azimuth_expansion(double angle, int m);
        bool converged(double I, double epsilon);

        double I() const { return value; }
    };

    template <int NSTOKES, int CNSTR = -1> class InputDerivatives;

} // namespace sasktran_disco
