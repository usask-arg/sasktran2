#include <sasktran2/internal_common.h>

#include <sasktran2/external/faddeeva.h>

#ifdef SKTRAN_OMP_SUPPORT
#include <omp.h>
#endif

namespace sasktran2::math::spectroscopy {
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    const double sqrt_2pi = std::sqrt(2.0 * EIGEN_PI);
    const double sqrt_pi = std::sqrt(EIGEN_PI);

    // Quadratic interpolation helper

    /**
     * Interpolates to w, given f(a), f(b), f(c) using quadratic interpolation
     *
     * @param w value to interpolate two
     * @param a
     * @param b
     * @param c
     * @param fa
     * @param fb
     * @param fc
     * @return double
     */
    inline double quadratic_interpolate(double w, double a, double b, double c,
                                        double fa, double fb, double fc) {
        double denom_ab = (a - b);
        double denom_ac = (a - c);
        double denom_bc = (b - c);

        double term_a = fa * (w - b) * (w - c) / ((a - b) * (a - c));
        double term_b = fb * (w - a) * (w - c) / ((b - a) * (b - c));
        double term_c = fc * (w - a) * (w - b) / ((c - a) * (c - b));
        return term_a + term_b + term_c;
    }

    /**
     * @brief Humlicek region 1 approximation
     *
     * Valid when |x| + y > 15
     *
     * @param x
     * @param y
     * @return double
     */
    inline double humlicek_region1(double x, double y) {
        double a1 = 0.2820948 * y + 0.5641896 * y * y * y;
        double b1 = 0.5641896 * y;
        double a2 = 0.25 + y * y + y * y * y * y;
        double b2 = -1.0 + 2 * y * y;

        return (a1 + b1 * x * x) / (a2 + b2 * x * x + x * x * x * x);
    }

    /**
     * @brief Humlicek region 2 approximation
     *
     * Valid when 5.5 < |x| + y < 15
     *
     * @param x
     * @param y
     * @return double
     */
    inline double humlicek_region2(double x, double y) {
        double y2 = y * y;
        double y3 = y2 * y;
        double y4 = y3 * y;
        double y5 = y4 * y;
        double y6 = y5 * y;
        double y7 = y6 * y;
        double y8 = y7 * y;

        double a3 = 1.05786 * y + 4.65456 * y3 + 3.10304 * y5 + 0.56419 * y7;
        double b3 = 2.962 * y + 0.56419 * y3 + 1.69257 * y5;
        double c3 = 1.69257 * y - 2.53885 * y3;
        double d3 = 0.56419 * y;

        double a4 = 0.5625 + 4.5 * y2 + 10.5 * y4 + 6 * y6 + y8;
        double b4 = -4.5 + 9.0 * y2 + 6.0 * y4 + 4.0 * y6;
        double c4 = 10.5 - 6.0 * y2 + 6.0 * y4;
        double d4 = -6.0 + 4.0 * y2;

        double x2 = x * x;
        double x4 = x2 * x2;
        double x6 = x4 * x2;
        double x8 = x6 * x2;

        return (a3 + b3 * x2 + c3 * x4 + d3 * x6) /
               (a4 + b4 * x2 + c4 * x4 + d4 * x6 + x8);
    }

    /**
     * @brief Analytic approximation to the Voigt profile using the
     * Tepper-Garcia function
     *
     * Reference: https://arxiv.org/pdf/astro-ph/0602124
     *
     * Valid for gamma / doppler_width < 1e-4, has a special case at x=0
     *
     * @param x
     * @param gamma
     * @param doppler_width
     * @return double
     */
    inline double tepper_garcia(double x, double gamma, double doppler_width) {
        if (x == 0) {
            std::complex<double> z(0, gamma * inv_sqrt2 / doppler_width);

            std::complex<double> val = Faddeeva::w(z);
            return val.real();
        }

        double u = x * inv_sqrt2 / (doppler_width);
        double a = gamma * inv_sqrt2 / (doppler_width);

        double P = u * u;
        double Q = 1.5 / P;
        double R = std::exp(-P);

        double T = R - (a / (sqrt_pi * P)) *
                           (R * R * (4.0 * P * P + 7 * P + 4 + Q) - Q - 1);

        return T;
    }

    /**
     * @brief Lorentzian line shape
     *
     * @param x
     * @param gamma
     * @return double
     */
    inline double lorentzian(double x, double y) {
        return y / (sqrt_pi * (x * x + y * y));
    }

    /**
     * @brief Gaussian line shape
     *
     * @param x
     * @param sigma
     * @return double
     */
    inline double gaussian(double x, double y) { return std::exp(-x * x); }

    /**
     * @brief Voigt line shape
     *
     * Uses the Faddeeva function to calculate the voigt function.  Using other
     * approximations have been tried, but the branching made them slower than
     * just using the Faddeeva function.
     *
     * @param x
     * @param y
     * @return double
     */
    inline double voigt(double x, double y) {
        return Faddeeva::w(std::complex<double>(x, y)).real();

        // All of this branching just ended up being slower
        const double epsilon = 1e-4;
        if (x * x > 1.52 / epsilon - 2.84 * y * y) {
            return lorentzian(x, y);
        }

        if (abs(x) < 2.15 - 2.53 * y / epsilon) {
            return gaussian(x, y);
        }

        if (std::abs(x) + y > 15) {
            return humlicek_region1(x, y);
        } else if (std::abs(x) + y > 5.5 && false) {
            return humlicek_region2(x, y);
        } else {
            std::complex<double> val = Faddeeva::w(std::complex<double>(x, y));
            return val.real();
        }
    }

    /**
     * @brief Broadens a line using the Voigt function, and adds the result to
     * the result matrix
     *
     * @param result
     * @param wavenumbers
     * @param geo_index
     * @param startidx
     * @param endidx
     * @param doppler_width
     * @param gamma
     * @param line_center
     * @param line_intensity
     */
    inline void sum_line(Eigen::Ref<Eigen::MatrixXd>& result,
                         const Eigen::VectorXd& wavenumbers, int geo_index,
                         size_t startidx, size_t endidx, double doppler_width,
                         double gamma, double line_center,
                         double line_intensity) {
        double y = gamma * inv_sqrt2 / doppler_width;
        // Compute line shape
        double norm_factor = 1.0 / (sqrt_2pi * doppler_width);

        double normalized_intensity = line_intensity * norm_factor;

        for (size_t w = startidx; w < endidx; ++w) {
            double x =
                (wavenumbers[w] - line_center) * inv_sqrt2 / doppler_width;

            result(w, geo_index) += voigt(x, y) * normalized_intensity;
        }
    }

    /**
     * @brief Broadens a line using the Voigt function, and adds the result to
     * the result matrix by interpolating between points
     *
     * This function is mostly just experimental, it turned out to be slower
     * than just using the voigt function
     *
     * @param result
     * @param wavenumbers
     * @param geo_index
     * @param startidx
     * @param endidx
     * @param doppler_width
     * @param gamma
     * @param line_center
     * @param line_intensity
     * @param dx
     */
    inline void sum_line_interpolated(Eigen::Ref<Eigen::MatrixXd>& result,
                                      const Eigen::VectorXd& wavenumbers,
                                      int geo_index, size_t startidx,
                                      size_t endidx, double doppler_width,
                                      double gamma, double line_center,
                                      double line_intensity, double dx) {
        double y = gamma * inv_sqrt2 / doppler_width;
        // Compute line shape
        double norm_factor = 1.0 / (sqrt_2pi * doppler_width);

        double normalized_intensity = line_intensity * norm_factor;

        double a =
            (wavenumbers[startidx] - line_center) * inv_sqrt2 / doppler_width;
        ;
        double b = a + dx;
        double c = b + dx;

        double fa = voigt(a, y);
        double fb = voigt(b, y);
        double fc = voigt(c, y);

        for (size_t w = startidx; w < endidx; ++w) {
            double x =
                (wavenumbers[w] - line_center) * inv_sqrt2 / doppler_width;

            if (x > c) {
                a = b;
                b = c;
                c = b + dx;

                fa = fb;
                fb = fc;
                fc = voigt(c, y);
            }
            result(w, geo_index) +=
                quadratic_interpolate(x, a, b, c, fa, fb, fc) *
                normalized_intensity;
        }
    }

    /**
     * Broadens the input line data using the Voigt profile. All values are
     * assumed to be in cgs units unless otherwise specified. The result matrix
     * is filled with the broadened line data.
     *
     * @param line_center Central vacuum wavenumber of the line [line]
     * @param line_intensity Line intensity [line]
     * @param lower_energy Lower enregy level of the line [line]
     * @param gamma_air Lorentz broadening due to air [line]
     * @param gamma_self Lorentz broadening due to self [line]
     * @param delta_air Pressure shift due to air [line]
     * @param n_air [line]
     * @param iso_id Identifier for the isotopalog [line]
     * @param partitions Partition function ratios at the specified temperatures
     * for each isotopalog [ngeo, num_isotop]
     * @param molecular_mass Molecular mass of each isotopalog [num_isotop]
     * @param pressure Partial pressure at each geometry (Pressure in Pa /
     * 101325) [geometry]
     * @param pself Self partial pressure at each geometry [geometry]
     * @param temperature Temperature in K at each geometry [geometry]
     * @param wavenumber_grid Wavenumber grid to produce the result on
     * [wavenumber]
     * @param result Resulting cross sections [geometry, wavenumber]
     * @param line_contribution_width +/- from the line center to consider.
     * Defaults to 25 [cm^-1]
     * @param cull_factor Lines with total estimated vertical column relative
     * optical depth less than this value are culled. Defaults to 0.0
     * @param num_threads Number of threads to use if OMP is enabled. Defaults
     * to 1.
     */
    inline void voigt_broaden(
        Eigen::Ref<const Eigen::VectorXd> line_center,     // [line]
        Eigen::Ref<const Eigen::VectorXd> line_intensity,  // [line]
        Eigen::Ref<const Eigen::VectorXd> lower_energy,    // [line]
        Eigen::Ref<const Eigen::VectorXd> gamma_air,       // [line]
        Eigen::Ref<const Eigen::VectorXd> gamma_self,      // [line]
        Eigen::Ref<const Eigen::VectorXd> delta_air,       // [line]
        Eigen::Ref<const Eigen::VectorXd> n_air,           // [line]
        Eigen::Ref<const Eigen::VectorXi> iso_id,          // [line]
        Eigen::Ref<const Eigen::MatrixXd> partitions,      // [ngeo, num_isotop]
        Eigen::Ref<const Eigen::VectorXd> molecular_mass,  // [num_isotop]
        Eigen::Ref<const Eigen::VectorXd> pressure,        // [geometry]
        Eigen::Ref<const Eigen::VectorXd> pself,           // [geometry]
        Eigen::Ref<const Eigen::VectorXd> temperature,     // [geometry]
        Eigen::Ref<const Eigen::VectorXd> wavenumber_grid, // [wavenumber]
        Eigen::Ref<Eigen::MatrixXd> result, // [wavenumber, geometry]
        double line_contribution_width = 25.0, double cull_factor = 0.0,
        const int num_threads = 1, const double interpolation_delta = 0.0) {
        // Constants (cgs)
        const double c2 = 1.4387769;
        const double SPEED_OF_LIGHT = 2.99792458e10;
        const double NA = 6.02214179e23;
        const double K_B = 1.38064852e-16;

        const int n_line = static_cast<int>(line_center.size());
        const int n_geo = static_cast<int>(pressure.size());
        const int n_wavenumber = static_cast<int>(wavenumber_grid.size());

        // Initialize result to zero
        result.setZero();

        double max_p_self = pself.maxCoeff();
        if (max_p_self == 0) {
            max_p_self = 1.0;
        }

        for (int i = 0; i < n_line; ++i) {
            const double lc = line_center[i];

            if (lc > wavenumber_grid[n_wavenumber - 1] +
                         line_contribution_width ||
                lc < (wavenumber_grid[0] - line_contribution_width)) {
                continue;
            }

            if (line_intensity[i] * 101325 * max_p_self / (K_B * 1e-7 * 296) <
                cull_factor) {
                continue;
            }

            const double le = lower_energy[i];

            // Determine wavenumber bounds
            double start_val = lc - line_contribution_width;
            double end_val = lc + line_contribution_width;

            // Search for start and end indices
            const double* wbegin = wavenumber_grid.data();
            const double* wend = wavenumber_grid.data() + n_wavenumber;
            const double* start_it = std::lower_bound(wbegin, wend, start_val);
            const double* end_it = std::lower_bound(wbegin, wend, end_val);

            int start_idx = static_cast<int>(std::distance(wbegin, start_it));
            int end_idx = static_cast<int>(std::distance(wbegin, end_it));

            if (start_idx == end_idx) {
                continue;
            }

            // Common reference temperature normalization
            double denominator =
                (1.0 - std::exp(-c2 * lc / 296.0)) * std::exp(-c2 * le / 296.0);

#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (int g = 0; g < n_geo; ++g) {
                const double T = temperature[g];
                const double P = pressure[g];
                const double Pself = pself[g];

                double numerator =
                    std::exp(-c2 * le / T) * (1.0 - std::exp(-c2 * lc / T));

                double common_factor = numerator / denominator;
                int iso_index = iso_id[i] - 1;
                double partition_factor = 1.0 / partitions(g, iso_index);
                double adjusted_line_intensity =
                    line_intensity[i] * common_factor * partition_factor;

                double mol_mass = molecular_mass[iso_index];
                double doppler_width =
                    lc / SPEED_OF_LIGHT * std::sqrt(NA * K_B * T / mol_mass);

                // Compute gamma (Lorentz broadening)
                double g_air = gamma_air[i];
                double g_self = gamma_self[i];
                double da = delta_air[i];
                double na = n_air[i];

                double gamma_val = std::pow(296.0 / T, na) *
                                   (g_air * (P - Pself) + g_self * Pself);

                double shifted_center = lc + da * P;

                if (interpolation_delta == 0.0) {
                    sum_line(result, wavenumber_grid, g, start_idx, end_idx,
                             doppler_width, gamma_val, shifted_center,
                             adjusted_line_intensity);
                } else {
                    sum_line_interpolated(
                        result, wavenumber_grid, g, start_idx, end_idx,
                        doppler_width, gamma_val, shifted_center,
                        adjusted_line_intensity, interpolation_delta);
                }
            }
        }
    }

} // namespace sasktran2::math::spectroscopy
