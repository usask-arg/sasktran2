#include <sasktran2/internal_common.h>

#include <sasktran2/external/faddeeva.h>
#include <sasktran2/math/custom_faddeeva.h>

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
        // return Faddeeva::w(std::complex<double>(x, y)).real();

        // All of this branching just ended up being slower
        const double epsilon = 1e-4;
        if (x * x > 1.52 / epsilon - 2.84 * y * y) {
            return lorentzian(x, y);
        }

        if (abs(x) < 2.15 - 2.53 * y / epsilon) {
            return gaussian(x, y);
        }

        return CustomFaddeeva::faddeeva_w(x, y).real();
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
                         double line_intensity, bool subtract_pedastal) {
        double y = gamma * inv_sqrt2 / doppler_width;
        // Compute line shape
        double norm_factor = 1.0 / (sqrt_2pi * doppler_width);

        double normalized_intensity = line_intensity * norm_factor;

        for (size_t w = startidx; w < endidx; ++w) {
            double x =
                (wavenumbers[w] - line_center) * inv_sqrt2 / doppler_width;

            result(w, geo_index) += voigt(x, y) * normalized_intensity;
        }

        if (subtract_pedastal) {
            double x = 25 * inv_sqrt2 / doppler_width;
            result(Eigen::seq(startidx, endidx - 1), geo_index).array() -=
                voigt(x, y) * normalized_intensity;
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
    inline void voigt_broaden_uniform(
        Eigen::Ref<const Eigen::VectorXd> line_center,    // [line]
        Eigen::Ref<const Eigen::VectorXd> line_intensity, // [line]
        Eigen::Ref<const Eigen::VectorXd> lower_energy,   // [line]
        Eigen::Ref<const Eigen::VectorXd> gamma_air,      // [line]
        Eigen::Ref<const Eigen::VectorXd> gamma_self,     // [line]
        Eigen::Ref<const Eigen::VectorXd> delta_air,      // [line]
        Eigen::Ref<const Eigen::VectorXd> n_air,          // [line]
        Eigen::Ref<const Eigen::VectorXi> iso_id,         // [line]
        Eigen::Ref<const Eigen::MatrixXd> partitions,     // [ngeo, num_isotop]
        Eigen::Ref<const Eigen::VectorXd> molecular_mass, // [num_isotop]
        Eigen::Ref<const Eigen::VectorXd> pressure,       // [geometry]
        Eigen::Ref<const Eigen::VectorXd> pself,          // [geometry]
        Eigen::Ref<const Eigen::VectorXd> temperature,    // [geometry]
        double first_wavenumber, double wavenumber_spacing,
        Eigen::Ref<Eigen::MatrixXd> result, // [wavenumber, geometry]
        double line_contribution_width = 25.0, double cull_factor = 0.0,
        const int num_threads = 1, const double interpolation_delta = 0.0,
        bool subtract_pedastal = false) {
        // Constants (cgs)
        const double c2 = 1.4387769;
        const double SPEED_OF_LIGHT = 2.99792458e10;
        const double NA = 6.02214179e23;
        const double K_B = 1.38064852e-16;

        const int n_line = static_cast<int>(line_center.size());
        const int n_geo = static_cast<int>(pressure.size());
        const int n_wavenumber = static_cast<int>(result.rows());

        double last_wavenumber =
            first_wavenumber + wavenumber_spacing * (n_wavenumber - 1);

        const double epsilon = 1e-4;

        // Initialize result to zero
        result.setZero();

        double max_p_self = pself.maxCoeff();
        if (max_p_self == 0) {
            max_p_self = 1.0;
        }

        // Find the first line is that is > 25 cm^-1 from the first wavenumber
        const double* line_begin = line_center.data();
        const double* line_end = line_center.data() + n_line;
        const double* first_line_it = std::lower_bound(
            line_begin, line_end, first_wavenumber - line_contribution_width);

        const double* last_line_it = std::lower_bound(
            line_begin, line_end, last_wavenumber + line_contribution_width);

        int start_line_idx =
            static_cast<int>(std::distance(line_begin, first_line_it));
        int end_line_idx =
            static_cast<int>(std::distance(line_begin, last_line_it));

        int zero_pivot = 0;

        for (int i = start_line_idx; i < end_line_idx; ++i) {
            const double lc = line_center[i];

            if (line_intensity[i] * 101325 * max_p_self / (K_B * 1e-7 * 296) <
                cull_factor) {
                continue;
            }

            int start_wavenumber_idx =
                floor((lc - line_contribution_width - first_wavenumber) /
                      wavenumber_spacing);
            int end_wavenumber_idx =
                ceil((lc + line_contribution_width - first_wavenumber) /
                     wavenumber_spacing);

            if (start_wavenumber_idx < 0) {
                start_wavenumber_idx = 0;
            }
            if (end_wavenumber_idx > n_wavenumber) {
                end_wavenumber_idx = n_wavenumber;
            }

            while (zero_pivot < end_wavenumber_idx &&
                   zero_pivot * wavenumber_spacing + first_wavenumber < lc) {
                zero_pivot++;
            }

            if (start_wavenumber_idx == end_wavenumber_idx) {
                continue;
            }

            const double le = lower_energy[i];

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
                double doppler_width = lc / SPEED_OF_LIGHT *
                                       std::sqrt(NA * K_B * T / mol_mass) /
                                       inv_sqrt2;

                // Compute gamma (Lorentz broadening)
                double g_air = gamma_air[i];
                double g_self = gamma_self[i];
                double da = delta_air[i];
                double na = n_air[i];

                double gamma_val = std::pow(296.0 / T, na) *
                                   (g_air * (P - Pself) + g_self * Pself);

                double shifted_center = lc + da * P;

                double y = gamma_val / doppler_width;
                // Compute line shape
                double norm_factor =
                    1.0 / (sqrt_2pi * doppler_width * inv_sqrt2);

                double normalized_intensity =
                    adjusted_line_intensity * norm_factor;

                // condition for lorentzian x * x > 1.52 / epsilon - 2.84 * y *
                // y => x = sqrt(1.52 / epsilon - 2.84 * y * y)
                int lorentzian_low = 0;
                int lorentzian_high = 0;
                // condition for gaussian is abs(x) < 2.15 - 2.53 * y / epsilon

                if (2.84 * y * y > 1.52 / epsilon) {
                    // Entirely lorentzian, skip the split
                    for (size_t w = start_wavenumber_idx;
                         w < end_wavenumber_idx; ++w) {
                        double x = (w * wavenumber_spacing + first_wavenumber -
                                    shifted_center) /
                                   doppler_width;

                        result(w, g) += lorentzian(x, y) * normalized_intensity;
                    }
                } else {
                    double max_x = end_wavenumber_idx * wavenumber_spacing +
                                   first_wavenumber / doppler_width;
                    double min_x = start_wavenumber_idx * wavenumber_spacing +
                                   first_wavenumber / doppler_width;

                    double max_abs_x =
                        std::max(std::abs(max_x), std::abs(min_x));

                    if (max_abs_x < 2.15 - 2.53 * y / epsilon) {
                        // Entirely gaussian, skip the split
                        for (size_t w = start_wavenumber_idx;
                             w < end_wavenumber_idx; ++w) {
                            double x = (w * wavenumber_spacing +
                                        first_wavenumber - shifted_center) /
                                       doppler_width;

                            result(w, g) +=
                                gaussian(x, y) * normalized_intensity;
                        }
                    } else {
                        // Find the lorentzian split point
                        double split_x =
                            std::sqrt(1.52 / epsilon - 2.84 * y * y);

                        // Distance from the line center to the split point in
                        // index units
                        int lorentzian_split = floor((split_x * doppler_width) /
                                                     wavenumber_spacing);

                        if (zero_pivot + lorentzian_split >
                            end_wavenumber_idx) {
                            lorentzian_split = end_wavenumber_idx - zero_pivot;
                        }
                        if (zero_pivot - lorentzian_split <
                            start_wavenumber_idx) {
                            lorentzian_split =
                                zero_pivot - start_wavenumber_idx;
                        }

                        for (size_t w = start_wavenumber_idx;
                             w < zero_pivot - lorentzian_split; ++w) {
                            double x = (w * wavenumber_spacing +
                                        first_wavenumber - shifted_center) /
                                       doppler_width;

                            result(w, g) +=
                                lorentzian(x, y) * normalized_intensity;
                        }

                        for (size_t w = zero_pivot - lorentzian_split;
                             w < zero_pivot + lorentzian_split; ++w) {
                            double x = (w * wavenumber_spacing +
                                        first_wavenumber - shifted_center) /
                                       doppler_width;

                            result(w, g) +=
                                CustomFaddeeva::faddeeva_w(x, y).real() *
                                normalized_intensity;
                        }

                        for (size_t w = zero_pivot + lorentzian_split;
                             w < end_wavenumber_idx; ++w) {
                            double x = (w * wavenumber_spacing +
                                        first_wavenumber - shifted_center) /
                                       doppler_width;

                            result(w, g) +=
                                lorentzian(x, y) * normalized_intensity;
                        }
                    }
                    if (subtract_pedastal) {
                        double x = line_contribution_width / doppler_width;
                        result(Eigen::seq(start_wavenumber_idx,
                                          end_wavenumber_idx - 1),
                               g)
                            .array() -= voigt(x, y) * normalized_intensity;
                    }
                }
            }
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
        const int num_threads = 1, const double interpolation_delta = 0.0,
        bool subtract_pedastal = false) {
        // Constants (cgs)
        const double c2 = 1.4387769;
        const double SPEED_OF_LIGHT = 2.99792458e10;
        const double NA = 6.02214179e23;
        const double K_B = 1.38064852e-16;

        const int n_line = static_cast<int>(line_center.size());
        const int n_geo = static_cast<int>(pressure.size());
        const int n_wavenumber = static_cast<int>(wavenumber_grid.size());

        const double epsilon = 1e-4;

        // Initialize result to zero
        result.setZero();

        double max_p_self = pself.maxCoeff();
        if (max_p_self == 0) {
            max_p_self = 1.0;
        }

        double first_wavenumber = wavenumber_grid[0];
        double last_wavenumber = wavenumber_grid[n_wavenumber - 1];
        // Find the first line is that is > 25 cm^-1 from the first wavenumber
        const double* line_begin = line_center.data();
        const double* line_end = line_center.data() + n_line;
        const double* first_line_it = std::lower_bound(
            line_begin, line_end, first_wavenumber - line_contribution_width);

        const double* last_line_it = std::lower_bound(
            line_begin, line_end, last_wavenumber + line_contribution_width);

        int start_line_idx =
            static_cast<int>(std::distance(line_begin, first_line_it));
        int end_line_idx =
            static_cast<int>(std::distance(line_begin, last_line_it));

        // Now we have our lines that can potentially contribute, we want to
        // keep track of the lower and upper indicies on the wavenumber grid
        // that the line can contribute to
        int start_wavenumber_idx = 0;
        int end_wavenumber_idx = 0;
        int zero_pivot = 0;

        // For certian y values we can switch to lorentzian profiles,
        // if 2.84 * y *y > 1.52 / epsilon we can use lorentzian for all x

        // condition for lorentzian x * x > 1.52 / epsilon - 2.84 * y * y
        int lorentzian_low = 0;
        int lorentzian_high = 0;
        // condition for gaussian is abs(x) < 2.15 - 2.53 * y / epsilon

        for (int i = start_line_idx; i < end_line_idx; ++i) {
            const double lc = line_center[i];

            if (line_intensity[i] * 101325 * max_p_self / (K_B * 1e-7 * 296) <
                cull_factor) {
                continue;
            }

            while (wavenumber_grid[start_wavenumber_idx] <
                       lc - line_contribution_width &&
                   start_wavenumber_idx < n_wavenumber) {
                start_wavenumber_idx++;
            }
            while (wavenumber_grid[end_wavenumber_idx] <
                       lc + line_contribution_width &&
                   end_wavenumber_idx < n_wavenumber) {
                end_wavenumber_idx++;
            }

            while (wavenumber_grid[zero_pivot] < 0 &&
                   zero_pivot < end_wavenumber_idx) {
                zero_pivot++;
            }

            if (start_wavenumber_idx == end_wavenumber_idx) {
                continue;
            }

            const double le = lower_energy[i];

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
                double doppler_width = lc / SPEED_OF_LIGHT *
                                       std::sqrt(NA * K_B * T / mol_mass) /
                                       inv_sqrt2;

                // Compute gamma (Lorentz broadening)
                double g_air = gamma_air[i];
                double g_self = gamma_self[i];
                double da = delta_air[i];
                double na = n_air[i];

                double gamma_val = std::pow(296.0 / T, na) *
                                   (g_air * (P - Pself) + g_self * Pself);

                double shifted_center = lc + da * P;

                double y = gamma_val / doppler_width;
                // Compute line shape
                double norm_factor =
                    1.0 / (sqrt_2pi * doppler_width * inv_sqrt2);

                double normalized_intensity =
                    adjusted_line_intensity * norm_factor;

                for (size_t w = start_wavenumber_idx; w < end_wavenumber_idx;
                     ++w) {
                    double x =
                        (wavenumber_grid[w] - shifted_center) / doppler_width;

                    result(w, g) += voigt(x, y) * normalized_intensity;
                }
                if (subtract_pedastal) {
                    double x = line_contribution_width / doppler_width;
                    result(Eigen::seq(start_wavenumber_idx,
                                      end_wavenumber_idx - 1),
                           g)
                        .array() -= voigt(x, y) * normalized_intensity;
                }
            }
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
    inline void voigt_broaden_with_line_coupling(
        Eigen::Ref<const Eigen::VectorXd> line_center,     // [line]
        Eigen::Ref<const Eigen::VectorXd> line_intensity,  // [line]
        Eigen::Ref<const Eigen::VectorXd> lower_energy,    // [line]
        Eigen::Ref<const Eigen::VectorXd> gamma_air,       // [line]
        Eigen::Ref<const Eigen::VectorXd> gamma_self,      // [line]
        Eigen::Ref<const Eigen::VectorXd> delta_air,       // [line]
        Eigen::Ref<const Eigen::VectorXd> n_air,           // [line]
        Eigen::Ref<const Eigen::VectorXi> iso_id,          // [line]
        Eigen::Ref<const Eigen::MatrixXd> partitions,      // [ngeo, num_isotop]
        Eigen::Ref<const Eigen::MatrixXd> Y_coupling,      // [ngeo, line]
        Eigen::Ref<const Eigen::MatrixXd> G_coupling,      // [ngeo, line]
        Eigen::Ref<const Eigen::VectorXd> molecular_mass,  // [num_isotop]
        Eigen::Ref<const Eigen::VectorXd> pressure,        // [geometry]
        Eigen::Ref<const Eigen::VectorXd> pself,           // [geometry]
        Eigen::Ref<const Eigen::VectorXd> temperature,     // [geometry]
        Eigen::Ref<const Eigen::VectorXd> wavenumber_grid, // [wavenumber]
        Eigen::Ref<Eigen::MatrixXd> result, // [wavenumber, geometry]
        double line_contribution_width = 25.0, double cull_factor = 0.0,
        const int num_threads = 1, const double interpolation_delta = 0.0,
        bool subtract_pedastal = false) {
        // Constants (cgs)
        const double c2 = 1.4387769;
        const double SPEED_OF_LIGHT = 2.99792458e10;
        const double NA = 6.02214179e23;
        const double K_B = 1.38064852e-16;

        const int n_line = static_cast<int>(line_center.size());
        const int n_geo = static_cast<int>(pressure.size());
        const int n_wavenumber = static_cast<int>(wavenumber_grid.size());

        const double epsilon = 1e-4;

        // Initialize result to zero
        result.setZero();

        double max_p_self = pself.maxCoeff();
        if (max_p_self == 0) {
            max_p_self = 1.0;
        }

        double first_wavenumber = wavenumber_grid[0];
        double last_wavenumber = wavenumber_grid[n_wavenumber - 1];
        // Find the first line is that is > 25 cm^-1 from the first wavenumber
        const double* line_begin = line_center.data();
        const double* line_end = line_center.data() + n_line;
        const double* first_line_it = std::lower_bound(
            line_begin, line_end, first_wavenumber - line_contribution_width);

        const double* last_line_it = std::lower_bound(
            line_begin, line_end, last_wavenumber + line_contribution_width);

        int start_line_idx =
            static_cast<int>(std::distance(line_begin, first_line_it));
        int end_line_idx =
            static_cast<int>(std::distance(line_begin, last_line_it));

        // Now we have our lines that can potentially contribute, we want to
        // keep track of the lower and upper indicies on the wavenumber grid
        // that the line can contribute to
        int start_wavenumber_idx = 0;
        int end_wavenumber_idx = 0;
        int zero_pivot = 0;

        // For certian y values we can switch to lorentzian profiles,
        // if 2.84 * y *y > 1.52 / epsilon we can use lorentzian for all x

        // condition for lorentzian x * x > 1.52 / epsilon - 2.84 * y * y
        int lorentzian_low = 0;
        int lorentzian_high = 0;
        // condition for gaussian is abs(x) < 2.15 - 2.53 * y / epsilon

        for (int i = start_line_idx; i < end_line_idx; ++i) {
            const double lc = line_center[i];

            if (line_intensity[i] * 101325 * max_p_self / (K_B * 1e-7 * 296) <
                cull_factor) {
                continue;
            }

            while (wavenumber_grid[start_wavenumber_idx] <
                       lc - line_contribution_width &&
                   start_wavenumber_idx < n_wavenumber) {
                start_wavenumber_idx++;
            }
            while (wavenumber_grid[end_wavenumber_idx] <
                       lc + line_contribution_width &&
                   end_wavenumber_idx < n_wavenumber) {
                end_wavenumber_idx++;
            }

            while (wavenumber_grid[zero_pivot] < 0 &&
                   zero_pivot < end_wavenumber_idx) {
                zero_pivot++;
            }

            if (start_wavenumber_idx == end_wavenumber_idx) {
                continue;
            }

            const double le = lower_energy[i];

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
                double doppler_width = lc / SPEED_OF_LIGHT *
                                       std::sqrt(NA * K_B * T / mol_mass) /
                                       inv_sqrt2;

                // Compute gamma (Lorentz broadening)
                double g_air = gamma_air[i];
                double g_self = gamma_self[i];
                double da = delta_air[i];
                double na = n_air[i];

                double gamma_val = std::pow(296.0 / T, na) *
                                   (g_air * (P - Pself) + g_self * Pself);

                double shifted_center = lc + da * P;

                double y = gamma_val / doppler_width;
                // Compute line shape
                double norm_factor =
                    1.0 / (sqrt_2pi * doppler_width * inv_sqrt2);

                double normalized_intensity =
                    adjusted_line_intensity * norm_factor;

                if (Y_coupling(g, i) != 0.0 || G_coupling(g, i) != 0.0) {
                    for (size_t w = start_wavenumber_idx;
                         w < end_wavenumber_idx; ++w) {
                        double x = (wavenumber_grid[w] - shifted_center) /
                                   doppler_width;

                        result(w, g) += (CustomFaddeeva::faddeeva_w(x, y) *
                                         std::complex<double>(
                                             1 + pressure(g) * pressure(g) *
                                                     G_coupling(g, i),
                                             -pressure(g) * Y_coupling(g, i)))
                                            .real() *
                                        normalized_intensity;
                    }
                } else {
                    for (size_t w = start_wavenumber_idx;
                         w < end_wavenumber_idx; ++w) {
                        double x = (wavenumber_grid[w] - shifted_center) /
                                   doppler_width;

                        result(w, g) += voigt(x, y) * normalized_intensity;
                    }
                }

                if (subtract_pedastal) {
                    double x = line_contribution_width / doppler_width;
                    result(Eigen::seq(start_wavenumber_idx,
                                      end_wavenumber_idx - 1),
                           g)
                        .array() -= voigt(x, y) * normalized_intensity;
                }
            }
        }
    }

} // namespace sasktran2::math::spectroscopy
