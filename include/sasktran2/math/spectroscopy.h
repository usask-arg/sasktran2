#include <sasktran2/internal_common.h>

#include <sasktran2/external/faddeeva.h>

#ifdef SKTRAN_OMP_SUPPORT
#include <omp.h>
#endif

namespace sasktran2::math::spectroscopy {
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
        Eigen::Ref<Eigen::MatrixXd> result, // [geometry, wavenumber]
        double line_contribution_width = 25.0, double cull_factor = 0.0,
        const int num_threads = 1) {
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

        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        const double sqrt_2pi = std::sqrt(2.0 * EIGEN_PI);

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

                // Compute line shape
                double norm_factor = 1.0 / (sqrt_2pi * doppler_width);

                for (int w = start_idx; w < end_idx; ++w) {
                    double x = wavenumber_grid[w];
                    // z = (x - center)/(sqrt(2)*dw) +
                    // i*(gamma_val/(sqrt(2)*dw))
                    std::complex<double> z(
                        (x - shifted_center) * inv_sqrt2 / doppler_width,
                        gamma_val * inv_sqrt2 / doppler_width);

                    std::complex<double> val = Faddeeva::w(z);
                    double line_fn = norm_factor * val.real();

                    result(w, g) += line_fn * adjusted_line_intensity;
                }
            }
        }
    }

} // namespace sasktran2::math::spectroscopy
