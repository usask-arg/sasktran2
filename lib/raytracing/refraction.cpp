#include <sasktran2/refraction.h>
#include "sktran_disco/sktran_do.h"

namespace sasktran2::raytracing::refraction {
    std::pair<double, double>
    integrate_path(const sasktran2::Geometry1D& geometry, double rt, double nt,
                   double r1, double r2,
                   std::vector<std::pair<int, double>>& index_weights) {
        const double min_cell_length = 0.1;
        const int num_integration_points = 64;
        const double radius_dither =
            1e-6; // Amount to increase r1, r2 by if they are less than rt due
                  // to numerical issues

        // Always make sure r1 < r2
        if (r2 < r1) {
            std::swap(r1, r2);
        }

        // Make sure r1 and r2 are at least rt + radius_dither
        if (r1 <= (rt + radius_dither)) {
            r1 = rt + radius_dither;
        }

        if (r2 <= (rt + radius_dither)) {
            r2 = rt + radius_dither;
        }

        std::pair<double, double> result; // [path_length, path_angle]
        result.first = 0;
        result.second = 0;

        if (std::abs(r1 - r2) < min_cell_length) {
            result.first = sqrt(r2 * r2 - rt * rt) - sqrt(r1 * r1 - rt * rt);
            // Closed form expression for the deflection angle assuming n=nt=1
            result.second = acos(rt / r2) - acos(rt / r1);

            return result;
        }

        // Construct functions to integrate

        // Eq 21 from Thompson 1982, Ray tracing in a refracting spherically
        // symmetric atmosphere Integral to find the path length of the cell
        auto path_integrand = [&](double x, double n) {
            double sqf = sqrt(1.0 + (n - nt) / n * rt / (x * x)) *
                         sqrt(x * x + (n + nt) / n * rt);
            double F = sqrt(x * x + 2.0 * rt) + sqf;
            double G = sqrt(x * x + 2.0 * rt) * sqf;

            return 2.0 * rt * rt * (nt + n) / n * (nt - n) / n * (x * x + rt) /
                   (x * x * F * G);
        };

        // Also need to find the deflection angle, use a modified version of eq
        // 12
        auto integrand_angle = [&](double x, double n) {
            double r = x * x + rt;

            double t1 = sqrt((n - nt) / n * r / (x * x) + nt / n);
            double t2 = sqrt(r + rt * nt / n);

            return 2.0 * nt * rt / r * (1.0 / (n * t1 * t2));
        };

        // Now we do integration

        // Perform Gaussian quadrature
        double x_low = sqrt(r1 - rt);
        double x_high = sqrt(r2 - rt);

        const double* mu =
            sasktran_disco::getQuadratureAbscissae(num_integration_points);
        const double* wt =
            sasktran_disco::getQuadratureWeights(num_integration_points);

        double diff_d2 = (x_high - x_low) / 2.0;
        double sum_d2 = (x_high + x_low) / 2.0;

        double sum = 0;
        for (int i = 0; i < num_integration_points / 2; ++i) {
            double a1 = 0.5 * mu[i] + 0.5;
            double a2 = -0.5 * mu[i] + 0.5;
            double a3 = 0.5 * mu[i] - 0.5;
            double a4 = -0.5 * mu[i] - 0.5;
            double w = 0.5 * wt[i];

            double x1 = diff_d2 * a1 + sum_d2;
            double x2 = diff_d2 * a2 + sum_d2;
            double x3 = diff_d2 * a3 + sum_d2;
            double x4 = diff_d2 * a4 + sum_d2;

            double n1 = refractive_index_at_altitude(
                geometry, x1 * x1 + rt - geometry.coordinates().earth_radius(),
                index_weights);
            double n2 = refractive_index_at_altitude(
                geometry, x2 * x2 + rt - geometry.coordinates().earth_radius(),
                index_weights);
            double n3 = refractive_index_at_altitude(
                geometry, x3 * x3 + rt - geometry.coordinates().earth_radius(),
                index_weights);
            double n4 = refractive_index_at_altitude(
                geometry, x4 * x4 + rt - geometry.coordinates().earth_radius(),
                index_weights);

            result.first += w * path_integrand(x1, n1);
            result.first += w * path_integrand(x2, n2);
            result.first += w * path_integrand(x3, n3);
            result.first += w * path_integrand(x4, n4);

            result.second += w * integrand_angle(x1, n1);
            result.second += w * integrand_angle(x2, n2);
            result.second += w * integrand_angle(x3, n3);
            result.second += w * integrand_angle(x4, n4);
        }
        result.first *= diff_d2;
        result.second *= diff_d2;

        // The integral only gives the extra curvature path length, so we add on
        // the staright line path length
        result.first += sqrt(r2 * r2 - rt * rt) - sqrt(r1 * r1 - rt * rt);

        return result;
    }

}; // namespace sasktran2::raytracing::refraction
