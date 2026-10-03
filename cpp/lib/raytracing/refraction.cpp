#include <algorithm>
#include <array>
#include <cmath>
#include <sasktran2/refraction.h>
#include "sktran_disco/sktran_do.h"

namespace {
    // Eq 21 from Thompson 1982, Ray tracing in a refracting spherically
    // symmetric atmosphere. Integrand for the extra path length of a cell
    // relative to the straight line, in the variable x = sqrt(r - rt).
    inline double path_integrand(double x, double n, double rt, double nt) {
        double sqf = sqrt(1.0 + (n - nt) / n * rt / (x * x)) *
                     sqrt(x * x + (n + nt) / n * rt);
        double F = sqrt(x * x + 2.0 * rt) + sqf;
        double G = sqrt(x * x + 2.0 * rt) * sqf;

        return 2.0 * rt * rt * (nt + n) / n * (nt - n) / n * (x * x + rt) /
               (x * x * F * G);
    }

    // Modified version of eq 12 from Thompson 1982, integrand for the central
    // angle swept by the ray in the variable x = sqrt(r - rt).
    inline double angle_integrand(double x, double n, double rt, double nt) {
        double r = x * x + rt;

        double t1 = sqrt((n - nt) / n * r / (x * x) + nt / n);
        double t2 = sqrt(r + rt * nt / n);

        return 2.0 * nt * rt / r * (1.0 / (n * t1 * t2));
    }

    // Straight line distance between radii rt <= r1 <= r2 for a ray with
    // tangent radius rt. Factored so the argument is exactly zero at the
    // tangent point, r * r - rt * rt can round negative when contracted to an
    // FMA.
    inline double straight_path(double rt, double r1, double r2) {
        return std::sqrt((r2 - rt) * (r2 + rt)) -
               std::sqrt((r1 - rt) * (r1 + rt));
    }

    // Eight point Gauss-Legendre rule on [-1, 1], positive nodes only
    constexpr std::array<double, 4> gl8_nodes = {
        0.1834346424956498049394761, 0.5255324099163289858177390,
        0.7966664774136267395915539, 0.9602898564975362316835609};
    constexpr std::array<double, 4> gl8_weights = {
        0.3626837833783619829651504, 0.3137066458778872873379622,
        0.2223810344533744705443560, 0.1012285362903762591525314};

    /**
     * Central angle swept by a refracted ray between radii r1 and r2 within
     * a single altitude grid cell. Uses the same integrand as integrate_path
     * with a lower order quadrature, which is sufficient for the smooth
     * integrand inside one cell and is used where many evaluations are
     * required.
     */
    double integrate_angle(const sasktran2::Geometry1D& geometry, double rt,
                           double nt, double r1, double r2,
                           std::vector<std::pair<int, double>>& index_weights) {
        const double min_cell_length = 0.1;

        if (r2 < r1) {
            std::swap(r1, r2);
        }
        // The integrand is regular at the tangent point in x = sqrt(r - rt),
        // so integrate all the way down to it. Skipping even a micrometre
        // next to the tangent point would miss sqrt(2e-6 / rt) of angle.
        r1 = std::max(r1, rt);
        r2 = std::max(r2, rt);

        if (r2 - r1 < min_cell_length) {
            return std::acos(std::min(1.0, rt / r2)) -
                   std::acos(std::min(1.0, rt / r1));
        }

        const double x_low = sqrt(r1 - rt);
        const double x_high = sqrt(r2 - rt);
        const double half_width = (x_high - x_low) / 2.0;
        const double center = (x_high + x_low) / 2.0;
        const double earth_radius = geometry.coordinates().earth_radius();

        double result = 0.0;
        for (std::size_t i = 0; i < gl8_nodes.size(); ++i) {
            for (const double x : {center + half_width * gl8_nodes[i],
                                   center - half_width * gl8_nodes[i]}) {
                const double n = sasktran2::raytracing::refraction::
                    refractive_index_at_altitude(
                        geometry, x * x + rt - earth_radius, index_weights);
                result += gl8_weights[i] * angle_integrand(x, n, rt, nt);
            }
        }
        return result * half_width;
    }

    /**
     * Central angle swept by a refracted ray between radii r_low < r_high,
     * split at the altitude grid boundaries in the same way as the ray
     * tracer splits a ray into layers.
     */
    double integrate_angle_over_grid(
        const sasktran2::Geometry1D& geometry, double rt, double nt,
        double r_low, double r_high,
        std::vector<std::pair<int, double>>& index_weights) {
        const auto& grid = geometry.altitude_grid().grid();
        const double earth_radius = geometry.coordinates().earth_radius();

        double result = 0.0;
        double lower = r_low;
        auto boundary =
            std::upper_bound(grid.begin(), grid.end(), r_low - earth_radius);
        while (lower < r_high) {
            const double upper =
                boundary == grid.end()
                    ? r_high
                    : std::min(r_high, earth_radius + *boundary);
            result +=
                integrate_angle(geometry, rt, nt, lower, upper, index_weights);
            lower = upper;
            if (boundary == grid.end()) {
                break;
            }
            ++boundary;
        }
        return result;
    }

    /**
     * Angle between the start point and the asymptotic direction above the
     * atmosphere of a refracted ray that starts at the given radius with the
     * given local zenith angle, measured in the plane of the ray.
     *
     * @return false if the ray intersects the surface
     */
    bool asymptotic_direction_angle(
        const sasktran2::Geometry1D& geometry, double radius,
        double start_refractive_index, double zenith, double& angle,
        std::vector<std::pair<int, double>>& index_weights) {
        namespace refraction = sasktran2::raytracing::refraction;
        const auto& grid = geometry.altitude_grid().grid();
        const double earth_radius = geometry.coordinates().earth_radius();
        const double top_radius =
            earth_radius + grid(Eigen::placeholders::last);

        // Conserved ray invariant n r sin(zenith)
        const double invariant =
            start_refractive_index * radius * std::sin(zenith);
        const double rt =
            refraction::tangent_radius(geometry, invariant, index_weights);
        const double nt = refraction::refractive_index_at_altitude(
            geometry, rt - earth_radius, index_weights);

        double swept;
        if (std::cos(zenith) >= 0.0) {
            swept = integrate_angle_over_grid(geometry, rt, nt, radius,
                                              top_radius, index_weights);
        } else {
            if (rt - earth_radius <= grid(0)) {
                return false;
            }
            // Down to the tangent point, back up through the start radius,
            // then on to the top of the atmosphere
            swept = 2.0 * integrate_angle_over_grid(geometry, rt, nt, rt,
                                                    radius, index_weights) +
                    integrate_angle_over_grid(geometry, rt, nt, radius,
                                              top_radius, index_weights);
        }

        // Above the atmosphere the ray is straight with the same invariant
        angle = swept + std::asin(std::min(1.0, invariant / top_radius));
        return true;
    }

    /**
     * Starting zenith angle of the refracted ray whose tangent point lies on
     * the surface. Rays with larger zenith angles intersect the surface.
     */
    double grazing_zenith(const sasktran2::Geometry1D& geometry, double radius,
                          double start_refractive_index,
                          std::vector<std::pair<int, double>>& index_weights) {
        const auto& grid = geometry.altitude_grid().grid();
        const double ground_invariant =
            sasktran2::raytracing::refraction::refractive_index_at_altitude(
                geometry, grid(0), index_weights) *
            (geometry.coordinates().earth_radius() + grid(0));
        const double start_invariant = start_refractive_index * radius;
        if (ground_invariant >= start_invariant) {
            return EIGEN_PI / 2.0;
        }
        return EIGEN_PI - std::asin(ground_invariant / start_invariant);
    }
} // namespace

namespace sasktran2::raytracing::refraction {
    bool refracted_direction_to_sun(
        const sasktran2::Geometry1D& geometry, const Eigen::Vector3d& position,
        Eigen::Vector3d& direction_to_sun,
        std::vector<std::pair<int, double>>& index_weights) {
        const Eigen::Vector3d& sun = geometry.coordinates().sun_unit();
        direction_to_sun = sun;

        const double radius = position.norm();
        const Eigen::Vector3d up = position / radius;
        const double cos_sun = std::clamp(up.dot(sun), -1.0, 1.0);
        if (cos_sun > NADIR_VIEWING_CUTOFF) {
            // The ray tracer does not refract near-radial rays
            return true;
        }

        Eigen::Vector3d horizontal = sun - cos_sun * up;
        const double sin_sun = horizontal.norm();
        if (sin_sun == 0.0) {
            // Sun directly below the point
            return false;
        }
        horizontal /= sin_sun;
        const double geometric_zenith = std::atan2(sin_sun, cos_sun);

        const double start_refractive_index =
            observer_refractive_index(geometry, radius, index_weights);
        const double max_zenith = grazing_zenith(
            geometry, radius, start_refractive_index, index_weights);

        // residual(zenith) = asymptotic angle - geometric zenith, which
        // increases with the starting zenith angle. A radial ray is not bent,
        // so residual(0) = -geometric zenith.
        const auto residual = [&](double zenith, double& value) {
            double angle;
            if (!asymptotic_direction_angle(geometry, radius,
                                            start_refractive_index, zenith,
                                            angle, index_weights)) {
                return false;
            }
            value = angle - geometric_zenith;
            return true;
        };
        // Directly at the grazing limit the tangent point can round below the
        // surface, back off until the ray clears it
        const auto clamped_residual = [&](double& zenith, double& value) {
            double backoff = 1e-12;
            while (!residual(zenith, value)) {
                zenith -= backoff;
                backoff *= 2.0;
                if (zenith <= 0.0) {
                    return false;
                }
            }
            return true;
        };

        constexpr double angle_tolerance = 1e-11;
        double lower = 0.0;
        double lower_residual = -geometric_zenith;
        bool at_grazing_limit = geometric_zenith >= max_zenith;
        double upper = std::min(geometric_zenith, max_zenith);
        double upper_residual;
        if (!clamped_residual(upper, upper_residual)) {
            return false;
        }

        // Refraction normally bends the ray towards the surface, so the
        // apparent zenith angle is smaller than the geometric one. A
        // refractive index that increases with altitude bends rays the other
        // way, in which case the bracket is extended upwards.
        while (upper_residual < -angle_tolerance) {
            if (at_grazing_limit) {
                // Even the grazing ray does not reach the sun
                return false;
            }
            lower = upper;
            lower_residual = upper_residual;
            upper -= 2.0 * upper_residual;
            if (upper >= max_zenith) {
                upper = max_zenith;
                at_grazing_limit = true;
            }
            if (!clamped_residual(upper, upper_residual)) {
                return false;
            }
        }

        double zenith = upper;
        if (upper_residual > angle_tolerance) {
            // Illinois variant of regula falsi on the bracket [lower, upper]
            int last_side = 0;

            // The bending varies slowly with the zenith angle, so subtracting
            // the bending at the upper bound is an accurate first guess
            zenith = std::clamp(upper - upper_residual, lower, upper);
            for (int iteration = 0; iteration < 100; ++iteration) {
                double value;
                if (!residual(zenith, value)) {
                    // Only possible from roundoff at the grazing limit
                    value = std::numeric_limits<double>::infinity();
                }
                if (std::abs(value) <= angle_tolerance) {
                    break;
                }
                if (value < 0.0) {
                    lower = zenith;
                    lower_residual = value;
                    if (last_side < 0) {
                        upper_residual /= 2.0;
                    }
                    last_side = -1;
                } else {
                    upper = zenith;
                    upper_residual = std::isfinite(value) ? value : 1.0;
                    if (last_side > 0) {
                        lower_residual /= 2.0;
                    }
                    last_side = 1;
                }
                if (upper - lower <=
                    4.0 * std::numeric_limits<double>::epsilon()) {
                    break;
                }
                zenith = (lower * upper_residual - upper * lower_residual) /
                         (upper_residual - lower_residual);
            }
        }

        direction_to_sun =
            std::cos(zenith) * up + std::sin(zenith) * horizontal;
        return true;
    }

    double refracted_visibility_limit(
        const sasktran2::Geometry1D& geometry, double radius_m,
        std::vector<std::pair<int, double>>& index_weights) {
        const double start_refractive_index =
            observer_refractive_index(geometry, radius_m, index_weights);
        double zenith = grazing_zenith(geometry, radius_m,
                                       start_refractive_index, index_weights);

        // Directly at the grazing limit the tangent point can round below the
        // surface, back off until the ray clears it. Upward rays always do.
        double angle;
        double backoff = 1e-12;
        while (!asymptotic_direction_angle(geometry, radius_m,
                                           start_refractive_index, zenith,
                                           angle, index_weights)) {
            zenith -= backoff;
            backoff *= 2.0;
        }
        return angle;
    }

    std::pair<double, double>
    integrate_path(const sasktran2::Geometry1D& geometry, double rt, double nt,
                   double r1, double r2,
                   std::vector<std::pair<int, double>>& index_weights) {
        const double min_cell_length = 0.1;
        const int num_integration_points = 64;

        // Always make sure r1 < r2
        if (r2 < r1) {
            std::swap(r1, r2);
        }

        // The integrands are regular at the tangent point in x = sqrt(r - rt),
        // so integrate all the way down to it. Roundoff can place r1 or r2
        // marginally below rt.
        r1 = std::max(r1, rt);
        r2 = std::max(r2, rt);

        std::pair<double, double> result; // [path_length, path_angle]
        result.first = 0;
        result.second = 0;

        if (std::abs(r1 - r2) < min_cell_length) {
            result.first = straight_path(rt, r1, r2);
            // Closed form expression for the deflection angle assuming n=nt=1
            result.second =
                acos(std::min(1.0, rt / r2)) - acos(std::min(1.0, rt / r1));

            return result;
        }

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

            result.first += w * path_integrand(x1, n1, rt, nt);
            result.first += w * path_integrand(x2, n2, rt, nt);
            result.first += w * path_integrand(x3, n3, rt, nt);
            result.first += w * path_integrand(x4, n4, rt, nt);

            result.second += w * angle_integrand(x1, n1, rt, nt);
            result.second += w * angle_integrand(x2, n2, rt, nt);
            result.second += w * angle_integrand(x3, n3, rt, nt);
            result.second += w * angle_integrand(x4, n4, rt, nt);
        }
        result.first *= diff_d2;
        result.second *= diff_d2;

        // The integral only gives the extra curvature path length, so we add on
        // the staright line path length
        result.first += straight_path(rt, r1, r2);

        if (std::isnan(result.first) || std::isnan(result.second)) {
            spdlog::warn("NaN encountered in refraction integrals: r1={}, "
                         "r2={}, rt={}, nt={}",
                         r1, r2, rt, nt);
            result.first = 0;
            result.second = 0;
        }

        return result;
    }

}; // namespace sasktran2::raytracing::refraction
