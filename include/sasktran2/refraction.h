#pragma once
#include "sasktran2/grids.h"
#include <sasktran2/internal_common.h>
#include <sasktran2/math/scattering.h>

#include <sasktran2/geometry.h>
#include <sasktran2/viewinggeometry.h>

namespace sasktran2::raytracing::refraction {

    /**
     * @brief Calculates the refractive index at a given altitude by
     * interpolating log(refractive_index) linearly in altitude.  This
     * implicitly assumes a refractive index that varies only in altitude
     *
     * @param geometry The geometry object
     * @param altitude_m Altitude in [m] to calculate the refractive index at
     * @param index_weights Workspace memory
     * @return double The refractive index at altitude altitude_m
     */
    inline double refractive_index_at_altitude(
        const sasktran2::Geometry1D& geometry, double altitude_m,
        std::vector<std::pair<int, double>>& index_weights) {

        // Create a location object with the required altitude
        sasktran2::Location loc;
        loc.position = geometry.coordinates().reference_point(altitude_m);

        // Get interpolation weights
        geometry.assign_interpolation_weights(loc, index_weights);

        // Interpolate log of refractive index
        double log_n = 0;
        for (const auto& [index, weight] : index_weights) {
            log_n += weight * log(geometry.refractive_index()[index]);
        }

        return exp(log_n);
    }

    /**
     * Calculates the tangent radius of a ray taking into account refraction.
     * This assumes that the refractive index varies only in altitude.
     *
     * @param geometry The geometry object
     * @param straight_line_tangent_radius_m The tangent radius of the ray if
     * there were no refraction
     * @param index_weights Workspace memory
     * @return double The tangent radius of the ray taking into account
     * refraction
     */
    inline double
    tangent_radius(const sasktran2::Geometry1D& geometry,
                   double straight_line_tangent_radius_m,
                   std::vector<std::pair<int, double>>& index_weights) {
        const size_t maxiter = 500;
        const double tolerance = 1e-6;

        // There is probably a real formula for this, but this iterative
        // approach is good enough and copied from SASKTRAN1

        // Essentially we want to solve rt = rt / n(rt) where n(rt) is the
        // refractive index at the tangent radius
        size_t currentiter = 0;
        double currentrt = straight_line_tangent_radius_m;
        double nextrt = 0.0;

        while (currentiter < maxiter) {
            double n = refractive_index_at_altitude(
                geometry, currentrt - geometry.coordinates().earth_radius(),
                index_weights);

            nextrt = straight_line_tangent_radius_m / n;

            if (fabs(nextrt - currentrt) < tolerance) {
                break;
            }
            ++currentiter;
            if (currentiter != maxiter) {
                currentrt = nextrt;
            }
        }

        if (currentiter == maxiter) {
            if (fabs(nextrt - currentrt) < 1) {
                // Usually okay
                spdlog::info("Poor convergence of tangent radius");
                currentrt = (nextrt + currentrt) / 2.0;
            } else {
                spdlog::warn("Refractive tangent radius failed to converge");
            }
        }

        return currentrt;
    }

    /**
     * Performs the path and deflection angle integrals found in Thompson 1982
     *
     * @param geometry The base geometry object, necessary to get earth radius,
     * refractive index, etc.
     * @param rt The REFRACTED tangent radius of the ray
     * @param nt The index of refraction at radius rt
     * @param r1 The start radius of the integration
     * @param r2 The end radius of the integration
     * @param index_weights Workspace memory
     * @return std::pair<double, double> (path_length in m, path_angle in m)
     */
    std::pair<double, double>
    integrate_path(const sasktran2::Geometry1D& geometry, double rt, double nt,
                   double r1, double r2,
                   std::vector<std::pair<int, double>>& index_weights);
} // namespace sasktran2::raytracing::refraction
