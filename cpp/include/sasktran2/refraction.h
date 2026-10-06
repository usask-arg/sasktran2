#pragma once
#include "sasktran2/grids.h"
#include <array>
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

        // Same weights as Geometry1D::assign_interpolation_weights for a
        // location at this altitude, without constructing the location. This
        // is evaluated at every quadrature node of the refraction integrals.
        std::array<int, 2> index;
        std::array<double, 2> weight;
        int num_contributing;
        geometry.altitude_grid().calculate_interpolation_weights(
            altitude_m, index, weight, num_contributing);
        index_weights.resize(num_contributing);
        for (int i = 0; i < num_contributing; ++i) {
            index_weights[i] = {index[i], weight[i]};
        }

        // Interpolate log of refractive index
        double log_n = 0;
        for (const auto& [index, weight] : index_weights) {
            log_n += weight * log(geometry.refractive_index()[index]);
        }

        return exp(log_n);
    }

    /**
     * Refractive index that enters the ray invariant n r sin(zenith) for a
     * ray starting at the given radius. Points at or above the top of the
     * atmosphere are treated as being in vacuum.
     *
     * @param geometry The geometry object
     * @param radius_m Radius of the ray start point in [m]
     * @param index_weights Workspace memory
     * @return double The refractive index at the start point
     */
    inline double observer_refractive_index(
        const sasktran2::Geometry1D& geometry, double radius_m,
        std::vector<std::pair<int, double>>& index_weights) {
        const double altitude_m =
            radius_m - geometry.coordinates().earth_radius();
        if (altitude_m >=
            geometry.altitude_grid().grid()(Eigen::placeholders::last)) {
            return 1.0;
        }
        return refractive_index_at_altitude(geometry, altitude_m,
                                            index_weights);
    }

    /**
     * Calculates the tangent radius of a ray taking into account refraction.
     * This assumes that the refractive index varies only in altitude.
     *
     * @param geometry The geometry object
     * @param straight_line_tangent_radius_m The ray invariant n r sin(zenith)
     * evaluated at the start of the ray. For a ray starting in vacuum this is
     * the tangent radius the ray would have without refraction.
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

    /**
     * Finds the apparent direction of the sun at a point in a spherically
     * symmetric refracting atmosphere.
     *
     * The returned direction is the local tangent of the refracted ray that
     * leaves the point and, after bending through the atmosphere, travels
     * parallel to the geometric sun direction above the top of the
     * atmosphere. Its negative is the local propagation direction of the
     * direct solar beam. Tracing a refracted ray from the point along the
     * returned direction reproduces the solar path.
     *
     * Points whose geometric sun direction is within the near-zenith cutoff
     * used by the ray tracer return the geometric sun direction unchanged.
     *
     * @param geometry The geometry object, must be spherical
     * @param position Location of the point
     * @param direction_to_sun Output apparent unit direction to the sun
     * @param index_weights Workspace memory
     * @return false if every ray from the point that reaches the sun
     * intersects the surface, i.e. the point is in the refracted shadow of
     * the Earth. direction_to_sun is set to the geometric sun direction in
     * that case.
     */
    bool refracted_direction_to_sun(
        const sasktran2::Geometry1D& geometry, const Eigen::Vector3d& position,
        Eigen::Vector3d& direction_to_sun,
        std::vector<std::pair<int, double>>& index_weights);

    /**
     * Largest geometric solar zenith angle for which a point at the given
     * radius is illuminated through the refracting atmosphere. This is the
     * asymptotic direction of the ray that leaves the point and grazes the
     * surface.
     *
     * @param geometry The geometry object, must be spherical
     * @param radius_m Radius of the point in [m]
     * @param index_weights Workspace memory
     * @return double Geometric solar zenith angle in [rad]
     */
    double refracted_visibility_limit(
        const sasktran2::Geometry1D& geometry, double radius_m,
        std::vector<std::pair<int, double>>& index_weights);
} // namespace sasktran2::raytracing::refraction
