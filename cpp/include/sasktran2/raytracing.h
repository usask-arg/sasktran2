#pragma once
#include "sasktran2/grids.h"
#include <sasktran2/internal_common.h>
#include <sasktran2/math/scattering.h>

#include <sasktran2/geometry.h>
#include <sasktran2/viewinggeometry.h>

namespace sasktran2::raytracing {
    /** Defines the type of layer.  Knowing the layer type makes some ray
     * tracing and interpolation easier later on.
     *
     * 'complete' Indicates that the end points of the layer are directly on
     * grid point boundaries
     *
     * 'partial' Means the layer either begins or ends or both somewhere inside
     * of a layer
     *
     * 'tangent' Means that either the start or end point of the layer is the
     * tangent point.
     *
     */
    enum LayerType { complete, partial, tangent };

    /** Heap-free interpolation weights on one structured 1D cell. */
    struct InterpolationStencil1D {
        int cell_index = -1;
        double upper_weight = 0.0;

        bool valid() const { return cell_index >= 0; }
        bool empty() const { return !valid(); }

        std::size_t size() const {
            if (!valid()) {
                return 0;
            }
            return upper_weight == 0.0 || upper_weight == 1.0 ? 1 : 2;
        }

        std::pair<int, double> operator[](std::size_t index) const {
            assert(index < size());
            if (size() == 1) {
                return upper_weight == 0.0
                           ? std::make_pair(cell_index, 1.0)
                           : std::make_pair(cell_index + 1, 1.0);
            }
            return index == 0 ? std::make_pair(cell_index, 1.0 - upper_weight)
                              : std::make_pair(cell_index + 1, upper_weight);
        }

        void reset() {
            cell_index = -1;
            upper_weight = 0.0;
        }

        void assign(const std::vector<std::pair<int, double>>& weights,
                    int grid_size) {
            reset();
            if (grid_size < 2) {
                return;
            }
            int first_index = grid_size;
            int last_index = -1;
            for (const auto& weight : weights) {
                if (weight.second != 0.0) {
                    first_index = std::min(first_index, weight.first);
                    last_index = std::max(last_index, weight.first);
                }
            }
            if (last_index < 0 || first_index < 0 || first_index >= grid_size) {
                return;
            }

            cell_index = std::min(first_index, grid_size - 2);
            if (last_index > cell_index + 1) {
                reset();
                return;
            }

            double lower_weight = 0.0;
            upper_weight = 0.0;
            for (const auto& weight : weights) {
                if (weight.first == cell_index) {
                    lower_weight += weight.second;
                } else if (weight.first == cell_index + 1) {
                    upper_weight += weight.second;
                }
            }
            if (lower_weight + upper_weight != 1.0) {
                reset();
            }
        }
    };

    /** Integrated optical-depth weights for one structured 1D cell.
     *
     * A traced 1D layer is confined to one altitude cell, so its optical depth
     * can be written directly as weights on the cell's lower and upper grid
     * points.  Keeping these combined weights avoids repeating endpoint
     * interpolation when constructing geometry matrices.
     */
    struct IntegratedCellPath1D {
        int cell_index = -1;
        std::array<double, 2> weights = {0.0, 0.0};

        bool valid() const { return cell_index >= 0; }

        void reset() {
            cell_index = -1;
            weights = {0.0, 0.0};
        }
    };

    /** Cell-local interpolation coordinates for a structured 2D endpoint. */
    struct InterpolationCoordinates2D {
        double altitude_upper_weight = 0.0;
        double horizontal_upper_weight = 0.0;

        std::array<double, 4> weights() const {
            const double altitude_lower = 1.0 - altitude_upper_weight;
            const double horizontal_lower = 1.0 - horizontal_upper_weight;
            return {horizontal_lower * altitude_lower,
                    horizontal_lower * altitude_upper_weight,
                    horizontal_upper_weight * altitude_lower,
                    horizontal_upper_weight * altitude_upper_weight};
        }
    };

    struct IntegratedCellPath2D {
        std::array<double, 4> weights = {0.0, 0.0, 0.0, 0.0};
    };

    // Values needed to integrate a layer.  Rather than interpolating anything
    // we directly store weights/indices to the optical table to speed up
    // calculations over multiple wavelengths and allow for easier calculation
    // of derivatives.  Everything here is wavelength independent

    /** Quantities that define a layer along a ray as well as any necessary
     * quantities to perform integrations over the layer.  Rather than
     * interpolating directly we store weights/indices to the atmospher tables
     * to speed up calculations over multiple wavelengths as well as to be able
     * to calculate derivatives easier.  Everything in this structure is
     * considered to be wavelength independent and only a function of geometry.
     */
    struct LayerGeometry {
        Location entrance; /**< The boundary location of the layer closer to the
                              observer */
        Location exit; /**< The boundary location of the layer farther from the
                          observer */

        double r_entrance; /**< Radius of the entrance point in [m] */
        double r_exit;     /**< Radius of the exit point in [m] */

        Eigen::Vector3d
            average_look_away; /**< Look vector away from the observer */

        double layer_distance;   /**< Total distance of the ray within the layer
                                    in [m], assuming not refracted */
        double curvature_factor; /**< 1 for unrefracted, the effective path
                                    length in the layer is layer_distance *
                                    curvature_factor */

        // Quadrature parameters
        double od_quad_start; /**< The OD in the layer is od_quad_start *
                                 k_entrance + od_quad_end * k_exit */
        double od_quad_end;   /**< The OD in the layer is od_quad_start *
                                 k_entrance + od_quad_end * k_exit */

        double od_quad_start_fraction; /**< The fraction the start matters */
        double od_quad_end_fraction;   /**< The fraction the end matters */

        double saz_entrance; /**< Relative solar azimuth angle at entrance in
                                radians, unrefracted */
        double saz_exit; /**< Relative solar azimuth angle at exit in radians,
                            unrefracted */
        double cos_sza_entrance; /**< Cosine solar zenith at entrance,
                                    unrefracted */
        double cos_sza_exit; /**< Cosine solar zenith at exit, unrefracted */

        LayerType type; /**< The type of layer, complete, partial, or tangent */
    };

    struct SphericalLayer : public LayerGeometry {
        InterpolationStencil1D entrance_interpolation_weights;
        InterpolationStencil1D exit_interpolation_weights;

        IntegratedCellPath1D integrated_od; /**< Combined grid-point weights
                                               for optical depth in 1D. */
    };

    /** The geometry information of a fully traced ray.
     */
    struct TracedRay {
        sasktran2::viewinggeometry::ViewingRay
            observer_and_look; /**< The observer location and look direction */

        bool is_straight;   /**< True if the ray is straight, i.e., the look
                               vector does not change along the ray */
        bool ground_is_hit; /**< True if the ground is hit by the ray */

        std::vector<SphericalLayer>
            layers; /**< Set of traced ray layers, starting from the end of the
                       ray moving towards the observer */

        double tangent_radius; /**< Radius of the tangent point in [m] INCLUDING
                                  refractive effects. Could be negative if the
                                  ray hits the ground. */

        std::vector<std::pair<int, double>> interpolation_index_weights; /**<
            Workspace memory to be used by the raytracer.
         */

        /** Resets the storage for the TracedRay so it can be traced again
         */
        void reset() {
            ground_is_hit = false;
            layers.resize(0);
            tangent_radius = std::numeric_limits<double>::quiet_NaN();
        }
    };

    /** A layer traced through one structured 2D cell. */
    struct StructuredLayer2D : public LayerGeometry {
        int altitude_cell = -1; /**< Altitude-cell index for this segment. */
        int horizontal_cell =
            -1; /**< Horizontal-cell index for this segment. */
        InterpolationCoordinates2D entrance_interpolation;
        InterpolationCoordinates2D exit_interpolation;
        IntegratedCellPath2D integrated_od;
    };

    /** Standalone traced-ray output for Geometry2D.
     *
     * This is intentionally separate from TracedRay so the 2D topology is not
     * coupled into the existing 1D engine path.
     */
    struct TracedRay2D {
        sasktran2::viewinggeometry::ViewingRay observer_and_look;
        bool is_straight = true;
        bool ground_is_hit = false;
        std::vector<StructuredLayer2D> layers;
        double tangent_radius = std::numeric_limits<double>::quiet_NaN();

        void reset() {
            ground_is_hit = false;
            layers.resize(0);
            tangent_radius = std::numeric_limits<double>::quiet_NaN();
        }
    };

    /** Calculates cosine of solar zenith angle and relative solar azimuth angle
     * for a given location, look vector, and sun position
     *
     * @param sun_unit Sun position
     * @param loc Location
     * @param look_away look vector away from loc
     * @param csz Returned cosine solar zenith angle
     * @param saa Returned relative solar azimuth angle in radians
     */
    inline void calculate_csz_saz(
        const Eigen::Vector3d& sun_unit, const Location& loc,
        const Eigen::Vector3d& look_away, double& csz, double& saa,
        sasktran2::geometrytype geo = sasktran2::geometrytype::spherical) {
        Eigen::Vector3d local_up;

        if (geo == sasktran2::geometrytype::spherical) {
            local_up = loc.position.normalized();
        } else if (geo == sasktran2::geometrytype::planeparallel ||
                   geo == sasktran2::geometrytype::pseudospherical) {
            // TODO: Is this correct?
            local_up = Eigen::Vector3d(0, 0, 1);
        } else {
            spdlog::error(
                "calculate_csz_saz does not support this geometry type");
        }

        csz = local_up.dot(sun_unit);

        Eigen::Vector3d los_projected =
            (look_away - local_up * (look_away.dot(local_up))).normalized();
        Eigen::Vector3d sun_projceted =
            (sun_unit - local_up * (sun_unit.dot(local_up))).normalized();

        // Take sun to by the x axis, then the y axis is up cross sun
        Eigen::Vector3d y_axis = local_up.cross(sun_projceted);

        // put -1 since these are look away
        saa =
            atan2(y_axis.dot(los_projected), sun_projceted.dot(los_projected));
    }

    /** Populates a layer with the solar parameters
     *
     * @param sun_unit Unit vector to the sun
     * @param layer layer to populate
     */
    inline void add_solar_parameters(
        const Eigen::Vector3d& sun_unit, SphericalLayer& layer,
        sasktran2::geometrytype geo = sasktran2::geometrytype::spherical) {
        // Set the entrance quantities
        calculate_csz_saz(sun_unit, layer.entrance, layer.average_look_away,
                          layer.cos_sza_entrance, layer.saz_entrance, geo);

        // Set the exit quantities
        calculate_csz_saz(sun_unit, layer.exit, layer.average_look_away,
                          layer.cos_sza_exit, layer.saz_exit, geo);
    }

    inline void
    add_interpolation_weights(SphericalLayer& layer,
                              const sasktran2::Geometry1D& geometry,
                              std::vector<std::pair<int, double>>& workspace) {
        geometry.assign_interpolation_weights(layer.exit, workspace);
        layer.exit_interpolation_weights.assign(workspace, geometry.size());
        geometry.assign_interpolation_weights(layer.entrance, workspace);
        layer.entrance_interpolation_weights.assign(workspace, geometry.size());
    }

    /** Combines endpoint interpolation and quadrature into one 1D cell
     * stencil.  If a layer does not satisfy the structured-cell invariant, the
     * stencil remains invalid and matrix construction falls back to the legacy
     * endpoint calculation.
     */
    inline void
    add_integrated_od_weights(SphericalLayer& layer,
                              const sasktran2::Geometry1D& geometry) {
        layer.integrated_od.reset();

        int first_index = geometry.size();
        int last_index = -1;
        const auto add_support =
            [&first_index, &last_index](const auto& interpolation_weights,
                                        double quadrature_weight) {
                for (std::size_t index = 0;
                     index < interpolation_weights.size(); ++index) {
                    const auto iw = interpolation_weights[index];
                    if (iw.second * quadrature_weight != 0.0) {
                        first_index = std::min(first_index, iw.first);
                        last_index = std::max(last_index, iw.first);
                    }
                }
            };
        add_support(layer.entrance_interpolation_weights, layer.od_quad_start);
        add_support(layer.exit_interpolation_weights, layer.od_quad_end);

        if (last_index < 0) {
            first_index = std::min(layer.entrance.lower_alt_index,
                                   layer.exit.lower_alt_index);
            if (first_index < 0) {
                return;
            }
            first_index = std::min(first_index, geometry.size() - 2);
            last_index = first_index;
        }
        if (first_index < 0 || first_index >= geometry.size()) {
            return;
        }

        const int cell_index = std::min(first_index, geometry.size() - 2);
        if (last_index > cell_index + 1) {
            return;
        }

        layer.integrated_od.cell_index = cell_index;
        for (std::size_t index = 0;
             index < layer.entrance_interpolation_weights.size(); ++index) {
            const auto iw = layer.entrance_interpolation_weights[index];
            const double weight = iw.second * layer.od_quad_start;
            if (weight != 0.0) {
                layer.integrated_od.weights[iw.first - cell_index] += weight;
            }
        }
        for (std::size_t index = 0;
             index < layer.exit_interpolation_weights.size(); ++index) {
            const auto iw = layer.exit_interpolation_weights[index];
            const double weight = iw.second * layer.od_quad_end;
            if (weight != 0.0) {
                layer.integrated_od.weights[iw.first - cell_index] += weight;
            }
        }
    }

    /** Populates a layer with the optical depth quadrature parameters
     *
     * @param layer layer to populate
     */
    inline void add_od_quadrature(
        SphericalLayer& layer,
        sasktran2::geometrytype geometry_type =
            sasktran2::geometrytype::spherical,
        grids::interpolation interpolation = grids::interpolation::linear) {
        double r0 = layer.entrance.radius();
        double r1 = layer.exit.radius();
        double dr = r1 - r0;

        layer.average_look_away =
            (layer.exit.position - layer.entrance.position).normalized();

        if (interpolation == grids::interpolation::lower) {
            // For lower interpolation we just use the point with the lower
            // altitude
            if (r0 < r1) {
                layer.od_quad_start =
                    layer.layer_distance * layer.curvature_factor;
                layer.od_quad_end = 0.0;

                layer.od_quad_start_fraction = 0.5;
                layer.od_quad_end_fraction = 0.5;
            } else {
                layer.od_quad_start = 0.0;
                layer.od_quad_end =
                    layer.layer_distance * layer.curvature_factor;

                layer.od_quad_start_fraction = 0.5;
                layer.od_quad_end_fraction = 0.5;
            }

            return;
        }

        if (abs(dr) < 0.001 ||
            geometry_type == sasktran2::geometrytype::planeparallel ||
            interpolation == grids::interpolation::shell) {
            // Tiny layer, we can just use the average of the extinctions
            // Or we are doing shell interpolation
            layer.od_quad_start =
                layer.layer_distance / 2 * layer.curvature_factor;
            layer.od_quad_end =
                layer.layer_distance / 2 * layer.curvature_factor;

            layer.od_quad_start_fraction = 0.5;
            layer.od_quad_end_fraction = 0.5;
        } else {
            double costheta0 =
                layer.entrance.cos_zenith_angle(layer.average_look_away);
            double costheta1 =
                layer.exit.cos_zenith_angle(layer.average_look_away);

            double t0 = r0 * costheta0;
            double t1 = r1 * costheta1;
            double rt = r0 * sqrt(1.0 - costheta0 * costheta0);

            double dt1;
            double dt2;

            if (t1 >= t0) {
                dt1 = t1 - t0;
                if (abs(rt) < 10) {
                    dt2 = 0.5 * ((r1 * t1 - r0 * t0));
                } else {
                    dt2 = 0.5 * ((r1 * t1 - r0 * t0) +
                                 rt * rt * log((r1 + t1) / (r0 + t0)));
                }
            } else {
                dt1 = t0 - t1;
                if (abs(rt) < 10) {
                    dt2 = 0.5 * ((r0 + t0 - r1 * t1));
                } else {
                    dt2 = 0.5 * ((r0 * t0 - r1 * t1) +
                                 rt * rt * log((r0 + t0) / (r1 + t1)));
                }
            }
            layer.od_quad_start =
                (r1 * dt1 - dt2) / dr * layer.curvature_factor;
            layer.od_quad_end =
                -1 * (r0 * dt1 - dt2) / dr * layer.curvature_factor;

            layer.od_quad_start_fraction =
                layer.od_quad_start / (layer.od_quad_start + layer.od_quad_end);
            layer.od_quad_end_fraction =
                layer.od_quad_end / (layer.od_quad_start + layer.od_quad_end);
        }

#ifdef SASKTRAN_DEBUG_ASSERTS
        if (std::isnan(layer.od_quad_start) || std::isnan(layer.od_quad_end) ||
            std::isnan(layer.od_quad_end_fraction) ||
            std::isnan(layer.od_quad_end_fraction)) {
            spdlog::error("One of layer quadrature parameters was nan");

            spdlog::error("od_quad_start: {}, od_quad_end: {}",
                          layer.od_quad_start, layer.od_quad_end);
        }
#endif
    }

    /** Takes a list of traced rays and returns the min and max COS sza of any
     * layer in the traced rays.
     *
     *  Usually this function would be used in constructing grids for various
     * source function discretizations. For example it can be used in
     * constructing a solar transmission table to determine what SZA's need to
     * be used. Or in construction of the diffuse HR source.
     *
     * @param rays
     * @return (min_cos_sza, max_cos_sza)
     */
    inline std::pair<double, double>
    min_max_cos_sza_of_all_rays(const std::vector<TracedRay>& rays) {
        double min_cos_sza = 1;
        double max_cos_sza = -1;

        for (const auto& ray : rays) {
            for (const auto& layer : ray.layers) {
                min_cos_sza = std::min(layer.cos_sza_entrance, min_cos_sza);
                min_cos_sza = std::min(layer.cos_sza_exit, min_cos_sza);

                max_cos_sza = std::max(layer.cos_sza_entrance, max_cos_sza);
                max_cos_sza = std::max(layer.cos_sza_exit, max_cos_sza);
            }
        }

        return std::make_pair(min_cos_sza, max_cos_sza);
    }

    /** A traced_ray has \f$N_l\f$ layers.  The vector of optical depths for
     * these layers can always we written as \f$M k\f$ where \f$M\f$ is a sparse
     * matrix, and \f$k\f$ is the vector of extinction coefficients.  This
     * method calculates the sparse matrix \f$M\f$.
     *
     *  The primary advantage of using this matrix formalism is that the matrix
     * only needs to be calculated once and is purely geometry dependent.
     * Knowing the matrix also allows for easy calculation of derivatives of
     *  layer optical depths.
     *
     * @param traced_ray A traced ray
     * @param geometry An object defining the model geometry, needed for
     * interpolation
     * @param result The resulting sparse matrix \f$M\f$
     */
    inline void
    construct_od_matrix(const TracedRay& traced_ray, const Geometry& geometry,
                        Eigen::SparseMatrix<double, Eigen::RowMajor>& result) {
        // We want to construct a matrix such that A * extinction = layer OD
        // The matrix should have size (nlayer, nextinction)
        result.resize((long)traced_ray.layers.size(), geometry.size());

        typedef Eigen::Triplet<double> T;
        std::vector<std::pair<int, double>> index_weights;

        std::vector<T> tripletList;
        tripletList.reserve(traced_ray.layers.size() * 2);

        for (int i = (int)traced_ray.layers.size() - 1; i >= 0; --i) {
            const auto& layer = traced_ray.layers[i];

            if (layer.integrated_od.valid()) {
                for (int j = 0; j < 2; ++j) {
                    if (layer.integrated_od.weights[j] != 0.0) {
                        tripletList.emplace_back(
                            i, layer.integrated_od.cell_index + j,
                            layer.integrated_od.weights[j]);
                    }
                }
                continue;
            }

            geometry.assign_interpolation_weights(layer.entrance,
                                                  index_weights);

            for (const auto& iw : index_weights) {
                tripletList.push_back(
                    T(i, iw.first, iw.second * layer.od_quad_start));
            }

            geometry.assign_interpolation_weights(layer.exit, index_weights);

            for (const auto& iw : index_weights) {
                tripletList.push_back(
                    T(i, iw.first, iw.second * layer.od_quad_end));
            }
        }
        result.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    /** Constructs the four-node-per-layer optical-depth matrix for a
     * standalone structured 2D traced ray. */
    inline void
    construct_od_matrix(const TracedRay2D& traced_ray,
                        const Geometry2D& geometry,
                        Eigen::SparseMatrix<double, Eigen::RowMajor>& result) {
        result.resize(static_cast<long>(traced_ray.layers.size()),
                      geometry.size());
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(traced_ray.layers.size() * 4);

        for (int layer_index = 0; layer_index < traced_ray.layers.size();
             ++layer_index) {
            const auto& layer = traced_ray.layers[layer_index];
            const std::array<int, 4> location_indices = {
                geometry.location_index(layer.altitude_cell,
                                        layer.horizontal_cell),
                geometry.location_index(layer.altitude_cell + 1,
                                        layer.horizontal_cell),
                geometry.location_index(layer.altitude_cell,
                                        layer.horizontal_cell + 1),
                geometry.location_index(layer.altitude_cell + 1,
                                        layer.horizontal_cell + 1)};
            for (int local_index = 0; local_index < 4; ++local_index) {
                if (layer.integrated_od.weights[local_index] != 0.0) {
                    triplets.emplace_back(
                        layer_index, location_indices[local_index],
                        layer.integrated_od.weights[local_index]);
                }
            }
        }
        result.setFromTriplets(triplets.begin(), triplets.end());
    }

    /** Abstract base interface class for a ray tracer

     */
    class RayTracerBase {
      public:
        virtual ~RayTracerBase(){};

        /** Traces a ray
         *
         * @param ray ViewingRay to trace
         * @param tracedray Result
         * @param include_refraction If true, the ray will be traced with
         * refractive effects enabled
         */
        virtual void
        trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray,
                  TracedRay& tracedray,
                  bool include_refraction = false) const = 0;
    };

    /** A ray tracer that traces straight (unrefracted) rays in a spherical
     * shell (1d) atmosphere.
     */
    class SphericalShellRayTracer : public RayTracerBase {
      public:
        /** Constructs the spherical shell ray tracer.  Note that we always
         * trace rays on the same grid that the extinction profile is specified
         * on, i.e., the global model geometry grid.
         *
         * @param geometry Geometry object that is used to obtain the altitude
         * grid as well as earth radius
         */
        SphericalShellRayTracer(const sasktran2::Geometry1D& geometry)
            : m_alt_grid(geometry.altitude_grid()),
              m_earth_radius(geometry.coordinates().earth_radius()),
              m_geometry(geometry) {}

        void trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray,
                       TracedRay& tracedray,
                       bool include_refraction = false) const override;

      private:
        const sasktran2::grids::AltitudeGrid&
            m_alt_grid; /**< Internal altitude grid */
        const sasktran2::Geometry1D&
            m_geometry;              /**< Internal geometry reference */
        const double m_earth_radius; /**< earth radius in [m] */

        /** Enum which defines if locally the look vector is looking up or down
         */
        enum ViewingDirection { up = -1, down = 1 };

        /** Enum which defines if we are on the far side or near side of the
         * tangent point
         */
        enum TangentSide { farside = -1, nearside = 1 };

        /** Traces a ray where the observes is outside of the atmosphere looking
         * at the ground
         *
         * @param ray
         * @param tracedray
         */
        void trace_ray_observer_outside_ground_viewing(
            const sasktran2::viewinggeometry::ViewingRay& ray,
            TracedRay& tracedray) const;

        /** Traces a ray where the observer is inside the atmosphere, looking
         * upwards (no tangent layer)
         *
         * @param ray
         * @param tracedray
         */
        void trace_ray_observer_inside_looking_up(
            const sasktran2::viewinggeometry::ViewingRay& ray,
            TracedRay& tracedray) const;

        /** Traces a ray where the observer is inside the atmosphere looking at
         * the ground
         *
         * @param ray
         * @param tracedray
         */
        void trace_ray_observer_inside_looking_ground(
            const sasktran2::viewinggeometry::ViewingRay& ray,
            TracedRay& tracedray) const;

        /** Traces a ray where the observer is inside the atmosphere, looking
         * through the limb, towards space
         *
         * @param ray
         * @param tracedray
         */
        void trace_ray_observer_inside_looking_limb(
            const sasktran2::viewinggeometry::ViewingRay& ray,
            TracedRay& tracedray) const;

        /** Traces a ray where the observer is outside the atmosphere, limb
         * viewing.
         *
         * @param ray
         * @param tracedray
         */
        void trace_ray_observer_outside_limb_viewing(
            const sasktran2::viewinggeometry::ViewingRay& ray,
            TracedRay& tracedray) const;

        /** Internal method to construct a 'complete' layer. This is the most
         * common case
         *
         * @param layer layer to construct
         * @param ray LOS ray
         * @param exit_index index of the exit point
         * @param direction direction we are viewing
         * @param side Side of the tangent point we are on
         */
        void complete_layer(SphericalLayer& layer,
                            const sasktran2::viewinggeometry::ViewingRay& ray,
                            size_t exit_index, ViewingDirection direction,
                            TangentSide side) const;

        /** Internal method to construct a 'partial' layer.  This case usually
         * happens when the observer is inside the atmosphere
         *
         * @param layer layer to construct
         * @param ray LOS ray
         * @param start_index index of the start point
         * @param direction directino we are viewing
         * @param side Side of the tangent point we are on
         */
        void partial_layer(SphericalLayer& layer,
                           const sasktran2::viewinggeometry::ViewingRay& ray,
                           size_t start_index, ViewingDirection direction,
                           TangentSide side) const;

        /** Internal method to construct a complete tangent layer.  This
         * involves a layer that spans from one boundary point to the tangent
         * point.
         *
         * @param layer layer to construct
         * @param ray LOS ray
         * @param upper_index index to the point above the tangent layer
         * @param tangent_altitude tangent altitude in [m]
         * @param direction direction we are viewing
         * @param side side of the tangent point we are on
         */
        void tangent_layer(SphericalLayer& layer,
                           const sasktran2::viewinggeometry::ViewingRay& ray,
                           size_t upper_index, double tangent_altitude,
                           ViewingDirection direction, TangentSide side) const;

        /** Internal method to construct a partial tangent layer, this happens
         * when the observer is within the tangent layer.
         *
         * @param layer layer to construct
         * @param ray LOS ray
         * @param start_index ?
         * @param tangent_altitude Tangent altitude in [m]
         * @param direction direction we are viewing
         * @param side side of the tangent point we are on
         */
        void
        partial_tangent_layer(SphericalLayer& layer,
                              const sasktran2::viewinggeometry::ViewingRay& ray,
                              size_t start_index, double tangent_altitude,
                              ViewingDirection direction,
                              TangentSide side) const;

        /**
         * Finalizes the ray geometry by calculating the intersection points for
         * each layer, the quadrature parameters, and the solar angles
         *
         * @param tracedray
         * @param rt
         */
        void finalize_ray_geometry(TracedRay& tracedray) const;

        /** Converts altitude to distance along the ray
         *
         * @param ray viewing ray
         * @param altitude altitude
         * @param direction viewing direction
         * @param side tangent side we are on
         * @return distance along ray to the specified altitude
         */
        double
        distance_to_altitude(const sasktran2::viewinggeometry::ViewingRay& ray,
                             double altitude, ViewingDirection direction,
                             TangentSide side) const {
            double cos_zenith =
                abs(ray.observer.cos_zenith_angle(ray.look_away));
            double ro = ray.observer.radius();
            double re = m_earth_radius + altitude;

            double rtsq = ro * ro * (1 - cos_zenith * cos_zenith);

            // Distance to the tangent point
            double tangent_distance = side * direction * ro * cos_zenith;

            // Distance from the tangent point to the requested altitude
            double dist_from_tangent;
            if (rtsq > re * re) {
                // Either a rounding error or a problem, TODO: This needs to be
                // figured out?
                if (abs(rtsq - re * re) < 100) {
                    dist_from_tangent = 0.0;
                } else {
                    throw("Error, requesting distance to a shell that does not "
                          "exist");
                }
            } else {
                dist_from_tangent =
                    side * direction * sqrt(abs(re * re - rtsq));
            }

#ifdef SASKTRAN_DEBUG_ASSERTS
            if (std::isnan(tangent_distance) || std::isnan(dist_from_tangent)) {
                spdlog::error("Error computing tangent distances");
            }
#endif

            if (side == TangentSide::nearside) {
                return tangent_distance - dist_from_tangent;
            } else {
                return tangent_distance + dist_from_tangent;
            }
        }
    };

    /** A ray tracer that traces straight (unrefracted) rays in a Plane Parallel
     * (1d) atmosphere
     */
    class PlaneParallelRayTracer : public RayTracerBase {
      public:
        /** Constructs the plane parallel ray tracer.  Note that we always trace
         * rays on the same grid that the extinction profile is specified on,
         * i.e., the global model geometry grid.

         *
         * @param geometry Geometry object that is used to obtain the altitude
         * grid as well as earth radius
         */
        PlaneParallelRayTracer(const sasktran2::Geometry1D& geometry)
            : m_alt_grid(geometry.altitude_grid()),
              m_earth_radius(geometry.coordinates().earth_radius()),
              m_geometry(geometry) {}

        void trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray,
                       TracedRay& tracedray,
                       bool include_refraction = false) const;

      private:
        const sasktran2::grids::AltitudeGrid&
            m_alt_grid; /**< Internal altitude grid */
        const sasktran2::Geometry1D&
            m_geometry;              /**< Internal geometry reference */
        const double m_earth_radius; /**< earth radius in [m] */

        /** Enum which defines if locally the look vector is looking up or down
         */
        enum ViewingDirection { up = -1, down = 1 };

        /** Traces a ray where the observes is outside of the atmosphere looking
         * at the ground
         *
         * @param ray
         * @param tracedray
         */
        void trace_ray_observer_outside_looking_ground(
            const sasktran2::viewinggeometry::ViewingRay& ray,
            TracedRay& tracedray) const;

        /** Traces a ray where the observer is inside the atmosphere, looking
         * upwards
         *
         * @param ray
         * @param tracedray
         */
        void trace_ray_observer_inside_looking_up(
            const sasktran2::viewinggeometry::ViewingRay& ray,
            TracedRay& tracedray) const;

        /** Traces a ray where the observer is inside the atmosphere looking at
         * the ground
         *
         * @param ray
         * @param tracedray
         */
        void trace_ray_observer_inside_looking_ground(
            const sasktran2::viewinggeometry::ViewingRay& ray,
            TracedRay& tracedray) const;

        /** Internal method to construct a 'complete' layer. This is the most
         * common case
         *
         * @param layer layer to construct
         * @param ray LOS ray
         * @param exit_index index of the exit point
         * @param direction direction we are viewing
         * @param side Side of the tangent point we are on
         */
        void complete_layer(SphericalLayer& layer,
                            const sasktran2::viewinggeometry::ViewingRay& ray,
                            size_t exit_index,
                            ViewingDirection direction) const;

        /** Internal method to construct a 'partial' layer.  This case usually
         * happens when the observer is inside the atmosphere
         *
         * @param layer layer to construct
         * @param ray LOS ray
         * @param start_index index of the start point
         * @param direction directino we are viewing
         * @param side Side of the tangent point we are on
         */
        void partial_layer(SphericalLayer& layer,
                           const sasktran2::viewinggeometry::ViewingRay& ray,
                           size_t start_index,
                           ViewingDirection direction) const;

        /** Converts altitude to distance along the ray
         *
         * @param ray viewing ray
         * @param altitude altitude
         * @param direction viewing direction
         * @return distance along ray to the specified altitude
         */
        double
        distance_to_altitude(const sasktran2::viewinggeometry::ViewingRay& ray,
                             double altitude, ViewingDirection direction) const;
    };

#ifdef SKTRAN_RUST_SUPPORT
    class RustRayTracerImpl;
    class RustRayTracer2DImpl;

    /** A ray tracer backed by the Rust vertical-grid ray tracing core.
     */
    class RustRayTracer : public RayTracerBase {
      public:
        RustRayTracer(const sasktran2::Geometry1D& geometry);
        ~RustRayTracer() override;

        void trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray,
                       TracedRay& tracedray,
                       bool include_refraction = false) const override;

      private:
        const sasktran2::Geometry1D& m_geometry;
        std::unique_ptr<RustRayTracerImpl> m_impl;
    };

    /** Standalone Rust structured-2D ray tracer.
     *
     * The refractive-index overload accepts one altitude-only profile for this
     * ray. Profiles are not stored on Geometry2D and may differ between calls.
     */
    class RustRayTracer2D {
      public:
        explicit RustRayTracer2D(const sasktran2::Geometry2D& geometry);
        ~RustRayTracer2D();

        void trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray,
                       TracedRay2D& tracedray) const;

        void trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray,
                       const Eigen::VectorXd& refractive_index,
                       TracedRay2D& tracedray) const;

      private:
        void trace_ray_impl(const sasktran2::viewinggeometry::ViewingRay& ray,
                            const Eigen::VectorXd* refractive_index,
                            TracedRay2D& tracedray) const;

        const sasktran2::Geometry2D& m_geometry;
        std::unique_ptr<RustRayTracer2DImpl> m_impl;
    };
#endif
} // namespace sasktran2::raytracing
