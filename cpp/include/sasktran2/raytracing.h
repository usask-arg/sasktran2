#pragma once
#include "sasktran2/grids.h"
#include <sasktran2/internal_common.h>
#include <sasktran2/math/scattering.h>

#include <cstdint>

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

    /** A non-owning view over grid interpolation or path-integration weights.
     *
     * The backing storage is owned by TracedRay. Keeping the grid indices
     * separate from the weights allows the same cell-node list to be reused
     * for entrance interpolation, exit interpolation, and integrated optical
     * depth without reserving a maximum-dimensional stencil in every layer.
     */
    class GridWeightStencilView {
      private:
        const int* m_indices = nullptr;
        const double* m_weights = nullptr;
        std::size_t m_size = 0;

      public:
        GridWeightStencilView() = default;
        GridWeightStencilView(const int* indices, const double* weights,
                              std::size_t size)
            : m_indices(indices), m_weights(weights), m_size(size) {}

        std::size_t size() const { return m_size; }
        bool empty() const { return m_size == 0; }

        std::pair<int, double> operator[](std::size_t index) const {
            assert(index < m_size);
            return {m_indices[index], m_weights[index]};
        }
    };

    /** One dimension-independent segment of a traced ray. */
    struct TracedLayer : public LayerGeometry {
        std::uint32_t grid_weight_offset = 0;
        std::uint8_t grid_weight_count = 0;
    };

    // Retain the established 1D name for downstream code. The layer itself is
    // now independent of atmosphere dimensionality.
    using SphericalLayer = TracedLayer;

    /** The geometry information of a fully traced ray.
     */
    struct TracedRay {
        sasktran2::viewinggeometry::ViewingRay
            observer_and_look; /**< The observer location and look direction */

        bool is_straight = true;    /**< True if the ray is straight, i.e., the
                                       look vector does not change along the ray */
        bool ground_is_hit = false; /**< True if the ground is hit by the ray */

        std::vector<TracedLayer>
            layers; /**< Set of traced ray layers, starting from the end of the
                       ray moving towards the observer */

        double tangent_radius =
            std::numeric_limits<double>::quiet_NaN(); /**< Radius of the tangent
                                                         point in [m] INCLUDING
                                                         refractive effects.
                                                         Could be negative if
                                                         the ray hits the
                                                         ground. */

        std::vector<std::pair<int, double>> interpolation_index_weights; /**<
            Workspace memory to be used by the raytracer.
         */

      private:
        std::vector<int> grid_weight_indices;
        std::vector<double> entrance_grid_weights;
        std::vector<double> exit_grid_weights;
        std::vector<double> integrated_od_weights;

      public:
        /** Reserves storage for compiled grid weights without changing the
         * ray. Views obtained before a mutation of the ray must not be kept. */
        void reserve_grid_weights(std::size_t count) {
            grid_weight_indices.reserve(count);
            entrance_grid_weights.reserve(count);
            exit_grid_weights.reserve(count);
            integrated_od_weights.reserve(count);
        }

        std::size_t num_grid_weights() const {
            return grid_weight_indices.size();
        }

        template <typename IndexContainer, typename EntranceContainer,
                  typename ExitContainer, typename ODContainer>
        void set_layer_weights(std::size_t layer_index,
                               const IndexContainer& indices,
                               const EntranceContainer& entrance_weights,
                               const ExitContainer& exit_weights,
                               const ODContainer& od_weights) {
            const std::size_t count = indices.size();
            if (entrance_weights.size() != count ||
                exit_weights.size() != count || od_weights.size() != count) {
                throw std::invalid_argument(
                    "Traced-ray layer weight arrays must have equal sizes");
            }
            set_layer_weights(layer_index, indices.data(),
                              entrance_weights.data(), exit_weights.data(),
                              od_weights.data(), count);
        }

        void set_layer_weights(std::size_t layer_index, const int* indices,
                               const double* entrance_weights,
                               const double* exit_weights,
                               const double* od_weights, std::size_t count) {
            if (count > std::numeric_limits<std::uint8_t>::max() ||
                grid_weight_indices.size() >
                    std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(
                    "Traced-ray grid-weight storage exceeds compact index "
                    "limits");
            }

            auto& layer = layers.at(layer_index);
            layer.grid_weight_offset = grid_weight_indices.size();
            layer.grid_weight_count = count;
            grid_weight_indices.insert(grid_weight_indices.end(), indices,
                                       indices + count);
            entrance_grid_weights.insert(entrance_grid_weights.end(),
                                         entrance_weights,
                                         entrance_weights + count);
            exit_grid_weights.insert(exit_grid_weights.end(), exit_weights,
                                     exit_weights + count);
            integrated_od_weights.insert(integrated_od_weights.end(),
                                         od_weights, od_weights + count);
        }

        GridWeightStencilView entrance_weights(std::size_t layer_index) const {
            const auto& layer = layers.at(layer_index);
            if (layer.grid_weight_count == 0) {
                return {};
            }
            return {grid_weight_indices.data() + layer.grid_weight_offset,
                    entrance_grid_weights.data() + layer.grid_weight_offset,
                    layer.grid_weight_count};
        }

        GridWeightStencilView exit_weights(std::size_t layer_index) const {
            const auto& layer = layers.at(layer_index);
            if (layer.grid_weight_count == 0) {
                return {};
            }
            return {grid_weight_indices.data() + layer.grid_weight_offset,
                    exit_grid_weights.data() + layer.grid_weight_offset,
                    layer.grid_weight_count};
        }

        GridWeightStencilView
        optical_depth_weights(std::size_t layer_index) const {
            const auto& layer = layers.at(layer_index);
            if (layer.grid_weight_count == 0) {
                return {};
            }
            return {grid_weight_indices.data() + layer.grid_weight_offset,
                    integrated_od_weights.data() + layer.grid_weight_offset,
                    layer.grid_weight_count};
        }

        /** Resets the storage for the TracedRay so it can be traced again
         */
        void reset() {
            ground_is_hit = false;
            layers.resize(0);
            grid_weight_indices.clear();
            entrance_grid_weights.clear();
            exit_grid_weights.clear();
            integrated_od_weights.clear();
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

    /** Resolves a 1D layer to flattened atmosphere-grid weights.
     *
     * All geometry-dependent interpolation is completed during tracing. The
     * source and optical-depth paths therefore consume the same representation
     * regardless of atmosphere dimensionality.
     */
    inline void
    add_interpolation_weights(TracedRay& ray, std::size_t layer_index,
                              const sasktran2::Geometry1D& geometry,
                              std::vector<std::pair<int, double>>& workspace) {
        auto& layer = ray.layers.at(layer_index);
        std::array<std::pair<int, double>, 2> entrance;
        std::array<std::pair<int, double>, 2> exit;
        std::size_t entrance_count = 0;
        std::size_t exit_count = 0;
        const auto assign = [&](const Location& location, auto& destination,
                                std::size_t& count) {
            geometry.assign_interpolation_weights(location, workspace);
            if (workspace.size() > destination.size()) {
                throw std::runtime_error(
                    "Traced 1D endpoint spans more than one grid cell");
            }
            count = workspace.size();
            std::copy(workspace.begin(), workspace.end(), destination.begin());
        };
        assign(layer.entrance, entrance, entrance_count);
        assign(layer.exit, exit, exit_count);

        std::array<int, 4> indices;
        std::size_t index_count = 0;
        const auto add_indices = [&indices, &index_count](const auto& stencil,
                                                          std::size_t count) {
            for (std::size_t index = 0; index < count; ++index) {
                if (stencil[index].second == 0.0) {
                    continue;
                }
                const int grid_index = stencil[index].first;
                if (std::find(indices.begin(), indices.begin() + index_count,
                              grid_index) == indices.begin() + index_count) {
                    if (index_count == indices.size()) {
                        throw std::runtime_error(
                            "Traced 1D layer touches more than four grid "
                            "points");
                    }
                    indices[index_count++] = grid_index;
                }
            }
        };
        add_indices(entrance, entrance_count);
        add_indices(exit, exit_count);
        if (index_count == 0) {
            throw std::runtime_error(
                "Unable to resolve traced 1D layer interpolation weights");
        }
        std::sort(indices.begin(), indices.begin() + index_count);

        std::array<double, 4> entrance_weights = {0.0, 0.0, 0.0, 0.0};
        std::array<double, 4> exit_weights = {0.0, 0.0, 0.0, 0.0};
        const auto accumulate = [&indices, index_count](const auto& stencil,
                                                        std::size_t count,
                                                        auto& weights) {
            for (std::size_t index = 0; index < count; ++index) {
                const auto weight = stencil[index];
                if (weight.second == 0.0) {
                    continue;
                }
                const auto position =
                    std::find(indices.begin(), indices.begin() + index_count,
                              weight.first);
                if (position == indices.begin() + index_count) {
                    throw std::runtime_error(
                        "Unable to map traced 1D interpolation weight");
                }
                weights[std::distance(indices.begin(), position)] +=
                    weight.second;
            }
        };
        accumulate(entrance, entrance_count, entrance_weights);
        accumulate(exit, exit_count, exit_weights);

        std::array<double, 4> od_weights = {0.0, 0.0, 0.0, 0.0};
        for (std::size_t index = 0; index < index_count; ++index) {
            od_weights[index] = entrance_weights[index] * layer.od_quad_start +
                                exit_weights[index] * layer.od_quad_end;
        }
        ray.set_layer_weights(layer_index, indices.data(),
                              entrance_weights.data(), exit_weights.data(),
                              od_weights.data(), index_count);
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
        std::vector<T> tripletList;
        tripletList.reserve(traced_ray.num_grid_weights());

        for (int i = (int)traced_ray.layers.size() - 1; i >= 0; --i) {
            const auto weights = traced_ray.optical_depth_weights(i);
            for (std::size_t j = 0; j < weights.size(); ++j) {
                const auto weight = weights[j];
                if (weight.second != 0.0) {
                    tripletList.emplace_back(i, weight.first, weight.second);
                }
            }
        }
        result.setFromTriplets(tripletList.begin(), tripletList.end());
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
                       TracedRay& tracedray) const;

        void trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray,
                       const Eigen::VectorXd& refractive_index,
                       TracedRay& tracedray) const;

      private:
        void trace_ray_impl(const sasktran2::viewinggeometry::ViewingRay& ray,
                            const Eigen::VectorXd* refractive_index,
                            TracedRay& tracedray) const;

        const sasktran2::Geometry2D& m_geometry;
        std::unique_ptr<RustRayTracer2DImpl> m_impl;
    };
#endif
} // namespace sasktran2::raytracing
