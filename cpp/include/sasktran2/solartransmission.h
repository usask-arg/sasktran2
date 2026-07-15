#pragma once

#include "sasktran2/atmosphere/grid_storage.h"
#include "sasktran2/geometry.h"
#include <sasktran2/source_interface.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/config.h>
#include <sasktran2/dual.h>

namespace sasktran2::solartransmission {
    /** Generates one row of what is known as the geometry matrix.  The optical
     * depth at every necessary integration point can be written as \f$M k\f$
     * where \f$k\f$ is the extinction coefficient vector.  In 1D geometry, the
     *  matrix \f$M\f$ is densish, but it may be sparse in higher dimensions.
     * This function constructs one row of the matrix \f$M\f$, which corresponds
     * to one traced ray.
     *
     * @param row row index
     * @param traced_ray ray traced to the sun
     * @param result matrix to store the result in
     */
    inline void assign_dense_matrix_column(
        int row, const sasktran2::raytracing::TracedRay& traced_ray,
        Eigen::MatrixXd& result) {
        for (int i = 0; i < traced_ray.layers.size(); ++i) {
            const auto weights = traced_ray.optical_depth_weights(i);
            for (std::size_t j = 0; j < weights.size(); ++j) {
                const auto weight = weights[j];
                result(row, weight.first) += weight.second;
            }
        }
    }

    class SolarTransmissionBase {
      protected:
        const Geometry& m_geometry;
        const Geometry1D* m_geometry_1d = nullptr;
        const sasktran2::raytracing::RayTracerBase* m_raytracer = nullptr;
#ifdef SKTRAN_RUST_SUPPORT
        const Geometry2D* m_geometry_2d = nullptr;
        const sasktran2::raytracing::RustRayTracer2D* m_raytracer_2d = nullptr;
#endif

      public:
        SolarTransmissionBase(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : m_geometry(geometry), m_geometry_1d(&geometry),
              m_raytracer(&raytracer) {}

#ifdef SKTRAN_RUST_SUPPORT
        SolarTransmissionBase(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : m_geometry(geometry), m_geometry_2d(&geometry),
              m_raytracer_2d(&raytracer) {}
#endif
    };

    class SolarTransmissionExact : public SolarTransmissionBase {
      private:
      public:
        SolarTransmissionExact(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : SolarTransmissionBase(geometry, raytracer) {}

        virtual void initialize_config(const sasktran2::Config& config){};

        virtual void
        initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>&
                                integration_rays){};

        void generate_geometry_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            Eigen::MatrixXd& od_matrix,
            std::vector<bool>& ground_hit_flag) const;

#ifdef SKTRAN_RUST_SUPPORT
        SolarTransmissionExact(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : SolarTransmissionBase(geometry, raytracer) {}

        void generate_geometry_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            Eigen::SparseMatrix<double, Eigen::RowMajor>& od_matrix,
            std::vector<bool>& ground_hit_flag) const;
#endif
    };

    class SolarTransmissionTable : public SolarTransmissionExact {
      private:
        std::unique_ptr<sasktran2::grids::SourceLocationInterpolator>
            m_location_interpolator;
        const sasktran2::Config* m_config;
        Eigen::MatrixXd m_geometry_matrix;

        std::vector<bool> m_ground_hit_flag;

      public:
        SolarTransmissionTable(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : SolarTransmissionExact(geometry, raytracer) {}

        void initialize_config(const sasktran2::Config& config) override {
            m_config = &config;
        };

        void
        initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>&
                                integration_rays) override;

        void generate_interpolation_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator,
            std::vector<bool>& ground_hit_flag) const;

        const Eigen::MatrixXd& geometry_matrix() const {
            return m_geometry_matrix;
        }
    };

    /**
     * The PhaseHandler is responsible for constructing the phase function for
     * the single scatter source term. The single scatter source term is needed
     * on a set of (wavelength, cos_angle, geometry) points.
     *
     * Construction can be done in one of two ways, depending on the
     * configuration and user input. The default is to construct the phase
     * function through the user input Legendre coefficients.
     *
     * The other option is for the user to directly input the phase function.
     * This is useful for cases where the phase function requires many terms in
     * the Legendre series.
     *
     * Storage for the phase function is handled in a semi-complicated hard to
     * understand manner. For each thread, we store the phase function (stokes,
     * internal_index) where internal_index represents a single scattering angle
     * at a single geomtry level in the atmosphere.
     *
     * The internal_index is mapped back to the geometry index through the
     * m_internal_to_geometry, and the scattering angle is determined by the
     * m_internal_to_cos_scatter.  This combination lets us calculate the phase
     * function.  But to actually use it, we need to map the entrance and exit
     * points of the ray to the internal indices.  This is done through the
     * m_geometry_entrance_to_internal and m_geometry_exit_to_internal
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class PhaseHandler {
      private:
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;
        const sasktran2::Config* m_config;
        const sasktran2::Geometry& m_geometry;

        std::vector<std::array<double, NSTOKES>>
            m_scatter_angles; /** Full list of  scattering angles that we need.
                                 for NSTOKES =3 this is (cos_scatter, C1, C2) */

        Eigen::MatrixXd m_wigner_d00; /** Wigner D matrix for the phase function
                                         (legendre_order, scatter_angle) */
        Eigen::MatrixXd m_wigner_d02; /** Wigner D matrix for the phase function
                                         (legendre_order, scatter_angle) */

        // Internal phase functions and derivatives on the actual grid
        Eigen::Tensor<double, 3>
            m_phase; /** (stokes eq, internal_index, thread) **/
        Eigen::Tensor<double, 4>
            m_d_phase; /** (stokes eq, internal_index, deriv, thread) **/

        std::vector<std::vector<std::vector<int>>>
            m_geometry_entrance_to_internal; /** Mapping from layer entrances to
                                                internal,
                                                [los][layer][interp_index] */
        std::vector<std::vector<std::vector<int>>>
            m_geometry_exit_to_internal;         /** Mapping from layer exits to
                                                    internal, [los][layer][interp_index]
                                                  */
        std::vector<int> m_internal_to_geometry; /** Maps the internal index
                                                      to the geometry index */
        std::vector<int>
            m_internal_to_cos_scatter; /** Determines what scattering angle to
                                          use for each internal index */

      public:
        PhaseHandler(const Geometry& geometry) : m_geometry(geometry) {}

        /**
         * Initializes the phase handler with the configuration object
         */
        void initialize_config(const sasktran2::Config& config) {
            m_config = &config;
        }

        /**
         *  Initializes the phase handler with the atmosphere object
         *
         *  @param atmosphere The atmosphere object
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        /**
         * Initializes the phase handler with the geometry object
         *
         * @param los_rays The traced line of sight rays
         * @param index_map The index map
         */
        void initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
            const std::vector<std::vector<int>>& index_map);

        /**
         *   Calculates the phase function from the legendre coefficients at the
         * necessary scatter angles
         *
         *   @param threadidx The thread index
         *   @param wavelidx The wavelength index
         */
        void calculate(int wavelidx, int threadidx);

        /**
         * Calculates the phase function at a given point and puts it into
         * source
         *
         * @param wavelidx The wavelength index
         * @param losidx The line of sight index
         * @param layeridx The layer index
         * @param index_weights The interpolation weights
         * @param is_entrance  True if we are at the entrance to a layer, false
         * if we are at the exit to a layer
         * @param source The source term
         */
        void scatter(int wavelidx, int losidx, int layeridx,
                     const raytracing::GridWeightStencilView& index_weights,
                     bool is_entrance,
                     sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& source) const;

      private:
        void initialize_geometry_impl(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
            const std::vector<std::vector<int>>& index_map);

        void
        scatter_impl(int wavelidx, int losidx, int layeridx,
                     const raytracing::GridWeightStencilView& index_weights,
                     bool is_entrance,
                     sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& source) const;
    };

    template <bool USE_ACTIVE_DERIVATIVES, int NSTOKES>
    inline void scattering_source_impl(
        const PhaseHandler<NSTOKES>& phase_handler, int threadidx, int losidx,
        int layeridx, int wavelidx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, double solar_trans,
        const atmosphere::Atmosphere<NSTOKES>& atmosphere,
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
            solar_trans_iter,
        bool calculate_derivatives, std::vector<int>& active_derivative_indices,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            source) {
        const auto& storage = atmosphere.storage();
        source.value.setZero();

        if (calculate_derivatives) {
            if constexpr (USE_ACTIVE_DERIVATIVES) {
                active_derivative_indices.clear();
                for (auto it = solar_trans_iter; it; ++it) {
                    active_derivative_indices.push_back(it.index());
                }
                for (std::size_t index = 0; index < index_weights.size();
                     ++index) {
                    const auto ele = index_weights[index];
                    if (ele.second == 0.0) {
                        continue;
                    }
                    active_derivative_indices.push_back(ele.first);
                    active_derivative_indices.push_back(
                        atmosphere.ssa_deriv_start_index() + ele.first);
                    for (int derivative_group = 0;
                         derivative_group <
                         atmosphere.num_scattering_deriv_groups();
                         ++derivative_group) {
                        active_derivative_indices.push_back(
                            atmosphere.scat_deriv_start_index() +
                            derivative_group *
                                atmosphere.storage().total_extinction.rows() +
                            ele.first);
                    }
                }
                std::sort(active_derivative_indices.begin(),
                          active_derivative_indices.end());
                active_derivative_indices.erase(
                    std::unique(active_derivative_indices.begin(),
                                active_derivative_indices.end()),
                    active_derivative_indices.end());
                for (const int derivative_index : active_derivative_indices) {
                    source.deriv.col(derivative_index).setZero();
                }
            } else {
                source.deriv.setZero();
            }
        }

        double ssa = 0;
        double k = 0;
        if (index_weights.size() == 2) {
            // Structured 1D layers always have two backing nodes. Keep the
            // common ray representation while retaining the direct endpoint
            // interpolation used by the previous 1D-specific stencil.
            const auto lower = index_weights[0];
            const auto upper = index_weights[1];
            if (lower.second == 0.0) {
                ssa = storage.ssa(upper.first, wavelidx) * upper.second;
                k = storage.total_extinction(upper.first, wavelidx) *
                    upper.second;
            } else if (upper.second == 0.0) {
                ssa = storage.ssa(lower.first, wavelidx) * lower.second;
                k = storage.total_extinction(lower.first, wavelidx) *
                    lower.second;
            } else {
                ssa = storage.ssa(lower.first, wavelidx) * lower.second +
                      storage.ssa(upper.first, wavelidx) * upper.second;
                k = storage.total_extinction(lower.first, wavelidx) *
                        lower.second +
                    storage.total_extinction(upper.first, wavelidx) *
                        upper.second;
            }
        } else {
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto ele = index_weights[index];
                if (ele.second == 0.0) {
                    continue;
                }
                ssa += storage.ssa(ele.first, wavelidx) * ele.second;
                k += storage.total_extinction(ele.first, wavelidx) * ele.second;
            }
        }

        const double source_amplitude = k * ssa * solar_trans / (EIGEN_PI * 4);
        const bool use_zero_safe_derivative_path =
            calculate_derivatives && (k == 0.0 || ssa == 0.0);
        if (!use_zero_safe_derivative_path) {
            // The common path lets the phase handler work directly with the
            // scaled source, avoiding a second dense derivative pass.
            source.value(0) = source_amplitude;
            phase_handler.scatter(threadidx, losidx, layeridx, index_weights,
                                  is_entrance, source);

            if (!calculate_derivatives) {
                return;
            }
            for (auto it = solar_trans_iter; it; ++it) {
                source.deriv(Eigen::placeholders::all, it.index()) -=
                    it.value() * source.value;
            }
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto ele = index_weights[index];
                if (ele.second == 0.0) {
                    continue;
                }
                source.deriv(Eigen::placeholders::all,
                             atmosphere.ssa_deriv_start_index() + ele.first) +=
                    ele.second * source.value / ssa;
                source.deriv(Eigen::placeholders::all, ele.first) +=
                    ele.second * source.value / k;
            }
            return;
        }

        // For derivatives at k == 0 or SSA == 0, evaluate the unit-amplitude
        // phase first so the nonzero boundary derivative remains well-defined.
        source.value(0) = 1.0;

        phase_handler.scatter(threadidx, losidx, layeridx, index_weights,
                              is_entrance, source);

        const Eigen::Vector<double, NSTOKES> phase = source.value;
        source.value *= source_amplitude;

        if (!calculate_derivatives) {
            return;
        }
        if constexpr (USE_ACTIVE_DERIVATIVES) {
            for (const int derivative_index : active_derivative_indices) {
                source.deriv.col(derivative_index) *= source_amplitude;
            }
        } else {
            source.deriv *= source_amplitude;
        }
        // Solar transmission derivative factors
        for (auto it = solar_trans_iter; it; ++it) {
            source.deriv(Eigen::placeholders::all, it.index()) -=
                it.value() * source.value;
        }

        // And SSA/k derivative factors
        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto ele = index_weights[index];
            if (ele.second == 0.0) {
                continue;
            }
            source.deriv(Eigen::placeholders::all,
                         atmosphere.ssa_deriv_start_index() + ele.first) +=
                ele.second * k * solar_trans / (EIGEN_PI * 4) * phase;

            source.deriv(Eigen::placeholders::all, ele.first) +=
                ele.second * ssa * solar_trans / (EIGEN_PI * 4) * phase;
        }
    }

    template <int NSTOKES>
    inline void scattering_source(
        const PhaseHandler<NSTOKES>& phase_handler, int threadidx, int losidx,
        int layeridx, int wavelidx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, double solar_trans,
        const atmosphere::Atmosphere<NSTOKES>& atmosphere,
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
            solar_trans_iter,
        bool calculate_derivatives, bool use_active_derivatives,
        std::vector<int>& active_derivative_indices,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            source) {
        if (use_active_derivatives) {
            scattering_source_impl<true>(
                phase_handler, threadidx, losidx, layeridx, wavelidx,
                index_weights, is_entrance, solar_trans, atmosphere,
                solar_trans_iter, calculate_derivatives,
                active_derivative_indices, source);
        } else {
            scattering_source_impl<false>(
                phase_handler, threadidx, losidx, layeridx, wavelidx,
                index_weights, is_entrance, solar_trans, atmosphere,
                solar_trans_iter, calculate_derivatives,
                active_derivative_indices, source);
        }
    }

    template <typename S, int NSTOKES>
    class SingleScatterSource : public SourceTermInterface<NSTOKES> {
      private:
        S m_solar_transmission;
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;

        Eigen::MatrixXd m_geometry_matrix;
        Eigen::SparseMatrix<double, Eigen::RowMajor> m_geometry_sparse;
        std::vector<bool> m_ground_hit_flag;

        std::vector<Eigen::VectorXd> m_solar_trans;
        std::vector<std::vector<int>> m_index_map;

        PhaseHandler<NSTOKES> m_phase_handler;

        mutable std::vector<std::vector<int>> m_start_active_derivative_indices;
        mutable std::vector<std::vector<int>> m_end_active_derivative_indices;

        mutable std::vector<
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
            m_start_source_cache;
        mutable std::vector<
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
            m_end_source_cache;

        const Geometry& m_geometry;
        const Geometry1D* m_geometry_1d = nullptr;
        const Geometry2D* m_geometry_2d = nullptr;
        const sasktran2::Config* m_config;

        std::vector<bool> m_los_ground_is_hit;
        std::vector<sasktran2::raytracing::LayerGeometry> m_los_end_layers;

        void integrated_source_constant(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::LayerGeometry& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction) const;

      public:
        SingleScatterSource(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : m_solar_transmission(geometry, raytracer), m_geometry(geometry),
              m_geometry_1d(&geometry), m_phase_handler(geometry){};

#ifdef SKTRAN_RUST_SUPPORT
        template <typename T = S,
                  std::enable_if_t<std::is_same_v<T, SolarTransmissionExact>,
                                   int> = 0>
        SingleScatterSource(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : m_solar_transmission(geometry, raytracer), m_geometry(geometry),
              m_geometry_2d(&geometry), m_phase_handler(geometry){};
#endif

        void initialize_config(const sasktran2::Config& config) override;

        /** Here the single scatter source term initializes the internal solar
         * transmission object, usually this involves tracing the required rays
         * and setting up any internal matrices.
         *
         *  @param internal_viewing Information on the internal viewing
         * geometry, los_rays and flux observers
         */
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;

        /**
         *
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        /** Triggers an internal calculation of the source term.  This method is
         * called at the beginning of each 'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        void calculate(int wavelidx, int threadidx) override;

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly
         * passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
        void integrated_source(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction =
                    SourceTermInterface<NSTOKES>::IntegrationDirection::none)
            const override;

        bool supports_geometry_dimension(int dimension) const override {
            return dimension == 1 ||
                   (dimension == 2 && m_geometry_2d != nullptr);
        }

        /** Calculates the source term at the end of the ray.  Common examples
         * of this are ground scattering, ground emission, or the solar radiance
         * if looking directly at the sun.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param surface The surface object
         * @param source The returned source term
         */
        void end_of_ray_source(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const override;

        /**
         * @brief Not used for the Single Scatter source.
         *
         * @param wavelidx
         * @param losidx
         * @param wavel_threadidx
         * @param threadidx
         * @param source
         */
        virtual void start_of_ray_source(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const override{};
    };

    template <int NSTOKES>
    class OccultationSource : public SourceTermInterface<NSTOKES> {
      private:
        std::vector<bool> m_ground_is_hit;

      public:
        void initialize_config(const sasktran2::Config& config) override;

        /** Initializes any geometry information that is required for
         * calculating the source term.  This method is called after the line of
         * sight rays ar traced.
         *
         * @param internal_viewing Information on the internal viewing geometry,
         * los_rays and flux observers
         */
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;

        /**
         *
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        /** Triggers an internal calculation of the source term.  This method is
         * called at the beginning of each 'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        void calculate(int wavelidx, int threadidx) override;

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly
         * passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
        void integrated_source(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction =
                    SourceTermInterface<NSTOKES>::IntegrationDirection::none)
            const override;

        /** Calculates the source term at the end of the ray.  Common examples
         * of this are ground scattering, ground emission, or the solar radiance
         * if looking directly at the sun.
         *
         * @param wavelidx Raw index for the wavelength we are calculating
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param source The returned source term
         */
        void end_of_ray_source(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const override;

        /**
         * @brief Not used for the occultation source.
         *
         * @param wavelidx
         * @param losidx
         * @param wavel_threadidx
         * @param threadidx
         * @param source
         */
        virtual void start_of_ray_source(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const override{};

        bool has_interior_source() const override { return false; }
    };

} // namespace sasktran2::solartransmission
