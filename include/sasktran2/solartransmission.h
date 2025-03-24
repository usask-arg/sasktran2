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
     * @param geometry global geometry
     * @param result matrix to store the result in
     * @param index_weights buffer to avoid allocs
     */
    inline void assign_dense_matrix_column(
        int row, const sasktran2::raytracing::TracedRay& traced_ray,
        const Geometry& geometry, Eigen::MatrixXd& result,
        std::vector<std::pair<int, double>>& index_weights) {
        for (int i = 0; i < traced_ray.layers.size(); ++i) {
            const auto& layer = traced_ray.layers[i];

            geometry.assign_interpolation_weights(layer.entrance,
                                                  index_weights);

            for (const auto& iw : index_weights) {
                result(row, iw.first) += iw.second * layer.od_quad_start;
            }

            geometry.assign_interpolation_weights(layer.exit, index_weights);

            for (const auto& iw : index_weights) {
                result(row, iw.first) += iw.second * layer.od_quad_end;
            }
        }
    }

    class SolarTransmissionBase {
      protected:
        const Geometry1D& m_geometry;
        const sasktran2::raytracing::RayTracerBase& m_raytracer;

      public:
        SolarTransmissionBase(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : m_geometry(geometry), m_raytracer(raytracer) {}
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
        const sasktran2::Geometry1D& m_geometry;

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
        PhaseHandler(const Geometry1D& geometry) : m_geometry(geometry) {}

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
                     const std::vector<std::pair<int, double>>& index_weights,
                     bool is_entrance,
                     sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& source) const;
    };

    template <int NSTOKES>
    inline void scattering_source(
        const PhaseHandler<NSTOKES>& phase_handler, int threadidx, int losidx,
        int layeridx, int wavelidx,
        const std::vector<std::pair<int, double>>& index_weights,
        bool is_entrance, double solar_trans,
        const atmosphere::Atmosphere<NSTOKES>& atmosphere,
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
            solar_trans_iter,
        bool calculate_derivatives,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source,
        double& ssa, double& k) {
        const auto& storage = atmosphere.storage();
        source.value.setZero();

        if (calculate_derivatives) {
            source.deriv.setZero();
        }

        ssa = 0;
        k = 0;
        for (auto& ele : index_weights) {
            ssa += storage.ssa(ele.first, wavelidx) * ele.second;

            k += storage.total_extinction(ele.first, wavelidx) * ele.second;
        }

        source.value(0) = k * ssa * solar_trans / (EIGEN_PI * 4);

        phase_handler.scatter(threadidx, losidx, layeridx, index_weights,
                              is_entrance, source);

        if (!calculate_derivatives) {
            return;
        }
        // Solar transmission derivative factors
        for (auto it = solar_trans_iter; it; ++it) {
            source.deriv(Eigen::all, it.index()) -= it.value() * source.value;
        }

        // And SSA/k derivative factors
        for (auto& ele : index_weights) {
            source.deriv(Eigen::all,
                         atmosphere.ssa_deriv_start_index() + ele.first) +=
                ele.second * source.value / ssa;

            source.deriv(Eigen::all, ele.first) +=
                ele.second * source.value / k;
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

        sasktran2::Dual<double> m_precomputed_sources;

        mutable std::vector<std::vector<std::pair<int, double>>>
            m_thread_index_cache_one;
        mutable std::vector<std::vector<std::pair<int, double>>>
            m_thread_index_cache_two;

        mutable std::vector<
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
            m_start_source_cache;
        mutable std::vector<
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
            m_end_source_cache;

        const Geometry1D& m_geometry;
        const sasktran2::Config* m_config;

        const std::vector<sasktran2::raytracing::TracedRay>* m_los_rays;

        int m_num_cells;

        void integrated_source_quadrature(
            int wavelidx, int losidx, int layeridx, int threadidx,
            const sasktran2::raytracing::SphericalLayer& layer,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const;

        void integrated_source_constant(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction) const;

        void integrated_source_linear(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const;

      public:
        SingleScatterSource(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : m_solar_transmission(geometry, raytracer), m_geometry(geometry),
              m_phase_handler(geometry){};

        void initialize_config(const sasktran2::Config& config) override;

        /** Here the single scatter source term initializes the internal solar
         * transmission object, usually this involves tracing the required rays
         * and setting up any internal matrices.
         *
         *  @param los_rays The traced line of sight rays
         */
        void initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays)
            override;

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
            int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
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
      public:
        void initialize_config(const sasktran2::Config& config) override;

        /** Initializes any geometry information that is required for
         * calculating the source term.  This method is called after the line of
         * sight rays ar traced.
         *
         * @param los_rays The traced line of sight rays
         */
        void initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays)
            override;

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
            int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
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
    };

} // namespace sasktran2::solartransmission
