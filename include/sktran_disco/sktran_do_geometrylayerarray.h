#pragma once
#include "sktran_disco/sktran_do.h"
#include <sasktran2/geometry.h>

namespace sasktran_disco {
    template <int NSTOKES, int CNSTR = -1>
    using GeometryLayerArrayROP =
        ReadOnlyProperties<BasicProperties<NSTOKES>, SolarProperties<NSTOKES>,
                           UserSpecProperties>;

    /** Class which calculates and stores the wavelength independent aspects of
     *   constructing the OpticalLayerArray
     *
     *   @tparam NSTOKES Number of stokes components
     *   @tparam CNSTR Number of constraints
     */
    template <int NSTOKES, int CNSTR = -1>
    class GeometryLayerArray : public GeometryLayerArrayROP<NSTOKES> {
      protected:
        const PersistentConfiguration<NSTOKES, CNSTR>&
            m_config;                      /** config object */
        Eigen::MatrixXd m_chapman_factors; /** Chapman factors, determines solar
                                              transmission in layers */
        Eigen::MatrixXd
            m_optical_interpolator;  /** Interpolating matrix from the
                                        sk2.Atmosphere object to layers */
        Eigen::VectorXd m_floor_h;   /** Floor heights of the layers */
        Eigen::VectorXd m_ceiling_h; /** Ceiling eights of the layers */

        bool m_no_interp; /** flag to indicate if interpolation from the
                             atmosphere grid is needed or not */

        GeometryLayerArray(
            const PersistentConfiguration<NSTOKES, CNSTR>& config)
            : GeometryLayerArrayROP<NSTOKES>(config), m_config(config){};

        /** Calculates the Chapman factors for the layers
         *
         *   @param earth_rad Earth radius
         */
        void calculate_chapman_factors(double earth_rad,
                                       const sasktran2::Geometry1D& geometry);

      public:
        GeometryLayerArray(
            const PersistentConfiguration<NSTOKES, CNSTR>& config,
            const sasktran2::Geometry1D& geometry);

        const Eigen::MatrixXd& chapman_factors() const {
            return m_chapman_factors;
        }
        const Eigen::MatrixXd& interpolating_matrix() const {
            return m_optical_interpolator;
        }

        bool no_interp() const { return m_no_interp; }

        const Eigen::VectorXd& layer_floor() const { return m_floor_h; }
        const Eigen::VectorXd& layer_ceiling() const { return m_ceiling_h; }
    };
} // namespace sasktran_disco
