#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/dual.h>
#include "../geometry.h"

namespace sasktran2::atmosphere {
    /** Base abstract storage container for specifying the atmospheric
     * constituent parameters on the full radiative transfer geometry grid,
     * which may be 1dimensional, 2dimensional, or 3dimensional. The full
     * required information is total extinction, scattering extinction, and
     * phase information for each geometry grid point.
     */
    class AtmosphereGridStorage {};

    template <int NSTOKES>
    class AtmosphereGridStorageFull : public AtmosphereGridStorage {
      public:
        Eigen::MatrixXd ssa;              // location, wavel
        Eigen::MatrixXd total_extinction; // location, wavel

        // Scattering parameters
        Eigen::Matrix<sasktran2::types::leg_coeff, -1, -1>
            f; // location, wavel, (Delta scaling factor)

        int applied_f_order;    // Order of the delta_m scaling
        int applied_f_location; // Index to the phase moment that defines the
                                // legendre scaling f

        int scatderivstart;
        Eigen::Tensor<double, 3>
            leg_coeff; // legendre order (polarized stacked), location, wavel
        Eigen::Tensor<double, 4>
            d_leg_coeff; // (legendre order, location, wavel, deriv)
        Eigen::Tensor<double, 3> d_f; // (location, wavel, deriv)

        int numscatderiv;

      public:
        AtmosphereGridStorageFull(int nwavel, int nlocation, int numlegendre) {
            ssa.resize(nlocation, nwavel);
            total_extinction.resize(nlocation, nwavel);
            f.resize(nlocation, nwavel);

            if constexpr (NSTOKES == 1) {
                leg_coeff.resize(numlegendre, nlocation, nwavel);
            } else {
                leg_coeff.resize(numlegendre * 4, nlocation, nwavel);
            }

            ssa.setZero();
            total_extinction.setZero();
            leg_coeff.setZero();
            f.setZero();

            applied_f_location = -1;
            applied_f_order = -1;

            numscatderiv = 0;
            scatderivstart = 0;
        }

        void resize_derivatives(int numderiv) {
            int legendre = leg_coeff.dimension(0);
            int numgeo = leg_coeff.dimension(1);
            int nwavel = leg_coeff.dimension(2);
            scatderivstart = 2 * numgeo;

            d_leg_coeff.resize(legendre, numgeo, nwavel, numderiv);
            d_f.resize(numgeo, nwavel, numderiv);

            d_leg_coeff.setZero();
            d_f.setZero();

            numscatderiv = numderiv;
        }

        int max_stored_legendre() const {
            const auto& d = leg_coeff.dimensions();

            if constexpr (NSTOKES == 1) {
                return (int)d[0];
            } else if constexpr (NSTOKES == 3) {
                return (int)(d[0] / 4);
            }
        }
    };

    /** Class which performs interpolation over the phase matrix.
     *  This can either mean direct interpolation of phase matrix values, or
     * calculation of the phase matrix from the Legendre/greek coefficients
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES, bool ssonly = false> class PhaseInterpolator {
        using ScatWeightType = typename std::conditional<
            ssonly, Eigen::Matrix<sasktran2::types::leg_coeff, NSTOKES, -1>,
            Eigen::Matrix<sasktran2::types::leg_coeff, NSTOKES * NSTOKES,
                          -1>>::type;

      private:
        ScatWeightType m_scattering_weights;
        bool m_geometry_loaded;

      public:
        PhaseInterpolator();

        void load_scattering_angle(int num_legendre,
                                   const Eigen::Vector3d& incoming_ray,
                                   const Eigen::Vector3d& outgoing_ray,
                                   bool outgoing_facing_away = true);

        template <sasktran2::dualstorage S>
        void scatter(const AtmosphereGridStorageFull<NSTOKES>& phase_storage,
                     int wavelidx,
                     const std::vector<std::pair<int, double>>& index_weights,
                     sasktran2::Dual<double, S, NSTOKES>& source) const;
    };
} // namespace sasktran2::atmosphere
