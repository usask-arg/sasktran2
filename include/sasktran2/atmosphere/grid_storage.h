#pragma once

#include <map>
#include <sasktran2/internal_common.h>
#include <sasktran2/dual.h>
#include <sasktran2/derivative_mapping.h>
#include <string>

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
      private:
        std::map<std::string, DerivativeMapping>
            m_derivative_mappings; /** Derivatives
                                       for the
                                       atmosphere */

      public:
        Eigen::MatrixXd ssa;              // location, wavel
        Eigen::MatrixXd total_extinction; // location, wavel
        Eigen::MatrixXd emission_source;  // location, wavel

        // Unscaled quantities, necessary for temporary storage during
        // derivative propagations
        Eigen::MatrixXd unscaled_ssa;              // location, wavel
        Eigen::MatrixXd unscaled_total_extinction; // location, wavel

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

        Eigen::MatrixXi
            max_order; // Maximum order of the phase function for each location

        std::vector<Eigen::MatrixXi>
            d_max_order; // Maximum order of the derivative of the phase
                         // function for each location

        int numscatderiv;

        Eigen::VectorXd solar_irradiance; // wavel
      public:
        AtmosphereGridStorageFull(int nwavel, int nlocation, int numlegendre) {
            ssa.resize(nlocation, nwavel);
            total_extinction.resize(nlocation, nwavel);
            unscaled_total_extinction.resize(nlocation, nwavel);
            unscaled_ssa.resize(nlocation, nwavel);
            emission_source.resize(nlocation, nwavel);
            f.resize(nlocation, nwavel);
            max_order.resize(nlocation, nwavel);
            solar_irradiance.resize(nwavel);

            if constexpr (NSTOKES == 1) {
                leg_coeff.resize(numlegendre, nlocation, nwavel);
            } else {
                leg_coeff.resize(numlegendre * 4, nlocation, nwavel);
            }

            ssa.setZero();
            total_extinction.setZero();
            emission_source.setZero();
            leg_coeff.setZero();
            f.setZero();
            solar_irradiance.setConstant(1.0);

            applied_f_location = -1;
            applied_f_order = -1;

            numscatderiv = 0;
            scatderivstart = 0;
        }

        /**
         *  Allocates the necessary storage for the given number of scattering
         * derivatives
         *
         * @param numderiv the number of scattering derivatives to allocate
         */
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

            d_max_order.resize(numderiv);
            for (int i = 0; i < numderiv; ++i) {
                d_max_order[i].resize(numgeo, nwavel);
                d_max_order[i].setZero();
            }
        }

        /*
         *  Determines the maximum stored legendre coefficient index.
         */
        int max_stored_legendre() const {
            const auto& d = leg_coeff.dimensions();

            if constexpr (NSTOKES == 1) {
                return (int)d[0];
            } else if constexpr (NSTOKES == 3) {
                return (int)(d[0] / 4);
            }
        }

        /**
         *  Determines the maximum order of the phase function for each
         * location.  This is done by finding the highest stored legendre
         * coefficient index that is non-zero.
         */
        void determine_maximum_order() {
            max_order.setConstant(1);

            for (int i = 0; i < max_order.rows(); ++i) {
                for (int j = 0; j < max_order.cols(); ++j) {

                    if constexpr (NSTOKES == 1) {
                        for (int k = 0; k < max_stored_legendre(); ++k) {
                            if (leg_coeff(k, i, j) != 0) {
                                max_order(i, j) = k + 1;
                            }
                        }
                    } else {
                        for (int k = 0; k < max_stored_legendre(); ++k) {
                            if (leg_coeff(4 * k, i, j) != 0) {
                                max_order(i, j) = k + 1;
                            }
                        }
                    }
                }
            }

            for (int d = 0; d < d_max_order.size(); ++d) {
                auto& max_order_deriv = d_max_order[d];

                for (int i = 0; i < max_order_deriv.rows(); ++i) {
                    for (int j = 0; j < max_order_deriv.cols(); ++j) {

                        if constexpr (NSTOKES == 1) {
                            for (int k = 0; k < max_stored_legendre(); ++k) {
                                if (d_leg_coeff(k, i, j, d) != 0) {
                                    max_order_deriv(i, j) = k + 1;
                                }
                            }
                        } else {
                            for (int k = 0; k < max_stored_legendre(); ++k) {
                                if (d_leg_coeff(4 * k, i, j, d) != 0) {
                                    max_order_deriv(i, j) = k + 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        std::map<std::string, DerivativeMapping>& derivative_mappings() {
            return m_derivative_mappings;
        }

        const std::map<std::string, DerivativeMapping>&
        derivative_mappings_const() const {
            return m_derivative_mappings;
        }

        DerivativeMapping& get_derivative_mapping(const std::string& name) {
            if (auto it = m_derivative_mappings.find(name);
                it != m_derivative_mappings.end()) {
                // Key exists; just return reference
                return it->second;
            } else {
                // Key does not exist; create it in-place without default
                // constructor
                auto [new_it, inserted] = m_derivative_mappings.emplace(
                    std::piecewise_construct, std::forward_as_tuple(name),
                    std::forward_as_tuple(ssa.cols(), ssa.rows(),
                                          leg_coeff.dimension(0)));
                return new_it->second;
            }
        }
    };

    /** Class which performs interpolation over the phase matrix.
     *  This can either mean direct interpolation of phase matrix values, or
     * calculation of the phase matrix from the Legendre/greek coefficients
     *
     * Note: This class isn't used anymore, we keep it around because it may be
     * used in future source functions that allocate less memory?
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
