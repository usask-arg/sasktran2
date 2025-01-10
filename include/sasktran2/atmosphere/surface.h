#pragma once

#include "sasktran2/derivative_mapping.h"
#include "sasktran2/grids.h"
#include <limits>
#include <sasktran2/internal_common.h>
#include <stdexcept>

namespace sasktran2::atmosphere {
    // Forward declaration
    template <int NSTOKES> struct Surface;

    namespace brdf {

        /**
         * Base class that defines the interface for a Biderctional Reflectance
         * Distribution Function (BRDF)
         *
         * @tparam NSTOKES
         */
        template <int NSTOKES> struct BRDF {
            virtual ~BRDF() = default;

            /**
             * Interface for the BRDF function. Each BRDF function can take in a
             * variable amount of args.  The number of args is defined by the
             * num_args() function.
             *
             * @param mu_in Cosine of the incoming zenith angle
             * @param mu_out Cosine of the outgoing zenith angle
             * @param phi_diff Absolute difference in azimuth angle in
             * radiances, 0 corresponds to the "forward scattering direction"
             * @param args
             * @return Eigen::Matrix<double, NSTOKES, NSTOKES> BRDF matrix
             * normalized to integrate to 1
             */
            virtual Eigen::Matrix<double, NSTOKES, NSTOKES>
            brdf(double mu_in, double mu_out, double phi_diff,
                 Eigen::Ref<const Eigen::VectorXd> args) const = 0;

            // Some BRDFs may not have derivatives, in which case this is never
            // called
            /**
             * The derivative of the brdf function for a given deriv index
             *
             * @param deriv_index
             * @param mu_in
             * @param mu_out
             * @param phi_diff
             * @param args
             * @param d_args
             * @return Eigen::Matrix<double, NSTOKES, NSTOKES>
             */
            virtual Eigen::Matrix<double, NSTOKES, NSTOKES>
            d_brdf(int deriv_index, double mu_in, double mu_out,
                   double phi_diff, Eigen::Ref<const Eigen::VectorXd> args,
                   Eigen::Ref<const Eigen::VectorXd> d_args) const {
                spdlog::critical("Derivative not implemented for this BRDF");
                throw std::runtime_error(
                    "Derivative not implemented for this BRDF");
            };

            /**
             * @brief Number of derivatives that the BRDF can return
             *
             * @return int
             */
            virtual int num_deriv() const = 0;

            /**
             * @brief Number of arguments that the BRDF function takes in
             *
             * @return int
             */
            virtual int num_args() const = 0;

            /**
             * In the discrete ordinates model we expand the BRDF in terms of a
             * fourier series in azimuth, this is the maximumm order that the
             * BRDF can be expanded to.
             *
             * @return int
             */
            virtual int max_azimuthal_order() const {
                return std::numeric_limits<int>::max();
            }
        };

        template <int NSTOKES> struct WeightedBRDF : public BRDF<NSTOKES> {
          private:
          public:
        };

        /**
         * A lambertian BRDF function, the brdf is a constant value of albedo /
         * pi
         *
         * args(0) = albedo
         *
         * @tparam NSTOKES
         */
        template <int NSTOKES> struct Lambertian : BRDF<NSTOKES> {
            Eigen::Matrix<double, NSTOKES, NSTOKES>
            brdf(double mu_in, double mu_out, double phi_diff,
                 Eigen::Ref<const Eigen::VectorXd> args) const override {
                Eigen::Matrix<double, NSTOKES, NSTOKES> res;
                res.setZero();
                res(0, 0) = args(0) / EIGEN_PI;

                return res;
            }

            Eigen::Matrix<double, NSTOKES, NSTOKES>
            d_brdf(int deriv_index, double mu_in, double mu_out,
                   double phi_diff, Eigen::Ref<const Eigen::VectorXd> args,
                   Eigen::Ref<const Eigen::VectorXd> d_args) const override {
                Eigen::Matrix<double, NSTOKES, NSTOKES> res;
                res.setZero();
                res(0, 0) = d_args(0) / EIGEN_PI;

                return res;
            }

            int num_deriv() const override { return 1; }

            int num_args() const override { return 1; }

            int max_azimuthal_order() const override { return 1; }
        };

        /**
         * A BRDF reperesenting a snow surface using the Kokhanovsky model.
         *
         * MODIS BRDF/Albedo Product ATBD
         * https://modis.gsfc.nasa.gov/data/atbd/atbd_mod09.pdf
         *
         * args(0) = (chi + M) / wavelennm * L
         *
         * @tparam NSTOKES
         */
        template <int NSTOKES> struct SnowKokhanovsky : BRDF<NSTOKES> {
          private:
            double K0(double mu) const {
                return (3.0 / 7.0) * (1.0 + 2.0 * mu);
            }

            double R0(double mus, double muv, double theta) const {
                const double a = 1.247;
                const double b = 1.186;
                const double c = 5.157;

                return (a + b * (mus + muv) + c * mus * muv + p(theta)) /
                       (4.0 * (mus + muv));
            }

            double p(double thetadegrees) const {
                return (11.1 * exp(-0.087 * thetadegrees) +
                        1.1 * exp(-0.014 * thetadegrees));
            }

          public:
            Eigen::Matrix<double, NSTOKES, NSTOKES>
            brdf(double mu_in, double mu_out, double phi_diff,
                 Eigen::Ref<const Eigen::VectorXd> args) const override {
                Eigen::Matrix<double, NSTOKES, NSTOKES> res;
                res.setZero();

                // Input arguments
                // args(0) = (chi + M) / wavelennm * L

                double mus = mu_in;
                double muv = mu_out;
                double ss = sqrt(1 - mus * mus);
                double sv = sqrt(1 - muv * muv);
                double cost = std::max(
                    -1.0, std::min(1.0, -mus * muv + ss * sv * cos(phi_diff)));
                double theta = acos(cost) * 180.0 / EIGEN_PI;
                double alpha;
                double r0;

                alpha = sqrt(4 * EIGEN_PI * args(0));

                r0 = R0(mus, muv, theta);
                res(0, 0) =
                    r0 * (exp(-alpha * K0(mus) * K0(muv) / r0)) / EIGEN_PI;

                return res;
            }

            Eigen::Matrix<double, NSTOKES, NSTOKES>
            d_brdf(int deriv_index, double mu_in, double mu_out,
                   double phi_diff, Eigen::Ref<const Eigen::VectorXd> args,
                   Eigen::Ref<const Eigen::VectorXd> d_args) const override {
                Eigen::Matrix<double, NSTOKES, NSTOKES> res;
                res.setZero();

                // Input arguments
                // args(0) = (chi + M) / wavelennm * L

                double mus = mu_in;
                double muv = mu_out;
                double ss = sqrt(1 - mus * mus);
                double sv = sqrt(1 - muv * muv);
                double cost = std::max(
                    -1.0, std::min(1.0, -mus * muv + ss * sv * cos(phi_diff)));
                double theta = acos(cost) * 180.0 / EIGEN_PI;
                double alpha;
                double r0;

                alpha = sqrt(4 * EIGEN_PI * args(0));

                // d_alpha = 0.5 / alpha * d_args(0) * 4 * EIGEN_PI;

                r0 = R0(mus, muv, theta);
                double factor = -K0(mus) * K0(muv) / r0;
                res(0, 0) = r0 * (exp(alpha * factor)) *
                            (2 / alpha * d_args(0) * factor);

                return res;
            }

            int num_deriv() const override { return 1; }

            int num_args() const override { return 1; }
        };

        /**
         * A BRDF corresponding to the kernel-based model used by MODIS.
         *
         * args(0) = isotropic (Lambertian) component
         * args(1) = volumetric (Ross-thick) component
         * args(2) = geometric (Li-Sparse-R) component
         *
         * @tparam NSTOKES
         */
        template <int NSTOKES> struct MODIS : BRDF<NSTOKES> {
          public:
            Eigen::Matrix<double, NSTOKES, NSTOKES>
            brdf(double mu_in, double mu_out, double phi_diff,
                 Eigen::Ref<const Eigen::VectorXd> args) const override {
                Eigen::Matrix<double, NSTOKES, NSTOKES> res;
                res.setZero();

                double csza = mu_in;
                double cvza = mu_out;
                double ssza = sqrt(1 - csza * csza);
                double svza = sqrt(1 - cvza * cvza);
                double tsza = ssza / csza;
                double tvza = svza / cvza;
                double craa =
                    -cos(phi_diff); // negate b/c input defines raa = 0 as
                                    // forward, but the formula below defines
                                    // raa = 0 as backward
                double sraa = sin(phi_diff);
                double csa = std::max(
                    -1.0, std::min(1.0, csza * cvza + ssza * svza * craa));
                double sa = acos(csa);
                double ssa = sin(sa);

                // volumetric kernel (Ross-Thick, Roujean et al. 1992, Eqn 38)
                double k_vol =
                    ((0.5 * EIGEN_PI - sa) * csa + ssa) / (csza + cvza) -
                    0.25 * EIGEN_PI;

                // geometric kernel (Li-Sparse-R, Wanner et al. 1995, Eqns
                // 39-44)
                double d2 = tsza * tsza + tvza * tvza - 2 * tsza * tvza * craa;
                double ct = std::max(
                    -1.0, std::min(1.0, 2 *
                                            sqrt(d2 + tsza * tsza * tvza *
                                                          tvza * sraa * sraa) *
                                            csza * cvza / (csza + cvza)));
                double t = acos(ct);
                double st = sin(t);
                double o =
                    (t - st * ct) * (csza + cvza) / (EIGEN_PI * csza * cvza);
                double k_geo =
                    o - (csza + cvza - 0.5 * (1 + csa)) / (csza * cvza);

                res(0, 0) =
                    (args(0) + args(1) * k_vol + args(2) * k_geo) / EIGEN_PI;
                return res;
            }

            Eigen::Matrix<double, NSTOKES, NSTOKES>
            d_brdf(int deriv_index, double mu_in, double mu_out,
                   double phi_diff, Eigen::Ref<const Eigen::VectorXd> args,
                   Eigen::Ref<const Eigen::VectorXd> d_args) const override {
                Eigen::Matrix<double, NSTOKES, NSTOKES> res;
                res.setZero();

                if (deriv_index == 0) {
                    res(0, 0) = 1.0 / EIGEN_PI;
                    return res;
                }

                double csza = mu_in;
                double cvza = mu_out;
                double ssza = sqrt(1 - csza * csza);
                double svza = sqrt(1 - cvza * cvza);
                double tsza = ssza / csza;
                double tvza = svza / cvza;
                double craa =
                    -cos(phi_diff); // negate b/c input defines raa = 0 as
                                    // forward, but the formula below defines
                                    // raa = 0 as backward
                double sraa = sin(phi_diff);
                double csa = std::max(
                    -1.0, std::min(1.0, csza * cvza + ssza * svza * craa));
                double sa = acos(csa);
                double ssa = sin(sa);

                // volumetric kernel (Ross-Thick, Roujean et al. 1992, Eqn 38)
                double k_vol =
                    ((0.5 * EIGEN_PI - sa) * csa + ssa) / (csza + cvza) -
                    0.25 * EIGEN_PI;

                // geometric kernel (Li-Sparse-R, Wanner et al. 1995, Eqns
                // 39-44)
                double d2 = tsza * tsza + tvza * tvza - 2 * tsza * tvza * craa;
                double ct = std::max(
                    -1.0, std::min(1.0, 2 *
                                            sqrt(d2 + tsza * tsza * tvza *
                                                          tvza * sraa * sraa) *
                                            csza * cvza / (csza + cvza)));
                double t = acos(ct);
                double st = sin(t);
                double o =
                    (t - st * ct) * (csza + cvza) / (EIGEN_PI * csza * cvza);
                double k_geo =
                    o - (csza + cvza - 0.5 * (1 + csa)) / (csza * cvza);

                if (deriv_index == 1) {
                    res(0, 0) = k_vol / EIGEN_PI;
                } else if (deriv_index == 2) {
                    res(0, 0) = k_geo / EIGEN_PI;
                }

                return res;
            }

            int num_deriv() const override { return 3; }

            int num_args() const override { return 3; }
        };

    } // namespace brdf
    /**
     * The full surface representation inside the SASKTRAN2 model.  Currently
     * this includes the BRDF object and surface emissions.
     *
     * The BRDF is handled by the user setting a BRDF object, which defaults to
     * Lambertian, and then setting the brdf_args and d_brdf_args parameters.
     * The brdf_args is a matrix of size num_args x num_wavel, and the
     * d_brdf_args is a vector of matrices of size num_args x num_wavel, where
     * the vector size is the number of derivatives that the BRDF can return.
     * These arguments are passed to the BRDF object for every wavelenght in the
     * calculation.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> struct Surface {
      private:
        int m_num_wavel;
        std::shared_ptr<brdf::BRDF<NSTOKES>> m_brdf_object;
        Eigen::MatrixXd m_brdf_args;
        std::vector<Eigen::MatrixXd> m_d_brdf_args;
        Eigen::VectorXd m_emission;

        std::map<std::string, SurfaceDerivativeMapping>
            m_derivative_mappings; /** Derivatives
                                for the
                                atmosphere */

        void allocate(int num_wavel) {
            m_brdf_args.resize(m_brdf_object->num_args(), num_wavel);
            m_d_brdf_args.resize(m_brdf_object->num_deriv());
            for (int i = 0; i < m_brdf_object->num_deriv(); ++i) {
                m_d_brdf_args[i].resize(m_brdf_object->num_args(), num_wavel);

                m_d_brdf_args[i].setConstant(0.0);
                m_d_brdf_args[i](i, Eigen::all).setConstant(1.0);
            }
            m_brdf_args.setZero();
            m_emission.resize(num_wavel);
            m_emission.setZero();
        }

      public:
        /**
         * @brief Calls the BRDF function for a given wavelength index
         *
         * @param wavel_idx
         * @param mu_in
         * @param mu_out
         * @param phi_diff
         * @return Eigen::Matrix<double, NSTOKES, NSTOKES>
         */
        Eigen::Matrix<double, NSTOKES, NSTOKES> brdf(int wavel_idx,
                                                     double mu_in,
                                                     double mu_out,
                                                     double phi_diff) const {
            return m_brdf_object->brdf(mu_in, mu_out, phi_diff,
                                       m_brdf_args.col(wavel_idx));
        }

        /**
         * @brief The derivative of the BRDF function for a given derivative and
         * wavelength index
         *
         * @param wavel_idx
         * @param mu_in
         * @param mu_out
         * @param phi_diff
         * @param deriv_index
         * @return Eigen::Matrix<double, NSTOKES, NSTOKES>
         */
        Eigen::Matrix<double, NSTOKES, NSTOKES>
        d_brdf(int wavel_idx, double mu_in, double mu_out, double phi_diff,
               int deriv_index) const {
            return m_brdf_object->d_brdf(
                deriv_index, mu_in, mu_out, phi_diff,
                m_brdf_args.col(wavel_idx),
                m_d_brdf_args[deriv_index].col(wavel_idx));
        }

        /**
         * @brief The number of derivatives
         *
         * @return int
         */
        int num_deriv() const { return (int)m_brdf_object->num_deriv(); }

        /**
         * @brief Maximum azimuthal order in the BRDF expansion
         *
         * @return int
         */
        int max_azimuthal_order() const {
            return m_brdf_object->max_azimuthal_order();
        }

        /**
         * @brief Construct a new Surface object
         *
         * Upon construction the BRDF object is set to a Lambertian BRDF
         *
         * @param num_wavel
         */
        Surface(int num_wavel) {
            m_num_wavel = num_wavel;
            m_brdf_object = std::make_shared<brdf::Lambertian<NSTOKES>>();
            allocate(num_wavel);
        }

        /**
         * @brief The internal BRDF object
         *
         * @return const brdf::BRDF<NSTOKES>&
         */
        const brdf::BRDF<NSTOKES>& brdf_object() const {
            return *m_brdf_object;
        }

        /**
         * @brief Set the brdf object object
         *
         * @param brdf
         */
        void set_brdf_object(std::shared_ptr<brdf::BRDF<NSTOKES>> brdf) {
            m_brdf_object = std::move(brdf);
            allocate(m_num_wavel);
        }

        /**
         * @brief The arguments for the BRDF function
         *
         * @return Eigen::MatrixXd& shape (num_args, num_wavel)
         */
        Eigen::MatrixXd& brdf_args() { return m_brdf_args; }

        /**
         * @brief The arguments for the BRDF function
         *
         * @return const Eigen::MatrixXd& shape (num_args, num_wavel)
         */
        const Eigen::MatrixXd& brdf_args() const { return m_brdf_args; }

        /**
         * @brief The derivative arguments for the BRDF function
         *
         * @return std::vector<Eigen::MatrixXd>& shape [num_deriv](num_args,
         * num_wavel)
         */
        std::vector<Eigen::MatrixXd>& d_brdf_args() { return m_d_brdf_args; }

        SurfaceDerivativeMapping&
        get_derivative_mapping(const std::string& name) {
            if (auto it = m_derivative_mappings.find(name);
                it != m_derivative_mappings.end()) {
                // Key exists; just return reference
                return it->second;
            } else {
                // Key does not exist; create it in-place without default
                // constructor
                auto [new_it, inserted] = m_derivative_mappings.emplace(
                    std::piecewise_construct, std::forward_as_tuple(name),
                    std::forward_as_tuple(m_num_wavel,
                                          m_brdf_object->num_args()));
                return new_it->second;
            }
        }

        const std::map<std::string, SurfaceDerivativeMapping>&
        derivative_mappings() const {
            return m_derivative_mappings;
        }

        /**
         * @brief The surface emission.
         *
         * @return Eigen::VectorXd& shape (num_wavel)
         */
        Eigen::VectorXd& emission() { return m_emission; }

        /**
         * @brief The surface emission.
         *
         * @return const Eigen::VectorXd& shape (num_wavel)
         */
        const Eigen::VectorXd& emission() const { return m_emission; }
    };

} // namespace sasktran2::atmosphere
