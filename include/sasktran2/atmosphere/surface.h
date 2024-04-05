#pragma once

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

    } // namespace brdf
    /**
     * The full surface representation inside the SASKTRAN2 model.  Currently
     * this is just the BRDF object, but will be expanded in the future to
     * include other surface properties such as surface emissions.
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
        std::shared_ptr<brdf::BRDF<NSTOKES>> m_brdf_object;
        Eigen::MatrixXd m_brdf_args;
        std::vector<Eigen::MatrixXd> m_d_brdf_args;

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
            m_brdf_object = std::make_shared<brdf::Lambertian<NSTOKES>>();

            m_brdf_args.resize(m_brdf_object->num_args(), num_wavel);
            m_d_brdf_args.resize(m_brdf_object->num_deriv());
            for (int i = 0; i < m_brdf_object->num_deriv(); ++i) {
                m_d_brdf_args[i].resize(m_brdf_object->num_args(), num_wavel);

                m_d_brdf_args[i].setConstant(1.0);
            }
            m_brdf_args.setZero();
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
    };

} // namespace sasktran2::atmosphere
