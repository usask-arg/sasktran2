#pragma once

#include "geometry.h"
#include "scattering.h"

#include <sasktran2/atmosphere/atmosphere.h>

#include <Eigen/Core>

#include <cstddef>
#include <memory>
#include <vector>

namespace sasktran2::successive_orders {

    /** Geometry-only assembly of the scalar point scattering operator.
     *
     * Atmospheric points retain Legendre coefficients rather than dense
     * angular matrices. Ground BRDF geometry is precomputed once, while BRDF
     * values remain wavelength and atmosphere dependent. Native derivative
     * products use Atmosphere's standard scattering-group and surface slots.
     */
    class ScalarScatteringAssembler {
      public:
        ScalarScatteringAssembler(const SourceGeometry1D& geometry,
                                  int num_coefficients);

        const ScatteringBlockLayout& layout() const { return m_layout; }
        int num_coefficients() const {
            return m_angular_basis->num_coefficients();
        }
        int ground_value_size() const { return m_ground_value_offsets.back(); }

        ScatteringOperator<1> create_operator() const;

        void
        assemble_values(const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
                        int wavelength,
                        ScatteringOperator<1>& scattering) const;

        /** Maps a native atmospheric tangent into scattering parameters. */
        void
        assemble_jvp(const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
                     int wavelength,
                     Eigen::Ref<const Eigen::VectorXd> native_tangent,
                     Eigen::MatrixXd& atmospheric_coefficient_tangent,
                     Eigen::VectorXd& ground_value_tangent) const;

        /** Accumulates scattering-parameter cotangents into native slots. */
        void accumulate_vjp(
            const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
            int wavelength,
            Eigen::Ref<const Eigen::MatrixXd> atmospheric_coefficient_gradient,
            Eigen::Ref<const Eigen::VectorXd> ground_value_gradient,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            bool accumulate_scattering = true,
            bool accumulate_surface = true) const;

        std::size_t storage_bytes() const;

      private:
        struct GroundAngularGeometry {
            Eigen::VectorXd mu_in;
            Eigen::VectorXd weighted_mu_in;
            Eigen::VectorXd mu_out;
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::RowMajor>
                phi_difference;

            std::size_t storage_bytes() const {
                return static_cast<std::size_t>(
                           mu_in.size() + weighted_mu_in.size() +
                           mu_out.size() + phi_difference.size()) *
                       sizeof(double);
            }
        };

        void validate_atmosphere(
            const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
            int wavelength) const;
        void validate_operator(const ScatteringOperator<1>& scattering) const;
        void validate_native_product(
            const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
            Eigen::Index native_size) const;
        void
        assemble_ground(const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
                        int wavelength, Eigen::VectorXd& values) const;

        const SourceGeometry1D* m_geometry;
        ScatteringBlockLayout m_layout;
        std::shared_ptr<const ScalarAngularBasis> m_angular_basis;
        std::vector<int> m_ground_value_offsets;
        std::vector<GroundAngularGeometry> m_ground_geometry;
    };

    /** Geometry-only assembly of coefficient-space I/Q/U scattering. */
    class VectorScatteringAssembler {
      public:
        VectorScatteringAssembler(const SourceGeometry1D& geometry,
                                  int num_coefficients);

        const ScatteringBlockLayout& layout() const { return m_layout; }
        int num_coefficients() const {
            return m_angular_basis->num_coefficients();
        }
        int coefficient_value_size() const {
            return m_layout.atmospheric_blocks() * 4 * num_coefficients();
        }
        int ground_value_size() const { return m_ground_value_offsets.back(); }

        ScatteringOperator<3> create_operator() const;

        void
        assemble_values(const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
                        int wavelength,
                        ScatteringOperator<3>& scattering) const;

        /** Maps a native atmospheric tangent into dense scattering values. */
        void
        assemble_jvp(const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
                     int wavelength,
                     Eigen::Ref<const Eigen::VectorXd> native_tangent,
                     Eigen::MatrixXd& atmospheric_coefficient_tangent,
                     Eigen::VectorXd& ground_value_tangent) const;

        /** Accumulates dense scattering cotangents into native slots. */
        void accumulate_vjp(
            const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
            int wavelength,
            Eigen::Ref<const Eigen::MatrixXd> atmospheric_coefficient_gradient,
            Eigen::Ref<const Eigen::VectorXd> ground_value_gradient,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            bool accumulate_scattering = true,
            bool accumulate_surface = true) const;

        std::size_t storage_bytes() const;

      private:
        using RowMajorMatrix = ScatteringOperator<3>::RowMajorMatrix;

        struct GroundAngularGeometry {
            Eigen::VectorXd mu_in;
            Eigen::VectorXd weighted_mu_in;
            Eigen::VectorXd mu_out;
            RowMajorMatrix phi_difference;

            std::size_t storage_bytes() const {
                return static_cast<std::size_t>(
                           mu_in.size() + weighted_mu_in.size() +
                           mu_out.size() + phi_difference.size()) *
                       sizeof(double);
            }
        };

        void validate_atmosphere(
            const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
            int wavelength) const;
        void validate_operator(const ScatteringOperator<3>& scattering) const;
        void validate_native_product(
            const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
            Eigen::Index native_size) const;
        void
        assemble_ground(const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
                        int wavelength, Eigen::VectorXd& values) const;

        const SourceGeometry1D* m_geometry;
        ScatteringBlockLayout m_layout;
        std::shared_ptr<const VectorAngularBasis> m_angular_basis;
        std::vector<int> m_ground_value_offsets;
        std::vector<GroundAngularGeometry> m_ground_geometry;
    };

} // namespace sasktran2::successive_orders
