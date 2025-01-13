#include <sasktran2/derivative_mapping.h>

namespace sasktran2 {
    DerivativeMapping::DerivativeMapping(int nwavel, int ninternallocation,
                                         int nlegendre)
        : m_nwavel(nwavel), m_ninternallocation(ninternallocation),
          m_nlegendre(nlegendre), m_log_radiance_space(false) {}

    void DerivativeMapping::allocate_legendre_derivatives() {
        // If the legendre derivatives are already allocated, we don't need to
        // do anything
        if (m_native_mapping.d_legendre.has_value()) {
            return;
        }
        // Else create the value
        m_native_mapping.d_legendre = Eigen::Tensor<double, 3>(
            m_nlegendre, m_ninternallocation, m_nwavel);
        m_native_mapping.d_legendre.value().setZero();

        // Also allocate the scat factor
        m_native_mapping.scat_factor =
            Eigen::MatrixXd(m_ninternallocation, m_nwavel);
        m_native_mapping.scat_factor.value().setZero();
    }

    void DerivativeMapping::allocate_extinction_derivatives() {
        // If the od derivatives are already allocated, we don't need to do
        // anything
        if (m_native_mapping.d_extinction.has_value()) {
            return;
        }
        // Else create the value
        m_native_mapping.d_extinction =
            Eigen::MatrixXd(m_ninternallocation, m_nwavel);
        m_native_mapping.d_extinction.value().setZero();
    }

    void DerivativeMapping::allocate_ssa_derivatives() {
        // If the ssa derivatives are already allocated, we don't need to do
        // anything
        if (m_native_mapping.d_ssa.has_value()) {
            return;
        }
        // Else create the value
        m_native_mapping.d_ssa = Eigen::MatrixXd(m_ninternallocation, m_nwavel);
        m_native_mapping.d_ssa.value().setZero();
    }

    void DerivativeMapping::allocate_emission_derivatives() {
        // If the emission derivatives are already allocated, we don't need to
        // do anything
        if (m_native_mapping.d_emission.has_value()) {
            return;
        }
        // Else create the value
        m_native_mapping.d_emission =
            Eigen::MatrixXd(m_ninternallocation, m_nwavel);
        m_native_mapping.d_emission.value().setZero();
    }

    void DerivativeMapping::set_zero() {
        if (m_native_mapping.d_extinction.has_value()) {
            m_native_mapping.d_extinction.value().setZero();
        }
        if (m_native_mapping.d_ssa.has_value()) {
            m_native_mapping.d_ssa.value().setZero();
        }
        if (m_native_mapping.d_legendre.has_value()) {
            m_native_mapping.d_legendre.value().setZero();
        }
        if (m_native_mapping.scat_factor.has_value()) {
            m_native_mapping.scat_factor.value().setZero();
        }
        if (m_native_mapping.d_emission.has_value()) {
            m_native_mapping.d_emission.value().setZero();
        }
    }

    SurfaceDerivativeMapping::SurfaceDerivativeMapping(int nwavel,
                                                       int nbrdf_args)
        : m_nwavel(nwavel), m_nbrdf_args(nbrdf_args) {}

    void SurfaceDerivativeMapping::allocate_brdf_derivatives() {
        // If the brdf derivatives are already allocated, we don't need to do
        // anything
        if (m_native_surface_mapping.d_brdf.has_value()) {
            return;
        }
        // Else create the value
        m_native_surface_mapping.d_brdf =
            Eigen::MatrixXd(m_nwavel, m_nbrdf_args);
        m_native_surface_mapping.d_brdf.value().setZero();
    }

    void SurfaceDerivativeMapping::allocate_emission_derivatives() {
        // If the emission derivatives are already allocated, we don't need to
        // do anything
        if (m_native_surface_mapping.d_emission.has_value()) {
            return;
        }
        // Else create the value
        m_native_surface_mapping.d_emission = Eigen::MatrixXd(m_nwavel, 1);
        m_native_surface_mapping.d_emission.value().setZero();
    }

    void SurfaceDerivativeMapping::set_zero() {
        if (m_native_surface_mapping.d_brdf.has_value()) {
            m_native_surface_mapping.d_brdf.value().setZero();
        }
        if (m_native_surface_mapping.d_emission.has_value()) {
            m_native_surface_mapping.d_emission.value().setZero();
        }
    }
}; // namespace sasktran2
