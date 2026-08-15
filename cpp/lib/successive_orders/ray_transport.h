#pragma once

#include "interpolation.h"
#include "transport.h"

#include <sasktran2/atmosphere/atmosphere.h>

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace sasktran2::successive_orders {

    struct RayTransportWorkspace {
        Eigen::VectorXd optical_depth;
        Eigen::VectorXd albedo;
        Eigen::VectorXd transmission_before;
        Eigen::VectorXd source_fraction;
        Eigen::VectorXd factor_cotangent;

        void resize(int maximum_layers);
        std::size_t storage_bytes() const;
    };

    /** Atmosphere-dependent values on fixed ray/source interpolation geometry.
     *
     * One scalar CSR topology is shared by all Stokes components. Values are
     * rebuilt in place when the atmospheric state changes. Directional and
     * reverse derivatives operate directly on extinction and single-scatter
     * albedo native columns without materializing a derivative matrix.
     */
    class RayTransportMap {
      public:
        RayTransportMap(const std::vector<RayInterpolation>& rays,
                        int num_source_columns,
                        const std::vector<int>& row_offsets,
                        const std::vector<int>& column_indices);

        const TransportSparsity& sparsity() const { return m_sparsity; }
        int num_rays() const { return m_sparsity.rows(); }
        int maximum_layers() const { return m_maximum_layers; }
        std::size_t storage_bytes() const;

        template <int NSTOKES>
        void assemble_values(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            int wavelength, TransportOperator& transport) const;

        template <int NSTOKES>
        void assemble_jvp(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            int wavelength, Eigen::Ref<const Eigen::VectorXd> native_tangent,
            Eigen::Ref<Eigen::VectorXd> value_tangent) const;

        template <int NSTOKES>
        void accumulate_vjp(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            int wavelength, Eigen::Ref<const Eigen::VectorXd> value_gradient,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            RayTransportWorkspace& workspace) const;

      private:
        void validate_operator(const TransportOperator& transport) const;
        void validate_wavelength(int wavelength, int num_wavelengths) const;

        const std::vector<RayInterpolation>* m_rays;
        TransportSparsity m_sparsity;
        int m_maximum_layers = 0;
    };

    extern template void RayTransportMap::assemble_values<1>(
        const sasktran2::atmosphere::Atmosphere<1>&, int,
        TransportOperator&) const;
    extern template void RayTransportMap::assemble_values<3>(
        const sasktran2::atmosphere::Atmosphere<3>&, int,
        TransportOperator&) const;
    extern template void RayTransportMap::assemble_jvp<1>(
        const sasktran2::atmosphere::Atmosphere<1>&, int,
        Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>) const;
    extern template void RayTransportMap::assemble_jvp<3>(
        const sasktran2::atmosphere::Atmosphere<3>&, int,
        Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>) const;
    extern template void RayTransportMap::accumulate_vjp<1>(
        const sasktran2::atmosphere::Atmosphere<1>&, int,
        Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>,
        RayTransportWorkspace&) const;
    extern template void RayTransportMap::accumulate_vjp<3>(
        const sasktran2::atmosphere::Atmosphere<3>&, int,
        Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>,
        RayTransportWorkspace&) const;

} // namespace sasktran2::successive_orders
