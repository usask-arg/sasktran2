#include "ray_transport.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace sasktran2::successive_orders {
    namespace {
        template <typename Weights, typename Values>
        double interpolate(const Weights& weights, const Values& values,
                           int wavelength) {
            double result = 0.0;
            for (const auto& weight : weights) {
                result += weight.weight * values(weight.index, wavelength);
            }
            return result;
        }

        template <typename Weights>
        double interpolate_tangent(const Weights& weights, int tangent_offset,
                                   Eigen::Ref<const Eigen::VectorXd> tangent) {
            double result = 0.0;
            for (const auto& weight : weights) {
                result +=
                    weight.weight * tangent(tangent_offset + weight.index);
            }
            return result;
        }

        double optical_depth(const RayInterpolation& ray, int layer_index,
                             Eigen::Ref<const Eigen::VectorXd> extinction) {
            double result = 0.0;
            const auto weights = ray.optical_depth_for_layer(layer_index);
            for (std::size_t index = 0; index < weights.size(); ++index) {
                const auto [atmosphere_index, weight] = weights[index];
                result += weight * extinction(atmosphere_index);
            }
            return result;
        }

        double
        optical_depth_tangent(const RayInterpolation& ray, int layer_index,
                              Eigen::Ref<const Eigen::VectorXd> tangent) {
            double result = 0.0;
            const auto weights = ray.optical_depth_for_layer(layer_index);
            for (std::size_t index = 0; index < weights.size(); ++index) {
                const auto [atmosphere_index, weight] = weights[index];
                result += weight * tangent(atmosphere_index);
            }
            return result;
        }

        template <int NSTOKES>
        double ground_transport_albedo(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            int wavelength, const RayInterpolation& ray) {
            if (!atmosphere.surface().has_spatial_lambertian_albedo()) {
                return 1.0;
            }
            return atmosphere.surface().spatial_lambertian_albedo(
                wavelength, ray.ground_horizontal());
        }

        template <int NSTOKES>
        double ground_transport_albedo_tangent(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            const RayInterpolation& ray,
            Eigen::Ref<const Eigen::VectorXd> native_tangent) {
            if (!atmosphere.surface().has_spatial_lambertian_albedo()) {
                return 0.0;
            }
            double result = 0.0;
            for (const auto& [horizontal_index, weight] :
                 ray.ground_horizontal()) {
                result += weight * native_tangent(
                                       atmosphere.surface_deriv_start_index() +
                                       horizontal_index);
            }
            return result;
        }

        template <int NSTOKES>
        void accumulate_ground_transport_albedo_vjp(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            const RayInterpolation& ray, double albedo_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) {
            if (!atmosphere.surface().has_spatial_lambertian_albedo() ||
                albedo_cotangent == 0.0) {
                return;
            }
            for (const auto& [horizontal_index, weight] :
                 ray.ground_horizontal()) {
                native_gradient(atmosphere.surface_deriv_start_index() +
                                horizontal_index) += weight * albedo_cotangent;
            }
        }
    } // namespace

    void RayTransportWorkspace::resize(int maximum_layers) {
        optical_depth.resize(maximum_layers);
        albedo.resize(maximum_layers);
        transmission_before.resize(maximum_layers);
        source_fraction.resize(maximum_layers);
        factor_cotangent.resize(maximum_layers);
    }

    std::size_t RayTransportWorkspace::storage_bytes() const {
        return static_cast<std::size_t>(optical_depth.size() + albedo.size() +
                                        transmission_before.size() +
                                        source_fraction.size() +
                                        factor_cotangent.size()) *
               sizeof(double);
    }

    RayTransportMap::RayTransportMap(const std::vector<RayInterpolation>& rays,
                                     int num_source_columns,
                                     const std::vector<int>& row_offsets,
                                     const std::vector<int>& column_indices)
        : m_rays(&rays),
          m_sparsity(num_source_columns, row_offsets, column_indices) {
        if (m_sparsity.rows() != static_cast<int>(rays.size())) {
            throw std::invalid_argument(
                "successive-orders transport row count does not match rays");
        }
        for (std::size_t row = 0; row < rays.size(); ++row) {
            const auto& ray = rays[row];
            const bool owns_optical_depth =
                !ray.optical_depth_indices.empty() ||
                !ray.optical_depth_weights.empty();
            if ((ray.optical_depth_indices.size() !=
                 ray.optical_depth_weights.size()) ||
                (!ray.layers.empty() && !owns_optical_depth &&
                 ray.traced_ray == nullptr) ||
                (ray.traced_ray != nullptr &&
                 ray.layers.size() != ray.traced_ray->layers.size())) {
                throw std::invalid_argument(
                    "invalid successive-orders compiled ray interpolation");
            }
            const int row_start = row_offsets[row];
            const int row_nnz = row_offsets[row + 1] - row_start;
            if (ray.transport_value_offset !=
                    static_cast<std::size_t>(row_start) ||
                ray.transport_row_nnz != static_cast<std::uint32_t>(row_nnz)) {
                throw std::invalid_argument(
                    "successive-orders ray CSR metadata is inconsistent");
            }
            const auto validate_weights = [&](const auto& weights) {
                for (const auto& weight : weights) {
                    if (weight.row_inner_index >=
                            static_cast<std::uint32_t>(row_nnz) ||
                        column_indices[row_start + weight.row_inner_index] !=
                            weight.source_index) {
                        throw std::invalid_argument(
                            "successive-orders source interpolation does not "
                            "match its CSR row");
                    }
                }
            };
            validate_weights(ray.source_weights);
            validate_weights(ray.ground_weights);
            m_maximum_layers =
                std::max(m_maximum_layers, static_cast<int>(ray.layers.size()));
        }
    }

    std::size_t RayTransportMap::storage_bytes() const {
        return m_sparsity.storage_bytes();
    }

    void RayTransportMap::validate_operator(
        const TransportOperator& transport) const {
        if (&transport.sparsity() != &m_sparsity ||
            transport.values().size() != m_sparsity.nonzeros()) {
            throw std::invalid_argument(
                "successive-orders transport operator has the wrong topology");
        }
    }

    void RayTransportMap::validate_wavelength(int wavelength,
                                              int num_wavelengths) const {
        if (wavelength < 0 || wavelength >= num_wavelengths) {
            throw std::out_of_range(
                "successive-orders wavelength index is out of range");
        }
    }

    template <int NSTOKES>
    void RayTransportMap::assemble_values(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
        int wavelength, TransportOperator& transport) const {
        validate_operator(transport);
        validate_wavelength(wavelength, atmosphere.num_wavel());
        auto& values = transport.values();
        values.setZero();
        const auto extinction =
            atmosphere.storage().total_extinction.col(wavelength);
        const auto& ssa = atmosphere.storage().ssa;
        const auto& offsets = m_sparsity.row_offsets();

        for (int row = 0; row < num_rays(); ++row) {
            const auto& interpolation = (*m_rays)[row];
            double transmission_before = 1.0;
            for (int layer_index =
                     static_cast<int>(interpolation.layers.size()) - 1;
                 layer_index >= 0; --layer_index) {
                const auto atmosphere_weights =
                    interpolation.atmosphere_for_layer(layer_index);
                const auto source_weights =
                    interpolation.source_for_layer(layer_index);
                const double layer_optical_depth =
                    optical_depth(interpolation, layer_index, extinction);
                const double layer_transmission =
                    std::exp(-layer_optical_depth);
                const double albedo =
                    interpolate(atmosphere_weights, ssa, wavelength);
                const double source_fraction =
                    -std::expm1(-layer_optical_depth);
                const double factor =
                    transmission_before * albedo * source_fraction;
                for (const auto& source : source_weights) {
                    values(offsets[row] + source.row_inner_index) +=
                        source.weight * factor;
                }
                transmission_before *= layer_transmission;
            }

            if (interpolation.ground_is_hit()) {
                const double ground_albedo = ground_transport_albedo(
                    atmosphere, wavelength, interpolation);
                for (const auto& source : interpolation.ground()) {
                    values(offsets[row] + source.row_inner_index) +=
                        source.weight * transmission_before * ground_albedo;
                }
            }
        }
    }

    template <int NSTOKES>
    void RayTransportMap::assemble_jvp(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
        int wavelength, Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::Ref<Eigen::VectorXd> value_tangent) const {
        validate_wavelength(wavelength, atmosphere.num_wavel());
        if (native_tangent.size() != atmosphere.num_deriv() ||
            value_tangent.size() != m_sparsity.nonzeros()) {
            throw std::invalid_argument(
                "invalid successive-orders ray transport JVP dimensions");
        }
        value_tangent.setZero();
        const auto extinction =
            atmosphere.storage().total_extinction.col(wavelength);
        const auto& ssa = atmosphere.storage().ssa;
        const auto& offsets = m_sparsity.row_offsets();

        for (int row = 0; row < num_rays(); ++row) {
            const auto& interpolation = (*m_rays)[row];
            double transmission_before = 1.0;
            double cumulative_tangent = 0.0;
            for (int layer_index =
                     static_cast<int>(interpolation.layers.size()) - 1;
                 layer_index >= 0; --layer_index) {
                const auto atmosphere_weights =
                    interpolation.atmosphere_for_layer(layer_index);
                const auto source_weights =
                    interpolation.source_for_layer(layer_index);
                const double layer_optical_depth =
                    optical_depth(interpolation, layer_index, extinction);
                const double layer_tangent = optical_depth_tangent(
                    interpolation, layer_index, native_tangent);
                const double layer_transmission =
                    std::exp(-layer_optical_depth);
                const double albedo =
                    interpolate(atmosphere_weights, ssa, wavelength);
                const double albedo_tangent = interpolate_tangent(
                    atmosphere_weights, atmosphere.ssa_deriv_start_index(),
                    native_tangent);
                const double source_fraction =
                    -std::expm1(-layer_optical_depth);
                const double factor_tangent =
                    transmission_before *
                    (albedo_tangent * source_fraction +
                     albedo * layer_transmission * layer_tangent -
                     albedo * source_fraction * cumulative_tangent);
                for (const auto& source : source_weights) {
                    value_tangent(offsets[row] + source.row_inner_index) +=
                        source.weight * factor_tangent;
                }
                transmission_before *= layer_transmission;
                cumulative_tangent += layer_tangent;
            }

            if (interpolation.ground_is_hit()) {
                const double ground_albedo = ground_transport_albedo(
                    atmosphere, wavelength, interpolation);
                const double ground_albedo_tangent =
                    ground_transport_albedo_tangent(atmosphere, interpolation,
                                                    native_tangent);
                const double factor_tangent =
                    transmission_before * (ground_albedo_tangent -
                                           ground_albedo * cumulative_tangent);
                for (const auto& source : interpolation.ground()) {
                    value_tangent(offsets[row] + source.row_inner_index) +=
                        source.weight * factor_tangent;
                }
            }
        }
    }

    template <int NSTOKES>
    void RayTransportMap::accumulate_vjp(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
        int wavelength, Eigen::Ref<const Eigen::VectorXd> value_gradient,
        Eigen::Ref<Eigen::VectorXd> native_gradient,
        RayTransportWorkspace& workspace) const {
        validate_wavelength(wavelength, atmosphere.num_wavel());
        if (value_gradient.size() != m_sparsity.nonzeros() ||
            native_gradient.size() != atmosphere.num_deriv()) {
            throw std::invalid_argument(
                "invalid successive-orders ray transport VJP dimensions");
        }
        if (workspace.optical_depth.size() < m_maximum_layers) {
            workspace.resize(m_maximum_layers);
        }
        const auto extinction =
            atmosphere.storage().total_extinction.col(wavelength);
        const auto& ssa = atmosphere.storage().ssa;
        const auto& offsets = m_sparsity.row_offsets();

        for (int row = 0; row < num_rays(); ++row) {
            const auto& interpolation = (*m_rays)[row];
            const int num_layers =
                static_cast<int>(interpolation.layers.size());
            double transmission_before = 1.0;
            for (int layer_index = num_layers - 1; layer_index >= 0;
                 --layer_index) {
                const auto atmosphere_weights =
                    interpolation.atmosphere_for_layer(layer_index);
                const auto source_weights =
                    interpolation.source_for_layer(layer_index);
                workspace.optical_depth(layer_index) =
                    optical_depth(interpolation, layer_index, extinction);
                workspace.albedo(layer_index) =
                    interpolate(atmosphere_weights, ssa, wavelength);
                workspace.transmission_before(layer_index) =
                    transmission_before;
                const double layer_transmission =
                    std::exp(-workspace.optical_depth(layer_index));
                workspace.source_fraction(layer_index) =
                    -std::expm1(-workspace.optical_depth(layer_index));
                double factor_cotangent = 0.0;
                for (const auto& source : source_weights) {
                    factor_cotangent +=
                        source.weight *
                        value_gradient(offsets[row] + source.row_inner_index);
                }
                workspace.factor_cotangent(layer_index) = factor_cotangent;
                transmission_before *= layer_transmission;
            }

            double cumulative_cotangent = 0.0;
            if (interpolation.ground_is_hit()) {
                const double ground_albedo = ground_transport_albedo(
                    atmosphere, wavelength, interpolation);
                double ground_cotangent = 0.0;
                for (const auto& source : interpolation.ground()) {
                    ground_cotangent +=
                        source.weight *
                        value_gradient(offsets[row] + source.row_inner_index);
                }
                cumulative_cotangent =
                    -transmission_before * ground_albedo * ground_cotangent;
                accumulate_ground_transport_albedo_vjp(
                    atmosphere, interpolation,
                    transmission_before * ground_cotangent, native_gradient);
            }

            // Reverse the observer-to-end source traversal. A layer optical
            // depth affects every source closer to the observer and the
            // terminal ground term, while its own source fraction contributes
            // a direct derivative.
            for (int layer_index = 0; layer_index < num_layers; ++layer_index) {
                const auto atmosphere_weights =
                    interpolation.atmosphere_for_layer(layer_index);
                const double factor_cotangent =
                    workspace.factor_cotangent(layer_index);
                const double transmission_before =
                    workspace.transmission_before(layer_index);
                const double layer_transmission =
                    std::exp(-workspace.optical_depth(layer_index));
                const double albedo = workspace.albedo(layer_index);
                const double source_fraction =
                    workspace.source_fraction(layer_index);

                const double optical_depth_cotangent =
                    cumulative_cotangent + factor_cotangent *
                                               transmission_before * albedo *
                                               layer_transmission;
                const auto od_weights =
                    interpolation.optical_depth_for_layer(layer_index);
                for (std::size_t index = 0; index < od_weights.size();
                     ++index) {
                    const auto [atmosphere_index, weight] = od_weights[index];
                    native_gradient(atmosphere_index) +=
                        weight * optical_depth_cotangent;
                }

                const double albedo_cotangent =
                    factor_cotangent * transmission_before * source_fraction;
                for (const auto& weight : atmosphere_weights) {
                    native_gradient(atmosphere.ssa_deriv_start_index() +
                                    weight.index) +=
                        weight.weight * albedo_cotangent;
                }

                cumulative_cotangent -= factor_cotangent * transmission_before *
                                        albedo * source_fraction;
            }
        }
    }

    template void RayTransportMap::assemble_values<1>(
        const sasktran2::atmosphere::Atmosphere<1>&, int,
        TransportOperator&) const;
    template void RayTransportMap::assemble_values<3>(
        const sasktran2::atmosphere::Atmosphere<3>&, int,
        TransportOperator&) const;
    template void RayTransportMap::assemble_jvp<1>(
        const sasktran2::atmosphere::Atmosphere<1>&, int,
        Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>) const;
    template void RayTransportMap::assemble_jvp<3>(
        const sasktran2::atmosphere::Atmosphere<3>&, int,
        Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>) const;
    template void RayTransportMap::accumulate_vjp<1>(
        const sasktran2::atmosphere::Atmosphere<1>&, int,
        Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>,
        RayTransportWorkspace&) const;
    template void RayTransportMap::accumulate_vjp<3>(
        const sasktran2::atmosphere::Atmosphere<3>&, int,
        Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>,
        RayTransportWorkspace&) const;

} // namespace sasktran2::successive_orders
