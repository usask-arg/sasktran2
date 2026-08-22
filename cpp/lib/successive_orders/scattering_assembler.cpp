#include "scattering_assembler.h"

#include <sasktran2/math/scattering.h>
#include <sasktran2/math/wigner.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace sasktran2::successive_orders {
    namespace {
        ScatteringBlockLayout make_layout(const SourceGeometry1D& geometry) {
            if (geometry.num_interior_points() < 1 ||
                geometry.num_points() < 1 ||
                geometry.incoming_point_offsets().size() !=
                    static_cast<std::size_t>(geometry.num_points()) + 1 ||
                geometry.outgoing_point_offsets().size() !=
                    static_cast<std::size_t>(geometry.num_points()) + 1) {
                throw std::invalid_argument(
                    "successive-orders scattering assembly requires "
                    "initialized source geometry");
            }
            return {geometry.num_interior_points(), 1,
                    geometry.incoming_point_offsets(),
                    geometry.outgoing_point_offsets()};
        }

        std::shared_ptr<const ScalarAngularBasis>
        make_basis(const SourceGeometry1D& geometry, int num_coefficients) {
            if (num_coefficients < 1 || geometry.num_interior_points() < 1) {
                throw std::invalid_argument(
                    "invalid scalar successive-orders coefficient count");
            }
            const auto& point = geometry.source_point(0);
            return std::make_shared<const ScalarAngularBasis>(
                point.incoming_sphere(), point.outgoing_sphere(),
                num_coefficients);
        }

        int checked_dense_size(int rows, int columns) {
            const auto size = static_cast<long long>(rows) * columns;
            if (rows < 0 || columns < 0 ||
                size > std::numeric_limits<int>::max()) {
                throw std::length_error(
                    "successive-orders ground scattering size exceeds the "
                    "supported range");
            }
            return static_cast<int>(size);
        }

        Eigen::Vector3d
        deterministic_horizontal(const Eigen::Vector3d& normal) {
            const Eigen::Vector3d reference = std::abs(normal.z()) < 0.9
                                                  ? Eigen::Vector3d::UnitZ()
                                                  : Eigen::Vector3d::UnitX();
            return (reference - reference.dot(normal) * normal).normalized();
        }

        Eigen::Vector3d
        horizontal_direction(const Eigen::Vector3d& direction,
                             const Eigen::Vector3d& normal,
                             const Eigen::Vector3d& vertical_fallback) {
            const Eigen::Vector3d projection =
                direction - direction.dot(normal) * normal;
            const double norm = projection.norm();
            if (!std::isfinite(norm) ||
                norm <= 64.0 * std::numeric_limits<double>::epsilon()) {
                // Azimuth is physically undefined for a vertical direction,
                // but an anisotropic BRDF still requires a finite convention.
                return vertical_fallback;
            }
            return projection / norm;
        }

        double relative_azimuth(const Eigen::Vector3d& horizontal_in,
                                const Eigen::Vector3d& horizontal_out) {
            const double cosine =
                std::clamp(horizontal_in.dot(horizontal_out), -1.0, 1.0);
            // Preserve the legacy HR relative-azimuth convention.
            return EIGEN_PI - std::acos(cosine);
        }

        std::vector<int> stokes_offsets(const std::vector<int>& offsets,
                                        int stokes_components) {
            std::vector<int> result(offsets.size());
            for (std::size_t index = 0; index < offsets.size(); ++index) {
                if (offsets[index] < 0 ||
                    offsets[index] >
                        std::numeric_limits<int>::max() / stokes_components) {
                    throw std::length_error(
                        "successive-orders Stokes layout exceeds the "
                        "supported range");
                }
                result[index] = offsets[index] * stokes_components;
            }
            return result;
        }

        ScatteringBlockLayout
        make_vector_layout(const SourceGeometry1D& geometry) {
            // Reuse the scalar validation before applying the Stokes stride.
            static_cast<void>(make_layout(geometry));
            return {geometry.num_interior_points(), 3,
                    stokes_offsets(geometry.incoming_point_offsets(), 3),
                    stokes_offsets(geometry.outgoing_point_offsets(), 3)};
        }

        std::vector<int> make_dense_offsets(const ScatteringBlockLayout& layout,
                                            int first_block, int block_count) {
            std::vector<int> offsets(static_cast<std::size_t>(block_count) + 1,
                                     0);
            for (int local_block = 0; local_block < block_count;
                 ++local_block) {
                const int block = first_block + local_block;
                const int block_size =
                    checked_dense_size(layout.output_block_size(block),
                                       layout.input_block_size(block));
                if (offsets[local_block] >
                    std::numeric_limits<int>::max() - block_size) {
                    throw std::length_error(
                        "successive-orders dense scattering storage exceeds "
                        "the supported range");
                }
                offsets[local_block + 1] = offsets[local_block] + block_size;
            }
            return offsets;
        }

        template <typename Storage>
        double delta_m_coefficient_scale(const Storage& storage, int location,
                                         int wavelength) {
            if (storage.applied_f_order <= 0) {
                return 0.0;
            }
            const double f = storage.f(location, wavelength);
            return f / (1.0 - f);
        }

        template <typename Storage>
        double delta_m_derivative_scale(const Storage& storage, int location,
                                        int wavelength, int group) {
            if (storage.applied_f_order <= 0) {
                return 0.0;
            }
            const double f = storage.f(location, wavelength);
            return storage.d_f(location, wavelength, group) /
                   ((1.0 - f) * (1.0 - f));
        }

        bool delta_m_corrects_vector_component(int coefficient) {
            // The forward delta peak contributes to the three generalized
            // spherical-function "a" coefficients, but not to b1.
            return coefficient % 4 != 3;
        }
    } // namespace

    ScalarScatteringAssembler::ScalarScatteringAssembler(
        const SourceGeometry1D& geometry, int num_coefficients)
        : m_geometry(&geometry), m_layout(make_layout(geometry)),
          m_angular_basis(make_basis(geometry, num_coefficients)) {
        for (int point_index = 0; point_index < geometry.num_interior_points();
             ++point_index) {
            const auto& point = geometry.source_point(point_index);
            if (point.num_incoming() != m_angular_basis->input_size() ||
                point.num_outgoing() != m_angular_basis->output_size()) {
                throw std::invalid_argument(
                    "scalar successive-orders atmospheric angular grids "
                    "must share one basis");
            }
        }

        m_ground_value_offsets.assign(
            static_cast<std::size_t>(geometry.num_ground_points()) + 1, 0);
        m_ground_geometry.resize(geometry.num_ground_points());
        for (int ground_index = 0; ground_index < geometry.num_ground_points();
             ++ground_index) {
            const int point_index =
                geometry.num_interior_points() + ground_index;
            const auto& point = geometry.source_point(point_index);
            const int rows = point.num_outgoing();
            const int columns = point.num_incoming();
            const int block_size = checked_dense_size(rows, columns);
            if (m_ground_value_offsets[ground_index] >
                std::numeric_limits<int>::max() - block_size) {
                throw std::length_error(
                    "successive-orders ground scattering storage exceeds "
                    "the supported range");
            }
            m_ground_value_offsets[ground_index + 1] =
                m_ground_value_offsets[ground_index] + block_size;

            auto& angular = m_ground_geometry[ground_index];
            angular.mu_in.resize(columns);
            angular.weighted_mu_in.resize(columns);
            angular.mu_out.resize(rows);
            angular.phi_difference.resize(rows, columns);
            const Eigen::Vector3d normal =
                point.location().position.normalized();
            const Eigen::Vector3d vertical_fallback =
                deterministic_horizontal(normal);
            std::vector<Eigen::Vector3d> horizontal_in(columns);
            std::vector<Eigen::Vector3d> horizontal_out(rows);
            for (int input = 0; input < columns; ++input) {
                const Eigen::Vector3d direction =
                    point.incoming_sphere().get_quad_position(input);
                angular.mu_in(input) =
                    point.location().cos_zenith_angle(direction);
                angular.weighted_mu_in(input) =
                    4.0 * EIGEN_PI * angular.mu_in(input) *
                    point.incoming_sphere().quadrature_weight(input);
                horizontal_in[input] =
                    horizontal_direction(direction, normal, vertical_fallback);
            }
            for (int output = 0; output < rows; ++output) {
                const Eigen::Vector3d direction =
                    point.outgoing_sphere().get_quad_position(output);
                angular.mu_out(output) =
                    point.location().cos_zenith_angle(direction);
                horizontal_out[output] =
                    horizontal_direction(direction, normal, vertical_fallback);
            }
            for (int output = 0; output < rows; ++output) {
                for (int input = 0; input < columns; ++input) {
                    angular.phi_difference(output, input) = relative_azimuth(
                        horizontal_in[input], horizontal_out[output]);
                }
            }
        }
    }

    ScatteringOperator<1> ScalarScatteringAssembler::create_operator() const {
        return {m_layout, m_angular_basis};
    }

    void ScalarScatteringAssembler::validate_atmosphere(
        const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
        int wavelength) const {
        const auto& storage = atmosphere.storage();
        if (wavelength < 0 || wavelength >= atmosphere.num_wavel() ||
            m_angular_basis->num_coefficients() >
                storage.max_stored_legendre() ||
            atmosphere.surface().brdf_args().cols() != atmosphere.num_wavel()) {
            throw std::invalid_argument(
                "atmosphere does not match scalar successive-orders "
                "scattering assembly");
        }
        for (const auto& point : m_geometry->source_points()) {
            for (const auto& weight : point.atmosphere_weights()) {
                if (weight.index < 0 ||
                    weight.index >= storage.total_extinction.rows()) {
                    throw std::invalid_argument(
                        "successive-orders source interpolation exceeds the "
                        "atmosphere grid");
                }
            }
        }
    }

    void ScalarScatteringAssembler::validate_operator(
        const ScatteringOperator<1>& scattering) const {
        if (scattering.input_offsets() != m_layout.input_offsets() ||
            scattering.output_offsets() != m_layout.output_offsets() ||
            scattering.layout().atmospheric_blocks() !=
                m_layout.atmospheric_blocks() ||
            scattering.atmospheric_coefficients().rows() !=
                m_geometry->num_interior_points() ||
            scattering.atmospheric_coefficients().cols() !=
                m_angular_basis->num_coefficients() ||
            scattering.ground_value_size() != ground_value_size()) {
            throw std::invalid_argument(
                "scattering operator does not match its successive-orders "
                "assembler");
        }
    }

    void ScalarScatteringAssembler::validate_native_product(
        const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
        Eigen::Index native_size) const {
        if (native_size != atmosphere.num_deriv()) {
            throw std::invalid_argument(
                "invalid native tangent or gradient size for scalar "
                "successive-orders scattering");
        }
        const auto& storage = atmosphere.storage();
        if (atmosphere.num_scattering_deriv_groups() !=
                storage.d_leg_coeff.dimension(3) ||
            (atmosphere.num_scattering_deriv_groups() != 0 &&
             (storage.d_leg_coeff.dimension(0) <
                  m_angular_basis->num_coefficients() ||
              storage.d_leg_coeff.dimension(1) !=
                  storage.total_extinction.rows() ||
              storage.d_leg_coeff.dimension(2) != atmosphere.num_wavel()))) {
            throw std::invalid_argument(
                "invalid native Legendre derivative storage for scalar "
                "successive-orders scattering");
        }
    }

    void ScalarScatteringAssembler::assemble_ground(
        const sasktran2::atmosphere::Atmosphere<1>& atmosphere, int wavelength,
        Eigen::VectorXd& values) const {
        values.resize(ground_value_size());
        for (int ground_index = 0;
             ground_index < m_geometry->num_ground_points(); ++ground_index) {
            const auto& angular = m_ground_geometry[ground_index];
            const int rows = static_cast<int>(angular.mu_out.size());
            const int columns = static_cast<int>(angular.mu_in.size());
            const int offset = m_ground_value_offsets[ground_index];
            for (int output = 0; output < rows; ++output) {
                for (int input = 0; input < columns; ++input) {
                    double brdf_00 = 1.0 / EIGEN_PI;
                    if (!atmosphere.surface().has_spatial_lambertian_albedo()) {
                        brdf_00 = atmosphere.surface().brdf(
                            wavelength, angular.mu_in(input),
                            angular.mu_out(output),
                            angular.phi_difference(output, input),
                            m_geometry->ground_horizontal_weights(
                                ground_index))(0, 0);
                    }
                    values(offset + output * columns + input) =
                        angular.weighted_mu_in(input) * brdf_00;
                }
            }
        }
    }

    void ScalarScatteringAssembler::assemble_values(
        const sasktran2::atmosphere::Atmosphere<1>& atmosphere, int wavelength,
        ScatteringOperator<1>& scattering) const {
        validate_atmosphere(atmosphere, wavelength);
        validate_operator(scattering);
        auto& coefficients = scattering.atmospheric_coefficients();
        coefficients.setZero();
        const auto& storage = atmosphere.storage();
        for (int point_index = 0;
             point_index < m_geometry->num_interior_points(); ++point_index) {
            for (const auto& weight :
                 m_geometry->source_point(point_index).atmosphere_weights()) {
                const double delta_m_scale = delta_m_coefficient_scale(
                    storage, weight.index, wavelength);
                for (int degree = 0;
                     degree < m_angular_basis->num_coefficients(); ++degree) {
                    const double coefficient =
                        storage.leg_coeff(degree, weight.index, wavelength) -
                        (2.0 * degree + 1.0) * delta_m_scale;
                    coefficients(point_index, degree) +=
                        weight.weight * coefficient;
                }
            }
        }
        int active_coefficients = 1;
        for (int degree = coefficients.cols() - 1; degree >= 1; --degree) {
            if (!coefficients.col(degree).isZero(0.0)) {
                active_coefficients = degree + 1;
                break;
            }
        }
        scattering.set_active_coefficients(active_coefficients);
        assemble_ground(atmosphere, wavelength, scattering.ground_values());
    }

    void ScalarScatteringAssembler::assemble_jvp(
        const sasktran2::atmosphere::Atmosphere<1>& atmosphere, int wavelength,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::MatrixXd& atmospheric_coefficient_tangent,
        Eigen::VectorXd& ground_value_tangent) const {
        validate_atmosphere(atmosphere, wavelength);
        validate_native_product(atmosphere, native_tangent.size());
        atmospheric_coefficient_tangent.setZero(
            m_geometry->num_interior_points(),
            m_angular_basis->num_coefficients());
        ground_value_tangent.setZero(ground_value_size());
        if (atmosphere.num_deriv() == 0) {
            return;
        }

        const auto& storage = atmosphere.storage();
        const int locations = storage.total_extinction.rows();
        for (int point_index = 0;
             point_index < m_geometry->num_interior_points(); ++point_index) {
            for (const auto& weight :
                 m_geometry->source_point(point_index).atmosphere_weights()) {
                for (int group = 0;
                     group < atmosphere.num_scattering_deriv_groups();
                     ++group) {
                    const double direction =
                        native_tangent(atmosphere.scat_deriv_start_index() +
                                       group * locations + weight.index);
                    if (direction == 0.0) {
                        continue;
                    }
                    const double delta_m_scale = delta_m_derivative_scale(
                        storage, weight.index, wavelength, group);
                    for (int degree = 0;
                         degree < m_angular_basis->num_coefficients();
                         ++degree) {
                        const double coefficient_tangent =
                            storage.d_leg_coeff(degree, weight.index,
                                                wavelength, group) -
                            (2.0 * degree + 1.0) * delta_m_scale;
                        atmospheric_coefficient_tangent(point_index, degree) +=
                            weight.weight * direction * coefficient_tangent;
                    }
                }
            }
        }

        // Spatial Lambertian albedo is applied by each ray at its physical
        // ground intersection. The ground source points therefore represent
        // unit-albedo angular redistribution and have no albedo derivative.
        if (atmosphere.surface().has_spatial_lambertian_albedo()) {
            return;
        }

        for (int ground_index = 0;
             ground_index < m_geometry->num_ground_points(); ++ground_index) {
            const auto& angular = m_ground_geometry[ground_index];
            const int rows = static_cast<int>(angular.mu_out.size());
            const int columns = static_cast<int>(angular.mu_in.size());
            const int offset = m_ground_value_offsets[ground_index];
            for (int derivative = 0;
                 derivative < atmosphere.surface().num_deriv(); ++derivative) {
                const double direction = native_tangent(
                    atmosphere.surface_deriv_start_index() + derivative);
                if (direction == 0.0) {
                    continue;
                }
                for (int output = 0; output < rows; ++output) {
                    for (int input = 0; input < columns; ++input) {
                        const auto derivative_brdf =
                            atmosphere.surface().d_brdf(
                                wavelength, angular.mu_in(input),
                                angular.mu_out(output),
                                angular.phi_difference(output, input),
                                derivative,
                                m_geometry->ground_horizontal_weights(
                                    ground_index));
                        ground_value_tangent(offset + output * columns +
                                             input) +=
                            direction * angular.weighted_mu_in(input) *
                            derivative_brdf(0, 0);
                    }
                }
            }
        }
    }

    void ScalarScatteringAssembler::accumulate_vjp(
        const sasktran2::atmosphere::Atmosphere<1>& atmosphere, int wavelength,
        Eigen::Ref<const Eigen::MatrixXd> atmospheric_coefficient_gradient,
        Eigen::Ref<const Eigen::VectorXd> ground_value_gradient,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        validate_atmosphere(atmosphere, wavelength);
        validate_native_product(atmosphere, native_gradient.size());
        if (atmospheric_coefficient_gradient.rows() !=
                m_geometry->num_interior_points() ||
            atmospheric_coefficient_gradient.cols() !=
                m_angular_basis->num_coefficients() ||
            ground_value_gradient.size() != ground_value_size()) {
            throw std::invalid_argument(
                "invalid scattering cotangent dimensions for native VJP");
        }
        if (atmosphere.num_deriv() == 0) {
            return;
        }

        const auto& storage = atmosphere.storage();
        const int locations = storage.total_extinction.rows();
        for (int point_index = 0;
             point_index < m_geometry->num_interior_points(); ++point_index) {
            for (const auto& weight :
                 m_geometry->source_point(point_index).atmosphere_weights()) {
                for (int group = 0;
                     group < atmosphere.num_scattering_deriv_groups();
                     ++group) {
                    double value = 0.0;
                    const double delta_m_scale = delta_m_derivative_scale(
                        storage, weight.index, wavelength, group);
                    for (int degree = 0;
                         degree < m_angular_basis->num_coefficients();
                         ++degree) {
                        const double coefficient_derivative =
                            storage.d_leg_coeff(degree, weight.index,
                                                wavelength, group) -
                            (2.0 * degree + 1.0) * delta_m_scale;
                        value += atmospheric_coefficient_gradient(point_index,
                                                                  degree) *
                                 weight.weight * coefficient_derivative;
                    }
                    native_gradient(atmosphere.scat_deriv_start_index() +
                                    group * locations + weight.index) += value;
                }
            }
        }
        if (atmosphere.surface().has_spatial_lambertian_albedo()) {
            return;
        }
        for (int ground_index = 0;
             ground_index < m_geometry->num_ground_points(); ++ground_index) {
            const auto& angular = m_ground_geometry[ground_index];
            const int rows = static_cast<int>(angular.mu_out.size());
            const int columns = static_cast<int>(angular.mu_in.size());
            const int offset = m_ground_value_offsets[ground_index];
            for (int derivative = 0;
                 derivative < atmosphere.surface().num_deriv(); ++derivative) {
                double value = 0.0;
                for (int output = 0; output < rows; ++output) {
                    for (int input = 0; input < columns; ++input) {
                        const auto derivative_brdf =
                            atmosphere.surface().d_brdf(
                                wavelength, angular.mu_in(input),
                                angular.mu_out(output),
                                angular.phi_difference(output, input),
                                derivative,
                                m_geometry->ground_horizontal_weights(
                                    ground_index));
                        value += ground_value_gradient(
                                     offset + output * columns + input) *
                                 angular.weighted_mu_in(input) *
                                 derivative_brdf(0, 0);
                    }
                }
                native_gradient(atmosphere.surface_deriv_start_index() +
                                derivative) += value;
            }
        }
    }

    std::size_t ScalarScatteringAssembler::storage_bytes() const {
        std::size_t result =
            m_layout.storage_bytes() + m_angular_basis->storage_bytes() +
            m_ground_value_offsets.capacity() * sizeof(int) +
            m_ground_geometry.capacity() * sizeof(GroundAngularGeometry);
        for (const auto& geometry : m_ground_geometry) {
            result += geometry.storage_bytes();
        }
        return result;
    }

    VectorScatteringAssembler::VectorScatteringAssembler(
        const SourceGeometry1D& geometry, int num_coefficients)
        : m_geometry(&geometry), m_layout(make_vector_layout(geometry)),
          m_ground_value_offsets(
              make_dense_offsets(m_layout, geometry.num_interior_points(),
                                 geometry.num_ground_points())) {
        if (num_coefficients < 1) {
            throw std::invalid_argument(
                "invalid vector successive-orders coefficient count");
        }

        const auto& atmospheric_point = geometry.source_point(0);
        m_angular_basis = std::make_shared<const VectorAngularBasis>(
            atmospheric_point.incoming_sphere(),
            atmospheric_point.outgoing_sphere(), num_coefficients);

        for (int point_index = 0; point_index < geometry.num_interior_points();
             ++point_index) {
            const auto& point = geometry.source_point(point_index);
            if (point.num_incoming() != m_angular_basis->input_directions() ||
                point.num_outgoing() != m_angular_basis->output_directions()) {
                throw std::invalid_argument(
                    "vector successive-orders atmospheric angular grids "
                    "must share one basis");
            }
        }

        m_ground_geometry.resize(geometry.num_ground_points());
        for (int ground_index = 0; ground_index < geometry.num_ground_points();
             ++ground_index) {
            const int point_index =
                geometry.num_interior_points() + ground_index;
            const auto& point = geometry.source_point(point_index);
            const int rows = point.num_outgoing();
            const int columns = point.num_incoming();
            auto& ground = m_ground_geometry[ground_index];
            ground.mu_in.resize(columns);
            ground.weighted_mu_in.resize(columns);
            ground.mu_out.resize(rows);
            ground.phi_difference.resize(rows, columns);
            const Eigen::Vector3d normal =
                point.location().position.normalized();
            const Eigen::Vector3d vertical_fallback =
                deterministic_horizontal(normal);
            std::vector<Eigen::Vector3d> horizontal_in(columns);
            std::vector<Eigen::Vector3d> horizontal_out(rows);
            for (int input = 0; input < columns; ++input) {
                const Eigen::Vector3d direction =
                    point.incoming_sphere().get_quad_position(input);
                ground.mu_in(input) =
                    point.location().cos_zenith_angle(direction);
                ground.weighted_mu_in(input) =
                    4.0 * EIGEN_PI * ground.mu_in(input) *
                    point.incoming_sphere().quadrature_weight(input);
                horizontal_in[input] =
                    horizontal_direction(direction, normal, vertical_fallback);
            }
            for (int output = 0; output < rows; ++output) {
                const Eigen::Vector3d direction =
                    point.outgoing_sphere().get_quad_position(output);
                ground.mu_out(output) =
                    point.location().cos_zenith_angle(direction);
                horizontal_out[output] =
                    horizontal_direction(direction, normal, vertical_fallback);
            }
            for (int output = 0; output < rows; ++output) {
                for (int input = 0; input < columns; ++input) {
                    ground.phi_difference(output, input) = relative_azimuth(
                        horizontal_in[input], horizontal_out[output]);
                }
            }
        }
    }

    ScatteringOperator<3> VectorScatteringAssembler::create_operator() const {
        return ScatteringOperator<3>(m_layout, m_angular_basis);
    }

    void VectorScatteringAssembler::validate_atmosphere(
        const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
        int wavelength) const {
        const auto& storage = atmosphere.storage();
        if (wavelength < 0 || wavelength >= atmosphere.num_wavel() ||
            num_coefficients() > storage.max_stored_legendre() ||
            atmosphere.surface().brdf_args().cols() != atmosphere.num_wavel()) {
            throw std::invalid_argument(
                "atmosphere does not match vector successive-orders "
                "scattering assembly");
        }
        for (const auto& point : m_geometry->source_points()) {
            for (const auto& weight : point.atmosphere_weights()) {
                if (weight.index < 0 ||
                    weight.index >= storage.total_extinction.rows()) {
                    throw std::invalid_argument(
                        "successive-orders source interpolation exceeds the "
                        "atmosphere grid");
                }
            }
        }
    }

    void VectorScatteringAssembler::validate_operator(
        const ScatteringOperator<3>& scattering) const {
        if (scattering.input_offsets() != m_layout.input_offsets() ||
            scattering.output_offsets() != m_layout.output_offsets() ||
            scattering.layout().atmospheric_blocks() !=
                m_layout.atmospheric_blocks() ||
            scattering.coefficient_value_size() != coefficient_value_size() ||
            scattering.ground_value_size() != ground_value_size()) {
            throw std::invalid_argument(
                "scattering operator does not match its vector "
                "successive-orders assembler");
        }
    }

    void VectorScatteringAssembler::validate_native_product(
        const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
        Eigen::Index native_size) const {
        if (native_size != atmosphere.num_deriv()) {
            throw std::invalid_argument(
                "invalid native tangent or gradient size for vector "
                "successive-orders scattering");
        }
        const auto& storage = atmosphere.storage();
        if (atmosphere.num_scattering_deriv_groups() !=
                storage.d_leg_coeff.dimension(3) ||
            (atmosphere.num_scattering_deriv_groups() != 0 &&
             (storage.d_leg_coeff.dimension(0) < 4 * num_coefficients() ||
              storage.d_leg_coeff.dimension(1) !=
                  storage.total_extinction.rows() ||
              storage.d_leg_coeff.dimension(2) != atmosphere.num_wavel()))) {
            throw std::invalid_argument(
                "invalid native Greek-coefficient derivative storage for "
                "vector successive-orders scattering");
        }
    }

    void VectorScatteringAssembler::assemble_ground(
        const sasktran2::atmosphere::Atmosphere<3>& atmosphere, int wavelength,
        Eigen::VectorXd& values) const {
        values.setZero(ground_value_size());
        for (int ground_index = 0;
             ground_index < m_geometry->num_ground_points(); ++ground_index) {
            const auto& angular = m_ground_geometry[ground_index];
            const int rows = static_cast<int>(angular.mu_out.size());
            const int columns = static_cast<int>(angular.mu_in.size());
            Eigen::Map<RowMajorMatrix> block(
                values.data() + m_ground_value_offsets[ground_index], 3 * rows,
                3 * columns);
            for (int output = 0; output < rows; ++output) {
                for (int input = 0; input < columns; ++input) {
                    double brdf_00 = 1.0 / EIGEN_PI;
                    if (!atmosphere.surface().has_spatial_lambertian_albedo()) {
                        brdf_00 = atmosphere.surface().brdf(
                            wavelength, angular.mu_in(input),
                            angular.mu_out(output),
                            angular.phi_difference(output, input),
                            m_geometry->ground_horizontal_weights(
                                ground_index))(0, 0);
                    }
                    // Match the legacy HR ground convention: the BRDF only
                    // contributes to I <- I for vector radiance.
                    block(3 * output, 3 * input) =
                        angular.weighted_mu_in(input) * brdf_00;
                }
            }
        }
    }

    void VectorScatteringAssembler::assemble_values(
        const sasktran2::atmosphere::Atmosphere<3>& atmosphere, int wavelength,
        ScatteringOperator<3>& scattering) const {
        validate_atmosphere(atmosphere, wavelength);
        validate_operator(scattering);
        auto& coefficients = scattering.atmospheric_coefficients();
        coefficients.setZero();
        const auto& storage = atmosphere.storage();
        for (int point_index = 0;
             point_index < m_geometry->num_interior_points(); ++point_index) {
            for (const auto& weight :
                 m_geometry->source_point(point_index).atmosphere_weights()) {
                const double delta_m_scale = delta_m_coefficient_scale(
                    storage, weight.index, wavelength);
                for (int coefficient = 0; coefficient < coefficients.cols();
                     ++coefficient) {
                    const int degree = coefficient / 4;
                    const double delta_m_correction =
                        delta_m_corrects_vector_component(coefficient)
                            ? (2.0 * degree + 1.0) * delta_m_scale
                            : 0.0;
                    coefficients(point_index, coefficient) +=
                        weight.weight *
                        (storage.leg_coeff(coefficient, weight.index,
                                           wavelength) -
                         delta_m_correction);
                }
            }
        }
        int active_coefficients = 1;
        for (int degree = num_coefficients() - 1; degree >= 1; --degree) {
            if (!coefficients.middleCols(4 * degree, 4).isZero(0.0)) {
                active_coefficients = degree + 1;
                break;
            }
        }
        scattering.set_active_coefficients(active_coefficients);
        assemble_ground(atmosphere, wavelength, scattering.ground_values());
    }

    void VectorScatteringAssembler::assemble_jvp(
        const sasktran2::atmosphere::Atmosphere<3>& atmosphere, int wavelength,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::MatrixXd& atmospheric_coefficient_tangent,
        Eigen::VectorXd& ground_value_tangent) const {
        validate_atmosphere(atmosphere, wavelength);
        validate_native_product(atmosphere, native_tangent.size());
        atmospheric_coefficient_tangent.setZero(
            m_geometry->num_interior_points(), 4 * num_coefficients());
        ground_value_tangent.setZero(ground_value_size());
        if (atmosphere.num_deriv() == 0) {
            return;
        }

        const auto& storage = atmosphere.storage();
        const int locations = storage.total_extinction.rows();
        for (int point_index = 0;
             point_index < m_geometry->num_interior_points(); ++point_index) {
            for (const auto& weight :
                 m_geometry->source_point(point_index).atmosphere_weights()) {
                for (int group = 0;
                     group < atmosphere.num_scattering_deriv_groups();
                     ++group) {
                    const double direction =
                        native_tangent(atmosphere.scat_deriv_start_index() +
                                       group * locations + weight.index);
                    if (direction == 0.0) {
                        continue;
                    }
                    const double delta_m_scale = delta_m_derivative_scale(
                        storage, weight.index, wavelength, group);
                    for (int coefficient = 0;
                         coefficient < atmospheric_coefficient_tangent.cols();
                         ++coefficient) {
                        const int degree = coefficient / 4;
                        const double delta_m_correction =
                            delta_m_corrects_vector_component(coefficient)
                                ? (2.0 * degree + 1.0) * delta_m_scale
                                : 0.0;
                        atmospheric_coefficient_tangent(point_index,
                                                        coefficient) +=
                            weight.weight * direction *
                            (storage.d_leg_coeff(coefficient, weight.index,
                                                 wavelength, group) -
                             delta_m_correction);
                    }
                }
            }
        }

        if (atmosphere.surface().has_spatial_lambertian_albedo()) {
            return;
        }

        for (int ground_index = 0;
             ground_index < m_geometry->num_ground_points(); ++ground_index) {
            const auto& angular = m_ground_geometry[ground_index];
            const int rows = static_cast<int>(angular.mu_out.size());
            const int columns = static_cast<int>(angular.mu_in.size());
            Eigen::Map<RowMajorMatrix> block(
                ground_value_tangent.data() +
                    m_ground_value_offsets[ground_index],
                3 * rows, 3 * columns);
            for (int derivative = 0;
                 derivative < atmosphere.surface().num_deriv(); ++derivative) {
                const double direction = native_tangent(
                    atmosphere.surface_deriv_start_index() + derivative);
                if (direction == 0.0) {
                    continue;
                }
                for (int output = 0; output < rows; ++output) {
                    for (int input = 0; input < columns; ++input) {
                        const auto derivative_brdf =
                            atmosphere.surface().d_brdf(
                                wavelength, angular.mu_in(input),
                                angular.mu_out(output),
                                angular.phi_difference(output, input),
                                derivative,
                                m_geometry->ground_horizontal_weights(
                                    ground_index));
                        block(3 * output, 3 * input) +=
                            direction * angular.weighted_mu_in(input) *
                            derivative_brdf(0, 0);
                    }
                }
            }
        }
    }

    void VectorScatteringAssembler::accumulate_vjp(
        const sasktran2::atmosphere::Atmosphere<3>& atmosphere, int wavelength,
        Eigen::Ref<const Eigen::MatrixXd> atmospheric_coefficient_gradient,
        Eigen::Ref<const Eigen::VectorXd> ground_value_gradient,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        validate_atmosphere(atmosphere, wavelength);
        validate_native_product(atmosphere, native_gradient.size());
        if (atmospheric_coefficient_gradient.rows() !=
                m_geometry->num_interior_points() ||
            atmospheric_coefficient_gradient.cols() != 4 * num_coefficients() ||
            ground_value_gradient.size() != ground_value_size()) {
            throw std::invalid_argument(
                "invalid vector scattering cotangent dimensions for native "
                "VJP");
        }
        if (atmosphere.num_deriv() == 0) {
            return;
        }

        const auto& storage = atmosphere.storage();
        const int locations = storage.total_extinction.rows();
        for (int point_index = 0;
             point_index < m_geometry->num_interior_points(); ++point_index) {
            for (const auto& weight :
                 m_geometry->source_point(point_index).atmosphere_weights()) {
                for (int group = 0;
                     group < atmosphere.num_scattering_deriv_groups();
                     ++group) {
                    double value = 0.0;
                    const double delta_m_scale = delta_m_derivative_scale(
                        storage, weight.index, wavelength, group);
                    for (int coefficient = 0;
                         coefficient < atmospheric_coefficient_gradient.cols();
                         ++coefficient) {
                        const int degree = coefficient / 4;
                        const double delta_m_correction =
                            delta_m_corrects_vector_component(coefficient)
                                ? (2.0 * degree + 1.0) * delta_m_scale
                                : 0.0;
                        value += atmospheric_coefficient_gradient(point_index,
                                                                  coefficient) *
                                 weight.weight *
                                 (storage.d_leg_coeff(coefficient, weight.index,
                                                      wavelength, group) -
                                  delta_m_correction);
                    }
                    native_gradient(atmosphere.scat_deriv_start_index() +
                                    group * locations + weight.index) += value;
                }
            }
        }
        if (atmosphere.surface().has_spatial_lambertian_albedo()) {
            return;
        }
        for (int ground_index = 0;
             ground_index < m_geometry->num_ground_points(); ++ground_index) {
            const auto& angular = m_ground_geometry[ground_index];
            const int rows = static_cast<int>(angular.mu_out.size());
            const int columns = static_cast<int>(angular.mu_in.size());
            Eigen::Map<const RowMajorMatrix> block_gradient(
                ground_value_gradient.data() +
                    m_ground_value_offsets[ground_index],
                3 * rows, 3 * columns);
            for (int derivative = 0;
                 derivative < atmosphere.surface().num_deriv(); ++derivative) {
                double value = 0.0;
                for (int output = 0; output < rows; ++output) {
                    for (int input = 0; input < columns; ++input) {
                        const auto derivative_brdf =
                            atmosphere.surface().d_brdf(
                                wavelength, angular.mu_in(input),
                                angular.mu_out(output),
                                angular.phi_difference(output, input),
                                derivative,
                                m_geometry->ground_horizontal_weights(
                                    ground_index));
                        value += block_gradient(3 * output, 3 * input) *
                                 angular.weighted_mu_in(input) *
                                 derivative_brdf(0, 0);
                    }
                }
                native_gradient(atmosphere.surface_deriv_start_index() +
                                derivative) += value;
            }
        }
    }

    std::size_t VectorScatteringAssembler::storage_bytes() const {
        std::size_t result =
            m_layout.storage_bytes() + m_angular_basis->storage_bytes() +
            m_ground_value_offsets.capacity() * sizeof(int) +
            m_ground_geometry.capacity() * sizeof(GroundAngularGeometry);
        for (const auto& geometry : m_ground_geometry) {
            result += geometry.storage_bytes();
        }
        return result;
    }

} // namespace sasktran2::successive_orders
