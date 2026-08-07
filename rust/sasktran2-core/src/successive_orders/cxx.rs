use std::pin::Pin;
use std::sync::Arc;

use anyhow::{Result, anyhow};

use super::simd::{
    Workspace as SimdWorkspace,
    batch::{LANES, extract_lane, interleave_lane},
};
use super::{
    BlockDiagonalMatrix, CsrMatrix, FixedPointProblem, ScalarRayTransport, SolverConfig,
    SolverDiagnostics, SolverError, SuccessiveOrdersSolver,
};

use super::{
    ScalarCoefficientBasis, ScalarCoefficientScattering, ScatteringCoefficientInterpolator,
    SolarTransmissionOperator, SourceInterpolator, VectorCoefficientBasis,
    VectorCoefficientScattering, accumulate_atmospheric_coefficient_vjp,
    atmospheric_coefficient_jvp,
};

#[cxx::bridge(namespace = "sasktran2::rust::successive_orders")]
pub mod ffi {
    extern "Rust" {
        type RustRayTransport;
        type RustScatteringCoefficientInterpolator;
        type RustSolarTransmissionOperator;
        type RustSourceInterpolator;
        type RustSuccessiveOrdersSolver;
        type RustTransportSparsity;

        fn new_transport_sparsity(
            rows: usize,
            columns: usize,
            row_offsets: &[i32],
            column_indices: &[i32],
        ) -> Result<Box<RustTransportSparsity>>;

        fn transport_sparsity_storage_bytes(sparsity: &RustTransportSparsity) -> usize;

        #[allow(clippy::too_many_arguments)]
        fn new_ray_transport(
            ray_layer_offsets: &[u32],
            layer_atmosphere_offsets: &[u32],
            atmosphere_indices: &[u32],
            optical_depth_weights: &[f64],
            albedo_weights: &[f64],
            entrance_weights: &[f64],
            exit_weights: &[f64],
            layer_distance: &[f64],
            layer_start_fraction: &[f64],
            layer_end_fraction: &[f64],
            ray_scattering_cosine: &[f64],
            ray_phase_q_factor: &[f64],
            ray_phase_u_factor: &[f64],
            num_phase_moments: usize,
            solar_transmission_on_atmosphere_grid: bool,
            layer_source_offsets: &[u32],
            ray_transport_value_offsets: &[u32],
            ray_transport_row_nnz: &[u32],
            num_stokes: usize,
            source_value_inner_indices: &[u16],
            source_weights: &[f64],
            ray_ground_offsets: &[u32],
            ground_value_inner_indices: &[u16],
            ground_weights: &[f64],
        ) -> Result<Box<RustRayTransport>>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_ray_transport(
            transport: &RustRayTransport,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            transport_values: &mut [f64],
            layer_optical_depth: &mut [f64],
            layer_attenuation: &mut [f64],
            layer_prefix_attenuation: &mut [f64],
            ray_end_attenuation: &mut [f64],
        ) -> Result<()>;

        fn assemble_solver_ray_transport(
            transport: &RustRayTransport,
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            retain_layer_scratch: bool,
        ) -> Result<()>;

        fn prepare_solver_ray_transport_vjp_attenuation(
            transport: &RustRayTransport,
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_solver_ray_transport_with_first_order(
            transport: &RustRayTransport,
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            end_of_ray_source: &[f64],
            retain_layer_scratch: bool,
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_solver_ray_transport_with_first_order_batch4(
            transport: &RustRayTransport,
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            end_of_ray_source: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_solver_ray_transport_with_first_order_jvp(
            transport: &RustRayTransport,
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            extinction_tangent: &[f64],
            single_scatter_albedo_tangent: &[f64],
            legendre_coefficient_tangent: &[f64],
            solar_transmission_tangent: &[f64],
            end_of_ray_source: &[f64],
            end_of_ray_source_tangent: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn stage_solver_ray_transport_jvp_batch4_lane(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            lane: usize,
            zero_tangent: bool,
            extinction_tangent: &[f64],
            single_scatter_albedo_tangent: &[f64],
            legendre_coefficient_tangent: &[f64],
            end_of_ray_source: &[f64],
            end_of_ray_source_tangent: &[f64],
            first_order_forcing: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_solver_ray_transport_with_first_order_jvp_batch4(
            transport: &RustRayTransport,
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            solar_transmission_tangent: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_solver_ray_transport_with_first_order_vjp(
            transport: &RustRayTransport,
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            end_of_ray_source: &[f64],
            extinction_gradient: &mut [f64],
            single_scatter_albedo_gradient: &mut [f64],
            legendre_coefficient_gradient: &mut [f64],
            solar_transmission_gradient: &mut [f64],
            end_of_ray_source_gradient: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_solver_ray_transport_with_first_order_vjp_batch4(
            transport: &RustRayTransport,
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            end_of_ray_source: &[f64],
            solar_transmission_gradient: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn select_solver_ray_transport_vjp_batch4_lane(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            lane: usize,
            extinction_gradient: &mut [f64],
            single_scatter_albedo_gradient: &mut [f64],
            legendre_coefficient_gradient: &mut [f64],
            end_of_ray_source_gradient: &mut [f64],
        ) -> Result<()>;

        fn ray_transport_num_rays(transport: &RustRayTransport) -> usize;
        fn ray_transport_num_layers(transport: &RustRayTransport) -> usize;
        fn ray_transport_num_atmosphere_weights(transport: &RustRayTransport) -> usize;
        fn ray_transport_num_source_weights(transport: &RustRayTransport) -> usize;
        fn ray_transport_num_ground_weights(transport: &RustRayTransport) -> usize;
        fn ray_transport_num_phase_geometries(transport: &RustRayTransport) -> usize;
        fn ray_transport_num_grouped_phase_geometries(
            transport: &RustRayTransport,
            num_atmosphere_points: usize,
        ) -> usize;
        fn ray_transport_storage_bytes(transport: &RustRayTransport) -> usize;

        fn new_scattering_coefficient_interpolator(
            point_offsets: &[u32],
            atmosphere_indices: &[u32],
            weights: &[f64],
            num_atmosphere_points: usize,
            num_output_coefficients: usize,
        ) -> Result<Box<RustScatteringCoefficientInterpolator>>;

        fn scattering_coefficient_interpolator_num_weights(
            interpolator: &RustScatteringCoefficientInterpolator,
        ) -> usize;
        fn scattering_coefficient_interpolator_storage_bytes(
            interpolator: &RustScatteringCoefficientInterpolator,
        ) -> usize;

        #[allow(clippy::too_many_arguments)]
        fn calculate_atmospheric_coefficient_jvp(
            atmosphere_coefficient_derivatives: &[f64],
            atmosphere_coefficient_stride: usize,
            num_atmosphere_points: usize,
            num_wavelengths: usize,
            wavelength_index: usize,
            derivative_group_tangent: &[f64],
            output: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn accumulate_atmospheric_coefficient_gradient(
            coefficient_gradient: &[f64],
            atmosphere_coefficient_derivatives: &[f64],
            atmosphere_coefficient_stride: usize,
            num_atmosphere_points: usize,
            num_wavelengths: usize,
            wavelength_index: usize,
            active_derivative_groups: &[i32],
            atmosphere_gradient: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn new_solar_transmission_operator(
            path_rows: usize,
            path_columns: usize,
            dense_path_values: &[f64],
            path_row_offsets: &[u32],
            path_column_indices: &[u32],
            sparse_path_values: &[f64],
            interpolation_rows: usize,
            interpolation_columns: usize,
            interpolation_row_offsets: &[u32],
            interpolation_column_indices: &[u32],
            interpolation_values: &[f64],
            ground_hit: &[u8],
        ) -> Result<Box<RustSolarTransmissionOperator>>;

        fn calculate_solar_transmission(
            solar_operator: &RustSolarTransmissionOperator,
            extinction: &[f64],
            solar_irradiance: f64,
            transmission: &mut [f64],
            scratch: &mut [f64],
        ) -> Result<()>;

        fn calculate_solar_transmission_batch4(
            solar_operator: &RustSolarTransmissionOperator,
            extinction: &[f64],
            solar_irradiance: &[f64],
            transmission: &mut [f64],
            scratch: &mut [f64],
        ) -> Result<()>;

        fn calculate_solar_transmission_jvp(
            solar_operator: &RustSolarTransmissionOperator,
            extinction_tangent: &[f64],
            transmission: &[f64],
            transmission_tangent: &mut [f64],
            scratch: &mut [f64],
        ) -> Result<()>;

        fn calculate_solar_transmission_jvp_batch4(
            solar_operator: &RustSolarTransmissionOperator,
            extinction_tangent: &[f64],
            transmission: &[f64],
            transmission_tangent: &mut [f64],
            scratch: &mut [f64],
        ) -> Result<()>;

        fn accumulate_solar_transmission_vjp(
            solar_operator: &RustSolarTransmissionOperator,
            transmission: &[f64],
            transmission_cotangent: &[f64],
            extinction_gradient: &mut [f64],
            scratch: &mut [f64],
        ) -> Result<()>;

        fn accumulate_solar_transmission_vjp_batch4(
            solar_operator: &RustSolarTransmissionOperator,
            transmission: &[f64],
            transmission_cotangent: &mut [f64],
            extinction_gradient: &mut [f64],
            scratch: &mut [f64],
        ) -> Result<()>;

        fn solar_transmission_input_size(solar_operator: &RustSolarTransmissionOperator) -> usize;
        fn solar_transmission_output_size(solar_operator: &RustSolarTransmissionOperator) -> usize;
        fn solar_transmission_scratch_size(solar_operator: &RustSolarTransmissionOperator)
        -> usize;
        fn solar_transmission_forward_scratch_size(
            solar_operator: &RustSolarTransmissionOperator,
        ) -> usize;
        fn solar_transmission_storage_bytes(
            solar_operator: &RustSolarTransmissionOperator,
        ) -> usize;

        #[allow(clippy::too_many_arguments)]
        fn new_source_interpolator(
            ray_layer_offsets: &[u32],
            layer_atmosphere_offsets: &[u32],
            atmosphere_indices: &[u32],
            atmosphere_weights: &[f64],
            optical_depth_weights: &[f64],
            layer_source_offsets: &[u32],
            source_indices: &[u32],
            source_weights: &[f64],
            ray_ground_offsets: &[u32],
            ground_source_indices: &[u32],
            ground_source_weights: &[f64],
            num_atmosphere_points: usize,
            num_solution_points: usize,
            num_stokes: usize,
        ) -> Result<Box<RustSourceInterpolator>>;

        #[allow(clippy::too_many_arguments)]
        fn accumulate_source_ray_value(
            interpolator: &RustSourceInterpolator,
            ray_index: usize,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            wavelength_count: usize,
            solution: &[f64],
            solution_stride: usize,
            source: &mut [f64],
            source_stride: usize,
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn accumulate_source_ray_jacobian(
            interpolator: &RustSourceInterpolator,
            ray_index: usize,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            wavelength_count: usize,
            solution: &[f64],
            solution_stride: usize,
            source: &mut [f64],
            source_stride: usize,
            source_derivative: &mut [f64],
            num_derivatives: usize,
            derivative_stride: usize,
            single_scatter_albedo_derivative_start: usize,
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn accumulate_source_ray_jvp(
            interpolator: &RustSourceInterpolator,
            ray_index: usize,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            extinction_tangent: &[f64],
            single_scatter_albedo_tangent: &[f64],
            solution: &[f64],
            solution_tangent: &[f64],
            source: &mut [f64],
            source_tangent: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn accumulate_source_ray_vjp(
            interpolator: &RustSourceInterpolator,
            ray_index: usize,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            solution: &[f64],
            cotangent: &[f64],
            solution_cotangent: &mut [f64],
            single_scatter_albedo_gradient: &mut [f64],
            extinction_gradient: &mut [f64],
        ) -> Result<()>;

        fn source_interpolator_num_rays(interpolator: &RustSourceInterpolator) -> usize;
        fn source_interpolator_num_layers(interpolator: &RustSourceInterpolator) -> usize;
        fn source_interpolator_storage_bytes(interpolator: &RustSourceInterpolator) -> usize;

        fn set_scattering_coefficients_from_atmosphere(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            interpolator: &RustScatteringCoefficientInterpolator,
            atmosphere_coefficients: &[f64],
            atmosphere_coefficient_stride: usize,
        ) -> Result<()>;

        fn set_scattering_coefficients_from_atmosphere_batch4(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            interpolator: &RustScatteringCoefficientInterpolator,
            atmosphere_coefficients: &[f64],
            atmosphere_coefficient_stride: usize,
        ) -> Result<()>;

        fn boundary_scattering_values_mut(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
        ) -> &mut [f64];

        fn batch4_boundary_scattering_values_mut(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
        ) -> &mut [f64];

        fn boundary_scattering_value_tangent_mut(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
        ) -> &mut [f64];

        #[allow(clippy::too_many_arguments)]
        fn set_scattering_coefficient_jvp_from_atmosphere(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            interpolator: &RustScatteringCoefficientInterpolator,
            atmosphere_coefficient_derivatives: &[f64],
            atmosphere_coefficient_stride: usize,
            num_wavelengths: usize,
            wavelength_index: usize,
            derivative_group_tangent: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn accumulate_solver_scattering_coefficient_vjp(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            interpolator: &RustScatteringCoefficientInterpolator,
            atmosphere_coefficient_derivatives: &[f64],
            atmosphere_coefficient_stride: usize,
            num_wavelengths: usize,
            wavelength_index: usize,
            active_derivative_groups: &[i32],
            atmosphere_gradient: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn new_successive_orders_solver(
            transport_rows: usize,
            transport_columns: usize,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            scattering_output_offsets: &[usize],
            scattering_input_offsets: &[usize],
            scattering_value_offsets: &[usize],
            max_iterations: usize,
            relative_tolerance: f64,
            absolute_tolerance: f64,
            anderson_depth: usize,
            damping: f64,
        ) -> Result<Box<RustSuccessiveOrdersSolver>>;

        #[allow(clippy::too_many_arguments)]
        fn new_scalar_coefficient_successive_orders_solver(
            transport_rows: usize,
            transport_columns: usize,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            scattering_output_offsets: &[usize],
            scattering_input_offsets: &[usize],
            coefficient_blocks: usize,
            num_coefficients: usize,
            incoming_directions: &[f64],
            incoming_weights: &[f64],
            outgoing_directions: &[f64],
            max_iterations: usize,
            relative_tolerance: f64,
            absolute_tolerance: f64,
            anderson_depth: usize,
            damping: f64,
        ) -> Result<Box<RustSuccessiveOrdersSolver>>;

        #[allow(clippy::too_many_arguments)]
        fn new_vector_coefficient_successive_orders_solver(
            transport_rows: usize,
            transport_columns: usize,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            scattering_output_offsets: &[usize],
            scattering_input_offsets: &[usize],
            coefficient_blocks: usize,
            num_coefficients: usize,
            incoming_directions: &[f64],
            incoming_weights: &[f64],
            outgoing_directions: &[f64],
            max_iterations: usize,
            relative_tolerance: f64,
            absolute_tolerance: f64,
            anderson_depth: usize,
            damping: f64,
        ) -> Result<Box<RustSuccessiveOrdersSolver>>;

        #[allow(clippy::too_many_arguments)]
        fn new_scalar_coefficient_successive_orders_solver_with_sparsity(
            transport_sparsity: &RustTransportSparsity,
            scattering_output_offsets: &[usize],
            scattering_input_offsets: &[usize],
            coefficient_blocks: usize,
            num_coefficients: usize,
            incoming_directions: &[f64],
            incoming_weights: &[f64],
            outgoing_directions: &[f64],
            max_iterations: usize,
            relative_tolerance: f64,
            absolute_tolerance: f64,
            anderson_depth: usize,
            damping: f64,
        ) -> Result<Box<RustSuccessiveOrdersSolver>>;

        #[allow(clippy::too_many_arguments)]
        fn new_vector_coefficient_successive_orders_solver_with_sparsity(
            transport_sparsity: &RustTransportSparsity,
            scattering_output_offsets: &[usize],
            scattering_input_offsets: &[usize],
            coefficient_blocks: usize,
            num_coefficients: usize,
            incoming_directions: &[f64],
            incoming_weights: &[f64],
            outgoing_directions: &[f64],
            max_iterations: usize,
            relative_tolerance: f64,
            absolute_tolerance: f64,
            anderson_depth: usize,
            damping: f64,
        ) -> Result<Box<RustSuccessiveOrdersSolver>>;

        fn solve(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            transport_values: &[f64],
            scattering_values: &[f64],
            first_order_forcing: &[f64],
            initial_guess: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn solve_scalar_coefficients(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            transport_values: &[f64],
            boundary_scattering_values: &[f64],
            first_order_forcing: &[f64],
            initial_guess: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn solve_vector_coefficients(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            transport_values: &[f64],
            boundary_scattering_values: &[f64],
            first_order_forcing: &[f64],
            initial_guess: &[f64],
        ) -> Result<()>;

        fn solve_coefficients_assembled(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            initial_guess: &[f64],
        ) -> Result<()>;

        fn begin_coefficients_batch4(solver: Pin<&mut RustSuccessiveOrdersSolver>) -> Result<()>;

        fn stage_coefficients_batch4_lane(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            lane: usize,
            initial_guess: &[f64],
        ) -> Result<()>;

        fn solve_coefficients_batch4(solver: Pin<&mut RustSuccessiveOrdersSolver>) -> Result<()>;
        fn restore_coefficients_batch4_solution(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            solution: &[f64],
        ) -> Result<()>;
        fn select_coefficients_batch4_lane(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            lane: usize,
        ) -> Result<()>;

        fn batch4_solution(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn batch4_first_order_forcing(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn batch4_converged(solver: &RustSuccessiveOrdersSolver, lane: usize) -> bool;
        fn batch4_iterations(solver: &RustSuccessiveOrdersSolver, lane: usize) -> usize;
        fn batch4_initial_residual(solver: &RustSuccessiveOrdersSolver, lane: usize) -> f64;
        fn batch4_final_residual(solver: &RustSuccessiveOrdersSolver, lane: usize) -> f64;
        fn batch4_workspace_bytes(solver: &RustSuccessiveOrdersSolver) -> usize;

        fn begin_linearize_coefficients_jvp_batch4(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
        ) -> Result<()>;
        fn stage_linearize_coefficients_jvp_batch4_lane(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            lane: usize,
            zero_tangent: bool,
            first_order_forcing: &[f64],
        ) -> Result<()>;
        fn linearize_coefficients_jvp_batch4(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
        ) -> Result<()>;
        fn select_linearize_coefficients_jvp_batch4_lane(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            lane: usize,
        ) -> Result<()>;

        fn begin_linearize_coefficients_vjp_batch4(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
        ) -> Result<()>;
        fn stage_linearize_coefficients_vjp_batch4_lane(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            lane: usize,
            solution_cotangent: &[f64],
        ) -> Result<()>;
        fn linearize_coefficients_vjp_batch4(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
        ) -> Result<()>;
        fn select_linearize_coefficients_vjp_batch4_lane(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            lane: usize,
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn linearize_coefficients_jvp(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            transport_values: &[f64],
            transport_value_tangent: &[f64],
            boundary_scattering_value_tangent: &[f64],
            first_order_forcing: &[f64],
            first_order_forcing_tangent: &[f64],
        ) -> Result<()>;

        fn linearize_coefficients_jvp_assembled(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn linearize_coefficients_vjp(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            transport_values: &[f64],
            first_order_forcing: &[f64],
            solution_cotangent: &[f64],
        ) -> Result<()>;

        fn linearize_coefficients_vjp_assembled(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            solution_cotangent: &[f64],
        ) -> Result<()>;

        fn restore_coefficient_solution(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            boundary_scattering_values: &[f64],
            solution: &[f64],
        ) -> Result<()>;

        fn restore_coefficient_solution_assembled(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            first_order_forcing: &[f64],
            solution: &[f64],
        ) -> Result<()>;

        fn restore_solution_only(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            solution: &[f64],
        ) -> Result<()>;

        fn solution(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn first_order_forcing(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn solution_jvp(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn clear_solution_jvp(solver: Pin<&mut RustSuccessiveOrdersSolver>);
        fn transport_value_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn boundary_scattering_value_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn first_order_forcing_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn residual_history(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn iterations(solver: &RustSuccessiveOrdersSolver) -> usize;
        fn converged(solver: &RustSuccessiveOrdersSolver) -> bool;
        fn initial_residual(solver: &RustSuccessiveOrdersSolver) -> f64;
        fn final_residual(solver: &RustSuccessiveOrdersSolver) -> f64;
        fn transport_workspace_bytes(solver: &RustSuccessiveOrdersSolver) -> usize;
    }
}

pub struct RustRayTransport {
    transport: ScalarRayTransport,
}

pub struct RustTransportSparsity {
    rows: usize,
    columns: usize,
    row_offsets: Arc<[i32]>,
    column_indices: Arc<[i32]>,
}

type SharedTransportSparsity = (Arc<[i32]>, Arc<[i32]>);

fn new_transport_sparsity(
    rows: usize,
    columns: usize,
    row_offsets: &[i32],
    column_indices: &[i32],
) -> Result<Box<RustTransportSparsity>> {
    CsrMatrix::new(rows, columns, row_offsets, column_indices)?;
    Ok(Box::new(RustTransportSparsity {
        rows,
        columns,
        row_offsets: Arc::from(row_offsets),
        column_indices: Arc::from(column_indices),
    }))
}

fn transport_sparsity_storage_bytes(sparsity: &RustTransportSparsity) -> usize {
    sparsity.row_offsets.len() * std::mem::size_of::<i32>()
        + sparsity.column_indices.len() * std::mem::size_of::<i32>()
}

pub struct RustScatteringCoefficientInterpolator {
    interpolator: ScatteringCoefficientInterpolator,
}

pub struct RustSolarTransmissionOperator {
    operator: SolarTransmissionOperator,
}

pub struct RustSourceInterpolator {
    interpolator: SourceInterpolator,
}

#[allow(clippy::too_many_arguments)]
fn new_ray_transport(
    ray_layer_offsets: &[u32],
    layer_atmosphere_offsets: &[u32],
    atmosphere_indices: &[u32],
    optical_depth_weights: &[f64],
    albedo_weights: &[f64],
    entrance_weights: &[f64],
    exit_weights: &[f64],
    layer_distance: &[f64],
    layer_start_fraction: &[f64],
    layer_end_fraction: &[f64],
    ray_scattering_cosine: &[f64],
    ray_phase_q_factor: &[f64],
    ray_phase_u_factor: &[f64],
    num_phase_moments: usize,
    solar_transmission_on_atmosphere_grid: bool,
    layer_source_offsets: &[u32],
    ray_transport_value_offsets: &[u32],
    ray_transport_row_nnz: &[u32],
    num_stokes: usize,
    source_value_inner_indices: &[u16],
    source_weights: &[f64],
    ray_ground_offsets: &[u32],
    ground_value_inner_indices: &[u16],
    ground_weights: &[f64],
) -> Result<Box<RustRayTransport>> {
    let transport = ScalarRayTransport::new(
        ray_layer_offsets,
        layer_atmosphere_offsets,
        atmosphere_indices,
        optical_depth_weights,
        albedo_weights,
        entrance_weights,
        exit_weights,
        layer_distance,
        layer_start_fraction,
        layer_end_fraction,
        ray_scattering_cosine,
        ray_phase_q_factor,
        ray_phase_u_factor,
        num_phase_moments,
        solar_transmission_on_atmosphere_grid,
        layer_source_offsets,
        ray_transport_value_offsets,
        source_value_inner_indices,
        source_weights,
        ray_ground_offsets,
        ground_value_inner_indices,
        ground_weights,
    )?
    .with_stokes_layout(num_stokes, ray_transport_row_nnz)?;
    Ok(Box::new(RustRayTransport { transport }))
}

fn assemble_solver_ray_transport(
    transport: &RustRayTransport,
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    retain_layer_scratch: bool,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.prepare_transport_workspace(&transport.transport, retain_layer_scratch);
    transport.transport.assemble(
        extinction,
        single_scatter_albedo,
        &mut this.transport_values,
        &mut this.layer_optical_depth,
        &mut this.layer_attenuation,
        &mut this.layer_prefix_attenuation,
        &mut this.ray_end_attenuation,
    )
}

fn prepare_solver_ray_transport_vjp_attenuation(
    transport: &RustRayTransport,
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.layer_optical_depth
        .resize(transport.transport.num_layers(), 0.0);
    this.layer_attenuation
        .resize(transport.transport.num_layers(), 0.0);
    this.layer_prefix_attenuation
        .resize(transport.transport.num_layers(), 0.0);
    this.ray_end_attenuation
        .resize(transport.transport.num_rays(), 0.0);
    transport.transport.assemble_vjp_attenuation(
        extinction,
        single_scatter_albedo,
        &mut this.layer_optical_depth,
        &mut this.layer_attenuation,
        &mut this.layer_prefix_attenuation,
        &mut this.ray_end_attenuation,
    )
}

#[allow(clippy::too_many_arguments)]
fn assemble_solver_ray_transport_with_first_order(
    transport: &RustRayTransport,
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    end_of_ray_source: &[f64],
    retain_layer_scratch: bool,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.prepare_transport_workspace(&transport.transport, retain_layer_scratch);
    this.first_order_forcing.resize(
        transport.transport.num_rays() * transport.transport.num_stokes(),
        0.0,
    );
    transport.transport.assemble_with_first_order(
        extinction,
        single_scatter_albedo,
        legendre_coefficients,
        maximum_order,
        solar_transmission,
        end_of_ray_source,
        &mut this.transport_values,
        &mut this.first_order_forcing,
        &mut this.layer_optical_depth,
        &mut this.layer_attenuation,
        &mut this.layer_prefix_attenuation,
        &mut this.ray_end_attenuation,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn assemble_solver_ray_transport_with_first_order_batch4(
    transport: &RustRayTransport,
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    end_of_ray_source: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.transport_row_offsets.is_none() || this.transport_column_indices.is_none() {
        return Err(anyhow!(
            "wavelength SIMD requires Rust-owned transport sparsity"
        ));
    }
    this.simd
        .batch4_transport_values
        .resize(transport.transport.transport_value_size() * LANES, 0.0);
    this.simd
        .batch4_first_order_forcing
        .resize(transport.transport.num_rays() * LANES, 0.0);
    this.simd
        .batch4_solution
        .resize(this.solver.solution().len() * LANES, 0.0);
    this.simd.batch4_initial_guess.clear();
    this.simd.batch4_has_initial_guess = false;
    this.simd.batch4_staged_lanes = LANES;
    this.simd.batch4_diagnostics = std::array::from_fn(|_| SolverDiagnostics::default());
    transport
        .transport
        .assemble_batch4_with_first_order_scalar(
            extinction,
            single_scatter_albedo,
            legendre_coefficients,
            maximum_order,
            solar_transmission,
            end_of_ray_source,
            &mut this.simd.batch4_transport_values,
            &mut this.simd.batch4_first_order_forcing,
        )
        .map_err(|error| anyhow!("assembling four-wavelength ray transport: {error}"))?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn assemble_solver_ray_transport_with_first_order_jvp(
    transport: &RustRayTransport,
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    extinction_tangent: &[f64],
    single_scatter_albedo_tangent: &[f64],
    legendre_coefficient_tangent: &[f64],
    solar_transmission_tangent: &[f64],
    end_of_ray_source: &[f64],
    end_of_ray_source_tangent: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.prepare_transport_workspace(&transport.transport, false);
    this.transport_value_tangent
        .resize(transport.transport.transport_value_size(), 0.0);
    this.ray_end_attenuation_tangent
        .resize(transport.transport.num_rays(), 0.0);
    this.first_order_forcing_tangent.resize(
        transport.transport.num_rays() * transport.transport.num_stokes(),
        0.0,
    );
    transport.transport.assemble_jvp_with_first_order(
        extinction,
        single_scatter_albedo,
        legendre_coefficients,
        maximum_order,
        solar_transmission,
        extinction_tangent,
        single_scatter_albedo_tangent,
        legendre_coefficient_tangent,
        solar_transmission_tangent,
        end_of_ray_source,
        end_of_ray_source_tangent,
        &mut this.transport_values,
        &mut this.transport_value_tangent,
        &mut this.first_order_forcing_tangent,
        &mut this.ray_end_attenuation,
        &mut this.ray_end_attenuation_tangent,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn stage_solver_ray_transport_jvp_batch4_lane(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    lane: usize,
    zero_tangent: bool,
    extinction_tangent: &[f64],
    single_scatter_albedo_tangent: &[f64],
    legendre_coefficient_tangent: &[f64],
    end_of_ray_source: &[f64],
    end_of_ray_source_tangent: &[f64],
    first_order_forcing: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if lane >= LANES || lane != this.simd.batch4_jvp_staged_lanes {
        return Err(anyhow!(
            "four-wavelength transport JVP lanes must be staged in order"
        ));
    }
    if extinction_tangent.len() != single_scatter_albedo_tangent.len()
        || end_of_ray_source_tangent.len() != end_of_ray_source.len()
        || first_order_forcing.len() * LANES != this.simd.batch4_first_order_forcing.len()
    {
        return Err(anyhow!(
            "four-wavelength transport JVP lane arrays have inconsistent dimensions"
        ));
    }
    if lane == 0 {
        this.simd
            .batch4_extinction_tangent
            .resize(extinction_tangent.len() * LANES, 0.0);
        this.simd
            .batch4_single_scatter_albedo_tangent
            .resize(single_scatter_albedo_tangent.len() * LANES, 0.0);
        this.simd
            .batch4_legendre_coefficient_tangent
            .resize(legendre_coefficient_tangent.len() * LANES, 0.0);
        this.simd
            .batch4_end_of_ray_source
            .resize(end_of_ray_source.len() * LANES, 0.0);
        this.simd
            .batch4_end_of_ray_source_tangent
            .resize(end_of_ray_source_tangent.len() * LANES, 0.0);
    }
    if this.simd.batch4_extinction_tangent.len() != extinction_tangent.len() * LANES
        || this.simd.batch4_single_scatter_albedo_tangent.len()
            != single_scatter_albedo_tangent.len() * LANES
        || this.simd.batch4_legendre_coefficient_tangent.len()
            != legendre_coefficient_tangent.len() * LANES
        || this.simd.batch4_end_of_ray_source.len() != end_of_ray_source.len() * LANES
        || this.simd.batch4_end_of_ray_source_tangent.len()
            != end_of_ray_source_tangent.len() * LANES
    {
        return Err(anyhow!(
            "four-wavelength transport JVP lane sizes changed within a block"
        ));
    }

    interleave_lane(
        extinction_tangent,
        lane,
        &mut this.simd.batch4_extinction_tangent,
    );
    interleave_lane(
        single_scatter_albedo_tangent,
        lane,
        &mut this.simd.batch4_single_scatter_albedo_tangent,
    );
    interleave_lane(
        legendre_coefficient_tangent,
        lane,
        &mut this.simd.batch4_legendre_coefficient_tangent,
    );
    interleave_lane(
        end_of_ray_source,
        lane,
        &mut this.simd.batch4_end_of_ray_source,
    );
    interleave_lane(
        end_of_ray_source_tangent,
        lane,
        &mut this.simd.batch4_end_of_ray_source_tangent,
    );
    interleave_lane(
        first_order_forcing,
        lane,
        &mut this.simd.batch4_first_order_forcing,
    );
    if zero_tangent {
        for target in [
            &mut this.simd.batch4_scattering_coefficient_tangent,
            &mut this.simd.batch4_boundary_scattering_value_tangent,
        ] {
            for element in 0..target.len() / LANES {
                target[element * LANES + lane] = 0.0;
            }
        }
    } else {
        if this.scattering_coefficient_tangent.len() * LANES
            != this.simd.batch4_scattering_coefficient_tangent.len()
            || this.boundary_scattering_value_tangent.len() * LANES
                != this.simd.batch4_boundary_scattering_value_tangent.len()
        {
            return Err(anyhow!(
                "four-wavelength scattering JVP storage has inconsistent dimensions"
            ));
        }
        interleave_lane(
            &this.scattering_coefficient_tangent,
            lane,
            &mut this.simd.batch4_scattering_coefficient_tangent,
        );
        interleave_lane(
            &this.boundary_scattering_value_tangent,
            lane,
            &mut this.simd.batch4_boundary_scattering_value_tangent,
        );
    }
    this.simd.batch4_jvp_staged_lanes += 1;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn assemble_solver_ray_transport_with_first_order_jvp_batch4(
    transport: &RustRayTransport,
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    solar_transmission_tangent: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.simd.batch4_jvp_staged_lanes != LANES {
        return Err(anyhow!("four-wavelength transport JVP batch is incomplete"));
    }
    transport
        .transport
        .assemble_batch4_jvp_with_first_order_scalar(
            extinction,
            single_scatter_albedo,
            legendre_coefficients,
            maximum_order,
            solar_transmission,
            &this.simd.batch4_extinction_tangent,
            &this.simd.batch4_single_scatter_albedo_tangent,
            &this.simd.batch4_legendre_coefficient_tangent,
            solar_transmission_tangent,
            &this.simd.batch4_end_of_ray_source,
            &this.simd.batch4_end_of_ray_source_tangent,
            &mut this.simd.batch4_transport_values,
            &mut this.simd.batch4_transport_value_tangent,
            &mut this.simd.batch4_first_order_forcing_tangent,
        )
        .map_err(|error| anyhow!("assembling four-wavelength transport JVP: {error}"))?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn assemble_solver_ray_transport_with_first_order_vjp(
    transport: &RustRayTransport,
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    end_of_ray_source: &[f64],
    extinction_gradient: &mut [f64],
    single_scatter_albedo_gradient: &mut [f64],
    legendre_coefficient_gradient: &mut [f64],
    solar_transmission_gradient: &mut [f64],
    end_of_ray_source_gradient: &mut [f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let (_, transport_column_indices) = this.transport_sparsity()?;
    this.ray_end_attenuation_cotangent
        .resize(transport.transport.num_rays(), 0.0);
    this.ray_end_attenuation_cotangent.fill(0.0);
    transport.transport.assemble_vjp_with_first_order(
        extinction,
        single_scatter_albedo,
        legendre_coefficients,
        maximum_order,
        solar_transmission,
        &this.transport_value_gradient,
        &transport_column_indices,
        this.solver.solution(),
        &this.first_order_forcing_gradient,
        &this.ray_end_attenuation,
        &this.ray_end_attenuation_cotangent,
        end_of_ray_source,
        &this.layer_optical_depth,
        &this.layer_attenuation,
        &this.layer_prefix_attenuation,
        extinction_gradient,
        single_scatter_albedo_gradient,
        legendre_coefficient_gradient,
        solar_transmission_gradient,
        end_of_ray_source_gradient,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn assemble_solver_ray_transport_with_first_order_vjp_batch4(
    transport: &RustRayTransport,
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    end_of_ray_source: &[f64],
    solar_transmission_gradient: &mut [f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.simd.batch4_vjp_staged_lanes != LANES {
        return Err(anyhow!(
            "four-wavelength transport VJP requires a completed packed implicit VJP"
        ));
    }
    let transport_column_indices = this
        .transport_column_indices
        .as_ref()
        .ok_or_else(|| anyhow!("four-wavelength transport VJP requires Rust-owned sparsity"))?;
    this.simd
        .batch4_extinction_gradient
        .resize(extinction.len(), 0.0);
    this.simd
        .batch4_single_scatter_albedo_gradient
        .resize(single_scatter_albedo.len(), 0.0);
    this.simd
        .batch4_legendre_coefficient_gradient
        .resize(legendre_coefficients.len(), 0.0);
    this.simd
        .batch4_end_of_ray_source_gradient
        .resize(end_of_ray_source.len(), 0.0);
    transport
        .transport
        .assemble_batch4_vjp_with_first_order_scalar(
            extinction,
            single_scatter_albedo,
            legendre_coefficients,
            maximum_order,
            solar_transmission,
            transport_column_indices,
            &this.simd.batch4_solution,
            &this.simd.batch4_first_order_forcing_gradient,
            end_of_ray_source,
            &mut this.simd.batch4_extinction_gradient,
            &mut this.simd.batch4_single_scatter_albedo_gradient,
            &mut this.simd.batch4_legendre_coefficient_gradient,
            solar_transmission_gradient,
            &mut this.simd.batch4_end_of_ray_source_gradient,
        )
        .map_err(|error| anyhow!("assembling four-wavelength transport VJP: {error}"))?;
    this.simd.batch4_ray_vjp_ready = true;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn select_solver_ray_transport_vjp_batch4_lane(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    lane: usize,
    extinction_gradient: &mut [f64],
    single_scatter_albedo_gradient: &mut [f64],
    legendre_coefficient_gradient: &mut [f64],
    end_of_ray_source_gradient: &mut [f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if lane >= LANES
        || !this.simd.batch4_ray_vjp_ready
        || this.simd.batch4_extinction_gradient.len() != extinction_gradient.len() * LANES
        || this.simd.batch4_single_scatter_albedo_gradient.len()
            != single_scatter_albedo_gradient.len() * LANES
        || this.simd.batch4_legendre_coefficient_gradient.len()
            != legendre_coefficient_gradient.len() * LANES
        || this.simd.batch4_end_of_ray_source_gradient.len()
            != end_of_ray_source_gradient.len() * LANES
    {
        return Err(anyhow!("invalid four-wavelength transport VJP lane"));
    }
    extract_lane(
        &this.simd.batch4_extinction_gradient,
        lane,
        extinction_gradient,
    );
    extract_lane(
        &this.simd.batch4_single_scatter_albedo_gradient,
        lane,
        single_scatter_albedo_gradient,
    );
    extract_lane(
        &this.simd.batch4_legendre_coefficient_gradient,
        lane,
        legendre_coefficient_gradient,
    );
    extract_lane(
        &this.simd.batch4_end_of_ray_source_gradient,
        lane,
        end_of_ray_source_gradient,
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn assemble_ray_transport(
    transport: &RustRayTransport,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    transport_values: &mut [f64],
    layer_optical_depth: &mut [f64],
    layer_attenuation: &mut [f64],
    layer_prefix_attenuation: &mut [f64],
    ray_end_attenuation: &mut [f64],
) -> Result<()> {
    transport.transport.assemble(
        extinction,
        single_scatter_albedo,
        transport_values,
        layer_optical_depth,
        layer_attenuation,
        layer_prefix_attenuation,
        ray_end_attenuation,
    )?;
    Ok(())
}

fn ray_transport_num_rays(transport: &RustRayTransport) -> usize {
    transport.transport.num_rays()
}

fn ray_transport_num_layers(transport: &RustRayTransport) -> usize {
    transport.transport.num_layers()
}

fn ray_transport_num_atmosphere_weights(transport: &RustRayTransport) -> usize {
    transport.transport.num_atmosphere_weights()
}

fn ray_transport_num_source_weights(transport: &RustRayTransport) -> usize {
    transport.transport.num_source_weights()
}

fn ray_transport_num_ground_weights(transport: &RustRayTransport) -> usize {
    transport.transport.num_ground_weights()
}

fn ray_transport_num_phase_geometries(transport: &RustRayTransport) -> usize {
    transport.transport.num_phase_geometries()
}

fn ray_transport_num_grouped_phase_geometries(
    transport: &RustRayTransport,
    num_atmosphere_points: usize,
) -> usize {
    transport
        .transport
        .num_grouped_phase_geometries(num_atmosphere_points)
}

fn ray_transport_storage_bytes(transport: &RustRayTransport) -> usize {
    transport.transport.storage_bytes()
}

fn new_scattering_coefficient_interpolator(
    point_offsets: &[u32],
    atmosphere_indices: &[u32],
    weights: &[f64],
    num_atmosphere_points: usize,
    num_output_coefficients: usize,
) -> Result<Box<RustScatteringCoefficientInterpolator>> {
    Ok(Box::new(RustScatteringCoefficientInterpolator {
        interpolator: ScatteringCoefficientInterpolator::new(
            point_offsets,
            atmosphere_indices,
            weights,
            num_atmosphere_points,
            num_output_coefficients,
        )?,
    }))
}

fn scattering_coefficient_interpolator_num_weights(
    interpolator: &RustScatteringCoefficientInterpolator,
) -> usize {
    interpolator.interpolator.num_weights()
}

fn scattering_coefficient_interpolator_storage_bytes(
    interpolator: &RustScatteringCoefficientInterpolator,
) -> usize {
    interpolator.interpolator.storage_bytes()
}

#[allow(clippy::too_many_arguments)]
fn calculate_atmospheric_coefficient_jvp(
    atmosphere_coefficient_derivatives: &[f64],
    atmosphere_coefficient_stride: usize,
    num_atmosphere_points: usize,
    num_wavelengths: usize,
    wavelength_index: usize,
    derivative_group_tangent: &[f64],
    output: &mut [f64],
) -> Result<()> {
    atmospheric_coefficient_jvp(
        atmosphere_coefficient_derivatives,
        atmosphere_coefficient_stride,
        num_atmosphere_points,
        num_wavelengths,
        wavelength_index,
        derivative_group_tangent,
        output,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn accumulate_atmospheric_coefficient_gradient(
    coefficient_gradient: &[f64],
    atmosphere_coefficient_derivatives: &[f64],
    atmosphere_coefficient_stride: usize,
    num_atmosphere_points: usize,
    num_wavelengths: usize,
    wavelength_index: usize,
    active_derivative_groups: &[i32],
    atmosphere_gradient: &mut [f64],
) -> Result<()> {
    accumulate_atmospheric_coefficient_vjp(
        coefficient_gradient,
        atmosphere_coefficient_derivatives,
        atmosphere_coefficient_stride,
        num_atmosphere_points,
        num_wavelengths,
        wavelength_index,
        active_derivative_groups,
        atmosphere_gradient,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn new_solar_transmission_operator(
    path_rows: usize,
    path_columns: usize,
    dense_path_values: &[f64],
    path_row_offsets: &[u32],
    path_column_indices: &[u32],
    sparse_path_values: &[f64],
    interpolation_rows: usize,
    interpolation_columns: usize,
    interpolation_row_offsets: &[u32],
    interpolation_column_indices: &[u32],
    interpolation_values: &[f64],
    ground_hit: &[u8],
) -> Result<Box<RustSolarTransmissionOperator>> {
    Ok(Box::new(RustSolarTransmissionOperator {
        operator: SolarTransmissionOperator::new(
            path_rows,
            path_columns,
            dense_path_values,
            path_row_offsets,
            path_column_indices,
            sparse_path_values,
            interpolation_rows,
            interpolation_columns,
            interpolation_row_offsets,
            interpolation_column_indices,
            interpolation_values,
            ground_hit,
        )?,
    }))
}

fn calculate_solar_transmission(
    solar_operator: &RustSolarTransmissionOperator,
    extinction: &[f64],
    solar_irradiance: f64,
    transmission: &mut [f64],
    scratch: &mut [f64],
) -> Result<()> {
    solar_operator
        .operator
        .calculate(extinction, solar_irradiance, transmission, scratch)?;
    Ok(())
}

fn calculate_solar_transmission_batch4(
    solar_operator: &RustSolarTransmissionOperator,
    extinction: &[f64],
    solar_irradiance: &[f64],
    transmission: &mut [f64],
    scratch: &mut [f64],
) -> Result<()> {
    solar_operator.operator.calculate_batch4(
        extinction,
        solar_irradiance,
        transmission,
        scratch,
    )?;
    Ok(())
}

fn calculate_solar_transmission_jvp(
    solar_operator: &RustSolarTransmissionOperator,
    extinction_tangent: &[f64],
    transmission: &[f64],
    transmission_tangent: &mut [f64],
    scratch: &mut [f64],
) -> Result<()> {
    solar_operator.operator.calculate_jvp(
        extinction_tangent,
        transmission,
        transmission_tangent,
        scratch,
    )?;
    Ok(())
}

fn calculate_solar_transmission_jvp_batch4(
    solar_operator: &RustSolarTransmissionOperator,
    extinction_tangent: &[f64],
    transmission: &[f64],
    transmission_tangent: &mut [f64],
    scratch: &mut [f64],
) -> Result<()> {
    solar_operator.operator.calculate_jvp_batch4(
        extinction_tangent,
        transmission,
        transmission_tangent,
        scratch,
    )?;
    Ok(())
}

fn accumulate_solar_transmission_vjp(
    solar_operator: &RustSolarTransmissionOperator,
    transmission: &[f64],
    transmission_cotangent: &[f64],
    extinction_gradient: &mut [f64],
    scratch: &mut [f64],
) -> Result<()> {
    solar_operator.operator.accumulate_vjp(
        transmission,
        transmission_cotangent,
        extinction_gradient,
        scratch,
    )?;
    Ok(())
}

fn accumulate_solar_transmission_vjp_batch4(
    solar_operator: &RustSolarTransmissionOperator,
    transmission: &[f64],
    transmission_cotangent: &mut [f64],
    extinction_gradient: &mut [f64],
    scratch: &mut [f64],
) -> Result<()> {
    solar_operator.operator.accumulate_vjp_batch4(
        transmission,
        transmission_cotangent,
        extinction_gradient,
        scratch,
    )?;
    Ok(())
}

fn solar_transmission_input_size(solar_operator: &RustSolarTransmissionOperator) -> usize {
    solar_operator.operator.input_size()
}

fn solar_transmission_output_size(solar_operator: &RustSolarTransmissionOperator) -> usize {
    solar_operator.operator.output_size()
}

fn solar_transmission_scratch_size(solar_operator: &RustSolarTransmissionOperator) -> usize {
    solar_operator.operator.scratch_size()
}

fn solar_transmission_forward_scratch_size(
    solar_operator: &RustSolarTransmissionOperator,
) -> usize {
    solar_operator.operator.forward_scratch_size()
}

fn solar_transmission_storage_bytes(solar_operator: &RustSolarTransmissionOperator) -> usize {
    solar_operator.operator.storage_bytes()
}

#[allow(clippy::too_many_arguments)]
fn new_source_interpolator(
    ray_layer_offsets: &[u32],
    layer_atmosphere_offsets: &[u32],
    atmosphere_indices: &[u32],
    atmosphere_weights: &[f64],
    optical_depth_weights: &[f64],
    layer_source_offsets: &[u32],
    source_indices: &[u32],
    source_weights: &[f64],
    ray_ground_offsets: &[u32],
    ground_source_indices: &[u32],
    ground_source_weights: &[f64],
    num_atmosphere_points: usize,
    num_solution_points: usize,
    num_stokes: usize,
) -> Result<Box<RustSourceInterpolator>> {
    Ok(Box::new(RustSourceInterpolator {
        interpolator: SourceInterpolator::new(
            ray_layer_offsets,
            layer_atmosphere_offsets,
            atmosphere_indices,
            atmosphere_weights,
            optical_depth_weights,
            layer_source_offsets,
            source_indices,
            source_weights,
            ray_ground_offsets,
            ground_source_indices,
            ground_source_weights,
            num_atmosphere_points,
            num_solution_points,
            num_stokes,
        )?,
    }))
}

#[allow(clippy::too_many_arguments)]
fn accumulate_source_ray_value(
    interpolator: &RustSourceInterpolator,
    ray_index: usize,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    wavelength_count: usize,
    solution: &[f64],
    solution_stride: usize,
    source: &mut [f64],
    source_stride: usize,
) -> Result<()> {
    interpolator.interpolator.accumulate_ray_value(
        ray_index,
        extinction,
        single_scatter_albedo,
        wavelength_count,
        solution,
        solution_stride,
        source,
        source_stride,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn accumulate_source_ray_jacobian(
    interpolator: &RustSourceInterpolator,
    ray_index: usize,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    wavelength_count: usize,
    solution: &[f64],
    solution_stride: usize,
    source: &mut [f64],
    source_stride: usize,
    source_derivative: &mut [f64],
    num_derivatives: usize,
    derivative_stride: usize,
    single_scatter_albedo_derivative_start: usize,
) -> Result<()> {
    interpolator.interpolator.accumulate_ray_jacobian(
        ray_index,
        extinction,
        single_scatter_albedo,
        wavelength_count,
        solution,
        solution_stride,
        source,
        source_stride,
        source_derivative,
        num_derivatives,
        derivative_stride,
        single_scatter_albedo_derivative_start,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn accumulate_source_ray_jvp(
    interpolator: &RustSourceInterpolator,
    ray_index: usize,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    extinction_tangent: &[f64],
    single_scatter_albedo_tangent: &[f64],
    solution: &[f64],
    solution_tangent: &[f64],
    source: &mut [f64],
    source_tangent: &mut [f64],
) -> Result<()> {
    interpolator.interpolator.accumulate_ray_jvp(
        ray_index,
        extinction,
        single_scatter_albedo,
        extinction_tangent,
        single_scatter_albedo_tangent,
        solution,
        solution_tangent,
        source,
        source_tangent,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn accumulate_source_ray_vjp(
    interpolator: &RustSourceInterpolator,
    ray_index: usize,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    solution: &[f64],
    cotangent: &[f64],
    solution_cotangent: &mut [f64],
    single_scatter_albedo_gradient: &mut [f64],
    extinction_gradient: &mut [f64],
) -> Result<()> {
    interpolator.interpolator.accumulate_ray_vjp(
        ray_index,
        extinction,
        single_scatter_albedo,
        solution,
        cotangent,
        solution_cotangent,
        single_scatter_albedo_gradient,
        extinction_gradient,
    )?;
    Ok(())
}

fn source_interpolator_num_rays(interpolator: &RustSourceInterpolator) -> usize {
    interpolator.interpolator.num_rays()
}

fn source_interpolator_num_layers(interpolator: &RustSourceInterpolator) -> usize {
    interpolator.interpolator.num_layers()
}

fn source_interpolator_storage_bytes(interpolator: &RustSourceInterpolator) -> usize {
    interpolator.interpolator.storage_bytes()
}

pub struct RustSuccessiveOrdersSolver {
    solver: SuccessiveOrdersSolver,
    transport_row_offsets: Option<Arc<[i32]>>,
    transport_column_indices: Option<Arc<[i32]>>,
    scattering_coefficients: Vec<f64>,
    scattering_coefficient_tangent: Vec<f64>,
    boundary_scattering_values: Vec<f64>,
    boundary_scattering_value_tangent: Vec<f64>,
    transport_values: Vec<f64>,
    transport_value_tangent: Vec<f64>,
    layer_optical_depth: Vec<f64>,
    layer_attenuation: Vec<f64>,
    layer_prefix_attenuation: Vec<f64>,
    ray_end_attenuation: Vec<f64>,
    ray_end_attenuation_tangent: Vec<f64>,
    ray_end_attenuation_cotangent: Vec<f64>,
    first_order_forcing: Vec<f64>,
    first_order_forcing_tangent: Vec<f64>,
    solution_jvp: Vec<f64>,
    transport_value_gradient: Vec<f64>,
    scattering_coefficient_gradient: Vec<f64>,
    boundary_scattering_value_gradient: Vec<f64>,
    first_order_forcing_gradient: Vec<f64>,
    simd: SimdWorkspace,
}

impl RustSuccessiveOrdersSolver {
    fn new(solver: SuccessiveOrdersSolver) -> Self {
        Self {
            solver,
            transport_row_offsets: None,
            transport_column_indices: None,
            scattering_coefficients: Vec::new(),
            scattering_coefficient_tangent: Vec::new(),
            boundary_scattering_values: Vec::new(),
            boundary_scattering_value_tangent: Vec::new(),
            transport_values: Vec::new(),
            transport_value_tangent: Vec::new(),
            layer_optical_depth: Vec::new(),
            layer_attenuation: Vec::new(),
            layer_prefix_attenuation: Vec::new(),
            ray_end_attenuation: Vec::new(),
            ray_end_attenuation_tangent: Vec::new(),
            ray_end_attenuation_cotangent: Vec::new(),
            first_order_forcing: Vec::new(),
            first_order_forcing_tangent: Vec::new(),
            solution_jvp: Vec::new(),
            transport_value_gradient: Vec::new(),
            scattering_coefficient_gradient: Vec::new(),
            boundary_scattering_value_gradient: Vec::new(),
            first_order_forcing_gradient: Vec::new(),
            simd: SimdWorkspace::default(),
        }
    }

    fn attach_transport_sparsity(&mut self, sparsity: &RustTransportSparsity) {
        self.transport_row_offsets = Some(Arc::clone(&sparsity.row_offsets));
        self.transport_column_indices = Some(Arc::clone(&sparsity.column_indices));
        self.boundary_scattering_values
            .resize(self.solver.problem_mut().dense_scattering_value_size(), 0.0);
    }

    fn transport_sparsity(&self) -> Result<SharedTransportSparsity> {
        Ok((
            Arc::clone(
                self.transport_row_offsets
                    .as_ref()
                    .ok_or_else(|| anyhow!("solver has no Rust-owned transport sparsity"))?,
            ),
            Arc::clone(
                self.transport_column_indices
                    .as_ref()
                    .ok_or_else(|| anyhow!("solver has no Rust-owned transport sparsity"))?,
            ),
        ))
    }

    fn prepare_transport_workspace(
        &mut self,
        transport: &ScalarRayTransport,
        retain_layer_scratch: bool,
    ) {
        self.transport_values
            .resize(transport.transport_value_size(), 0.0);
        self.ray_end_attenuation.resize(transport.num_rays(), 0.0);
        if retain_layer_scratch {
            self.layer_optical_depth.resize(transport.num_layers(), 0.0);
            self.layer_attenuation.resize(transport.num_layers(), 0.0);
            self.layer_prefix_attenuation
                .resize(transport.num_layers(), 0.0);
        } else {
            // Layer scratch is only needed while forming a VJP. Return it to
            // the allocator after a derivative calculation instead of making
            // every wavelength worker retain the high-water mark forever.
            self.layer_optical_depth = Vec::new();
            self.layer_attenuation = Vec::new();
            self.layer_prefix_attenuation = Vec::new();
        }
    }
}

fn set_scattering_coefficients_from_atmosphere(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    interpolator: &RustScatteringCoefficientInterpolator,
    atmosphere_coefficients: &[f64],
    atmosphere_coefficient_stride: usize,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let coefficient_size = this.solver.problem_mut().scattering_coefficient_size();
    this.scattering_coefficients.resize(coefficient_size, 0.0);
    interpolator.interpolator.interpolate(
        atmosphere_coefficients,
        atmosphere_coefficient_stride,
        &mut this.scattering_coefficients,
    )?;
    this.solver
        .problem_mut()
        .set_scattering_coefficients(&this.scattering_coefficients)?;
    Ok(())
}

fn set_scattering_coefficients_from_atmosphere_batch4(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    interpolator: &RustScatteringCoefficientInterpolator,
    atmosphere_coefficients: &[f64],
    atmosphere_coefficient_stride: usize,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let coefficient_size = this.solver.problem_mut().scattering_coefficient_size();
    this.simd
        .batch4_scattering_coefficients
        .resize(coefficient_size * LANES, 0.0);
    interpolator
        .interpolator
        .interpolate_batch4(
            atmosphere_coefficients,
            atmosphere_coefficient_stride,
            &mut this.simd.batch4_scattering_coefficients,
        )
        .map_err(|error| {
            anyhow!("interpolating four-wavelength scattering coefficients: {error}")
        })?;
    Ok(())
}

fn boundary_scattering_values_mut(solver: Pin<&mut RustSuccessiveOrdersSolver>) -> &mut [f64] {
    &mut solver.get_mut().boundary_scattering_values
}

fn batch4_boundary_scattering_values_mut(
    solver: Pin<&mut RustSuccessiveOrdersSolver>,
) -> &mut [f64] {
    let this = solver.get_mut();
    this.simd
        .batch4_boundary_scattering_values
        .resize(this.boundary_scattering_values.len() * LANES, 0.0);
    this.simd.batch4_boundary_scattering_values.fill(0.0);
    &mut this.simd.batch4_boundary_scattering_values
}

fn boundary_scattering_value_tangent_mut(
    solver: Pin<&mut RustSuccessiveOrdersSolver>,
) -> &mut [f64] {
    let this = solver.get_mut();
    this.boundary_scattering_value_tangent
        .resize(this.boundary_scattering_values.len(), 0.0);
    &mut this.boundary_scattering_value_tangent
}

#[allow(clippy::too_many_arguments)]
fn set_scattering_coefficient_jvp_from_atmosphere(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    interpolator: &RustScatteringCoefficientInterpolator,
    atmosphere_coefficient_derivatives: &[f64],
    atmosphere_coefficient_stride: usize,
    num_wavelengths: usize,
    wavelength_index: usize,
    derivative_group_tangent: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let coefficient_size = this.solver.problem_mut().scattering_coefficient_size();
    this.scattering_coefficient_tangent
        .resize(coefficient_size, 0.0);
    interpolator.interpolator.interpolate_jvp(
        atmosphere_coefficient_derivatives,
        atmosphere_coefficient_stride,
        num_wavelengths,
        wavelength_index,
        derivative_group_tangent,
        &mut this.scattering_coefficient_tangent,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn accumulate_solver_scattering_coefficient_vjp(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    interpolator: &RustScatteringCoefficientInterpolator,
    atmosphere_coefficient_derivatives: &[f64],
    atmosphere_coefficient_stride: usize,
    num_wavelengths: usize,
    wavelength_index: usize,
    active_derivative_groups: &[i32],
    atmosphere_gradient: &mut [f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    interpolator.interpolator.accumulate_vjp(
        &this.scattering_coefficient_gradient,
        atmosphere_coefficient_derivatives,
        atmosphere_coefficient_stride,
        num_wavelengths,
        wavelength_index,
        active_derivative_groups,
        atmosphere_gradient,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn new_successive_orders_solver(
    transport_rows: usize,
    transport_columns: usize,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    scattering_output_offsets: &[usize],
    scattering_input_offsets: &[usize],
    scattering_value_offsets: &[usize],
    max_iterations: usize,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    anderson_depth: usize,
    damping: f64,
) -> Result<Box<RustSuccessiveOrdersSolver>> {
    let transport = CsrMatrix::new(
        transport_rows,
        transport_columns,
        transport_row_offsets,
        transport_column_indices,
    )?;
    let scattering_value_size = scattering_value_offsets
        .last()
        .copied()
        .ok_or_else(|| anyhow!("successive-orders scattering offsets are empty"))?;
    let scattering = BlockDiagonalMatrix::new(
        scattering_output_offsets.to_vec(),
        scattering_input_offsets.to_vec(),
        scattering_value_offsets.to_vec(),
        vec![0.0; scattering_value_size],
    )?;
    let problem = FixedPointProblem::new(transport, scattering)?;
    let config = SolverConfig {
        max_iterations,
        relative_tolerance,
        absolute_tolerance,
        anderson_depth,
        damping,
    };
    Ok(Box::new(RustSuccessiveOrdersSolver::new(
        SuccessiveOrdersSolver::new(problem, config)?,
    )))
}

#[allow(clippy::too_many_arguments)]
fn new_scalar_coefficient_successive_orders_solver(
    transport_rows: usize,
    transport_columns: usize,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    scattering_output_offsets: &[usize],
    scattering_input_offsets: &[usize],
    coefficient_blocks: usize,
    num_coefficients: usize,
    incoming_directions: &[f64],
    incoming_weights: &[f64],
    outgoing_directions: &[f64],
    max_iterations: usize,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    anderson_depth: usize,
    damping: f64,
) -> Result<Box<RustSuccessiveOrdersSolver>> {
    let transport = CsrMatrix::new(
        transport_rows,
        transport_columns,
        transport_row_offsets,
        transport_column_indices,
    )?;
    let incoming_directions = unpack_directions(incoming_directions)?;
    let outgoing_directions = unpack_directions(outgoing_directions)?;
    let basis = ScalarCoefficientBasis::from_directions(
        &incoming_directions,
        incoming_weights,
        &outgoing_directions,
        num_coefficients,
    )?;
    let scattering = ScalarCoefficientScattering::new(
        scattering_output_offsets.to_vec(),
        scattering_input_offsets.to_vec(),
        coefficient_blocks,
        basis,
    )?;
    let problem = FixedPointProblem::new(transport, scattering)?;
    let config = SolverConfig {
        max_iterations,
        relative_tolerance,
        absolute_tolerance,
        anderson_depth,
        damping,
    };
    Ok(Box::new(RustSuccessiveOrdersSolver::new(
        SuccessiveOrdersSolver::new(problem, config)?,
    )))
}

#[allow(clippy::too_many_arguments)]
fn new_vector_coefficient_successive_orders_solver(
    transport_rows: usize,
    transport_columns: usize,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    scattering_output_offsets: &[usize],
    scattering_input_offsets: &[usize],
    coefficient_blocks: usize,
    num_coefficients: usize,
    incoming_directions: &[f64],
    incoming_weights: &[f64],
    outgoing_directions: &[f64],
    max_iterations: usize,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    anderson_depth: usize,
    damping: f64,
) -> Result<Box<RustSuccessiveOrdersSolver>> {
    let transport = CsrMatrix::new(
        transport_rows,
        transport_columns,
        transport_row_offsets,
        transport_column_indices,
    )?;
    let incoming_directions = unpack_directions(incoming_directions)?;
    let outgoing_directions = unpack_directions(outgoing_directions)?;
    let basis = VectorCoefficientBasis::from_directions(
        &incoming_directions,
        incoming_weights,
        &outgoing_directions,
        num_coefficients,
    )?;
    let scattering = VectorCoefficientScattering::new_with_intensity_only_dense(
        scattering_output_offsets.to_vec(),
        scattering_input_offsets.to_vec(),
        coefficient_blocks,
        basis,
    )?;
    let problem = FixedPointProblem::new(transport, scattering)?;
    let config = SolverConfig {
        max_iterations,
        relative_tolerance,
        absolute_tolerance,
        anderson_depth,
        damping,
    };
    Ok(Box::new(RustSuccessiveOrdersSolver::new(
        SuccessiveOrdersSolver::new(problem, config)?,
    )))
}

#[allow(clippy::too_many_arguments)]
fn new_scalar_coefficient_successive_orders_solver_with_sparsity(
    transport_sparsity: &RustTransportSparsity,
    scattering_output_offsets: &[usize],
    scattering_input_offsets: &[usize],
    coefficient_blocks: usize,
    num_coefficients: usize,
    incoming_directions: &[f64],
    incoming_weights: &[f64],
    outgoing_directions: &[f64],
    max_iterations: usize,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    anderson_depth: usize,
    damping: f64,
) -> Result<Box<RustSuccessiveOrdersSolver>> {
    let mut solver = new_scalar_coefficient_successive_orders_solver(
        transport_sparsity.rows,
        transport_sparsity.columns,
        &transport_sparsity.row_offsets,
        &transport_sparsity.column_indices,
        scattering_output_offsets,
        scattering_input_offsets,
        coefficient_blocks,
        num_coefficients,
        incoming_directions,
        incoming_weights,
        outgoing_directions,
        max_iterations,
        relative_tolerance,
        absolute_tolerance,
        anderson_depth,
        damping,
    )?;
    solver.attach_transport_sparsity(transport_sparsity);
    Ok(solver)
}

#[allow(clippy::too_many_arguments)]
fn new_vector_coefficient_successive_orders_solver_with_sparsity(
    transport_sparsity: &RustTransportSparsity,
    scattering_output_offsets: &[usize],
    scattering_input_offsets: &[usize],
    coefficient_blocks: usize,
    num_coefficients: usize,
    incoming_directions: &[f64],
    incoming_weights: &[f64],
    outgoing_directions: &[f64],
    max_iterations: usize,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    anderson_depth: usize,
    damping: f64,
) -> Result<Box<RustSuccessiveOrdersSolver>> {
    let mut solver = new_vector_coefficient_successive_orders_solver(
        transport_sparsity.rows,
        transport_sparsity.columns,
        &transport_sparsity.row_offsets,
        &transport_sparsity.column_indices,
        scattering_output_offsets,
        scattering_input_offsets,
        coefficient_blocks,
        num_coefficients,
        incoming_directions,
        incoming_weights,
        outgoing_directions,
        max_iterations,
        relative_tolerance,
        absolute_tolerance,
        anderson_depth,
        damping,
    )?;
    solver.attach_transport_sparsity(transport_sparsity);
    Ok(solver)
}

fn solve(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    transport_values: &[f64],
    scattering_values: &[f64],
    first_order_forcing: &[f64],
    initial_guess: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.solver
        .problem_mut()
        .set_scattering_values(scattering_values)?;
    let initial_guess = (!initial_guess.is_empty()).then_some(initial_guess);
    this.solver.solve(
        transport_row_offsets,
        transport_column_indices,
        transport_values,
        first_order_forcing,
        initial_guess,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn solve_scalar_coefficients(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    transport_values: &[f64],
    boundary_scattering_values: &[f64],
    first_order_forcing: &[f64],
    initial_guess: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.solver
        .problem_mut()
        .set_scattering_values(boundary_scattering_values)?;
    let initial_guess = (!initial_guess.is_empty()).then_some(initial_guess);
    this.solver.solve(
        transport_row_offsets,
        transport_column_indices,
        transport_values,
        first_order_forcing,
        initial_guess,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn solve_vector_coefficients(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    transport_values: &[f64],
    boundary_scattering_values: &[f64],
    first_order_forcing: &[f64],
    initial_guess: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.solver
        .problem_mut()
        .set_scattering_values(boundary_scattering_values)?;
    let initial_guess = (!initial_guess.is_empty()).then_some(initial_guess);
    this.solver.solve(
        transport_row_offsets,
        transport_column_indices,
        transport_values,
        first_order_forcing,
        initial_guess,
    )?;
    Ok(())
}

fn solve_coefficients_assembled(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    initial_guess: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let (transport_row_offsets, transport_column_indices) = this.transport_sparsity()?;
    this.solver
        .problem_mut()
        .set_scattering_values(&this.boundary_scattering_values)?;
    let initial_guess = (!initial_guess.is_empty()).then_some(initial_guess);
    this.solver.solve(
        &transport_row_offsets,
        &transport_column_indices,
        &this.transport_values,
        &this.first_order_forcing,
        initial_guess,
    )?;
    Ok(())
}

fn begin_coefficients_batch4(mut solver: Pin<&mut RustSuccessiveOrdersSolver>) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.transport_row_offsets.is_none() || this.transport_column_indices.is_none() {
        return Err(anyhow!(
            "wavelength SIMD requires Rust-owned transport sparsity"
        ));
    }
    this.simd.batch4_staged_lanes = 0;
    this.simd.batch4_has_initial_guess = false;
    this.simd.batch4_diagnostics = std::array::from_fn(|_| SolverDiagnostics::default());
    Ok(())
}

fn stage_coefficients_batch4_lane(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    lane: usize,
    initial_guess: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if lane >= LANES || lane != this.simd.batch4_staged_lanes {
        return Err(anyhow!(
            "successive-orders SIMD lanes must be staged in order"
        ));
    }
    if lane == 0 {
        this.simd
            .batch4_transport_values
            .resize(this.transport_values.len() * LANES, 0.0);
        this.simd
            .batch4_scattering_coefficients
            .resize(this.scattering_coefficients.len() * LANES, 0.0);
        this.simd
            .batch4_boundary_scattering_values
            .resize(this.boundary_scattering_values.len() * LANES, 0.0);
        this.simd
            .batch4_first_order_forcing
            .resize(this.first_order_forcing.len() * LANES, 0.0);
        this.simd
            .batch4_solution
            .resize(this.solver.solution().len() * LANES, 0.0);
        this.simd.batch4_has_initial_guess = !initial_guess.is_empty();
        if this.simd.batch4_has_initial_guess {
            if initial_guess.len() != this.solver.solution().len() {
                return Err(anyhow!(
                    "successive-orders SIMD initial guess has the wrong size"
                ));
            }
            this.simd
                .batch4_initial_guess
                .resize(initial_guess.len() * LANES, 0.0);
        } else {
            this.simd.batch4_initial_guess.clear();
        }
    } else if this.simd.batch4_transport_values.len() != this.transport_values.len() * LANES
        || this.simd.batch4_scattering_coefficients.len()
            != this.scattering_coefficients.len() * LANES
        || this.simd.batch4_boundary_scattering_values.len()
            != this.boundary_scattering_values.len() * LANES
        || this.simd.batch4_first_order_forcing.len() != this.first_order_forcing.len() * LANES
        || this.simd.batch4_has_initial_guess == initial_guess.is_empty()
    {
        return Err(anyhow!(
            "successive-orders SIMD wavelength storage changed between lanes"
        ));
    }
    if this.simd.batch4_has_initial_guess && initial_guess.len() != this.solver.solution().len() {
        return Err(anyhow!(
            "successive-orders SIMD initial guess has the wrong size"
        ));
    }

    interleave_lane(
        &this.transport_values,
        lane,
        &mut this.simd.batch4_transport_values,
    );
    interleave_lane(
        &this.scattering_coefficients,
        lane,
        &mut this.simd.batch4_scattering_coefficients,
    );
    interleave_lane(
        &this.boundary_scattering_values,
        lane,
        &mut this.simd.batch4_boundary_scattering_values,
    );
    interleave_lane(
        &this.first_order_forcing,
        lane,
        &mut this.simd.batch4_first_order_forcing,
    );
    if this.simd.batch4_has_initial_guess {
        interleave_lane(initial_guess, lane, &mut this.simd.batch4_initial_guess);
    }
    this.simd.batch4_staged_lanes += 1;
    Ok(())
}

fn solve_coefficients_batch4(mut solver: Pin<&mut RustSuccessiveOrdersSolver>) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.simd.batch4_staged_lanes != LANES {
        return Err(anyhow!("successive-orders SIMD batch is incomplete"));
    }
    let problem = this.solver.problem_mut();
    let expected_transport = problem.transport_value_size() * LANES;
    let expected_scattering = problem.scattering_coefficient_size() * LANES;
    let expected_boundary = problem.dense_scattering_value_size() * LANES;
    let expected_forcing = problem.incoming_size() * LANES;
    if this.simd.batch4_transport_values.len() != expected_transport
        || this.simd.batch4_scattering_coefficients.len() != expected_scattering
        || this.simd.batch4_boundary_scattering_values.len() != expected_boundary
        || this.simd.batch4_first_order_forcing.len() != expected_forcing
    {
        return Err(anyhow!(
            "successive-orders SIMD buffer sizes are transport {}/{}, scattering {}/{}, boundary {}/{}, forcing {}/{}",
            this.simd.batch4_transport_values.len(),
            expected_transport,
            this.simd.batch4_scattering_coefficients.len(),
            expected_scattering,
            this.simd.batch4_boundary_scattering_values.len(),
            expected_boundary,
            this.simd.batch4_first_order_forcing.len(),
            expected_forcing,
        ));
    }
    for (name, values) in [
        ("transport", this.simd.batch4_transport_values.as_slice()),
        (
            "scattering",
            this.simd.batch4_scattering_coefficients.as_slice(),
        ),
        (
            "boundary",
            this.simd.batch4_boundary_scattering_values.as_slice(),
        ),
        ("forcing", this.simd.batch4_first_order_forcing.as_slice()),
    ] {
        if let Some((index, value)) = values
            .iter()
            .copied()
            .enumerate()
            .find(|(_, value)| !value.is_finite())
        {
            return Err(anyhow!(
                "successive-orders SIMD {name} value {index} is non-finite ({value})"
            ));
        }
    }
    let (transport_row_offsets, transport_column_indices) = this.transport_sparsity()?;
    let initial_guess = this
        .simd
        .batch4_has_initial_guess
        .then_some(this.simd.batch4_initial_guess.as_slice());
    let packed_result = this.solver.solve_batch4(
        &transport_row_offsets,
        &transport_column_indices,
        &this.simd.batch4_transport_values,
        &this.simd.batch4_scattering_coefficients,
        &this.simd.batch4_boundary_scattering_values,
        &this.simd.batch4_first_order_forcing,
        initial_guess,
    );
    match packed_result {
        Ok(result) => {
            this.simd.batch4_solution = result.solution;
            this.simd.batch4_diagnostics = result.diagnostics;
        }
        Err(SolverError::ImplicitLinearSolveDidNotConverge)
        | Err(SolverError::Operator(super::OperatorError::UnsupportedOperator)) => {
            // Preserve the established convergence and Anderson behavior for
            // uncommon breakdowns and non-SIMD solver configurations.
            let mut transport_values = vec![0.0; this.simd.batch4_transport_values.len() / LANES];
            let mut scattering_coefficients =
                vec![0.0; this.simd.batch4_scattering_coefficients.len() / LANES];
            let mut boundary_values =
                vec![0.0; this.simd.batch4_boundary_scattering_values.len() / LANES];
            let mut forcing = vec![0.0; this.simd.batch4_first_order_forcing.len() / LANES];
            let mut initial = vec![0.0; this.solver.solution().len()];
            for lane in 0..LANES {
                extract_lane(
                    &this.simd.batch4_transport_values,
                    lane,
                    &mut transport_values,
                );
                extract_lane(
                    &this.simd.batch4_scattering_coefficients,
                    lane,
                    &mut scattering_coefficients,
                );
                extract_lane(
                    &this.simd.batch4_boundary_scattering_values,
                    lane,
                    &mut boundary_values,
                );
                extract_lane(&this.simd.batch4_first_order_forcing, lane, &mut forcing);
                this.solver
                    .problem_mut()
                    .set_scattering_coefficients(&scattering_coefficients)?;
                this.solver
                    .problem_mut()
                    .set_scattering_values(&boundary_values)?;
                let initial_guess = if this.simd.batch4_has_initial_guess {
                    extract_lane(&this.simd.batch4_initial_guess, lane, &mut initial);
                    Some(initial.as_slice())
                } else {
                    None
                };
                this.solver.solve(
                    &transport_row_offsets,
                    &transport_column_indices,
                    &transport_values,
                    &forcing,
                    initial_guess,
                )?;
                interleave_lane(this.solver.solution(), lane, &mut this.simd.batch4_solution);
                this.simd.batch4_diagnostics[lane] = this.solver.diagnostics().clone();
            }
        }
        Err(error) => return Err(anyhow!("four-wavelength solve failed: {error}")),
    }
    Ok(())
}

fn restore_coefficients_batch4_solution(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    solution: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let expected = this.solver.solution().len() * LANES;
    if solution.len() != expected {
        return Err(anyhow!(
            "restored four-wavelength solution has size {}; expected {expected}",
            solution.len()
        ));
    }
    if solution.iter().any(|value| !value.is_finite()) {
        return Err(anyhow!(
            "restored four-wavelength solution contains a non-finite value"
        ));
    }
    this.simd.batch4_solution.clear();
    this.simd.batch4_solution.extend_from_slice(solution);
    this.simd.batch4_staged_lanes = LANES;
    for diagnostics in &mut this.simd.batch4_diagnostics {
        *diagnostics = SolverDiagnostics {
            iterations: 0,
            converged: true,
            reason: super::ConvergenceReason::Tolerance,
            initial_residual: f64::NAN,
            final_residual: f64::NAN,
            residual_history: Vec::new(),
        };
    }
    Ok(())
}

fn batch4_solution(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.simd.batch4_solution
}

fn select_coefficients_batch4_lane(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    lane: usize,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if lane >= LANES || this.simd.batch4_staged_lanes != LANES {
        return Err(anyhow!("invalid four-wavelength primal lane"));
    }
    this.simd
        .batch4_lane_solution
        .resize(this.solver.solution().len(), 0.0);
    extract_lane(
        &this.simd.batch4_solution,
        lane,
        &mut this.simd.batch4_lane_solution,
    );
    this.solver
        .restore_converged_solution(&this.simd.batch4_lane_solution)?;
    Ok(())
}

fn batch4_first_order_forcing(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.simd.batch4_first_order_forcing
}

fn batch4_diagnostics(
    solver: &RustSuccessiveOrdersSolver,
    lane: usize,
) -> Option<&SolverDiagnostics> {
    (lane < LANES && solver.simd.batch4_staged_lanes == LANES)
        .then(|| &solver.simd.batch4_diagnostics[lane])
}

fn batch4_converged(solver: &RustSuccessiveOrdersSolver, lane: usize) -> bool {
    batch4_diagnostics(solver, lane).is_some_and(|diagnostics| diagnostics.converged)
}

fn batch4_iterations(solver: &RustSuccessiveOrdersSolver, lane: usize) -> usize {
    batch4_diagnostics(solver, lane).map_or(0, |diagnostics| diagnostics.iterations)
}

fn batch4_initial_residual(solver: &RustSuccessiveOrdersSolver, lane: usize) -> f64 {
    batch4_diagnostics(solver, lane).map_or(f64::NAN, |diagnostics| diagnostics.initial_residual)
}

fn batch4_final_residual(solver: &RustSuccessiveOrdersSolver, lane: usize) -> f64 {
    batch4_diagnostics(solver, lane).map_or(f64::NAN, |diagnostics| diagnostics.final_residual)
}

fn batch4_workspace_bytes(solver: &RustSuccessiveOrdersSolver) -> usize {
    (solver.simd.batch4_transport_values.capacity()
        + solver.simd.batch4_scattering_coefficients.capacity()
        + solver.simd.batch4_boundary_scattering_values.capacity()
        + solver.simd.batch4_first_order_forcing.capacity()
        + solver.simd.batch4_initial_guess.capacity()
        + solver.simd.batch4_solution.capacity()
        + solver.simd.batch4_transport_value_tangent.capacity()
        + solver.simd.batch4_scattering_coefficient_tangent.capacity()
        + solver
            .simd
            .batch4_boundary_scattering_value_tangent
            .capacity()
        + solver.simd.batch4_first_order_forcing_tangent.capacity()
        + solver.simd.batch4_solution_jvp.capacity()
        + solver.simd.batch4_extinction_tangent.capacity()
        + solver.simd.batch4_single_scatter_albedo_tangent.capacity()
        + solver.simd.batch4_legendre_coefficient_tangent.capacity()
        + solver.simd.batch4_end_of_ray_source.capacity()
        + solver.simd.batch4_end_of_ray_source_tangent.capacity()
        + solver.simd.batch4_solution_cotangent.capacity()
        + solver
            .simd
            .batch4_scattering_coefficient_gradient
            .capacity()
        + solver
            .simd
            .batch4_boundary_scattering_value_gradient
            .capacity()
        + solver.simd.batch4_first_order_forcing_gradient.capacity()
        + solver.simd.batch4_extinction_gradient.capacity()
        + solver.simd.batch4_single_scatter_albedo_gradient.capacity()
        + solver.simd.batch4_legendre_coefficient_gradient.capacity()
        + solver.simd.batch4_end_of_ray_source_gradient.capacity()
        + solver.simd.batch4_lane_solution.capacity())
        * size_of::<f64>()
}

fn begin_linearize_coefficients_jvp_batch4(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.simd.batch4_staged_lanes != LANES
        || this.simd.batch4_solution.len() != this.solver.solution().len() * LANES
    {
        return Err(anyhow!(
            "four-wavelength JVP requires a completed packed primal solve"
        ));
    }
    if this
        .simd
        .batch4_diagnostics
        .iter()
        .any(|state| !state.converged)
    {
        return Err(anyhow!(
            "four-wavelength JVP requires converged packed primal lanes"
        ));
    }
    this.simd
        .batch4_transport_value_tangent
        .resize(this.simd.batch4_transport_values.len(), 0.0);
    this.simd
        .batch4_scattering_coefficient_tangent
        .resize(this.simd.batch4_scattering_coefficients.len(), 0.0);
    this.simd
        .batch4_boundary_scattering_value_tangent
        .resize(this.simd.batch4_boundary_scattering_values.len(), 0.0);
    this.simd
        .batch4_first_order_forcing_tangent
        .resize(this.simd.batch4_first_order_forcing.len(), 0.0);
    this.simd
        .batch4_solution_jvp
        .resize(this.simd.batch4_solution.len(), 0.0);
    this.simd.batch4_jvp_staged_lanes = 0;
    Ok(())
}

fn stage_linearize_coefficients_jvp_batch4_lane(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    lane: usize,
    zero_tangent: bool,
    first_order_forcing: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if lane >= LANES || lane != this.simd.batch4_jvp_staged_lanes {
        return Err(anyhow!("four-wavelength JVP lanes must be staged in order"));
    }
    let expected_transport = this.simd.batch4_transport_values.len() / LANES;
    let expected_scattering = this.simd.batch4_scattering_coefficients.len() / LANES;
    let expected_boundary = this.simd.batch4_boundary_scattering_values.len() / LANES;
    let expected_forcing = this.simd.batch4_first_order_forcing.len() / LANES;
    if this.transport_values.len() != expected_transport
        || first_order_forcing.len() != expected_forcing
    {
        return Err(anyhow!(
            "four-wavelength JVP primal transport storage has inconsistent dimensions"
        ));
    }
    interleave_lane(
        &this.transport_values,
        lane,
        &mut this.simd.batch4_transport_values,
    );
    interleave_lane(
        first_order_forcing,
        lane,
        &mut this.simd.batch4_first_order_forcing,
    );
    if zero_tangent {
        for target in [
            &mut this.simd.batch4_transport_value_tangent,
            &mut this.simd.batch4_scattering_coefficient_tangent,
            &mut this.simd.batch4_boundary_scattering_value_tangent,
            &mut this.simd.batch4_first_order_forcing_tangent,
        ] {
            for element in 0..target.len() / LANES {
                target[element * LANES + lane] = 0.0;
            }
        }
    } else {
        if this.transport_value_tangent.len() != expected_transport
            || this.scattering_coefficient_tangent.len() != expected_scattering
            || this.boundary_scattering_value_tangent.len() != expected_boundary
            || this.first_order_forcing_tangent.len() != expected_forcing
        {
            return Err(anyhow!(
                "four-wavelength JVP tangent storage has inconsistent dimensions"
            ));
        }
        interleave_lane(
            &this.transport_value_tangent,
            lane,
            &mut this.simd.batch4_transport_value_tangent,
        );
        interleave_lane(
            &this.scattering_coefficient_tangent,
            lane,
            &mut this.simd.batch4_scattering_coefficient_tangent,
        );
        interleave_lane(
            &this.boundary_scattering_value_tangent,
            lane,
            &mut this.simd.batch4_boundary_scattering_value_tangent,
        );
        interleave_lane(
            &this.first_order_forcing_tangent,
            lane,
            &mut this.simd.batch4_first_order_forcing_tangent,
        );
    }
    this.simd.batch4_jvp_staged_lanes += 1;
    Ok(())
}

fn linearize_coefficients_jvp_batch4(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.simd.batch4_jvp_staged_lanes != LANES {
        return Err(anyhow!("four-wavelength JVP batch is incomplete"));
    }
    let (transport_row_offsets, transport_column_indices) = this.transport_sparsity()?;
    match this.solver.solve_jvp_batch4(
        &transport_row_offsets,
        &transport_column_indices,
        &this.simd.batch4_transport_values,
        &this.simd.batch4_transport_value_tangent,
        &this.simd.batch4_scattering_coefficients,
        &this.simd.batch4_scattering_coefficient_tangent,
        &this.simd.batch4_boundary_scattering_values,
        &this.simd.batch4_boundary_scattering_value_tangent,
        &this.simd.batch4_first_order_forcing,
        &this.simd.batch4_first_order_forcing_tangent,
        &this.simd.batch4_solution,
    ) {
        Ok(solution_jvp) => this.simd.batch4_solution_jvp = solution_jvp,
        Err(SolverError::ImplicitLinearSolveDidNotConverge)
        | Err(SolverError::Operator(super::OperatorError::UnsupportedOperator)) => {
            let state_size = this.solver.solution().len();
            this.simd.batch4_lane_solution.resize(state_size, 0.0);
            let mut transport = vec![0.0; this.simd.batch4_transport_values.len() / LANES];
            let mut transport_tangent =
                vec![0.0; this.simd.batch4_transport_value_tangent.len() / LANES];
            let mut coefficients =
                vec![0.0; this.simd.batch4_scattering_coefficients.len() / LANES];
            let mut coefficient_tangent =
                vec![0.0; this.simd.batch4_scattering_coefficient_tangent.len() / LANES];
            let mut boundary = vec![0.0; this.simd.batch4_boundary_scattering_values.len() / LANES];
            let mut boundary_tangent =
                vec![0.0; this.simd.batch4_boundary_scattering_value_tangent.len() / LANES];
            let mut forcing = vec![0.0; this.simd.batch4_first_order_forcing.len() / LANES];
            let mut forcing_tangent =
                vec![0.0; this.simd.batch4_first_order_forcing_tangent.len() / LANES];
            for lane in 0..LANES {
                extract_lane(&this.simd.batch4_transport_values, lane, &mut transport);
                extract_lane(
                    &this.simd.batch4_transport_value_tangent,
                    lane,
                    &mut transport_tangent,
                );
                extract_lane(
                    &this.simd.batch4_scattering_coefficients,
                    lane,
                    &mut coefficients,
                );
                extract_lane(
                    &this.simd.batch4_scattering_coefficient_tangent,
                    lane,
                    &mut coefficient_tangent,
                );
                extract_lane(
                    &this.simd.batch4_boundary_scattering_values,
                    lane,
                    &mut boundary,
                );
                extract_lane(
                    &this.simd.batch4_boundary_scattering_value_tangent,
                    lane,
                    &mut boundary_tangent,
                );
                extract_lane(&this.simd.batch4_first_order_forcing, lane, &mut forcing);
                extract_lane(
                    &this.simd.batch4_first_order_forcing_tangent,
                    lane,
                    &mut forcing_tangent,
                );
                extract_lane(
                    &this.simd.batch4_solution,
                    lane,
                    &mut this.simd.batch4_lane_solution,
                );
                this.solver
                    .problem_mut()
                    .set_scattering_coefficients(&coefficients)?;
                this.solver.problem_mut().set_scattering_values(&boundary)?;
                this.solver
                    .restore_converged_solution(&this.simd.batch4_lane_solution)?;
                let lane_jvp = this.solver.solve_jvp(
                    &transport_row_offsets,
                    &transport_column_indices,
                    &transport,
                    &transport_tangent,
                    &forcing,
                    &forcing_tangent,
                    &coefficient_tangent,
                    &boundary_tangent,
                )?;
                interleave_lane(&lane_jvp, lane, &mut this.simd.batch4_solution_jvp);
            }
        }
        Err(error) => return Err(anyhow!("four-wavelength JVP solve failed: {error}")),
    }
    Ok(())
}

fn select_linearize_coefficients_jvp_batch4_lane(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    lane: usize,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if lane >= LANES || this.simd.batch4_jvp_staged_lanes != LANES {
        return Err(anyhow!("invalid four-wavelength JVP lane"));
    }
    this.simd
        .batch4_lane_solution
        .resize(this.solver.solution().len(), 0.0);
    extract_lane(
        &this.simd.batch4_solution,
        lane,
        &mut this.simd.batch4_lane_solution,
    );
    this.solver
        .restore_converged_solution(&this.simd.batch4_lane_solution)?;
    this.solution_jvp.resize(this.solver.solution().len(), 0.0);
    extract_lane(&this.simd.batch4_solution_jvp, lane, &mut this.solution_jvp);
    Ok(())
}

fn begin_linearize_coefficients_vjp_batch4(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.simd.batch4_staged_lanes != LANES
        || this.simd.batch4_solution.len() != this.solver.solution().len() * LANES
    {
        return Err(anyhow!(
            "four-wavelength VJP requires a completed packed primal solve"
        ));
    }
    if this
        .simd
        .batch4_diagnostics
        .iter()
        .any(|state| !state.converged)
    {
        return Err(anyhow!(
            "four-wavelength VJP requires converged packed primal lanes"
        ));
    }
    this.simd
        .batch4_solution_cotangent
        .resize(this.simd.batch4_solution.len(), 0.0);
    this.simd
        .batch4_scattering_coefficient_gradient
        .resize(this.simd.batch4_scattering_coefficients.len(), 0.0);
    this.simd
        .batch4_boundary_scattering_value_gradient
        .resize(this.simd.batch4_boundary_scattering_values.len(), 0.0);
    this.simd
        .batch4_first_order_forcing_gradient
        .resize(this.simd.batch4_first_order_forcing.len(), 0.0);
    this.simd.batch4_vjp_staged_lanes = 0;
    this.simd.batch4_ray_vjp_ready = false;
    Ok(())
}

fn stage_linearize_coefficients_vjp_batch4_lane(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    lane: usize,
    solution_cotangent: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if lane >= LANES
        || lane != this.simd.batch4_vjp_staged_lanes
        || solution_cotangent.len() != this.solver.solution().len()
    {
        return Err(anyhow!("invalid four-wavelength VJP lane"));
    }
    interleave_lane(
        solution_cotangent,
        lane,
        &mut this.simd.batch4_solution_cotangent,
    );
    this.simd.batch4_vjp_staged_lanes += 1;
    Ok(())
}

fn linearize_coefficients_vjp_batch4(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if this.simd.batch4_vjp_staged_lanes != LANES {
        return Err(anyhow!("four-wavelength VJP batch is incomplete"));
    }
    let (transport_row_offsets, transport_column_indices) = this.transport_sparsity()?;
    match this.solver.solve_vjp_batch4(
        &transport_row_offsets,
        &transport_column_indices,
        &this.simd.batch4_transport_values,
        &this.simd.batch4_scattering_coefficients,
        &this.simd.batch4_boundary_scattering_values,
        &this.simd.batch4_first_order_forcing,
        &this.simd.batch4_solution,
        &this.simd.batch4_solution_cotangent,
    ) {
        Ok(gradient) => {
            this.simd.batch4_scattering_coefficient_gradient = gradient.scattering_coefficients;
            this.simd.batch4_boundary_scattering_value_gradient = gradient.dense_scattering_values;
            this.simd.batch4_first_order_forcing_gradient = gradient.forcing;
        }
        Err(SolverError::ImplicitLinearSolveDidNotConverge)
        | Err(SolverError::Operator(super::OperatorError::UnsupportedOperator)) => {
            let state_size = this.solver.solution().len();
            this.simd.batch4_lane_solution.resize(state_size, 0.0);
            let mut solution_cotangent = vec![0.0; state_size];
            let mut transport = vec![0.0; this.simd.batch4_transport_values.len() / LANES];
            let mut coefficients =
                vec![0.0; this.simd.batch4_scattering_coefficients.len() / LANES];
            let mut boundary = vec![0.0; this.simd.batch4_boundary_scattering_values.len() / LANES];
            let mut forcing = vec![0.0; this.simd.batch4_first_order_forcing.len() / LANES];
            for lane in 0..LANES {
                extract_lane(&this.simd.batch4_transport_values, lane, &mut transport);
                extract_lane(
                    &this.simd.batch4_scattering_coefficients,
                    lane,
                    &mut coefficients,
                );
                extract_lane(
                    &this.simd.batch4_boundary_scattering_values,
                    lane,
                    &mut boundary,
                );
                extract_lane(&this.simd.batch4_first_order_forcing, lane, &mut forcing);
                extract_lane(
                    &this.simd.batch4_solution,
                    lane,
                    &mut this.simd.batch4_lane_solution,
                );
                extract_lane(
                    &this.simd.batch4_solution_cotangent,
                    lane,
                    &mut solution_cotangent,
                );
                this.solver
                    .problem_mut()
                    .set_scattering_coefficients(&coefficients)?;
                this.solver.problem_mut().set_scattering_values(&boundary)?;
                this.solver
                    .restore_converged_solution(&this.simd.batch4_lane_solution)?;
                let gradient = this.solver.solve_vjp_compact(
                    &transport_row_offsets,
                    &transport_column_indices,
                    &transport,
                    &forcing,
                    &solution_cotangent,
                )?;
                interleave_lane(
                    &gradient.scattering_coefficients,
                    lane,
                    &mut this.simd.batch4_scattering_coefficient_gradient,
                );
                interleave_lane(
                    &gradient.dense_scattering_values,
                    lane,
                    &mut this.simd.batch4_boundary_scattering_value_gradient,
                );
                interleave_lane(
                    &gradient.forcing,
                    lane,
                    &mut this.simd.batch4_first_order_forcing_gradient,
                );
            }
        }
        Err(error) => return Err(anyhow!("four-wavelength VJP solve failed: {error}")),
    }
    Ok(())
}

fn select_linearize_coefficients_vjp_batch4_lane(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    lane: usize,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    if lane >= LANES || this.simd.batch4_vjp_staged_lanes != LANES {
        return Err(anyhow!("invalid four-wavelength VJP lane"));
    }
    this.simd
        .batch4_lane_solution
        .resize(this.solver.solution().len(), 0.0);
    extract_lane(
        &this.simd.batch4_solution,
        lane,
        &mut this.simd.batch4_lane_solution,
    );
    this.solver
        .restore_converged_solution(&this.simd.batch4_lane_solution)?;
    this.transport_value_gradient.clear();
    this.scattering_coefficient_gradient.resize(
        this.simd.batch4_scattering_coefficient_gradient.len() / LANES,
        0.0,
    );
    this.boundary_scattering_value_gradient.resize(
        this.simd.batch4_boundary_scattering_value_gradient.len() / LANES,
        0.0,
    );
    this.first_order_forcing_gradient.resize(
        this.simd.batch4_first_order_forcing_gradient.len() / LANES,
        0.0,
    );
    extract_lane(
        &this.simd.batch4_scattering_coefficient_gradient,
        lane,
        &mut this.scattering_coefficient_gradient,
    );
    extract_lane(
        &this.simd.batch4_boundary_scattering_value_gradient,
        lane,
        &mut this.boundary_scattering_value_gradient,
    );
    extract_lane(
        &this.simd.batch4_first_order_forcing_gradient,
        lane,
        &mut this.first_order_forcing_gradient,
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn linearize_coefficients_jvp(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    transport_values: &[f64],
    transport_value_tangent: &[f64],
    boundary_scattering_value_tangent: &[f64],
    first_order_forcing: &[f64],
    first_order_forcing_tangent: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.solution_jvp = this.solver.solve_jvp(
        transport_row_offsets,
        transport_column_indices,
        transport_values,
        transport_value_tangent,
        first_order_forcing,
        first_order_forcing_tangent,
        &this.scattering_coefficient_tangent,
        boundary_scattering_value_tangent,
    )?;
    Ok(())
}

fn linearize_coefficients_jvp_assembled(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let (transport_row_offsets, transport_column_indices) = this.transport_sparsity()?;
    this.solution_jvp = this.solver.solve_jvp(
        &transport_row_offsets,
        &transport_column_indices,
        &this.transport_values,
        &this.transport_value_tangent,
        &this.first_order_forcing,
        &this.first_order_forcing_tangent,
        &this.scattering_coefficient_tangent,
        &this.boundary_scattering_value_tangent,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn linearize_coefficients_vjp(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    transport_values: &[f64],
    first_order_forcing: &[f64],
    solution_cotangent: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    // A VJP is evaluated one wavelength at a time. Drop the previous
    // wavelength's potentially very large transport gradient before the
    // solver allocates its replacement.
    this.transport_value_gradient = Vec::new();
    this.scattering_coefficient_gradient = Vec::new();
    this.boundary_scattering_value_gradient = Vec::new();
    this.first_order_forcing_gradient = Vec::new();
    let gradient = this.solver.solve_vjp_compact(
        transport_row_offsets,
        transport_column_indices,
        transport_values,
        first_order_forcing,
        solution_cotangent,
    )?;
    this.transport_value_gradient = gradient.transport_values;
    this.scattering_coefficient_gradient = gradient.scattering_coefficients;
    this.boundary_scattering_value_gradient = gradient.dense_scattering_values;
    this.first_order_forcing_gradient = gradient.forcing;
    Ok(())
}

fn linearize_coefficients_vjp_assembled(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    solution_cotangent: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let (transport_row_offsets, transport_column_indices) = this.transport_sparsity()?;
    this.transport_value_gradient = Vec::new();
    this.scattering_coefficient_gradient = Vec::new();
    this.boundary_scattering_value_gradient = Vec::new();
    this.first_order_forcing_gradient = Vec::new();
    let gradient = this.solver.solve_vjp_compact(
        &transport_row_offsets,
        &transport_column_indices,
        &this.transport_values,
        &this.first_order_forcing,
        solution_cotangent,
    )?;
    this.transport_value_gradient = gradient.transport_values;
    this.scattering_coefficient_gradient = gradient.scattering_coefficients;
    this.boundary_scattering_value_gradient = gradient.dense_scattering_values;
    this.first_order_forcing_gradient = gradient.forcing;
    Ok(())
}

fn restore_coefficient_solution(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    boundary_scattering_values: &[f64],
    solution: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.solver
        .problem_mut()
        .set_scattering_values(boundary_scattering_values)?;
    this.solver.restore_converged_solution(solution)?;
    Ok(())
}

fn restore_coefficient_solution_assembled(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    first_order_forcing: &[f64],
    solution: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    let (transport_row_offsets, _) = this.transport_sparsity()?;
    if first_order_forcing.len() + 1 != transport_row_offsets.len() {
        return Err(anyhow!("restored first-order forcing has the wrong size"));
    }
    if first_order_forcing.iter().any(|value| !value.is_finite()) {
        return Err(anyhow!(
            "restored first-order forcing contains a non-finite value"
        ));
    }
    this.first_order_forcing.clear();
    this.first_order_forcing
        .extend_from_slice(first_order_forcing);
    this.solver
        .problem_mut()
        .set_scattering_values(&this.boundary_scattering_values)?;
    this.solver.restore_converged_solution(solution)?;
    Ok(())
}

fn restore_solution_only(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    solution: &[f64],
) -> Result<()> {
    solver
        .as_mut()
        .get_mut()
        .solver
        .restore_converged_solution(solution)?;
    Ok(())
}

fn unpack_directions(values: &[f64]) -> Result<Vec<[f64; 3]>> {
    if !values.len().is_multiple_of(3) {
        return Err(anyhow!(
            "successive-orders direction array has the wrong shape"
        ));
    }
    Ok(values
        .chunks_exact(3)
        .map(|direction| [direction[0], direction[1], direction[2]])
        .collect())
}

fn solution(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    solver.solver.solution()
}

fn first_order_forcing(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.first_order_forcing
}

fn solution_jvp(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.solution_jvp
}

fn clear_solution_jvp(mut solver: Pin<&mut RustSuccessiveOrdersSolver>) {
    let this = solver.as_mut().get_mut();
    this.solution_jvp.clear();
    this.solution_jvp.resize(this.solver.solution().len(), 0.0);
}

fn transport_value_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.transport_value_gradient
}

fn boundary_scattering_value_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.boundary_scattering_value_gradient
}

fn first_order_forcing_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.first_order_forcing_gradient
}

fn residual_history(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.solver.diagnostics().residual_history
}

fn iterations(solver: &RustSuccessiveOrdersSolver) -> usize {
    solver.solver.diagnostics().iterations
}

fn converged(solver: &RustSuccessiveOrdersSolver) -> bool {
    solver.solver.diagnostics().converged
}

fn initial_residual(solver: &RustSuccessiveOrdersSolver) -> f64 {
    solver.solver.diagnostics().initial_residual
}

fn final_residual(solver: &RustSuccessiveOrdersSolver) -> f64 {
    solver.solver.diagnostics().final_residual
}

fn transport_workspace_bytes(solver: &RustSuccessiveOrdersSolver) -> usize {
    let element_size = std::mem::size_of::<f64>();
    (solver.transport_values.capacity()
        + solver.transport_value_tangent.capacity()
        + solver.boundary_scattering_values.capacity()
        + solver.boundary_scattering_value_tangent.capacity()
        + solver.layer_optical_depth.capacity()
        + solver.layer_attenuation.capacity()
        + solver.layer_prefix_attenuation.capacity()
        + solver.ray_end_attenuation.capacity()
        + solver.ray_end_attenuation_tangent.capacity()
        + solver.ray_end_attenuation_cotangent.capacity()
        + solver.first_order_forcing.capacity()
        + solver.first_order_forcing_tangent.capacity())
        * element_size
}
