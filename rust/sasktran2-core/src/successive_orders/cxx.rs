use std::pin::Pin;
use std::sync::Arc;

use anyhow::{Result, anyhow};

use super::{
    BlockDiagonalMatrix, CsrMatrix, FixedPointProblem, ScalarRayTransport, SolverConfig,
    SuccessiveOrdersSolver,
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

        fn calculate_solar_transmission_jvp(
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

        fn boundary_scattering_values_mut(
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
    )
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

fn boundary_scattering_values_mut(solver: Pin<&mut RustSuccessiveOrdersSolver>) -> &mut [f64] {
    &mut solver.get_mut().boundary_scattering_values
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
