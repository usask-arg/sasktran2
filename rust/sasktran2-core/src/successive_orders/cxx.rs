use std::pin::Pin;

use anyhow::{Result, anyhow};

use super::{
    BlockDiagonalMatrix, CsrMatrix, FixedPointProblem, ScalarRayTransport, SolverConfig,
    SuccessiveOrdersSolver,
};
use super::{
    ScalarCoefficientBasis, ScalarCoefficientScattering, VectorCoefficientBasis,
    VectorCoefficientScattering,
};

#[cxx::bridge(namespace = "sasktran2::rust::successive_orders")]
pub mod ffi {
    extern "Rust" {
        type RustScalarRayTransport;
        type RustSuccessiveOrdersSolver;

        #[allow(clippy::too_many_arguments)]
        fn new_scalar_ray_transport(
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
            num_phase_moments: usize,
            solar_transmission_on_atmosphere_grid: bool,
            layer_source_offsets: &[u32],
            ray_transport_value_offsets: &[u32],
            source_value_inner_indices: &[u16],
            source_weights: &[f64],
            ray_ground_offsets: &[u32],
            ground_value_inner_indices: &[u16],
            ground_weights: &[f64],
        ) -> Result<Box<RustScalarRayTransport>>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_scalar_ray_transport(
            transport: &RustScalarRayTransport,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            transport_values: &mut [f64],
            layer_optical_depth: &mut [f64],
            layer_attenuation: &mut [f64],
            layer_prefix_attenuation: &mut [f64],
            ray_end_attenuation: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_scalar_ray_transport_with_first_order(
            transport: &RustScalarRayTransport,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            transport_values: &mut [f64],
            first_order_radiance: &mut [f64],
            layer_optical_depth: &mut [f64],
            layer_attenuation: &mut [f64],
            layer_prefix_attenuation: &mut [f64],
            ray_end_attenuation: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_scalar_ray_transport_jvp(
            transport: &RustScalarRayTransport,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            extinction_tangent: &[f64],
            single_scatter_albedo_tangent: &[f64],
            legendre_coefficient_tangent: &[f64],
            solar_transmission_tangent: &[f64],
            transport_values: &mut [f64],
            transport_value_tangent: &mut [f64],
            first_order_radiance_tangent: &mut [f64],
            ray_end_attenuation: &mut [f64],
            ray_end_attenuation_tangent: &mut [f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn assemble_scalar_ray_transport_vjp(
            transport: &RustScalarRayTransport,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            maximum_order: &[i32],
            solar_transmission: &[f64],
            transport_value_gradient: &[f64],
            transport_column_indices: &[i32],
            solution: &[f64],
            first_order_radiance_gradient: &[f64],
            ray_end_attenuation_cotangent: &[f64],
            layer_optical_depth: &[f64],
            layer_attenuation: &[f64],
            layer_prefix_attenuation: &[f64],
            extinction_gradient: &mut [f64],
            single_scatter_albedo_gradient: &mut [f64],
            legendre_coefficient_gradient: &mut [f64],
            solar_transmission_gradient: &mut [f64],
        ) -> Result<()>;

        fn scalar_ray_transport_num_rays(transport: &RustScalarRayTransport) -> usize;
        fn scalar_ray_transport_num_layers(transport: &RustScalarRayTransport) -> usize;
        fn scalar_ray_transport_num_atmosphere_weights(transport: &RustScalarRayTransport)
        -> usize;
        fn scalar_ray_transport_num_source_weights(transport: &RustScalarRayTransport) -> usize;
        fn scalar_ray_transport_num_ground_weights(transport: &RustScalarRayTransport) -> usize;
        fn scalar_ray_transport_storage_bytes(transport: &RustScalarRayTransport) -> usize;

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
            scattering_coefficients: &[f64],
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
            scattering_coefficients: &[f64],
            boundary_scattering_values: &[f64],
            first_order_forcing: &[f64],
            initial_guess: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn linearize_coefficients_jvp(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            transport_row_offsets: &[i32],
            transport_column_indices: &[i32],
            transport_values: &[f64],
            transport_value_tangent: &[f64],
            scattering_coefficient_tangent: &[f64],
            boundary_scattering_value_tangent: &[f64],
            first_order_forcing: &[f64],
            first_order_forcing_tangent: &[f64],
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

        fn restore_coefficient_solution(
            solver: Pin<&mut RustSuccessiveOrdersSolver>,
            scattering_coefficients: &[f64],
            boundary_scattering_values: &[f64],
            solution: &[f64],
        ) -> Result<()>;

        fn solution(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn solution_jvp(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn transport_value_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn scattering_coefficient_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn boundary_scattering_value_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn first_order_forcing_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn residual_history(solver: &RustSuccessiveOrdersSolver) -> &[f64];
        fn iterations(solver: &RustSuccessiveOrdersSolver) -> usize;
        fn converged(solver: &RustSuccessiveOrdersSolver) -> bool;
        fn initial_residual(solver: &RustSuccessiveOrdersSolver) -> f64;
        fn final_residual(solver: &RustSuccessiveOrdersSolver) -> f64;
    }
}

pub struct RustScalarRayTransport {
    transport: ScalarRayTransport,
}

#[allow(clippy::too_many_arguments)]
fn new_scalar_ray_transport(
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
    num_phase_moments: usize,
    solar_transmission_on_atmosphere_grid: bool,
    layer_source_offsets: &[u32],
    ray_transport_value_offsets: &[u32],
    source_value_inner_indices: &[u16],
    source_weights: &[f64],
    ray_ground_offsets: &[u32],
    ground_value_inner_indices: &[u16],
    ground_weights: &[f64],
) -> Result<Box<RustScalarRayTransport>> {
    Ok(Box::new(RustScalarRayTransport {
        transport: ScalarRayTransport::new(
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
            num_phase_moments,
            solar_transmission_on_atmosphere_grid,
            layer_source_offsets,
            ray_transport_value_offsets,
            source_value_inner_indices,
            source_weights,
            ray_ground_offsets,
            ground_value_inner_indices,
            ground_weights,
        )?,
    }))
}

#[allow(clippy::too_many_arguments)]
fn assemble_scalar_ray_transport_with_first_order(
    transport: &RustScalarRayTransport,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    transport_values: &mut [f64],
    first_order_radiance: &mut [f64],
    layer_optical_depth: &mut [f64],
    layer_attenuation: &mut [f64],
    layer_prefix_attenuation: &mut [f64],
    ray_end_attenuation: &mut [f64],
) -> Result<()> {
    transport.transport.assemble_with_first_order(
        extinction,
        single_scatter_albedo,
        legendre_coefficients,
        maximum_order,
        solar_transmission,
        transport_values,
        first_order_radiance,
        layer_optical_depth,
        layer_attenuation,
        layer_prefix_attenuation,
        ray_end_attenuation,
    )
}

#[allow(clippy::too_many_arguments)]
fn assemble_scalar_ray_transport_jvp(
    transport: &RustScalarRayTransport,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    extinction_tangent: &[f64],
    single_scatter_albedo_tangent: &[f64],
    legendre_coefficient_tangent: &[f64],
    solar_transmission_tangent: &[f64],
    transport_values: &mut [f64],
    transport_value_tangent: &mut [f64],
    first_order_radiance_tangent: &mut [f64],
    ray_end_attenuation: &mut [f64],
    ray_end_attenuation_tangent: &mut [f64],
) -> Result<()> {
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
        transport_values,
        transport_value_tangent,
        first_order_radiance_tangent,
        ray_end_attenuation,
        ray_end_attenuation_tangent,
    )
}

#[allow(clippy::too_many_arguments)]
fn assemble_scalar_ray_transport_vjp(
    transport: &RustScalarRayTransport,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    maximum_order: &[i32],
    solar_transmission: &[f64],
    transport_value_gradient: &[f64],
    transport_column_indices: &[i32],
    solution: &[f64],
    first_order_radiance_gradient: &[f64],
    ray_end_attenuation_cotangent: &[f64],
    layer_optical_depth: &[f64],
    layer_attenuation: &[f64],
    layer_prefix_attenuation: &[f64],
    extinction_gradient: &mut [f64],
    single_scatter_albedo_gradient: &mut [f64],
    legendre_coefficient_gradient: &mut [f64],
    solar_transmission_gradient: &mut [f64],
) -> Result<()> {
    transport.transport.assemble_vjp_with_first_order(
        extinction,
        single_scatter_albedo,
        legendre_coefficients,
        maximum_order,
        solar_transmission,
        transport_value_gradient,
        transport_column_indices,
        solution,
        first_order_radiance_gradient,
        ray_end_attenuation_cotangent,
        layer_optical_depth,
        layer_attenuation,
        layer_prefix_attenuation,
        extinction_gradient,
        single_scatter_albedo_gradient,
        legendre_coefficient_gradient,
        solar_transmission_gradient,
    )
}

#[allow(clippy::too_many_arguments)]
fn assemble_scalar_ray_transport(
    transport: &RustScalarRayTransport,
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

fn scalar_ray_transport_num_rays(transport: &RustScalarRayTransport) -> usize {
    transport.transport.num_rays()
}

fn scalar_ray_transport_num_layers(transport: &RustScalarRayTransport) -> usize {
    transport.transport.num_layers()
}

fn scalar_ray_transport_num_atmosphere_weights(transport: &RustScalarRayTransport) -> usize {
    transport.transport.num_atmosphere_weights()
}

fn scalar_ray_transport_num_source_weights(transport: &RustScalarRayTransport) -> usize {
    transport.transport.num_source_weights()
}

fn scalar_ray_transport_num_ground_weights(transport: &RustScalarRayTransport) -> usize {
    transport.transport.num_ground_weights()
}

fn scalar_ray_transport_storage_bytes(transport: &RustScalarRayTransport) -> usize {
    transport.transport.storage_bytes()
}

pub struct RustSuccessiveOrdersSolver {
    solver: SuccessiveOrdersSolver,
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
            solution_jvp: Vec::new(),
            transport_value_gradient: Vec::new(),
            scattering_coefficient_gradient: Vec::new(),
            boundary_scattering_value_gradient: Vec::new(),
            first_order_forcing_gradient: Vec::new(),
        }
    }
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
    let scattering = VectorCoefficientScattering::new(
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
    scattering_coefficients: &[f64],
    boundary_scattering_values: &[f64],
    first_order_forcing: &[f64],
    initial_guess: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.solver
        .problem_mut()
        .set_scattering_coefficients(scattering_coefficients)?;
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
    scattering_coefficients: &[f64],
    boundary_scattering_values: &[f64],
    first_order_forcing: &[f64],
    initial_guess: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.solver
        .problem_mut()
        .set_scattering_coefficients(scattering_coefficients)?;
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
fn linearize_coefficients_jvp(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    transport_row_offsets: &[i32],
    transport_column_indices: &[i32],
    transport_values: &[f64],
    transport_value_tangent: &[f64],
    scattering_coefficient_tangent: &[f64],
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
        scattering_coefficient_tangent,
        boundary_scattering_value_tangent,
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

fn restore_coefficient_solution(
    mut solver: Pin<&mut RustSuccessiveOrdersSolver>,
    scattering_coefficients: &[f64],
    boundary_scattering_values: &[f64],
    solution: &[f64],
) -> Result<()> {
    let this = solver.as_mut().get_mut();
    this.solver
        .problem_mut()
        .set_scattering_coefficients(scattering_coefficients)?;
    this.solver
        .problem_mut()
        .set_scattering_values(boundary_scattering_values)?;
    this.solver.restore_converged_solution(solution)?;
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

fn solution_jvp(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.solution_jvp
}

fn transport_value_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.transport_value_gradient
}

fn scattering_coefficient_gradient(solver: &RustSuccessiveOrdersSolver) -> &[f64] {
    &solver.scattering_coefficient_gradient
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
