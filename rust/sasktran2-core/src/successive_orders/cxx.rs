use std::pin::Pin;

use anyhow::{Result, anyhow};

use super::{
    BlockDiagonalMatrix, CsrMatrix, FixedPointProblem, SolverConfig, SuccessiveOrdersSolver,
};
use super::{
    ScalarCoefficientBasis, ScalarCoefficientScattering, VectorCoefficientBasis,
    VectorCoefficientScattering,
};

#[cxx::bridge(namespace = "sasktran2::rust::successive_orders")]
pub mod ffi {
    extern "Rust" {
        type RustSuccessiveOrdersSolver;

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
    let gradient = this.solver.solve_vjp(
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
