use std::fmt::{Display, Formatter};

use super::{
    FixedPointProblem, OperatorError,
    simd::batch::{
        Batch4, LANES, dot_product as batch_dot_product, infinity_norm as batch_infinity_norm,
    },
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverConfig {
    pub max_iterations: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
    /// Zero selects damped Picard iteration. Positive values select Anderson
    /// acceleration with this history depth.
    pub anderson_depth: usize,
    pub damping: f64,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            max_iterations: 50,
            relative_tolerance: 1.0e-6,
            absolute_tolerance: 1.0e-12,
            anderson_depth: 0,
            damping: 1.0,
        }
    }
}

impl SolverConfig {
    fn validate(self) -> Result<Self, SolverError> {
        if self.max_iterations == 0
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance < 0.0
            || !self.absolute_tolerance.is_finite()
            || self.absolute_tolerance < 0.0
            || !self.damping.is_finite()
            || self.damping <= 0.0
            || self.damping > 1.0
        {
            return Err(SolverError::InvalidConfiguration);
        }
        Ok(self)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvergenceReason {
    Tolerance,
    MaximumIterations,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SolverDiagnostics {
    pub iterations: usize,
    pub converged: bool,
    pub reason: ConvergenceReason,
    pub initial_residual: f64,
    pub final_residual: f64,
    pub residual_history: Vec<f64>,
}

/// Reusable type-II Anderson history. Difference vectors are retained directly,
/// and their Gram matrix is extended by one row and column for each iteration
/// instead of being rebuilt from the full state history.
#[derive(Debug, Clone)]
struct AndersonWorkspace {
    depth: usize,
    vector_size: usize,
    previous_state: Vec<f64>,
    previous_residual: Vec<f64>,
    has_previous: bool,
    delta_states: Vec<Vec<f64>>,
    delta_residuals: Vec<Vec<f64>>,
    gram: Vec<f64>,
    right: Vec<f64>,
    normal_scratch: Vec<f64>,
    solution_scratch: Vec<f64>,
}

impl AndersonWorkspace {
    fn new(depth: usize, vector_size: usize) -> Self {
        let retained_vector_size = if depth == 0 { 0 } else { vector_size };
        Self {
            depth,
            vector_size,
            previous_state: vec![0.0; retained_vector_size],
            previous_residual: vec![0.0; retained_vector_size],
            has_previous: false,
            delta_states: Vec::with_capacity(depth),
            delta_residuals: Vec::with_capacity(depth),
            gram: vec![0.0; depth * depth],
            right: vec![0.0; depth],
            normal_scratch: vec![0.0; depth * depth],
            solution_scratch: vec![0.0; depth],
        }
    }

    fn reset(&mut self) {
        self.has_previous = false;
        self.delta_states.clear();
        self.delta_residuals.clear();
        self.gram.fill(0.0);
        self.right.fill(0.0);
    }

    /// Records a new `(state, residual)` pair and returns whether at least one
    /// usable difference is available.
    fn record(&mut self, state: &[f64], residual: &[f64]) -> bool {
        debug_assert_eq!(state.len(), self.vector_size);
        debug_assert_eq!(residual.len(), self.vector_size);
        if self.depth == 0 {
            return false;
        }
        if !self.has_previous {
            self.previous_state.copy_from_slice(state);
            self.previous_residual.copy_from_slice(residual);
            self.has_previous = true;
            return false;
        }

        let (mut delta_state, mut delta_residual) = if self.delta_states.len() == self.depth {
            let delta_state = self.delta_states.remove(0);
            let delta_residual = self.delta_residuals.remove(0);
            for row in 0..self.depth - 1 {
                self.right[row] = self.right[row + 1];
                for column in 0..self.depth - 1 {
                    self.gram[row * self.depth + column] =
                        self.gram[(row + 1) * self.depth + column + 1];
                }
            }
            (delta_state, delta_residual)
        } else {
            (vec![0.0; self.vector_size], vec![0.0; self.vector_size])
        };
        for index in 0..self.vector_size {
            delta_state[index] = state[index] - self.previous_state[index];
            delta_residual[index] = residual[index] - self.previous_residual[index];
        }

        let new_index = self.delta_residuals.len();
        for index in 0..new_index {
            let cross = dot_product(&self.delta_residuals[index], &delta_residual);
            self.gram[index * self.depth + new_index] = cross;
            self.gram[new_index * self.depth + index] = cross;
            // The current residual is the previous residual plus the newly
            // appended residual difference.
            self.right[index] += cross;
        }
        self.gram[new_index * self.depth + new_index] =
            dot_product(&delta_residual, &delta_residual);
        self.right[new_index] = dot_product(&delta_residual, residual);
        self.delta_states.push(delta_state);
        self.delta_residuals.push(delta_residual);
        self.previous_state.copy_from_slice(state);
        self.previous_residual.copy_from_slice(residual);
        true
    }

    fn update(&mut self, state: &mut [f64], mapped: &mut [f64], damping: f64) -> bool {
        let differences = self.delta_states.len();
        if differences == 0 {
            return false;
        }
        for row in 0..differences {
            self.solution_scratch[row] = self.right[row];
            for column in 0..differences {
                self.normal_scratch[row * differences + column] =
                    self.gram[row * self.depth + column];
            }
            self.normal_scratch[row * differences + row] += 1.0e-14;
        }
        if !solve_dense(
            &mut self.normal_scratch[..differences * differences],
            &mut self.solution_scratch[..differences],
            differences,
        ) {
            return false;
        }

        for row in 0..self.vector_size {
            let mut candidate = mapped[row];
            for column in 0..differences {
                candidate -= self.solution_scratch[column]
                    * (self.delta_states[column][row] + self.delta_residuals[column][row]);
            }
            mapped[row] = (1.0 - damping) * state[row] + damping * candidate;
        }
        if mapped.iter().any(|value| !value.is_finite()) {
            return false;
        }
        state.copy_from_slice(mapped);
        true
    }
}

impl Default for SolverDiagnostics {
    fn default() -> Self {
        Self {
            iterations: 0,
            converged: false,
            reason: ConvergenceReason::MaximumIterations,
            initial_residual: f64::NAN,
            final_residual: f64::NAN,
            residual_history: Vec::new(),
        }
    }
}

#[derive(Debug)]
pub enum SolverError {
    InvalidConfiguration,
    InitialGuessSize,
    LinearizationBeforeSolve,
    ImplicitLinearSolveDidNotConverge,
    NonFiniteIteration,
    Operator(OperatorError),
}

impl Display for SolverError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidConfiguration => {
                formatter.write_str("invalid successive-orders solver configuration")
            }
            Self::InitialGuessSize => {
                formatter.write_str("successive-orders initial guess has the wrong size")
            }
            Self::LinearizationBeforeSolve => {
                formatter.write_str("successive-orders must be solved before linearization")
            }
            Self::ImplicitLinearSolveDidNotConverge => {
                formatter.write_str("implicit successive-orders derivative solve did not converge")
            }
            Self::NonFiniteIteration => {
                formatter.write_str("successive-orders iteration produced a non-finite value")
            }
            Self::Operator(error) => Display::fmt(error, formatter),
        }
    }
}

impl std::error::Error for SolverError {}

impl From<OperatorError> for SolverError {
    fn from(value: OperatorError) -> Self {
        Self::Operator(value)
    }
}

#[derive(Debug, Clone)]
pub struct SuccessiveOrdersSolver {
    problem: FixedPointProblem,
    config: SolverConfig,
    solution: Vec<f64>,
    mapped: Vec<f64>,
    residual: Vec<f64>,
    incoming: Vec<f64>,
    initial_state: Vec<f64>,
    anderson: AndersonWorkspace,
    diagnostics: SolverDiagnostics,
}

#[derive(Debug, Clone)]
pub(crate) struct BatchSolveResult {
    pub(crate) solution: Vec<f64>,
    pub(crate) diagnostics: [SolverDiagnostics; LANES],
}

#[derive(Debug, Clone)]
pub(crate) struct BatchFixedPointGradient {
    pub(crate) scattering_coefficients: Vec<f64>,
    #[cfg_attr(test, allow(dead_code))]
    pub(crate) dense_scattering_values: Vec<f64>,
    pub(crate) forcing: Vec<f64>,
}

impl SuccessiveOrdersSolver {
    pub fn new(problem: FixedPointProblem, config: SolverConfig) -> Result<Self, SolverError> {
        let config = config.validate()?;
        let state_size = problem.state_size();
        let incoming_size = problem.incoming_size();
        Ok(Self {
            problem,
            config,
            solution: vec![0.0; state_size],
            mapped: vec![0.0; state_size],
            residual: vec![0.0; state_size],
            incoming: vec![0.0; incoming_size],
            initial_state: vec![0.0; state_size],
            anderson: AndersonWorkspace::new(config.anderson_depth, state_size),
            diagnostics: SolverDiagnostics::default(),
        })
    }

    #[inline]
    pub fn problem_mut(&mut self) -> &mut FixedPointProblem {
        &mut self.problem
    }

    #[inline]
    pub fn solution(&self) -> &[f64] {
        &self.solution
    }

    #[inline]
    pub fn diagnostics(&self) -> &SolverDiagnostics {
        &self.diagnostics
    }

    /// Restores the sufficient primal state for implicit differentiation of a
    /// converged fixed point. Anderson history is intentionally not retained:
    /// at convergence derivative products depend only on the fixed-point
    /// solution and operator.
    pub fn restore_converged_solution(&mut self, solution: &[f64]) -> Result<(), SolverError> {
        if solution.len() != self.solution.len() {
            return Err(SolverError::InitialGuessSize);
        }
        if solution.iter().any(|value| !value.is_finite()) {
            return Err(SolverError::NonFiniteIteration);
        }
        self.solution.copy_from_slice(solution);
        self.initial_state.copy_from_slice(solution);
        self.anderson.reset();
        self.diagnostics = SolverDiagnostics {
            iterations: 1,
            converged: true,
            reason: ConvergenceReason::Tolerance,
            ..SolverDiagnostics::default()
        };
        Ok(())
    }

    pub fn solve(
        &mut self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
        initial_guess: Option<&[f64]>,
    ) -> Result<&[f64], SolverError> {
        self.problem.validate_iteration_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
        )?;
        if let Some(initial_guess) = initial_guess {
            if initial_guess.len() != self.solution.len() {
                return Err(SolverError::InitialGuessSize);
            }
            self.solution.copy_from_slice(initial_guess);
        } else {
            self.solution.fill(0.0);
        }
        if self.solution.iter().any(|value| !value.is_finite()) {
            return Err(SolverError::NonFiniteIteration);
        }
        self.initial_state.copy_from_slice(&self.solution);

        self.anderson.reset();
        self.diagnostics = SolverDiagnostics::default();
        self.diagnostics
            .residual_history
            .reserve(self.config.max_iterations);

        // The fixed-point map is affine, so a converged solve is equivalently
        // `(I - S T) x = S b`. Try the bounded-memory Krylov solve first for
        // accelerated configurations, but accept it only after evaluating the
        // true fixed-point residual. Any breakdown or residual mismatch falls
        // back to the established Anderson path below.
        if self.config.anderson_depth > 0
            && self.config.max_iterations >= 2
            && (self.config.relative_tolerance > 0.0 || self.config.absolute_tolerance > 0.0)
            && self.try_solve_krylov(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                forcing,
            )?
        {
            return Ok(&self.solution);
        }
        self.diagnostics = SolverDiagnostics::default();
        self.diagnostics
            .residual_history
            .reserve(self.config.max_iterations);

        for iteration in 0..self.config.max_iterations {
            self.problem.apply(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                forcing,
                &self.solution,
                &mut self.mapped,
                &mut self.incoming,
            )?;
            for index in 0..self.solution.len() {
                self.residual[index] = self.mapped[index] - self.solution[index];
            }
            if self.mapped.iter().any(|value| !value.is_finite())
                || self.residual.iter().any(|value| !value.is_finite())
            {
                return Err(SolverError::NonFiniteIteration);
            }

            let residual_norm = infinity_norm(&self.residual);
            let scale = infinity_norm(&self.mapped).max(infinity_norm(&self.solution));
            let threshold = self.config.absolute_tolerance + self.config.relative_tolerance * scale;
            if iteration == 0 {
                self.diagnostics.initial_residual = residual_norm;
            }
            self.diagnostics.residual_history.push(residual_norm);
            self.diagnostics.iterations = iteration + 1;
            self.diagnostics.final_residual = residual_norm;

            if residual_norm <= threshold {
                self.solution.copy_from_slice(&self.mapped);
                self.diagnostics.converged = true;
                self.diagnostics.reason = ConvergenceReason::Tolerance;
                return Ok(&self.solution);
            }

            let accelerated = if self.config.anderson_depth > 0 {
                self.anderson.record(&self.solution, &self.residual)
                    && self.anderson.update(
                        &mut self.solution,
                        &mut self.mapped,
                        self.config.damping,
                    )
            } else {
                false
            };
            if !accelerated {
                damped_update(&mut self.solution, &self.residual, self.config.damping);
            }
        }

        self.diagnostics.reason = ConvergenceReason::MaximumIterations;
        Ok(&self.solution)
    }

    /// Solves four independent scalar wavelength systems in lockstep. The
    /// immutable transport sparsity and angular scattering basis are shared;
    /// every wavelength-dependent array is interleaved with four contiguous
    /// lanes. A caller should fall back to the scalar solver if this bounded
    /// Krylov path reports a breakdown.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_batch4(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        scattering_coefficients: &[f64],
        dense_scattering_values: &[f64],
        forcing: &[f64],
        initial_guess: Option<&[f64]>,
    ) -> Result<BatchSolveResult, SolverError> {
        if self.config.anderson_depth == 0
            || self.config.max_iterations < 2
            || (self.config.relative_tolerance == 0.0 && self.config.absolute_tolerance == 0.0)
        {
            return Err(SolverError::ImplicitLinearSolveDidNotConverge);
        }
        let active_num_coefficients = self.problem.validate_batch4_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            scattering_coefficients,
            dense_scattering_values,
            forcing,
        )?;
        let state_size = self.problem.state_size();
        let packed_state_size = state_size * LANES;
        if initial_guess.is_some_and(|values| values.len() != packed_state_size) {
            return Err(SolverError::InitialGuessSize);
        }

        let mut solution = initial_guess
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| vec![0.0; packed_state_size]);
        if solution.iter().any(|value| !value.is_finite()) {
            return Err(SolverError::NonFiniteIteration);
        }
        let mut mapped = vec![0.0; packed_state_size];
        let mut residual = vec![0.0; packed_state_size];
        let mut incoming = vec![0.0; self.problem.incoming_size() * LANES];
        let mut scattering_scratch = vec![
            0.0;
            self.problem.batch4_scattering_scratch_size(
                active_num_coefficients
            )?
        ];
        self.problem.apply_batch4(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            scattering_coefficients,
            dense_scattering_values,
            forcing,
            active_num_coefficients,
            &solution,
            &mut mapped,
            &mut incoming,
            &mut scattering_scratch,
        )?;
        for element in 0..state_size {
            (Batch4::load(&mapped, element) - Batch4::load(&solution, element))
                .store(&mut residual, element);
        }
        if mapped.iter().any(|value| !value.is_finite())
            || residual.iter().any(|value| !value.is_finite())
        {
            return Err(SolverError::NonFiniteIteration);
        }

        let initial_residual_norm = batch_infinity_norm(&residual);
        let initial_scale = batch_infinity_norm(&mapped);
        let solution_scale = batch_infinity_norm(&solution);
        let initial_threshold: [f64; LANES] = std::array::from_fn(|lane| {
            self.config.absolute_tolerance
                + self.config.relative_tolerance * initial_scale[lane].max(solution_scale[lane])
        });
        let mut active =
            std::array::from_fn(|lane| initial_residual_norm[lane] > initial_threshold[lane]);
        for (lane, &is_active) in active.iter().enumerate() {
            if !is_active {
                copy_batch_lane(&mapped, &mut solution, lane);
            }
        }
        let mut residual_history: [Vec<f64>; LANES] =
            std::array::from_fn(|lane| vec![initial_residual_norm[lane]]);
        if !active.iter().any(|&value| value) {
            return Ok(completed_batch_result(
                solution,
                initial_residual_norm,
                residual_history,
            ));
        }
        zero_inactive_lanes(&mut residual, active);

        let shadow_residual = residual.clone();
        let mut direction = vec![0.0; packed_state_size];
        let mut direction_image = vec![0.0; packed_state_size];
        let mut intermediate_residual = vec![0.0; packed_state_size];
        let mut intermediate_image = vec![0.0; packed_state_size];
        let mut previous_rho = Batch4::splat(1.0);
        let mut alpha = Batch4::splat(1.0);
        let mut omega = Batch4::splat(1.0);

        for iteration in 0..self.config.max_iterations {
            let rho = batch_dot_product(&shadow_residual, &residual);
            validate_batch_divisor(rho, active)?;
            if iteration == 0 {
                direction.copy_from_slice(&residual);
            } else {
                validate_batch_divisor(previous_rho, active)?;
                validate_batch_divisor(omega, active)?;
                let beta = safe_batch_divide(rho, previous_rho, active)
                    * safe_batch_divide(alpha, omega, active);
                for element in 0..state_size {
                    let value = Batch4::load(&residual, element)
                        + beta
                            * (Batch4::load(&direction, element)
                                - omega * Batch4::load(&direction_image, element));
                    Batch4::select(active, value, Batch4::splat(0.0))
                        .store(&mut direction, element);
                }
            }

            self.problem.apply_linear_system_batch4(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                scattering_coefficients,
                dense_scattering_values,
                active_num_coefficients,
                &direction,
                &mut direction_image,
                &mut incoming,
                &mut scattering_scratch,
            )?;
            let alpha_denominator = batch_dot_product(&shadow_residual, &direction_image);
            validate_batch_divisor(alpha_denominator, active)?;
            alpha = safe_batch_divide(rho, alpha_denominator, active);
            for element in 0..state_size {
                (Batch4::load(&residual, element)
                    - alpha * Batch4::load(&direction_image, element))
                .store(&mut intermediate_residual, element);
            }

            let intermediate_norm = batch_infinity_norm(&intermediate_residual);
            let solution_norm = batch_infinity_norm(&solution);
            let threshold: [f64; LANES] = std::array::from_fn(|lane| {
                self.config.absolute_tolerance
                    + self.config.relative_tolerance
                        * initial_residual_norm[lane].max(solution_norm[lane])
            });
            let intermediate_converged = std::array::from_fn(|lane| {
                active[lane] && intermediate_norm[lane] <= threshold[lane]
            });
            if intermediate_converged.iter().any(|&value| value) {
                for element in 0..state_size {
                    let candidate = Batch4::load(&solution, element)
                        + alpha * Batch4::load(&direction, element);
                    Batch4::select(
                        intermediate_converged,
                        candidate,
                        Batch4::load(&solution, element),
                    )
                    .store(&mut solution, element);
                }
                for lane in 0..LANES {
                    if intermediate_converged[lane] {
                        residual_history[lane].push(intermediate_norm[lane]);
                        active[lane] = false;
                    }
                }
                if !active.iter().any(|&value| value) {
                    break;
                }
            }
            zero_inactive_lanes(&mut intermediate_residual, active);

            self.problem.apply_linear_system_batch4(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                scattering_coefficients,
                dense_scattering_values,
                active_num_coefficients,
                &intermediate_residual,
                &mut intermediate_image,
                &mut incoming,
                &mut scattering_scratch,
            )?;
            let omega_denominator = batch_dot_product(&intermediate_image, &intermediate_image);
            validate_batch_positive_divisor(omega_denominator, active)?;
            omega = safe_batch_divide(
                batch_dot_product(&intermediate_image, &intermediate_residual),
                omega_denominator,
                active,
            );
            validate_batch_divisor(omega, active)?;
            for element in 0..state_size {
                let updated_solution = Batch4::load(&solution, element)
                    + alpha * Batch4::load(&direction, element)
                    + omega * Batch4::load(&intermediate_residual, element);
                let updated_residual = Batch4::load(&intermediate_residual, element)
                    - omega * Batch4::load(&intermediate_image, element);
                Batch4::select(active, updated_solution, Batch4::load(&solution, element))
                    .store(&mut solution, element);
                Batch4::select(active, updated_residual, Batch4::splat(0.0))
                    .store(&mut residual, element);
            }
            if solution.iter().any(|value| !value.is_finite())
                || residual.iter().any(|value| !value.is_finite())
            {
                return Err(SolverError::NonFiniteIteration);
            }
            let residual_norm = batch_infinity_norm(&residual);
            let solution_norm = batch_infinity_norm(&solution);
            for lane in 0..LANES {
                if !active[lane] {
                    continue;
                }
                residual_history[lane].push(residual_norm[lane]);
                let threshold = self.config.absolute_tolerance
                    + self.config.relative_tolerance
                        * initial_residual_norm[lane].max(solution_norm[lane]);
                if residual_norm[lane] <= threshold {
                    active[lane] = false;
                }
            }
            if !active.iter().any(|&value| value) {
                break;
            }
            previous_rho = rho;
        }
        if active.iter().any(|&value| value) {
            return Err(SolverError::ImplicitLinearSolveDidNotConverge);
        }

        // Match the scalar path: accept a recurrence solution only after an
        // explicit affine-map residual check, and retain the mapped value.
        self.problem.apply_batch4(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            scattering_coefficients,
            dense_scattering_values,
            forcing,
            active_num_coefficients,
            &solution,
            &mut mapped,
            &mut incoming,
            &mut scattering_scratch,
        )?;
        for element in 0..state_size {
            (Batch4::load(&mapped, element) - Batch4::load(&solution, element))
                .store(&mut residual, element);
        }
        let final_residual = batch_infinity_norm(&residual);
        let mapped_norm = batch_infinity_norm(&mapped);
        let solution_norm = batch_infinity_norm(&solution);
        for lane in 0..LANES {
            let threshold = self.config.absolute_tolerance
                + self.config.relative_tolerance * mapped_norm[lane].max(solution_norm[lane]);
            if !final_residual[lane].is_finite() || final_residual[lane] > threshold {
                return Err(SolverError::ImplicitLinearSolveDidNotConverge);
            }
            residual_history[lane].push(final_residual[lane]);
        }
        Ok(completed_batch_result(
            mapped,
            initial_residual_norm,
            residual_history,
        ))
    }

    /// Applies four independent implicit JVP systems in lockstep. Parameter
    /// tangents are assembled one lane at a time because they originate in
    /// the engine's native derivative layout; the repeated Krylov operator,
    /// which dominates the coefficient solve, remains fully packed.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_jvp_batch4(
        &mut self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        transport_value_tangent: &[f64],
        scattering_coefficients: &[f64],
        scattering_coefficient_tangent: &[f64],
        dense_scattering_values: &[f64],
        dense_scattering_value_tangent: &[f64],
        forcing: &[f64],
        forcing_tangent: &[f64],
        solution: &[f64],
    ) -> Result<Vec<f64>, SolverError> {
        let active_num_coefficients = self.problem.validate_batch4_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            scattering_coefficients,
            dense_scattering_values,
            forcing,
        )?;
        if transport_value_tangent.len() != self.problem.transport_value_size() * LANES
            || scattering_coefficient_tangent.len()
                != self.problem.scattering_coefficient_size() * LANES
            || dense_scattering_value_tangent.len()
                != self.problem.dense_scattering_value_size() * LANES
            || forcing_tangent.len() != self.problem.incoming_size() * LANES
            || solution.len() != self.problem.state_size() * LANES
        {
            return Err(SolverError::Operator(OperatorError::DimensionMismatch));
        }

        let transport_size = self.problem.transport_value_size();
        let coefficient_size = self.problem.scattering_coefficient_size();
        let dense_size = self.problem.dense_scattering_value_size();
        let forcing_size = self.problem.incoming_size();
        let state_size = self.problem.state_size();
        let mut lane_transport = vec![0.0; transport_size];
        let mut lane_transport_tangent = vec![0.0; transport_size];
        let mut lane_coefficients = vec![0.0; coefficient_size];
        let mut lane_coefficient_tangent = vec![0.0; coefficient_size];
        let mut lane_dense = vec![0.0; dense_size];
        let mut lane_dense_tangent = vec![0.0; dense_size];
        let mut lane_forcing = vec![0.0; forcing_size];
        let mut lane_forcing_tangent = vec![0.0; forcing_size];
        let mut lane_solution = vec![0.0; state_size];
        let mut lane_right_hand_side = vec![0.0; state_size];
        let mut incoming = vec![0.0; forcing_size];
        let mut incoming_tangent = vec![0.0; forcing_size];
        let mut right_hand_side = vec![0.0; state_size * LANES];
        for lane in 0..LANES {
            extract_batch_lane(transport_values, lane, &mut lane_transport);
            extract_batch_lane(transport_value_tangent, lane, &mut lane_transport_tangent);
            extract_batch_lane(scattering_coefficients, lane, &mut lane_coefficients);
            extract_batch_lane(
                scattering_coefficient_tangent,
                lane,
                &mut lane_coefficient_tangent,
            );
            extract_batch_lane(dense_scattering_values, lane, &mut lane_dense);
            extract_batch_lane(
                dense_scattering_value_tangent,
                lane,
                &mut lane_dense_tangent,
            );
            extract_batch_lane(forcing, lane, &mut lane_forcing);
            extract_batch_lane(forcing_tangent, lane, &mut lane_forcing_tangent);
            extract_batch_lane(solution, lane, &mut lane_solution);
            self.problem
                .set_scattering_coefficients(&lane_coefficients)?;
            self.problem.set_scattering_values(&lane_dense)?;
            self.problem.apply_parameter_jvp(
                transport_row_offsets,
                transport_column_indices,
                &lane_transport,
                &lane_transport_tangent,
                &lane_forcing,
                &lane_forcing_tangent,
                &lane_coefficient_tangent,
                &lane_dense_tangent,
                &lane_solution,
                &mut lane_right_hand_side,
                &mut incoming,
                &mut incoming_tangent,
            )?;
            interleave_batch_lane(&lane_right_hand_side, lane, &mut right_hand_side);
        }

        let mut packed_incoming = vec![0.0; forcing_size * LANES];
        let mut scattering_scratch = vec![
            0.0;
            self.problem.batch4_scattering_scratch_size(
                active_num_coefficients
            )?
        ];
        let problem = &self.problem;
        let config = self.config;
        solve_bicgstab_batch4(&right_hand_side, config, |input, output| {
            problem.apply_linear_system_batch4(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                scattering_coefficients,
                dense_scattering_values,
                active_num_coefficients,
                input,
                output,
                &mut packed_incoming,
                &mut scattering_scratch,
            )?;
            Ok(())
        })
    }

    /// Applies four independent converged fixed-point adjoints in lockstep.
    /// The compact gradient intentionally omits transport values; the ray
    /// geometry pullback reconstructs those factors from the packed solution
    /// and forcing cotangent, matching the scalar VJP path.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_vjp_batch4(
        &mut self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        scattering_coefficients: &[f64],
        dense_scattering_values: &[f64],
        forcing: &[f64],
        solution: &[f64],
        solution_cotangent: &[f64],
    ) -> Result<BatchFixedPointGradient, SolverError> {
        let active_num_coefficients = self.problem.validate_batch4_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            scattering_coefficients,
            dense_scattering_values,
            forcing,
        )?;
        if solution.len() != self.problem.state_size() * LANES
            || solution_cotangent.len() != self.problem.state_size() * LANES
        {
            return Err(SolverError::Operator(OperatorError::DimensionMismatch));
        }
        let mut packed_incoming = vec![0.0; self.problem.incoming_size() * LANES];
        let mut scattering_scratch = vec![
            0.0;
            self.problem.batch4_scattering_scratch_size(
                active_num_coefficients
            )?
        ];
        let problem = &self.problem;
        let config = self.config;
        let bicgstab_result = solve_bicgstab_batch4(solution_cotangent, config, |input, output| {
            problem.apply_linear_system_transpose_batch4(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                scattering_coefficients,
                dense_scattering_values,
                active_num_coefficients,
                input,
                output,
                &mut packed_incoming,
                &mut scattering_scratch,
            )?;
            Ok(())
        });
        let adjoint = match bicgstab_result {
            Ok(adjoint) => adjoint,
            Err(SolverError::ImplicitLinearSolveDidNotConverge) => {
                solve_affine_fixed_point_batch4(solution_cotangent, config, |input, output| {
                    problem.apply_linear_transpose_batch4(
                        transport_row_offsets,
                        transport_column_indices,
                        transport_values,
                        scattering_coefficients,
                        dense_scattering_values,
                        active_num_coefficients,
                        input,
                        output,
                        &mut packed_incoming,
                        &mut scattering_scratch,
                    )?;
                    Ok(())
                })?
            }
            Err(error) => return Err(error),
        };

        let transport_size = self.problem.transport_value_size();
        let coefficient_size = self.problem.scattering_coefficient_size();
        let dense_size = self.problem.dense_scattering_value_size();
        let forcing_size = self.problem.incoming_size();
        let state_size = self.problem.state_size();
        let mut lane_transport = vec![0.0; transport_size];
        let mut lane_coefficients = vec![0.0; coefficient_size];
        let mut lane_dense = vec![0.0; dense_size];
        let mut lane_forcing = vec![0.0; forcing_size];
        let mut lane_solution = vec![0.0; state_size];
        let mut lane_adjoint = vec![0.0; state_size];
        let mut coefficient_gradient = vec![0.0; coefficient_size * LANES];
        let mut dense_gradient = vec![0.0; dense_size * LANES];
        let mut forcing_gradient = vec![0.0; forcing_size * LANES];
        let mut lane_coefficient_gradient = vec![0.0; coefficient_size];
        let mut lane_dense_gradient = vec![0.0; dense_size];
        let mut lane_forcing_gradient = vec![0.0; forcing_size];
        let mut incoming = vec![0.0; forcing_size];
        let mut incoming_cotangent = vec![0.0; forcing_size];
        for lane in 0..LANES {
            extract_batch_lane(transport_values, lane, &mut lane_transport);
            extract_batch_lane(scattering_coefficients, lane, &mut lane_coefficients);
            extract_batch_lane(dense_scattering_values, lane, &mut lane_dense);
            extract_batch_lane(forcing, lane, &mut lane_forcing);
            extract_batch_lane(solution, lane, &mut lane_solution);
            extract_batch_lane(&adjoint, lane, &mut lane_adjoint);
            self.problem
                .set_scattering_coefficients(&lane_coefficients)?;
            self.problem.set_scattering_values(&lane_dense)?;
            self.problem.apply_parameter_vjp(
                transport_row_offsets,
                transport_column_indices,
                &lane_transport,
                &lane_forcing,
                &lane_solution,
                &lane_adjoint,
                &mut [],
                &mut lane_coefficient_gradient,
                &mut lane_dense_gradient,
                &mut lane_forcing_gradient,
                &mut incoming,
                &mut incoming_cotangent,
            )?;
            interleave_batch_lane(&lane_coefficient_gradient, lane, &mut coefficient_gradient);
            interleave_batch_lane(&lane_dense_gradient, lane, &mut dense_gradient);
            interleave_batch_lane(&lane_forcing_gradient, lane, &mut forcing_gradient);
        }
        Ok(BatchFixedPointGradient {
            scattering_coefficients: coefficient_gradient,
            dense_scattering_values: dense_gradient,
            forcing: forcing_gradient,
        })
    }

    fn try_solve_krylov(
        &mut self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
    ) -> Result<bool, SolverError> {
        self.problem.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
            &self.solution,
            &mut self.mapped,
            &mut self.incoming,
        )?;
        for index in 0..self.solution.len() {
            self.residual[index] = self.mapped[index] - self.solution[index];
        }
        if self.mapped.iter().any(|value| !value.is_finite())
            || self.residual.iter().any(|value| !value.is_finite())
        {
            return Err(SolverError::NonFiniteIteration);
        }

        let initial_residual_norm = infinity_norm(&self.residual);
        let initial_scale = infinity_norm(&self.mapped).max(infinity_norm(&self.solution));
        let initial_threshold =
            self.config.absolute_tolerance + self.config.relative_tolerance * initial_scale;
        if initial_residual_norm <= initial_threshold {
            self.solution.copy_from_slice(&self.mapped);
            self.diagnostics = SolverDiagnostics {
                iterations: 1,
                converged: true,
                reason: ConvergenceReason::Tolerance,
                initial_residual: initial_residual_norm,
                final_residual: initial_residual_norm,
                residual_history: vec![initial_residual_norm],
            };
            return Ok(true);
        }

        let initial_solution = &self.solution;
        let initial_residual = &self.residual;
        let krylov = {
            let problem = &self.problem;
            let incoming = &mut self.incoming;
            let mut apply_system = |input: &[f64], output: &mut [f64]| {
                problem.apply_linear(
                    transport_row_offsets,
                    transport_column_indices,
                    transport_values,
                    input,
                    output,
                    incoming,
                )?;
                for index in 0..input.len() {
                    output[index] = input[index] - output[index];
                }
                Ok(())
            };
            solve_bicgstab_from_initial(
                initial_solution,
                initial_residual,
                self.config,
                &mut apply_system,
            )
        };
        let (candidate, mut residual_history) = match krylov {
            Ok(result) => result,
            Err(SolverError::ImplicitLinearSolveDidNotConverge) => return Ok(false),
            Err(error) => return Err(error),
        };

        // Recurrence residuals can drift in finite precision. Re-evaluate the
        // actual affine map before changing the retained primal solution.
        self.problem.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
            &candidate,
            &mut self.mapped,
            &mut self.incoming,
        )?;
        for ((residual, mapped), candidate) in
            self.residual.iter_mut().zip(&self.mapped).zip(&candidate)
        {
            *residual = mapped - candidate;
        }
        if self.mapped.iter().any(|value| !value.is_finite())
            || self.residual.iter().any(|value| !value.is_finite())
        {
            return Err(SolverError::NonFiniteIteration);
        }
        let final_residual = infinity_norm(&self.residual);
        let scale = infinity_norm(&self.mapped).max(infinity_norm(&candidate));
        let threshold = self.config.absolute_tolerance + self.config.relative_tolerance * scale;
        if final_residual > threshold {
            return Ok(false);
        }

        residual_history.push(final_residual);
        self.solution.copy_from_slice(&self.mapped);
        self.diagnostics = SolverDiagnostics {
            iterations: residual_history.len(),
            converged: true,
            reason: ConvergenceReason::Tolerance,
            initial_residual: initial_residual_norm,
            final_residual,
            residual_history,
        };
        Ok(true)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn solve_jvp(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        transport_value_tangent: &[f64],
        forcing: &[f64],
        forcing_tangent: &[f64],
        scattering_coefficient_tangent: &[f64],
        dense_scattering_value_tangent: &[f64],
    ) -> Result<Vec<f64>, SolverError> {
        self.validate_linearization_configuration()?;
        if !self.diagnostics.converged && self.config.anderson_depth == 0 {
            return self.solve_jvp_replay(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                transport_value_tangent,
                forcing,
                forcing_tangent,
                scattering_coefficient_tangent,
                dense_scattering_value_tangent,
            );
        }

        let mut right_hand_side = vec![0.0; self.solution.len()];
        let mut incoming = vec![0.0; self.problem.incoming_size()];
        let mut incoming_tangent = vec![0.0; self.problem.incoming_size()];
        self.problem.apply_parameter_jvp(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            transport_value_tangent,
            forcing,
            forcing_tangent,
            scattering_coefficient_tangent,
            dense_scattering_value_tangent,
            &self.solution,
            &mut right_hand_side,
            &mut incoming,
            &mut incoming_tangent,
        )?;
        let mut linear_incoming = vec![0.0; self.problem.incoming_size()];
        let mut apply_linear = |state: &[f64], output: &mut [f64], system: bool| {
            self.problem
                .apply_linear(
                    transport_row_offsets,
                    transport_column_indices,
                    transport_values,
                    state,
                    output,
                    &mut linear_incoming,
                )
                .map_err(SolverError::from)?;
            if system {
                for index in 0..state.len() {
                    output[index] = state[index] - output[index];
                }
            }
            Ok(())
        };
        if self.diagnostics.converged {
            solve_converged_linear_system(&right_hand_side, self.config, apply_linear)
        } else {
            solve_affine_fixed_point(&right_hand_side, self.config, true, |input, output| {
                apply_linear(input, output, false)
            })
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_jvp_replay(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        transport_value_tangent: &[f64],
        forcing: &[f64],
        forcing_tangent: &[f64],
        scattering_coefficient_tangent: &[f64],
        dense_scattering_value_tangent: &[f64],
    ) -> Result<Vec<f64>, SolverError> {
        let iterations = self.diagnostics.iterations;
        let mut state = self.initial_state.clone();
        let mut state_tangent = vec![0.0; self.solution.len()];
        let mut mapped = vec![0.0; self.solution.len()];
        let mut mapped_tangent = vec![0.0; self.solution.len()];
        let mut incoming = vec![0.0; self.problem.incoming_size()];
        let mut incoming_tangent = vec![0.0; self.problem.incoming_size()];
        for iteration in 0..iterations {
            self.problem.apply_jvp(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                transport_value_tangent,
                forcing,
                forcing_tangent,
                scattering_coefficient_tangent,
                dense_scattering_value_tangent,
                &state,
                &state_tangent,
                &mut mapped,
                &mut mapped_tangent,
                &mut incoming,
                &mut incoming_tangent,
            )?;
            let damping = self.linearization_damping(iteration);
            for index in 0..state.len() {
                state[index] = (1.0 - damping) * state[index] + damping * mapped[index];
                state_tangent[index] =
                    (1.0 - damping) * state_tangent[index] + damping * mapped_tangent[index];
            }
        }
        Ok(state_tangent)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn solve_vjp(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
        solution_cotangent: &[f64],
    ) -> Result<FixedPointGradient, SolverError> {
        self.solve_vjp_impl(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
            solution_cotangent,
            true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn solve_vjp_compact(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
        solution_cotangent: &[f64],
    ) -> Result<FixedPointGradient, SolverError> {
        self.solve_vjp_impl(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
            solution_cotangent,
            false,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_vjp_impl(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
        solution_cotangent: &[f64],
        materialize_transport_gradient: bool,
    ) -> Result<FixedPointGradient, SolverError> {
        self.validate_linearization_configuration()?;
        self.problem.validate_iteration_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
        )?;
        if solution_cotangent.len() != self.solution.len() {
            return Err(SolverError::InitialGuessSize);
        }
        if !self.diagnostics.converged && self.config.anderson_depth == 0 {
            return self.solve_vjp_replay(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                forcing,
                solution_cotangent,
            );
        }

        let mut linear_incoming = vec![0.0; self.problem.incoming_size()];
        let mut scattering_scratch = vec![0.0; self.problem.scattering_transpose_scratch_size()];
        let mut apply_linear_transpose = |state: &[f64], output: &mut [f64], system: bool| {
            let result = if system {
                self.problem.apply_linear_system_transpose(
                    transport_row_offsets,
                    transport_column_indices,
                    transport_values,
                    state,
                    output,
                    &mut linear_incoming,
                    &mut scattering_scratch,
                )
            } else {
                self.problem.apply_linear_transpose(
                    transport_row_offsets,
                    transport_column_indices,
                    transport_values,
                    state,
                    output,
                    &mut linear_incoming,
                    &mut scattering_scratch,
                )
            };
            result.map_err(SolverError::from)
        };
        let adjoint = if self.diagnostics.converged {
            solve_converged_linear_system(solution_cotangent, self.config, apply_linear_transpose)?
        } else {
            solve_affine_fixed_point(solution_cotangent, self.config, true, |input, output| {
                apply_linear_transpose(input, output, false)
            })?
        };
        let mut gradient = if materialize_transport_gradient {
            FixedPointGradient::zeros(&self.problem)
        } else {
            // At convergence the transport-value VJP factorizes exactly as
            // incoming_cotangent[row] * solution[column]. Keep that form
            // compact for the geometry pullback.
            FixedPointGradient::zeros_without_transport(&self.problem)
        };
        let mut incoming_scratch = vec![0.0; self.problem.incoming_size()];
        let mut incoming_cotangent = vec![0.0; self.problem.incoming_size()];
        self.problem.apply_parameter_vjp(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
            &self.solution,
            &adjoint,
            &mut gradient.transport_values,
            &mut gradient.scattering_coefficients,
            &mut gradient.dense_scattering_values,
            &mut gradient.forcing,
            &mut incoming_scratch,
            &mut incoming_cotangent,
        )?;
        Ok(gradient)
    }

    fn solve_vjp_replay(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
        solution_cotangent: &[f64],
    ) -> Result<FixedPointGradient, SolverError> {
        let iterations = self.diagnostics.iterations;
        let mut states = Vec::with_capacity(iterations);
        let mut state = self.initial_state.clone();
        let mut mapped = vec![0.0; self.solution.len()];
        let mut incoming = vec![0.0; self.problem.incoming_size()];
        for iteration in 0..iterations {
            states.push(state.clone());
            self.problem.apply(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                forcing,
                &state,
                &mut mapped,
                &mut incoming,
            )?;
            let damping = self.linearization_damping(iteration);
            for index in 0..state.len() {
                state[index] = (1.0 - damping) * state[index] + damping * mapped[index];
            }
        }

        let mut gradient = FixedPointGradient::zeros(&self.problem);
        let mut state_cotangent = solution_cotangent.to_vec();
        let mut mapped_cotangent = vec![0.0; self.solution.len()];
        let mut previous_cotangent = vec![0.0; self.solution.len()];
        let mut incoming_scratch = vec![0.0; self.problem.incoming_size()];
        let mut incoming_cotangent = vec![0.0; self.problem.incoming_size()];
        for iteration in (0..iterations).rev() {
            let damping = self.linearization_damping(iteration);
            for index in 0..state_cotangent.len() {
                mapped_cotangent[index] = damping * state_cotangent[index];
                previous_cotangent[index] = (1.0 - damping) * state_cotangent[index];
            }
            let mut operator_state_cotangent = vec![0.0; self.solution.len()];
            self.problem.apply_vjp(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                forcing,
                &states[iteration],
                &mapped_cotangent,
                &mut operator_state_cotangent,
                &mut gradient.transport_values,
                &mut gradient.scattering_coefficients,
                &mut gradient.dense_scattering_values,
                &mut gradient.forcing,
                &mut incoming_scratch,
                &mut incoming_cotangent,
            )?;
            for index in 0..state_cotangent.len() {
                previous_cotangent[index] += operator_state_cotangent[index];
            }
            std::mem::swap(&mut state_cotangent, &mut previous_cotangent);
        }
        Ok(gradient)
    }

    fn validate_linearization_configuration(&self) -> Result<(), SolverError> {
        if self.diagnostics.iterations == 0 {
            return Err(SolverError::LinearizationBeforeSolve);
        }
        Ok(())
    }

    #[inline]
    fn linearization_damping(&self, iteration: usize) -> f64 {
        if self.diagnostics.converged && iteration + 1 == self.diagnostics.iterations {
            1.0
        } else {
            self.config.damping
        }
    }
}

fn completed_batch_result(
    solution: Vec<f64>,
    initial_residual: [f64; LANES],
    residual_history: [Vec<f64>; LANES],
) -> BatchSolveResult {
    let diagnostics = std::array::from_fn(|lane| SolverDiagnostics {
        iterations: residual_history[lane].len(),
        converged: true,
        reason: ConvergenceReason::Tolerance,
        initial_residual: initial_residual[lane],
        final_residual: residual_history[lane].last().copied().unwrap_or(0.0),
        residual_history: residual_history[lane].clone(),
    });
    BatchSolveResult {
        solution,
        diagnostics,
    }
}

#[inline]
fn copy_batch_lane(source: &[f64], target: &mut [f64], lane: usize) {
    debug_assert_eq!(source.len(), target.len());
    for element in 0..source.len() / LANES {
        target[element * LANES + lane] = source[element * LANES + lane];
    }
}

#[inline]
fn extract_batch_lane(source: &[f64], lane: usize, target: &mut [f64]) {
    debug_assert!(lane < LANES);
    debug_assert_eq!(source.len(), target.len() * LANES);
    for (element, value) in target.iter_mut().enumerate() {
        *value = source[element * LANES + lane];
    }
}

#[inline]
fn interleave_batch_lane(source: &[f64], lane: usize, target: &mut [f64]) {
    debug_assert!(lane < LANES);
    debug_assert_eq!(target.len(), source.len() * LANES);
    for (element, &value) in source.iter().enumerate() {
        target[element * LANES + lane] = value;
    }
}

#[inline]
fn zero_inactive_lanes(values: &mut [f64], active: [bool; LANES]) {
    for element in 0..values.len() / LANES {
        Batch4::select(active, Batch4::load(values, element), Batch4::splat(0.0))
            .store(values, element);
    }
}

#[inline]
fn safe_batch_divide(numerator: Batch4, denominator: Batch4, active: [bool; LANES]) -> Batch4 {
    let numerator = Batch4::select(active, numerator, Batch4::splat(0.0));
    let denominator = Batch4::select(active, denominator, Batch4::splat(1.0));
    numerator / denominator
}

fn validate_batch_divisor(value: Batch4, active: [bool; LANES]) -> Result<(), SolverError> {
    if value
        .to_array()
        .into_iter()
        .zip(active)
        .any(|(value, active)| active && (!value.is_finite() || value.abs() < 1.0e-30))
    {
        return Err(SolverError::ImplicitLinearSolveDidNotConverge);
    }
    Ok(())
}

fn validate_batch_positive_divisor(
    value: Batch4,
    active: [bool; LANES],
) -> Result<(), SolverError> {
    if value
        .to_array()
        .into_iter()
        .zip(active)
        .any(|(value, active)| active && (!value.is_finite() || value < 1.0e-30))
    {
        return Err(SolverError::ImplicitLinearSolveDidNotConverge);
    }
    Ok(())
}

/// Bounded-memory BiCGSTAB for four interleaved linear systems. Lanes stop
/// independently; completed lanes are masked out without changing the packed
/// storage layout used by the operator kernels.
fn solve_bicgstab_batch4<F>(
    right_hand_side: &[f64],
    config: SolverConfig,
    mut apply_system: F,
) -> Result<Vec<f64>, SolverError>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<(), SolverError>,
{
    if !right_hand_side.len().is_multiple_of(LANES) {
        return Err(SolverError::Operator(OperatorError::DimensionMismatch));
    }
    let state_size = right_hand_side.len() / LANES;
    let mut solution = vec![0.0; right_hand_side.len()];
    let mut residual = right_hand_side.to_vec();
    if residual.iter().any(|value| !value.is_finite()) {
        return Err(SolverError::NonFiniteIteration);
    }
    let initial_residual_norm = batch_infinity_norm(&residual);
    let mut active = std::array::from_fn(|lane| {
        initial_residual_norm[lane]
            > config.absolute_tolerance + config.relative_tolerance * initial_residual_norm[lane]
    });
    if !active.iter().any(|&value| value) {
        return Ok(solution);
    }
    zero_inactive_lanes(&mut residual, active);

    let shadow_residual = residual.clone();
    let mut direction = vec![0.0; right_hand_side.len()];
    let mut direction_image = vec![0.0; right_hand_side.len()];
    let mut intermediate_residual = vec![0.0; right_hand_side.len()];
    let mut intermediate_image = vec![0.0; right_hand_side.len()];
    let mut previous_rho = Batch4::splat(1.0);
    let mut alpha = Batch4::splat(1.0);
    let mut omega = Batch4::splat(1.0);

    for iteration in 0..config.max_iterations {
        let rho = batch_dot_product(&shadow_residual, &residual);
        validate_batch_divisor(rho, active)?;
        if iteration == 0 {
            direction.copy_from_slice(&residual);
        } else {
            validate_batch_divisor(previous_rho, active)?;
            validate_batch_divisor(omega, active)?;
            let beta = safe_batch_divide(rho, previous_rho, active)
                * safe_batch_divide(alpha, omega, active);
            for element in 0..state_size {
                let value = Batch4::load(&residual, element)
                    + beta
                        * (Batch4::load(&direction, element)
                            - omega * Batch4::load(&direction_image, element));
                Batch4::select(active, value, Batch4::splat(0.0)).store(&mut direction, element);
            }
        }

        apply_system(&direction, &mut direction_image)?;
        let alpha_denominator = batch_dot_product(&shadow_residual, &direction_image);
        validate_batch_divisor(alpha_denominator, active)?;
        alpha = safe_batch_divide(rho, alpha_denominator, active);
        for element in 0..state_size {
            (Batch4::load(&residual, element) - alpha * Batch4::load(&direction_image, element))
                .store(&mut intermediate_residual, element);
        }

        let intermediate_norm = batch_infinity_norm(&intermediate_residual);
        let solution_norm = batch_infinity_norm(&solution);
        let intermediate_converged = std::array::from_fn(|lane| {
            active[lane]
                && intermediate_norm[lane]
                    <= config.absolute_tolerance
                        + config.relative_tolerance
                            * initial_residual_norm[lane].max(solution_norm[lane])
        });
        if intermediate_converged.iter().any(|&value| value) {
            for element in 0..state_size {
                let candidate =
                    Batch4::load(&solution, element) + alpha * Batch4::load(&direction, element);
                Batch4::select(
                    intermediate_converged,
                    candidate,
                    Batch4::load(&solution, element),
                )
                .store(&mut solution, element);
            }
            for lane in 0..LANES {
                active[lane] &= !intermediate_converged[lane];
            }
            if !active.iter().any(|&value| value) {
                return Ok(solution);
            }
        }
        zero_inactive_lanes(&mut intermediate_residual, active);

        apply_system(&intermediate_residual, &mut intermediate_image)?;
        let omega_denominator = batch_dot_product(&intermediate_image, &intermediate_image);
        validate_batch_positive_divisor(omega_denominator, active)?;
        omega = safe_batch_divide(
            batch_dot_product(&intermediate_image, &intermediate_residual),
            omega_denominator,
            active,
        );
        validate_batch_divisor(omega, active)?;
        for element in 0..state_size {
            let updated_solution = Batch4::load(&solution, element)
                + alpha * Batch4::load(&direction, element)
                + omega * Batch4::load(&intermediate_residual, element);
            let updated_residual = Batch4::load(&intermediate_residual, element)
                - omega * Batch4::load(&intermediate_image, element);
            Batch4::select(active, updated_solution, Batch4::load(&solution, element))
                .store(&mut solution, element);
            Batch4::select(active, updated_residual, Batch4::splat(0.0))
                .store(&mut residual, element);
        }
        if solution.iter().any(|value| !value.is_finite())
            || residual.iter().any(|value| !value.is_finite())
        {
            return Err(SolverError::NonFiniteIteration);
        }
        let residual_norm = batch_infinity_norm(&residual);
        let solution_norm = batch_infinity_norm(&solution);
        for lane in 0..LANES {
            if active[lane]
                && residual_norm[lane]
                    <= config.absolute_tolerance
                        + config.relative_tolerance
                            * initial_residual_norm[lane].max(solution_norm[lane])
            {
                active[lane] = false;
            }
        }
        if !active.iter().any(|&value| value) {
            return Ok(solution);
        }
        previous_rho = rho;
    }
    Err(SolverError::ImplicitLinearSolveDidNotConverge)
}

fn solve_affine_fixed_point_batch4<F>(
    right_hand_side: &[f64],
    config: SolverConfig,
    mut apply_linear: F,
) -> Result<Vec<f64>, SolverError>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<(), SolverError>,
{
    if !right_hand_side.len().is_multiple_of(LANES) {
        return Err(SolverError::Operator(OperatorError::DimensionMismatch));
    }
    let state_size = right_hand_side.len() / LANES;
    let mut state = vec![0.0; right_hand_side.len()];
    let mut mapped = vec![0.0; right_hand_side.len()];
    let mut residual = vec![0.0; right_hand_side.len()];
    let right_hand_side_norm = batch_infinity_norm(right_hand_side);
    let mut active = std::array::from_fn(|lane| {
        right_hand_side_norm[lane]
            > config.absolute_tolerance + config.relative_tolerance * right_hand_side_norm[lane]
    });
    if !active.iter().any(|&value| value) {
        return Ok(state);
    }

    for _ in 0..config.max_iterations {
        apply_linear(&state, &mut mapped)?;
        for element in 0..state_size {
            let candidate = Batch4::load(&mapped, element) + Batch4::load(right_hand_side, element);
            (candidate - Batch4::load(&state, element)).store(&mut residual, element);
            let updated = Batch4::load(&state, element)
                + Batch4::splat(config.damping) * Batch4::load(&residual, element);
            Batch4::select(active, updated, Batch4::load(&state, element))
                .store(&mut state, element);
        }
        if state.iter().any(|value| !value.is_finite())
            || residual.iter().any(|value| !value.is_finite())
        {
            return Err(SolverError::NonFiniteIteration);
        }
        let residual_norm = batch_infinity_norm(&residual);
        let mapped_norm = batch_infinity_norm(&mapped);
        let state_norm = batch_infinity_norm(&state);
        for lane in 0..LANES {
            if active[lane]
                && residual_norm[lane]
                    <= config.absolute_tolerance
                        + config.relative_tolerance * mapped_norm[lane].max(state_norm[lane])
            {
                active[lane] = false;
            }
        }
        if !active.iter().any(|&value| value) {
            return Ok(state);
        }
    }
    Err(SolverError::ImplicitLinearSolveDidNotConverge)
}

#[derive(Debug, Clone, PartialEq)]
pub struct FixedPointGradient {
    pub transport_values: Vec<f64>,
    pub scattering_coefficients: Vec<f64>,
    pub dense_scattering_values: Vec<f64>,
    pub forcing: Vec<f64>,
}

impl FixedPointGradient {
    fn zeros(problem: &FixedPointProblem) -> Self {
        Self {
            transport_values: vec![0.0; problem.transport_value_size()],
            scattering_coefficients: vec![0.0; problem.scattering_coefficient_size()],
            dense_scattering_values: vec![0.0; problem.dense_scattering_value_size()],
            forcing: vec![0.0; problem.incoming_size()],
        }
    }

    fn zeros_without_transport(problem: &FixedPointProblem) -> Self {
        Self {
            transport_values: Vec::new(),
            scattering_coefficients: vec![0.0; problem.scattering_coefficient_size()],
            dense_scattering_values: vec![0.0; problem.dense_scattering_value_size()],
            forcing: vec![0.0; problem.incoming_size()],
        }
    }
}

fn solve_affine_fixed_point<F>(
    right_hand_side: &[f64],
    config: SolverConfig,
    return_last_iterate: bool,
    mut apply_linear: F,
) -> Result<Vec<f64>, SolverError>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<(), SolverError>,
{
    let mut state = vec![0.0; right_hand_side.len()];
    let mut mapped = vec![0.0; right_hand_side.len()];
    let mut residual = vec![0.0; right_hand_side.len()];
    let mut anderson = AndersonWorkspace::new(config.anderson_depth, right_hand_side.len());

    for _ in 0..config.max_iterations {
        apply_linear(&state, &mut mapped)?;
        for index in 0..state.len() {
            mapped[index] += right_hand_side[index];
            residual[index] = mapped[index] - state[index];
        }
        if mapped.iter().any(|value| !value.is_finite())
            || residual.iter().any(|value| !value.is_finite())
        {
            return Err(SolverError::NonFiniteIteration);
        }

        let residual_norm = infinity_norm(&residual);
        let scale = infinity_norm(&mapped).max(infinity_norm(&state));
        let threshold = config.absolute_tolerance + config.relative_tolerance * scale;
        if residual_norm <= threshold {
            return Ok(mapped);
        }

        let accelerated = if config.anderson_depth > 0 {
            anderson.record(&state, &residual)
                && anderson.update(&mut state, &mut mapped, config.damping)
        } else {
            false
        };
        if !accelerated {
            damped_update(&mut state, &residual, config.damping);
        }
    }

    if return_last_iterate {
        Ok(state)
    } else {
        Err(SolverError::ImplicitLinearSolveDidNotConverge)
    }
}

fn solve_converged_linear_system<F>(
    right_hand_side: &[f64],
    config: SolverConfig,
    mut apply_operator: F,
) -> Result<Vec<f64>, SolverError>
where
    F: FnMut(&[f64], &mut [f64], bool) -> Result<(), SolverError>,
{
    let bicgstab_result = {
        let mut apply_system =
            |input: &[f64], output: &mut [f64]| apply_operator(input, output, true);
        solve_bicgstab(right_hand_side, config, &mut apply_system)
    };
    match bicgstab_result {
        Ok(solution) => Ok(solution),
        Err(SolverError::ImplicitLinearSolveDidNotConverge) => {
            solve_affine_fixed_point(right_hand_side, config, false, |input, output| {
                apply_operator(input, output, false)
            })
        }
        Err(error) => Err(error),
    }
}

/// Solves `B x = b`, where `apply_system` evaluates `B x`.
///
/// BiCGSTAB is used here because the converged successive-orders derivative
/// system is linear and generally nonsymmetric. Its storage is bounded to
/// seven state-sized work vectors, less than depth-three Anderson iteration.
fn solve_bicgstab<F>(
    right_hand_side: &[f64],
    config: SolverConfig,
    apply_system: &mut F,
) -> Result<Vec<f64>, SolverError>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<(), SolverError>,
{
    let initial_solution = vec![0.0; right_hand_side.len()];
    solve_bicgstab_from_initial(&initial_solution, right_hand_side, config, apply_system)
        .map(|(solution, _)| solution)
}

/// BiCGSTAB correction solve starting from a supplied state and exact initial
/// residual `b - B x0`. The residual history is returned so primal callers can
/// retain useful convergence diagnostics and independently verify the final
/// fixed-point residual before accepting the Krylov result.
fn solve_bicgstab_from_initial<F>(
    initial_solution: &[f64],
    initial_residual: &[f64],
    config: SolverConfig,
    apply_system: &mut F,
) -> Result<(Vec<f64>, Vec<f64>), SolverError>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<(), SolverError>,
{
    let size = initial_residual.len();
    debug_assert_eq!(initial_solution.len(), size);
    let mut solution = initial_solution.to_vec();
    let mut residual = initial_residual.to_vec();
    let shadow_residual = residual.clone();
    let mut direction = vec![0.0; size];
    let mut direction_image = vec![0.0; size];
    let mut intermediate_residual = vec![0.0; size];
    let mut intermediate_image = vec![0.0; size];
    let initial_residual_norm = infinity_norm(initial_residual);
    let threshold = |state: &[f64]| {
        config.absolute_tolerance
            + config.relative_tolerance * initial_residual_norm.max(infinity_norm(state))
    };
    let mut residual_history = Vec::with_capacity(config.max_iterations + 1);
    residual_history.push(initial_residual_norm);
    if infinity_norm(&residual) <= threshold(&solution) {
        return Ok((solution, residual_history));
    }

    let mut previous_rho = 1.0_f64;
    let mut alpha = 1.0_f64;
    let mut omega = 1.0_f64;
    for iteration in 0..config.max_iterations {
        let rho = dot_product(&shadow_residual, &residual);
        if !rho.is_finite() || rho.abs() < 1.0e-30 {
            return Err(SolverError::ImplicitLinearSolveDidNotConverge);
        }
        if iteration == 0 {
            direction.copy_from_slice(&residual);
        } else {
            if !omega.is_finite() || omega.abs() < 1.0e-30 {
                return Err(SolverError::ImplicitLinearSolveDidNotConverge);
            }
            let beta = (rho / previous_rho) * (alpha / omega);
            for index in 0..size {
                direction[index] =
                    residual[index] + beta * (direction[index] - omega * direction_image[index]);
            }
        }

        apply_system(&direction, &mut direction_image)?;
        let alpha_denominator = dot_product(&shadow_residual, &direction_image);
        if !alpha_denominator.is_finite() || alpha_denominator.abs() < 1.0e-30 {
            return Err(SolverError::ImplicitLinearSolveDidNotConverge);
        }
        alpha = rho / alpha_denominator;
        for index in 0..size {
            intermediate_residual[index] = residual[index] - alpha * direction_image[index];
        }
        if infinity_norm(&intermediate_residual) <= threshold(&solution) {
            for index in 0..size {
                solution[index] += alpha * direction[index];
            }
            residual_history.push(infinity_norm(&intermediate_residual));
            return Ok((solution, residual_history));
        }

        apply_system(&intermediate_residual, &mut intermediate_image)?;
        let omega_denominator = dot_product(&intermediate_image, &intermediate_image);
        if !omega_denominator.is_finite() || omega_denominator < 1.0e-30 {
            return Err(SolverError::ImplicitLinearSolveDidNotConverge);
        }
        omega = dot_product(&intermediate_image, &intermediate_residual) / omega_denominator;
        if !omega.is_finite() || omega.abs() < 1.0e-30 {
            return Err(SolverError::ImplicitLinearSolveDidNotConverge);
        }
        for index in 0..size {
            solution[index] += alpha * direction[index] + omega * intermediate_residual[index];
            residual[index] = intermediate_residual[index] - omega * intermediate_image[index];
        }
        if solution.iter().any(|value| !value.is_finite())
            || residual.iter().any(|value| !value.is_finite())
        {
            return Err(SolverError::NonFiniteIteration);
        }
        let residual_norm = infinity_norm(&residual);
        residual_history.push(residual_norm);
        if residual_norm <= threshold(&solution) {
            return Ok((solution, residual_history));
        }
        previous_rho = rho;
    }
    Err(SolverError::ImplicitLinearSolveDidNotConverge)
}

fn infinity_norm(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |norm, value| norm.max(value.abs()))
}

fn dot_product(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| left * right)
        .sum()
}

fn damped_update(state: &mut [f64], residual: &[f64], damping: f64) {
    for (value, correction) in state.iter_mut().zip(residual) {
        *value += damping * correction;
    }
}

fn solve_dense(matrix: &mut [f64], right: &mut [f64], size: usize) -> bool {
    for pivot in 0..size {
        let mut largest = pivot;
        for row in pivot + 1..size {
            if matrix[row * size + pivot].abs() > matrix[largest * size + pivot].abs() {
                largest = row;
            }
        }
        if !matrix[largest * size + pivot].is_finite()
            || matrix[largest * size + pivot].abs() < 1.0e-20
        {
            return false;
        }
        if largest != pivot {
            for column in pivot..size {
                matrix.swap(pivot * size + column, largest * size + column);
            }
            right.swap(pivot, largest);
        }
        let diagonal = matrix[pivot * size + pivot];
        for row in pivot + 1..size {
            let factor = matrix[row * size + pivot] / diagonal;
            matrix[row * size + pivot] = 0.0;
            for column in pivot + 1..size {
                matrix[row * size + column] -= factor * matrix[pivot * size + column];
            }
            right[row] -= factor * right[pivot];
        }
    }
    for row in (0..size).rev() {
        let mut value = right[row];
        for column in row + 1..size {
            value -= matrix[row * size + column] * right[column];
        }
        right[row] = value / matrix[row * size + row];
    }
    right.iter().all(|value| value.is_finite())
}
