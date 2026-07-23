use std::fmt::{Display, Formatter};

use super::{FixedPointProblem, OperatorError};

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
    UnsupportedLinearizationConfiguration,
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
            Self::UnsupportedLinearizationConfiguration => formatter.write_str(
                "Anderson-accelerated successive-orders products require a converged primal solve",
            ),
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
    state_history: Vec<Vec<f64>>,
    residual_vector_history: Vec<Vec<f64>>,
    diagnostics: SolverDiagnostics,
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
            state_history: Vec::with_capacity(config.anderson_depth + 1),
            residual_vector_history: Vec::with_capacity(config.anderson_depth + 1),
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

        self.state_history.clear();
        self.residual_vector_history.clear();
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
                self.record_history();
                self.state_history.len() >= 2
                    && anderson_update(
                        &mut self.solution,
                        &mut self.mapped,
                        &self.state_history,
                        &self.residual_vector_history,
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
        if !self.diagnostics.converged {
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

        let mut mapped = vec![0.0; self.solution.len()];
        let mut right_hand_side = vec![0.0; self.solution.len()];
        let zero_state_tangent = vec![0.0; self.solution.len()];
        let mut incoming = vec![0.0; self.problem.incoming_size()];
        let mut incoming_tangent = vec![0.0; self.problem.incoming_size()];
        self.problem.apply_jvp(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            transport_value_tangent,
            forcing,
            forcing_tangent,
            scattering_coefficient_tangent,
            dense_scattering_value_tangent,
            &self.solution,
            &zero_state_tangent,
            &mut mapped,
            &mut right_hand_side,
            &mut incoming,
            &mut incoming_tangent,
        )?;

        let mut linear_incoming = vec![0.0; self.problem.incoming_size()];
        solve_affine_fixed_point(&right_hand_side, self.config, |state, output| {
            self.problem
                .apply_linear(
                    transport_row_offsets,
                    transport_column_indices,
                    transport_values,
                    state,
                    output,
                    &mut linear_incoming,
                )
                .map_err(SolverError::from)
        })
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
        if !self.diagnostics.converged {
            return self.solve_vjp_replay(
                transport_row_offsets,
                transport_column_indices,
                transport_values,
                forcing,
                solution_cotangent,
            );
        }

        let mut linear_incoming = vec![0.0; self.problem.incoming_size()];
        let adjoint =
            solve_affine_fixed_point(solution_cotangent, self.config, |state, output| {
                self.problem
                    .apply_linear_transpose(
                        transport_row_offsets,
                        transport_column_indices,
                        transport_values,
                        state,
                        output,
                        &mut linear_incoming,
                    )
                    .map_err(SolverError::from)
            })?;

        let mut gradient = if materialize_transport_gradient {
            FixedPointGradient::zeros(&self.problem)
        } else {
            // At convergence the transport-value VJP factorizes exactly as
            // incoming_cotangent[row] * solution[column]. Keep that form
            // compact for the geometry pullback.
            FixedPointGradient::zeros_without_transport(&self.problem)
        };
        let mut state_cotangent = vec![0.0; self.solution.len()];
        let mut incoming_scratch = vec![0.0; self.problem.incoming_size()];
        let mut incoming_cotangent = vec![0.0; self.problem.incoming_size()];
        self.problem.apply_vjp(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
            &self.solution,
            &adjoint,
            &mut state_cotangent,
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
        if !self.diagnostics.converged && self.config.anderson_depth != 0 {
            return Err(SolverError::UnsupportedLinearizationConfiguration);
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

    fn record_history(&mut self) {
        let capacity = self.config.anderson_depth + 1;
        if self.state_history.len() == capacity {
            self.state_history.remove(0);
            self.residual_vector_history.remove(0);
        }
        self.state_history.push(self.solution.clone());
        self.residual_vector_history.push(self.residual.clone());
    }
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
    mut apply_linear: F,
) -> Result<Vec<f64>, SolverError>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<(), SolverError>,
{
    let mut state = vec![0.0; right_hand_side.len()];
    let mut mapped = vec![0.0; right_hand_side.len()];
    let mut residual = vec![0.0; right_hand_side.len()];
    let mut state_history = Vec::with_capacity(config.anderson_depth + 1);
    let mut residual_history = Vec::with_capacity(config.anderson_depth + 1);

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
            let capacity = config.anderson_depth + 1;
            if state_history.len() == capacity {
                state_history.remove(0);
                residual_history.remove(0);
            }
            state_history.push(state.clone());
            residual_history.push(residual.clone());
            state_history.len() >= 2
                && anderson_update(
                    &mut state,
                    &mut mapped,
                    &state_history,
                    &residual_history,
                    config.damping,
                )
        } else {
            false
        };
        if !accelerated {
            damped_update(&mut state, &residual, config.damping);
        }
    }

    Err(SolverError::ImplicitLinearSolveDidNotConverge)
}

fn infinity_norm(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |norm, value| norm.max(value.abs()))
}

fn damped_update(state: &mut [f64], residual: &[f64], damping: f64) {
    for (value, correction) in state.iter_mut().zip(residual) {
        *value += damping * correction;
    }
}

/// Type-II Anderson update. Returns false when the small least-squares system
/// is singular, allowing the caller to fall back to damped Picard iteration.
fn anderson_update(
    state: &mut [f64],
    mapped: &mut [f64],
    states: &[Vec<f64>],
    residuals: &[Vec<f64>],
    damping: f64,
) -> bool {
    let differences = states.len() - 1;
    let vector_size = state.len();
    let current_residual = residuals.last().unwrap();
    let mut normal = vec![0.0; differences * differences];
    let mut right = vec![0.0; differences];
    for row in 0..differences {
        right[row] = dot_difference(&residuals[row + 1], &residuals[row], current_residual);
        for column in 0..differences {
            normal[row * differences + column] = dot_differences(
                &residuals[row + 1],
                &residuals[row],
                &residuals[column + 1],
                &residuals[column],
            );
        }
        normal[row * differences + row] += 1.0e-14;
    }
    if !solve_dense(&mut normal, &mut right, differences) {
        return false;
    }

    for row in 0..vector_size {
        let mut candidate = mapped[row];
        for column in 0..differences {
            let delta_state = states[column + 1][row] - states[column][row];
            let delta_residual = residuals[column + 1][row] - residuals[column][row];
            candidate -= right[column] * (delta_state + delta_residual);
        }
        mapped[row] = (1.0 - damping) * state[row] + damping * candidate;
    }
    if mapped.iter().any(|value| !value.is_finite()) {
        return false;
    }
    state.copy_from_slice(mapped);
    true
}

fn dot_difference(upper: &[f64], lower: &[f64], right: &[f64]) -> f64 {
    upper
        .iter()
        .zip(lower)
        .zip(right)
        .map(|((&upper, &lower), &right)| (upper - lower) * right)
        .sum()
}

fn dot_differences(
    left_upper: &[f64],
    left_lower: &[f64],
    right_upper: &[f64],
    right_lower: &[f64],
) -> f64 {
    left_upper
        .iter()
        .zip(left_lower)
        .zip(right_upper.iter().zip(right_lower))
        .map(
            |((&left_upper, &left_lower), (&right_upper, &right_lower))| {
                (left_upper - left_lower) * (right_upper - right_lower)
            },
        )
        .sum()
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
