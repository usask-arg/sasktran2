use crate::math::wigner::WignerDCalculator;

use super::OperatorError;

/// Analysis and synthesis matrices for a scalar real spherical-harmonic
/// expansion.
///
/// Modes belonging to degree `l` occupy `[l*l, (l + 1)*(l + 1))`. The basis
/// normalization is chosen so that the dot product of the degree-`l` modes at
/// two directions is exactly `P_l(incoming.dot(outgoing))`. Consequently the
/// atmospheric Legendre coefficients can multiply modes directly without any
/// additional normalization.
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarCoefficientBasis {
    input_size: usize,
    output_size: usize,
    num_coefficients: usize,
    mode_degrees: Vec<usize>,
    analysis: Vec<f64>,
    synthesis: Vec<f64>,
}

impl ScalarCoefficientBasis {
    pub fn from_directions(
        incoming_directions: &[[f64; 3]],
        incoming_weights: &[f64],
        outgoing_directions: &[[f64; 3]],
        num_coefficients: usize,
    ) -> Result<Self, OperatorError> {
        if num_coefficients == 0
            || incoming_directions.len() != incoming_weights.len()
            || incoming_directions.is_empty()
            || outgoing_directions.is_empty()
            || incoming_weights.iter().any(|weight| !weight.is_finite())
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let incoming_directions = normalized_directions(incoming_directions)?;
        let outgoing_directions = normalized_directions(outgoing_directions)?;
        let num_modes = num_coefficients
            .checked_mul(num_coefficients)
            .ok_or(OperatorError::DimensionMismatch)?;
        let mut mode_degrees = vec![0; num_modes];
        let mut incoming_basis = vec![0.0; num_modes * incoming_directions.len()];
        let mut synthesis = vec![0.0; outgoing_directions.len() * num_modes];

        fill_basis(
            &incoming_directions,
            num_coefficients,
            &mut incoming_basis,
            &mut mode_degrees,
        );
        let mut outgoing_mode_degrees = vec![0; num_modes];
        fill_basis(
            &outgoing_directions,
            num_coefficients,
            &mut synthesis,
            &mut outgoing_mode_degrees,
        );
        if mode_degrees != outgoing_mode_degrees {
            return Err(OperatorError::DimensionMismatch);
        }

        // fill_basis stores [direction, mode]. Transpose and apply quadrature
        // weights to form the analysis matrix [mode, incoming direction].
        let mut analysis = vec![0.0; num_modes * incoming_directions.len()];
        for mode_index in 0..num_modes {
            for input_index in 0..incoming_directions.len() {
                analysis[mode_index * incoming_directions.len() + input_index] = incoming_basis
                    [input_index * num_modes + mode_index]
                    * incoming_weights[input_index];
            }
        }

        Ok(Self {
            input_size: incoming_directions.len(),
            output_size: outgoing_directions.len(),
            num_coefficients,
            mode_degrees,
            analysis,
            synthesis,
        })
    }

    #[inline]
    pub fn input_size(&self) -> usize {
        self.input_size
    }

    #[inline]
    pub fn output_size(&self) -> usize {
        self.output_size
    }

    #[inline]
    pub fn num_coefficients(&self) -> usize {
        self.num_coefficients
    }

    #[inline]
    pub fn num_modes(&self) -> usize {
        self.mode_degrees.len()
    }

    fn apply(&self, input: &[f64], coefficients: &[f64], output: &mut [f64], moments: &mut [f64]) {
        debug_assert_eq!(input.len(), self.input_size);
        debug_assert_eq!(coefficients.len(), self.num_coefficients);
        debug_assert_eq!(output.len(), self.output_size);
        debug_assert_eq!(moments.len(), self.num_modes());

        for (mode_index, moment) in moments.iter_mut().enumerate() {
            let row_start = mode_index * self.input_size;
            *moment = self.analysis[row_start..row_start + self.input_size]
                .iter()
                .zip(input)
                .map(|(&basis, &radiance)| basis * radiance)
                .sum::<f64>()
                * coefficients[self.mode_degrees[mode_index]];
        }
        for (output_index, result) in output.iter_mut().enumerate() {
            let row_start = output_index * self.num_modes();
            *result = self.synthesis[row_start..row_start + self.num_modes()]
                .iter()
                .zip(moments.iter())
                .map(|(&basis, &moment)| basis * moment)
                .sum();
        }
    }
}

/// A block-diagonal scattering operator with coefficient-space atmospheric
/// blocks followed by optional dense boundary blocks.
///
/// All atmospheric blocks share an angular basis but carry independent phase
/// coefficients. Spatial geometry is represented only by the block offsets,
/// so the operator is agnostic to whether those blocks came from a 1-D or 2-D
/// atmosphere.
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarCoefficientScattering {
    output_size: usize,
    input_size: usize,
    output_offsets: Vec<usize>,
    input_offsets: Vec<usize>,
    coefficient_blocks: usize,
    basis: ScalarCoefficientBasis,
    coefficients: Vec<f64>,
    dense_value_offsets: Vec<usize>,
    dense_values: Vec<f64>,
}

impl ScalarCoefficientScattering {
    pub fn new(
        output_offsets: Vec<usize>,
        input_offsets: Vec<usize>,
        coefficient_blocks: usize,
        basis: ScalarCoefficientBasis,
    ) -> Result<Self, OperatorError> {
        if output_offsets.len() < 2
            || output_offsets.len() != input_offsets.len()
            || output_offsets[0] != 0
            || input_offsets[0] != 0
            || output_offsets
                .windows(2)
                .any(|window| window[0] > window[1])
            || input_offsets.windows(2).any(|window| window[0] > window[1])
            || coefficient_blocks > output_offsets.len() - 1
        {
            return Err(OperatorError::InvalidBlockOffsets);
        }
        for block in 0..coefficient_blocks {
            if output_offsets[block + 1] - output_offsets[block] != basis.output_size()
                || input_offsets[block + 1] - input_offsets[block] != basis.input_size()
            {
                return Err(OperatorError::DimensionMismatch);
            }
        }

        let num_blocks = output_offsets.len() - 1;
        let mut dense_value_offsets: Vec<usize> = vec![0];
        for block in coefficient_blocks..num_blocks {
            let rows = output_offsets[block + 1] - output_offsets[block];
            let columns = input_offsets[block + 1] - input_offsets[block];
            dense_value_offsets.push(
                dense_value_offsets
                    .last()
                    .copied()
                    .unwrap()
                    .checked_add(
                        rows.checked_mul(columns)
                            .ok_or(OperatorError::DimensionMismatch)?,
                    )
                    .ok_or(OperatorError::DimensionMismatch)?,
            );
        }
        let dense_value_size = dense_value_offsets.last().copied().unwrap();
        let coefficient_size = coefficient_blocks
            .checked_mul(basis.num_coefficients())
            .ok_or(OperatorError::DimensionMismatch)?;
        Ok(Self {
            output_size: *output_offsets.last().unwrap(),
            input_size: *input_offsets.last().unwrap(),
            output_offsets,
            input_offsets,
            coefficient_blocks,
            basis,
            coefficients: vec![0.0; coefficient_size],
            dense_value_offsets,
            dense_values: vec![0.0; dense_value_size],
        })
    }

    #[inline]
    pub fn output_size(&self) -> usize {
        self.output_size
    }

    #[inline]
    pub fn input_size(&self) -> usize {
        self.input_size
    }

    #[inline]
    pub fn coefficient_value_size(&self) -> usize {
        self.coefficients.len()
    }

    #[inline]
    pub fn dense_value_size(&self) -> usize {
        self.dense_values.len()
    }

    pub fn set_coefficients(&mut self, coefficients: &[f64]) -> Result<(), OperatorError> {
        if coefficients.len() != self.coefficients.len() {
            return Err(OperatorError::DimensionMismatch);
        }
        if coefficients.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        self.coefficients.copy_from_slice(coefficients);
        Ok(())
    }

    pub fn set_dense_values(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        if values.len() != self.dense_values.len() {
            return Err(OperatorError::DimensionMismatch);
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        self.dense_values.copy_from_slice(values);
        Ok(())
    }

    pub fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        if input.len() != self.input_size || output.len() != self.output_size {
            return Err(OperatorError::DimensionMismatch);
        }
        let num_coefficients = self.basis.num_coefficients();
        let mut moments = vec![0.0; self.basis.num_modes()];
        for block in 0..self.coefficient_blocks {
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let coefficient_start = block * num_coefficients;
            self.basis.apply(
                &input[input_start..input_end],
                &self.coefficients[coefficient_start..coefficient_start + num_coefficients],
                &mut output[output_start..output_end],
                &mut moments,
            );
        }

        for block in self.coefficient_blocks..self.output_offsets.len() - 1 {
            let dense_block = block - self.coefficient_blocks;
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let columns = input_end - input_start;
            let value_start = self.dense_value_offsets[dense_block];
            for (local_row, result) in output[output_start..output_end].iter_mut().enumerate() {
                let row_start = value_start + local_row * columns;
                *result = self.dense_values[row_start..row_start + columns]
                    .iter()
                    .zip(&input[input_start..input_end])
                    .map(|(&value, &incoming)| value * incoming)
                    .sum();
            }
        }
        Ok(())
    }
}

fn normalized_directions(directions: &[[f64; 3]]) -> Result<Vec<[f64; 3]>, OperatorError> {
    directions
        .iter()
        .map(|direction| {
            let norm = direction
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt();
            if !norm.is_finite() || norm <= 0.0 {
                return Err(OperatorError::NonFiniteValue);
            }
            Ok([
                direction[0] / norm,
                direction[1] / norm,
                direction[2] / norm,
            ])
        })
        .collect()
}

fn fill_basis(
    directions: &[[f64; 3]],
    num_coefficients: usize,
    values: &mut [f64],
    mode_degrees: &mut [usize],
) {
    let num_modes = num_coefficients * num_coefficients;
    for (direction_index, direction) in directions.iter().enumerate() {
        let theta = direction[2].clamp(-1.0, 1.0).acos();
        let phi = direction[1].atan2(direction[0]);
        for m in 0..num_coefficients {
            let calculator = WignerDCalculator::new(m as i32, 0);
            let mut degrees = vec![0.0; num_coefficients];
            calculator.vector_d(theta, &mut degrees);
            for (degree, &degree_value) in degrees.iter().enumerate().skip(m) {
                let mode_start = degree * degree;
                mode_degrees[mode_start] = degree;
                if m == 0 {
                    values[direction_index * num_modes + mode_start] = degree_value;
                } else {
                    let scale = std::f64::consts::SQRT_2 * degree_value;
                    let cos_mode = mode_start + 2 * m - 1;
                    let sin_mode = mode_start + 2 * m;
                    mode_degrees[cos_mode] = degree;
                    mode_degrees[sin_mode] = degree;
                    values[direction_index * num_modes + cos_mode] = scale * (m as f64 * phi).cos();
                    values[direction_index * num_modes + sin_mode] = scale * (m as f64 * phi).sin();
                }
            }
        }
    }
}
