use crate::math::wigner::WignerDCalculator;
use num::complex::Complex64;

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

    #[allow(clippy::too_many_arguments)]
    fn apply_vjp(
        &self,
        input: &[f64],
        coefficients: &[f64],
        output_cotangent: &[f64],
        input_cotangent: &mut [f64],
        coefficient_gradient: &mut [f64],
        analyzed: &mut [f64],
        moment_cotangent: &mut [f64],
    ) {
        debug_assert_eq!(input.len(), self.input_size);
        debug_assert_eq!(coefficients.len(), self.num_coefficients);
        debug_assert_eq!(output_cotangent.len(), self.output_size);
        debug_assert_eq!(input_cotangent.len(), self.input_size);
        debug_assert_eq!(coefficient_gradient.len(), self.num_coefficients);
        for mode_index in 0..self.num_modes() {
            let analysis_start = mode_index * self.input_size;
            analyzed[mode_index] = self.analysis[analysis_start..analysis_start + self.input_size]
                .iter()
                .zip(input)
                .map(|(&basis, &radiance)| basis * radiance)
                .sum();
            moment_cotangent[mode_index] = (0..self.output_size)
                .map(|output_index| {
                    self.synthesis[output_index * self.num_modes() + mode_index]
                        * output_cotangent[output_index]
                })
                .sum();
            let degree = self.mode_degrees[mode_index];
            coefficient_gradient[degree] += analyzed[mode_index] * moment_cotangent[mode_index];
            let scaled_cotangent = coefficients[degree] * moment_cotangent[mode_index];
            for (input_index, result) in input_cotangent.iter_mut().enumerate() {
                *result += self.analysis[analysis_start + input_index] * scaled_cotangent;
            }
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

    #[allow(clippy::too_many_arguments)]
    pub fn apply_jvp(
        &self,
        input: &[f64],
        input_tangent: &[f64],
        coefficient_tangent: &[f64],
        dense_value_tangent: &[f64],
        output_tangent: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.input_size
            || input_tangent.len() != self.input_size
            || coefficient_tangent.len() != self.coefficients.len()
            || dense_value_tangent.len() != self.dense_values.len()
            || output_tangent.len() != self.output_size
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let num_coefficients = self.basis.num_coefficients();
        let mut moments = vec![0.0; self.basis.num_modes()];
        let mut parameter_output = vec![0.0; self.basis.output_size()];
        for block in 0..self.coefficient_blocks {
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let coefficient_start = block * num_coefficients;
            self.basis.apply(
                &input_tangent[input_start..input_end],
                &self.coefficients[coefficient_start..coefficient_start + num_coefficients],
                &mut output_tangent[output_start..output_end],
                &mut moments,
            );
            self.basis.apply(
                &input[input_start..input_end],
                &coefficient_tangent[coefficient_start..coefficient_start + num_coefficients],
                &mut parameter_output,
                &mut moments,
            );
            for (result, parameter) in output_tangent[output_start..output_end]
                .iter_mut()
                .zip(&parameter_output)
            {
                *result += parameter;
            }
        }
        for block in self.coefficient_blocks..self.output_offsets.len() - 1 {
            let dense_block = block - self.coefficient_blocks;
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let columns = input_end - input_start;
            let value_start = self.dense_value_offsets[dense_block];
            for (local_row, result) in output_tangent[output_start..output_end]
                .iter_mut()
                .enumerate()
            {
                let row_start = value_start + local_row * columns;
                *result = self.dense_values[row_start..row_start + columns]
                    .iter()
                    .zip(&input_tangent[input_start..input_end])
                    .zip(&dense_value_tangent[row_start..row_start + columns])
                    .zip(&input[input_start..input_end])
                    .map(
                        |(((&value, &incoming_tangent), &value_tangent), &incoming)| {
                            value * incoming_tangent + value_tangent * incoming
                        },
                    )
                    .sum();
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply_vjp(
        &self,
        input: &[f64],
        output_cotangent: &[f64],
        input_cotangent: &mut [f64],
        coefficient_gradient: &mut [f64],
        dense_value_gradient: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.input_size
            || output_cotangent.len() != self.output_size
            || input_cotangent.len() != self.input_size
            || coefficient_gradient.len() != self.coefficients.len()
            || dense_value_gradient.len() != self.dense_values.len()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        input_cotangent.fill(0.0);
        coefficient_gradient.fill(0.0);
        dense_value_gradient.fill(0.0);
        let num_coefficients = self.basis.num_coefficients();
        let mut analyzed = vec![0.0; self.basis.num_modes()];
        let mut moment_cotangent = vec![0.0; self.basis.num_modes()];
        for block in 0..self.coefficient_blocks {
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let coefficient_start = block * num_coefficients;
            self.basis.apply_vjp(
                &input[input_start..input_end],
                &self.coefficients[coefficient_start..coefficient_start + num_coefficients],
                &output_cotangent[output_start..output_end],
                &mut input_cotangent[input_start..input_end],
                &mut coefficient_gradient[coefficient_start..coefficient_start + num_coefficients],
                &mut analyzed,
                &mut moment_cotangent,
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
            for (local_row, &outgoing_cotangent) in output_cotangent[output_start..output_end]
                .iter()
                .enumerate()
            {
                let row_start = value_start + local_row * columns;
                for local_column in 0..columns {
                    input_cotangent[input_start + local_column] +=
                        self.dense_values[row_start + local_column] * outgoing_cotangent;
                    dense_value_gradient[row_start + local_column] +=
                        input[input_start + local_column] * outgoing_cotangent;
                }
            }
        }
        Ok(())
    }
}

/// Geometry-only spin-0/spin-2 basis for polarized I/Q/U scattering.
///
/// Each `(l, m)` mode stores three complex channels: intensity, `Q + iU`, and
/// `Q - iU`.  In this representation the four Greek coefficient families
/// (`a1`, `a2`, `a3`, `b1`) act through the same small real mixing matrix for
/// every `m`; all meridian-plane rotations are contained in the basis.
#[derive(Debug, Clone, PartialEq)]
pub struct VectorCoefficientBasis {
    input_directions: usize,
    output_directions: usize,
    num_coefficients: usize,
    mode_degrees: Vec<usize>,
    analysis: [Vec<Complex64>; 3],
    synthesis: [Vec<Complex64>; 3],
    legacy_frame_corrections: Vec<VectorKernelCorrection>,
}

/// The legacy C++ Stokes-frame convention suppresses rotations when either
/// direction is polar (or the two directions are parallel). Those exceptional
/// pairs do not obey the otherwise separable spin-harmonic addition theorem,
/// so retain only their small coefficient-space residual here.
#[derive(Debug, Clone, PartialEq)]
struct VectorKernelCorrection {
    input_index: usize,
    output_index: usize,
    // [degree, Greek coefficient family, output Stokes, input Stokes]
    values: Vec<f64>,
}

impl VectorCoefficientBasis {
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
        let legacy_incoming_directions = incoming_directions.to_vec();
        let legacy_outgoing_directions = outgoing_directions.to_vec();
        let incoming_directions = normalized_directions(incoming_directions)?;
        let outgoing_directions = normalized_directions(outgoing_directions)?;
        let num_modes = num_coefficients
            .checked_mul(num_coefficients)
            .ok_or(OperatorError::DimensionMismatch)?;
        let mut mode_degrees = vec![0; num_modes];
        let mut analysis: [Vec<Complex64>; 3] = std::array::from_fn(|_| {
            vec![Complex64::new(0.0, 0.0); num_modes * incoming_directions.len()]
        });
        let mut synthesis: [Vec<Complex64>; 3] = std::array::from_fn(|_| {
            vec![Complex64::new(0.0, 0.0); outgoing_directions.len() * num_modes]
        });

        fill_vector_basis(
            &incoming_directions,
            num_coefficients,
            &mut analysis,
            &mut mode_degrees,
            true,
            Some(incoming_weights),
        );
        let mut outgoing_mode_degrees = vec![0; num_modes];
        fill_vector_basis(
            &outgoing_directions,
            num_coefficients,
            &mut synthesis,
            &mut outgoing_mode_degrees,
            false,
            None,
        );
        if mode_degrees != outgoing_mode_degrees {
            return Err(OperatorError::DimensionMismatch);
        }
        let legacy_frame_corrections = build_legacy_frame_corrections(
            &legacy_incoming_directions,
            incoming_weights,
            &legacy_outgoing_directions,
            num_coefficients,
            &analysis,
            &synthesis,
        );

        Ok(Self {
            input_directions: incoming_directions.len(),
            output_directions: outgoing_directions.len(),
            num_coefficients,
            mode_degrees,
            analysis,
            synthesis,
            legacy_frame_corrections,
        })
    }

    #[inline]
    pub fn input_size(&self) -> usize {
        self.input_directions * 3
    }

    #[inline]
    pub fn output_size(&self) -> usize {
        self.output_directions * 3
    }

    #[inline]
    pub fn num_coefficients(&self) -> usize {
        self.num_coefficients
    }

    #[inline]
    pub fn num_modes(&self) -> usize {
        self.mode_degrees.len()
    }

    fn apply(
        &self,
        input: &[f64],
        coefficients: &[f64],
        output: &mut [f64],
        moments: &mut [[Complex64; 3]],
    ) {
        debug_assert_eq!(input.len(), self.input_size());
        debug_assert_eq!(coefficients.len(), self.num_coefficients * 4);
        debug_assert_eq!(output.len(), self.output_size());
        debug_assert_eq!(moments.len(), self.num_modes());

        for (mode_index, moment) in moments.iter_mut().enumerate() {
            let mut intensity = Complex64::new(0.0, 0.0);
            let mut plus = Complex64::new(0.0, 0.0);
            let mut minus = Complex64::new(0.0, 0.0);
            let row_start = mode_index * self.input_directions;
            for input_index in 0..self.input_directions {
                let stokes_start = input_index * 3;
                let q = input[stokes_start + 1];
                let u = input[stokes_start + 2];
                intensity += self.analysis[0][row_start + input_index] * input[stokes_start];
                plus += self.analysis[1][row_start + input_index] * Complex64::new(q, u);
                minus += self.analysis[2][row_start + input_index] * Complex64::new(q, -u);
            }

            let coefficient_start = self.mode_degrees[mode_index] * 4;
            let a1 = coefficients[coefficient_start];
            let a2 = coefficients[coefficient_start + 1];
            let a3 = coefficients[coefficient_start + 2];
            let b1 = coefficients[coefficient_start + 3];
            *moment = mix_vector_moment(intensity, plus, minus, a1, a2, a3, b1);
        }

        for output_index in 0..self.output_directions {
            let mut intensity = Complex64::new(0.0, 0.0);
            let mut plus = Complex64::new(0.0, 0.0);
            let mut minus = Complex64::new(0.0, 0.0);
            let row_start = output_index * self.num_modes();
            for (mode_index, moment) in moments.iter().enumerate() {
                intensity += self.synthesis[0][row_start + mode_index] * moment[0];
                plus += self.synthesis[1][row_start + mode_index] * moment[1];
                minus += self.synthesis[2][row_start + mode_index] * moment[2];
            }
            let stokes_start = output_index * 3;
            output[stokes_start] = intensity.re;
            output[stokes_start + 1] = 0.5 * (plus.re + minus.re);
            output[stokes_start + 2] = 0.5 * (plus.im - minus.im);
        }

        for correction in &self.legacy_frame_corrections {
            let input_start = correction.input_index * 3;
            let output_start = correction.output_index * 3;
            for degree in 0..self.num_coefficients {
                for family in 0..4 {
                    let coefficient = coefficients[degree * 4 + family];
                    if coefficient == 0.0 {
                        continue;
                    }
                    let matrix_start = (degree * 4 + family) * 9;
                    for row in 0..3 {
                        output[output_start + row] += coefficient
                            * correction.values
                                [matrix_start + row * 3..matrix_start + (row + 1) * 3]
                                .iter()
                                .zip(&input[input_start..input_start + 3])
                                .map(|(&value, &radiance)| value * radiance)
                                .sum::<f64>();
                    }
                }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn apply_vjp(
        &self,
        input: &[f64],
        coefficients: &[f64],
        output_cotangent: &[f64],
        input_cotangent: &mut [f64],
        coefficient_gradient: &mut [f64],
        analyzed: &mut [[Complex64; 3]],
        moment_cotangent: &mut [[Complex64; 3]],
    ) {
        debug_assert_eq!(input.len(), self.input_size());
        debug_assert_eq!(coefficients.len(), self.num_coefficients * 4);
        debug_assert_eq!(output_cotangent.len(), self.output_size());
        debug_assert_eq!(input_cotangent.len(), self.input_size());
        debug_assert_eq!(coefficient_gradient.len(), coefficients.len());

        for mode_index in 0..self.num_modes() {
            let analysis_start = mode_index * self.input_directions;
            let mut intensity = Complex64::new(0.0, 0.0);
            let mut plus = Complex64::new(0.0, 0.0);
            let mut minus = Complex64::new(0.0, 0.0);
            for input_index in 0..self.input_directions {
                let stokes_start = input_index * 3;
                let q = input[stokes_start + 1];
                let u = input[stokes_start + 2];
                intensity += self.analysis[0][analysis_start + input_index] * input[stokes_start];
                plus += self.analysis[1][analysis_start + input_index] * Complex64::new(q, u);
                minus += self.analysis[2][analysis_start + input_index] * Complex64::new(q, -u);
            }
            analyzed[mode_index] = [intensity, plus, minus];

            let mut outgoing_adjoint = [Complex64::new(0.0, 0.0); 3];
            for output_index in 0..self.output_directions {
                let stokes_start = output_index * 3;
                let intensity_adjoint = Complex64::new(output_cotangent[stokes_start], 0.0);
                let plus_adjoint = Complex64::new(
                    0.5 * output_cotangent[stokes_start + 1],
                    0.5 * output_cotangent[stokes_start + 2],
                );
                let minus_adjoint = Complex64::new(
                    0.5 * output_cotangent[stokes_start + 1],
                    -0.5 * output_cotangent[stokes_start + 2],
                );
                let synthesis_index = output_index * self.num_modes() + mode_index;
                outgoing_adjoint[0] +=
                    self.synthesis[0][synthesis_index].conj() * intensity_adjoint;
                outgoing_adjoint[1] += self.synthesis[1][synthesis_index].conj() * plus_adjoint;
                outgoing_adjoint[2] += self.synthesis[2][synthesis_index].conj() * minus_adjoint;
            }
            moment_cotangent[mode_index] = outgoing_adjoint;

            let [intensity, plus, minus] = analyzed[mode_index];
            let [intensity_adjoint, plus_adjoint, minus_adjoint] = outgoing_adjoint;
            let coefficient_start = self.mode_degrees[mode_index] * 4;
            coefficient_gradient[coefficient_start] += (intensity_adjoint.conj() * intensity).re;
            coefficient_gradient[coefficient_start + 1] += 0.5
                * ((plus_adjoint.conj() * (plus + minus)).re
                    + (minus_adjoint.conj() * (plus + minus)).re);
            coefficient_gradient[coefficient_start + 2] += 0.5
                * ((plus_adjoint.conj() * (plus - minus)).re
                    + (minus_adjoint.conj() * (minus - plus)).re);
            coefficient_gradient[coefficient_start + 3] +=
                (intensity_adjoint.conj() * (-0.5 * (plus + minus))).re
                    + (plus_adjoint.conj() * -intensity).re
                    + (minus_adjoint.conj() * -intensity).re;

            let a1 = coefficients[coefficient_start];
            let a2 = coefficients[coefficient_start + 1];
            let a3 = coefficients[coefficient_start + 2];
            let b1 = coefficients[coefficient_start + 3];
            let analyzed_adjoint = [
                a1 * intensity_adjoint - b1 * (plus_adjoint + minus_adjoint),
                -0.5 * b1 * intensity_adjoint
                    + 0.5 * ((a2 + a3) * plus_adjoint + (a2 - a3) * minus_adjoint),
                -0.5 * b1 * intensity_adjoint
                    + 0.5 * ((a2 - a3) * plus_adjoint + (a2 + a3) * minus_adjoint),
            ];
            for input_index in 0..self.input_directions {
                let stokes_start = input_index * 3;
                let intensity_basis = self.analysis[0][analysis_start + input_index];
                let plus_basis = self.analysis[1][analysis_start + input_index];
                let minus_basis = self.analysis[2][analysis_start + input_index];
                input_cotangent[stokes_start] += (analyzed_adjoint[0].conj() * intensity_basis).re;
                input_cotangent[stokes_start + 1] += (analyzed_adjoint[1].conj() * plus_basis).re
                    + (analyzed_adjoint[2].conj() * minus_basis).re;
                input_cotangent[stokes_start + 2] += (analyzed_adjoint[1].conj()
                    * (Complex64::new(0.0, 1.0) * plus_basis))
                    .re
                    + (analyzed_adjoint[2].conj() * (Complex64::new(0.0, -1.0) * minus_basis)).re;
            }
        }

        for correction in &self.legacy_frame_corrections {
            let input_start = correction.input_index * 3;
            let output_start = correction.output_index * 3;
            for degree in 0..self.num_coefficients {
                for family in 0..4 {
                    let coefficient_index = degree * 4 + family;
                    let matrix_start = coefficient_index * 9;
                    for row in 0..3 {
                        for column in 0..3 {
                            let value = correction.values[matrix_start + row * 3 + column];
                            input_cotangent[input_start + column] += coefficients
                                [coefficient_index]
                                * value
                                * output_cotangent[output_start + row];
                            coefficient_gradient[coefficient_index] += value
                                * input[input_start + column]
                                * output_cotangent[output_start + row];
                        }
                    }
                }
            }
        }
    }
}

/// Polarized coefficient-space atmospheric blocks followed by optional dense
/// boundary blocks.
#[derive(Debug, Clone, PartialEq)]
pub struct VectorCoefficientScattering {
    output_size: usize,
    input_size: usize,
    output_offsets: Vec<usize>,
    input_offsets: Vec<usize>,
    coefficient_blocks: usize,
    basis: VectorCoefficientBasis,
    coefficients: Vec<f64>,
    dense_value_offsets: Vec<usize>,
    dense_values: Vec<f64>,
}

impl VectorCoefficientScattering {
    pub fn new(
        output_offsets: Vec<usize>,
        input_offsets: Vec<usize>,
        coefficient_blocks: usize,
        basis: VectorCoefficientBasis,
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
            .checked_mul(basis.num_coefficients() * 4)
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
        let coefficient_stride = self.basis.num_coefficients() * 4;
        let mut moments = vec![[Complex64::new(0.0, 0.0); 3]; self.basis.num_modes()];
        for block in 0..self.coefficient_blocks {
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let coefficient_start = block * coefficient_stride;
            self.basis.apply(
                &input[input_start..input_end],
                &self.coefficients[coefficient_start..coefficient_start + coefficient_stride],
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

    #[allow(clippy::too_many_arguments)]
    pub fn apply_jvp(
        &self,
        input: &[f64],
        input_tangent: &[f64],
        coefficient_tangent: &[f64],
        dense_value_tangent: &[f64],
        output_tangent: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.input_size
            || input_tangent.len() != self.input_size
            || coefficient_tangent.len() != self.coefficients.len()
            || dense_value_tangent.len() != self.dense_values.len()
            || output_tangent.len() != self.output_size
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let coefficient_stride = self.basis.num_coefficients() * 4;
        let mut moments = vec![[Complex64::new(0.0, 0.0); 3]; self.basis.num_modes()];
        let mut parameter_output = vec![0.0; self.basis.output_size()];
        for block in 0..self.coefficient_blocks {
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let coefficient_start = block * coefficient_stride;
            self.basis.apply(
                &input_tangent[input_start..input_end],
                &self.coefficients[coefficient_start..coefficient_start + coefficient_stride],
                &mut output_tangent[output_start..output_end],
                &mut moments,
            );
            self.basis.apply(
                &input[input_start..input_end],
                &coefficient_tangent[coefficient_start..coefficient_start + coefficient_stride],
                &mut parameter_output,
                &mut moments,
            );
            for (result, parameter) in output_tangent[output_start..output_end]
                .iter_mut()
                .zip(&parameter_output)
            {
                *result += parameter;
            }
        }
        for block in self.coefficient_blocks..self.output_offsets.len() - 1 {
            let dense_block = block - self.coefficient_blocks;
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let columns = input_end - input_start;
            let value_start = self.dense_value_offsets[dense_block];
            for (local_row, result) in output_tangent[output_start..output_end]
                .iter_mut()
                .enumerate()
            {
                let row_start = value_start + local_row * columns;
                *result = self.dense_values[row_start..row_start + columns]
                    .iter()
                    .zip(&input_tangent[input_start..input_end])
                    .zip(&dense_value_tangent[row_start..row_start + columns])
                    .zip(&input[input_start..input_end])
                    .map(
                        |(((&value, &incoming_tangent), &value_tangent), &incoming)| {
                            value * incoming_tangent + value_tangent * incoming
                        },
                    )
                    .sum();
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply_vjp(
        &self,
        input: &[f64],
        output_cotangent: &[f64],
        input_cotangent: &mut [f64],
        coefficient_gradient: &mut [f64],
        dense_value_gradient: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.input_size
            || output_cotangent.len() != self.output_size
            || input_cotangent.len() != self.input_size
            || coefficient_gradient.len() != self.coefficients.len()
            || dense_value_gradient.len() != self.dense_values.len()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        input_cotangent.fill(0.0);
        coefficient_gradient.fill(0.0);
        dense_value_gradient.fill(0.0);
        let coefficient_stride = self.basis.num_coefficients() * 4;
        let mut analyzed = vec![[Complex64::new(0.0, 0.0); 3]; self.basis.num_modes()];
        let mut moment_cotangent = vec![[Complex64::new(0.0, 0.0); 3]; self.basis.num_modes()];
        for block in 0..self.coefficient_blocks {
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let coefficient_start = block * coefficient_stride;
            self.basis.apply_vjp(
                &input[input_start..input_end],
                &self.coefficients[coefficient_start..coefficient_start + coefficient_stride],
                &output_cotangent[output_start..output_end],
                &mut input_cotangent[input_start..input_end],
                &mut coefficient_gradient
                    [coefficient_start..coefficient_start + coefficient_stride],
                &mut analyzed,
                &mut moment_cotangent,
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
            for (local_row, &outgoing_cotangent) in output_cotangent[output_start..output_end]
                .iter()
                .enumerate()
            {
                let row_start = value_start + local_row * columns;
                for local_column in 0..columns {
                    input_cotangent[input_start + local_column] +=
                        self.dense_values[row_start + local_column] * outgoing_cotangent;
                    dense_value_gradient[row_start + local_column] +=
                        input[input_start + local_column] * outgoing_cotangent;
                }
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

fn fill_vector_basis(
    directions: &[[f64; 3]],
    num_coefficients: usize,
    values: &mut [Vec<Complex64>; 3],
    mode_degrees: &mut [usize],
    transpose_and_conjugate: bool,
    weights: Option<&[f64]>,
) {
    let num_modes = num_coefficients * num_coefficients;
    for (direction_index, direction) in directions.iter().enumerate() {
        let theta = direction[2].clamp(-1.0, 1.0).acos();
        let phi = direction[1].atan2(direction[0]);
        for degree in 0..num_coefficients {
            for m in -(degree as i32)..=degree as i32 {
                let mode_index = degree * degree + (m + degree as i32) as usize;
                mode_degrees[mode_index] = degree;
                let phase = Complex64::from_polar(1.0, m as f64 * phi);
                for (channel, spin) in [0, 2, -2].into_iter().enumerate() {
                    let value = WignerDCalculator::new(m, spin).d(theta, degree as i32) * phase;
                    if transpose_and_conjugate {
                        values[channel][mode_index * directions.len() + direction_index] =
                            value.conj() * weights.unwrap()[direction_index];
                    } else {
                        values[channel][direction_index * num_modes + mode_index] = value;
                    }
                }
            }
        }
    }
}

#[inline]
fn mix_vector_moment(
    intensity: Complex64,
    plus: Complex64,
    minus: Complex64,
    a1: f64,
    a2: f64,
    a3: f64,
    b1: f64,
) -> [Complex64; 3] {
    [
        a1 * intensity - 0.5 * b1 * (plus + minus),
        -b1 * intensity + 0.5 * ((a2 + a3) * plus + (a2 - a3) * minus),
        -b1 * intensity + 0.5 * ((a2 - a3) * plus + (a2 + a3) * minus),
    ]
}

fn build_legacy_frame_corrections(
    incoming_directions: &[[f64; 3]],
    incoming_weights: &[f64],
    outgoing_directions: &[[f64; 3]],
    num_coefficients: usize,
    analysis: &[Vec<Complex64>; 3],
    synthesis: &[Vec<Complex64>; 3],
) -> Vec<VectorKernelCorrection> {
    let num_modes = num_coefficients * num_coefficients;
    let mut corrections = Vec::new();
    for (output_index, &outgoing) in outgoing_directions.iter().enumerate() {
        for (input_index, &incoming) in incoming_directions.iter().enumerate() {
            if !uses_legacy_frame_special_case(incoming, outgoing) {
                continue;
            }
            let mut values = vec![0.0; num_coefficients * 4 * 9];
            for degree in 0..num_coefficients {
                let mode_range = degree * degree..(degree + 1) * (degree + 1);
                for family in 0..4 {
                    let mut coefficients = [0.0; 4];
                    coefficients[family] = 1.0;
                    let expected =
                        legacy_polarized_greek_matrix(incoming, outgoing, degree, coefficients);
                    let matrix_start = (degree * 4 + family) * 9;
                    for input_stokes in 0..3 {
                        let mut result = [Complex64::new(0.0, 0.0); 3];
                        for mode_index in mode_range.clone() {
                            let analysis_index =
                                mode_index * incoming_directions.len() + input_index;
                            let (intensity, plus, minus) = match input_stokes {
                                0 => (
                                    analysis[0][analysis_index],
                                    Complex64::new(0.0, 0.0),
                                    Complex64::new(0.0, 0.0),
                                ),
                                1 => (
                                    Complex64::new(0.0, 0.0),
                                    analysis[1][analysis_index],
                                    analysis[2][analysis_index],
                                ),
                                2 => (
                                    Complex64::new(0.0, 0.0),
                                    Complex64::new(0.0, 1.0) * analysis[1][analysis_index],
                                    Complex64::new(0.0, -1.0) * analysis[2][analysis_index],
                                ),
                                _ => unreachable!(),
                            };
                            let moment = mix_vector_moment(
                                intensity,
                                plus,
                                minus,
                                coefficients[0],
                                coefficients[1],
                                coefficients[2],
                                coefficients[3],
                            );
                            let synthesis_start = output_index * num_modes + mode_index;
                            for channel in 0..3 {
                                result[channel] +=
                                    synthesis[channel][synthesis_start] * moment[channel];
                            }
                        }
                        let factored = [
                            result[0].re,
                            0.5 * (result[1].re + result[2].re),
                            0.5 * (result[1].im - result[2].im),
                        ];
                        for output_stokes in 0..3 {
                            values[matrix_start + output_stokes * 3 + input_stokes] =
                                incoming_weights[input_index]
                                    * expected[output_stokes][input_stokes]
                                    - factored[output_stokes];
                        }
                    }
                }
            }
            if values.iter().any(|value| value.abs() > 1.0e-15) {
                corrections.push(VectorKernelCorrection {
                    input_index,
                    output_index,
                    values,
                });
            }
        }
    }
    corrections
}

#[inline]
fn uses_legacy_frame_special_case(incoming: [f64; 3], outgoing: [f64; 3]) -> bool {
    let cosine = incoming
        .iter()
        .zip(outgoing)
        .map(|(&left, right)| left * right)
        .sum::<f64>()
        .clamp(-1.0, 1.0);
    // Include a narrow guard band around the legacy 1e-8 branch. Its direct
    // angle formulas become ill-conditioned just outside the branch, whereas
    // the harmonic form remains well-conditioned.
    (1.0 - cosine * cosine).sqrt() < 1.0e-6
        || incoming[0].hypot(incoming[1]) < 1.0e-8
        || outgoing[0].hypot(outgoing[1]) < 1.0e-8
}

fn legacy_polarized_greek_matrix(
    incoming: [f64; 3],
    outgoing: [f64; 3],
    degree: usize,
    coefficients: [f64; 4],
) -> [[f64; 3]; 3] {
    let incoming = [-incoming[0], -incoming[1], -incoming[2]];
    let outgoing = [-outgoing[0], -outgoing[1], -outgoing[2]];
    let (theta, c1, c2, s1, s2) = legacy_stokes_scattering_factors(incoming, outgoing);
    let d00 = WignerDCalculator::new(0, 0).d(theta, degree as i32);
    let d22 = WignerDCalculator::new(2, 2).d(theta, degree as i32);
    let d02 = WignerDCalculator::new(0, 2).d(theta, degree as i32);
    let d2m2 = WignerDCalculator::new(2, -2).d(theta, degree as i32);
    let w_add = d22 + d2m2;
    let w_minus = d22 - d2m2;
    let [a1, a2, a3, b1] = coefficients;

    let mut result = [[0.0; 3]; 3];
    result[0][0] = a1 * d00;
    result[1][0] = -b1 * c2 * d02;
    result[2][0] = -b1 * s2 * d02;
    result[0][1] = -b1 * c1 * d02;
    result[0][2] = b1 * s1 * d02;
    result[1][1] = 0.5
        * (a2 * (c1 * c2 * w_add - s1 * s2 * w_minus) + a3 * (c1 * c2 * w_minus - s1 * s2 * w_add));
    result[2][1] = 0.5
        * (a2 * (c1 * s2 * w_add + s1 * c2 * w_minus) + a3 * (c1 * s2 * w_minus + s1 * c2 * w_add));
    result[1][2] = -0.5
        * (a2 * (s1 * c2 * w_add + c1 * s2 * w_minus) + a3 * (s1 * c2 * w_minus + c1 * s2 * w_add));
    result[2][2] = 0.5
        * (a2 * (-s1 * s2 * w_add + c1 * c2 * w_minus)
            + a3 * (-s1 * s2 * w_minus + c1 * c2 * w_add));
    result
}

fn legacy_stokes_scattering_factors(
    incoming: [f64; 3],
    outgoing: [f64; 3],
) -> (f64, f64, f64, f64, f64) {
    let cosine = incoming
        .iter()
        .zip(outgoing)
        .map(|(&left, right)| left * right)
        .sum::<f64>()
        .clamp(-1.0, 1.0);
    let theta = cosine.acos();
    let sin_scatter = theta.sin();
    if sin_scatter.abs() < 1.0e-8 {
        return (theta, 1.0, 1.0, 0.0, 0.0);
    }

    let cosine_incoming = incoming[2].clamp(-1.0, 1.0);
    let cosine_outgoing = outgoing[2].clamp(-1.0, 1.0);
    let sin_incoming = cosine_incoming.acos().sin();
    let sin_outgoing = cosine_outgoing.acos().sin();
    if sin_incoming.abs() < 1.0e-8 || sin_outgoing.abs() < 1.0e-8 {
        return (theta, 1.0, 1.0, 0.0, 0.0);
    }

    let cosine1 = ((cosine_outgoing - cosine_incoming * cosine) / (sin_incoming * sin_scatter))
        .clamp(-1.0, 1.0);
    let cosine2 = ((cosine_incoming - cosine_outgoing * cosine) / (sin_outgoing * sin_scatter))
        .clamp(-1.0, 1.0);
    let phi_incoming = incoming[1].atan2(incoming[0]);
    let phi_outgoing = outgoing[1].atan2(outgoing[0]);
    let mut sine1 = 2.0 * (1.0 - cosine1 * cosine1).sqrt() * cosine1;
    let mut sine2 = 2.0 * (1.0 - cosine2 * cosine2).sqrt() * cosine2;
    let mut phi_difference = phi_incoming - phi_outgoing;
    if phi_difference < 0.0 {
        phi_difference += 2.0 * std::f64::consts::PI;
    }
    if phi_difference < std::f64::consts::PI {
        sine1 *= -1.0;
        sine2 *= -1.0;
    }
    (
        theta,
        2.0 * cosine1 * cosine1 - 1.0,
        2.0 * cosine2 * cosine2 - 1.0,
        sine1,
        sine2,
    )
}
