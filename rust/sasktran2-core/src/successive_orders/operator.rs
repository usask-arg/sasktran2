use std::fmt::{Display, Formatter};

use super::{
    ScalarCoefficientScattering, VectorCoefficientScattering,
    simd::batch::{Batch4, LANES},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperatorError {
    InvalidRowOffsets,
    InvalidBlockOffsets,
    ColumnOutOfBounds,
    DimensionMismatch,
    NonFiniteValue,
    UnsupportedOperator,
}

impl Display for OperatorError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidRowOffsets => formatter.write_str("invalid compressed row offsets"),
            Self::InvalidBlockOffsets => formatter.write_str("invalid block diagonal offsets"),
            Self::ColumnOutOfBounds => {
                formatter.write_str("operator column index is out of bounds")
            }
            Self::DimensionMismatch => formatter.write_str("operator dimensions do not match"),
            Self::NonFiniteValue => formatter.write_str("operator contains a non-finite value"),
            Self::UnsupportedOperator => {
                formatter.write_str("operation is not supported by this operator")
            }
        }
    }
}

impl std::error::Error for OperatorError {}

/// Sparse transport operator metadata.
///
/// The row offsets, column indices, and values remain owned by the caller and
/// are borrowed for each solve. This keeps the C++/Rust integration from
/// retaining a second copy of the transport operator.
#[derive(Debug, Clone, PartialEq)]
pub struct CsrMatrix {
    rows: usize,
    columns: usize,
    num_nonzero: usize,
}

impl CsrMatrix {
    pub fn new(
        rows: usize,
        columns: usize,
        row_offsets: &[i32],
        column_indices: &[i32],
    ) -> Result<Self, OperatorError> {
        validate_csr_structure(rows, columns, row_offsets, column_indices)?;
        Ok(Self {
            rows,
            columns,
            num_nonzero: column_indices.len(),
        })
    }

    pub fn validate_data(
        &self,
        row_offsets: &[i32],
        column_indices: &[i32],
        values: &[f64],
    ) -> Result<(), OperatorError> {
        validate_csr_structure(self.rows, self.columns, row_offsets, column_indices)?;
        if column_indices.len() != self.num_nonzero || values.len() != self.num_nonzero {
            return Err(OperatorError::DimensionMismatch);
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        Ok(())
    }

    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    pub fn columns(&self) -> usize {
        self.columns
    }

    #[inline]
    pub fn num_nonzero(&self) -> usize {
        self.num_nonzero
    }

    pub fn apply(
        &self,
        row_offsets: &[i32],
        column_indices: &[i32],
        values: &[f64],
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.columns || output.len() != self.rows {
            return Err(OperatorError::DimensionMismatch);
        }
        for (row, result) in output.iter_mut().enumerate() {
            let start = row_offsets[row] as usize;
            let end = row_offsets[row + 1] as usize;
            *result = column_indices[start..end]
                .iter()
                .zip(&values[start..end])
                .map(|(&column, &value)| value * input[column as usize])
                .sum();
        }
        Ok(())
    }

    pub(crate) fn apply_batch4(
        &self,
        row_offsets: &[i32],
        column_indices: &[i32],
        values: &[f64],
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        if row_offsets.len() != self.rows + 1
            || column_indices.len() != self.num_nonzero
            || values.len() != self.num_nonzero * LANES
            || input.len() != self.columns * LANES
            || output.len() != self.rows * LANES
        {
            return Err(OperatorError::DimensionMismatch);
        }
        for row in 0..self.rows {
            let start = row_offsets[row] as usize;
            let end = row_offsets[row + 1] as usize;
            let mut result = Batch4::splat(0.0);
            for (index, &column) in column_indices.iter().enumerate().take(end).skip(start) {
                result =
                    result + Batch4::load(values, index) * Batch4::load(input, column as usize);
            }
            result.store(output, row);
        }
        Ok(())
    }

    pub fn apply_transpose(
        &self,
        row_offsets: &[i32],
        column_indices: &[i32],
        values: &[f64],
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.rows || output.len() != self.columns {
            return Err(OperatorError::DimensionMismatch);
        }
        debug_assert_eq!(row_offsets.len(), self.rows + 1);
        debug_assert_eq!(column_indices.len(), self.num_nonzero);
        debug_assert_eq!(values.len(), self.num_nonzero);
        output.fill(0.0);
        for (row, &incoming) in input.iter().enumerate() {
            let start = row_offsets[row] as usize;
            let end = row_offsets[row + 1] as usize;
            for index in start..end {
                output[column_indices[index] as usize] += values[index] * incoming;
            }
        }
        Ok(())
    }

    pub(crate) fn apply_transpose_batch4(
        &self,
        row_offsets: &[i32],
        column_indices: &[i32],
        values: &[f64],
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        if row_offsets.len() != self.rows + 1
            || column_indices.len() != self.num_nonzero
            || values.len() != self.num_nonzero * LANES
            || input.len() != self.rows * LANES
            || output.len() != self.columns * LANES
        {
            return Err(OperatorError::DimensionMismatch);
        }
        output.fill(0.0);
        for row in 0..self.rows {
            let incoming = Batch4::load(input, row);
            let start = row_offsets[row] as usize;
            let end = row_offsets[row + 1] as usize;
            for (index, &column) in column_indices.iter().enumerate().take(end).skip(start) {
                let column = column as usize;
                let updated = Batch4::load(output, column) + Batch4::load(values, index) * incoming;
                updated.store(output, column);
            }
        }
        Ok(())
    }

    pub(crate) fn apply_transpose_system_batch4(
        &self,
        row_offsets: &[i32],
        column_indices: &[i32],
        values: &[f64],
        fixed_point_input: &[f64],
        transpose_input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        if fixed_point_input.len() != self.columns * LANES
            || transpose_input.len() != self.rows * LANES
            || output.len() != self.columns * LANES
        {
            return Err(OperatorError::DimensionMismatch);
        }
        output.copy_from_slice(fixed_point_input);
        for row in 0..self.rows {
            let incoming = Batch4::load(transpose_input, row);
            let start = row_offsets[row] as usize;
            let end = row_offsets[row + 1] as usize;
            for (index, &column) in column_indices.iter().enumerate().take(end).skip(start) {
                let column = column as usize;
                let updated = Batch4::load(output, column) - Batch4::load(values, index) * incoming;
                updated.store(output, column);
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply_transpose_system(
        &self,
        row_offsets: &[i32],
        column_indices: &[i32],
        values: &[f64],
        fixed_point_input: &[f64],
        transpose_input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        if fixed_point_input.len() != self.columns
            || transpose_input.len() != self.rows
            || output.len() != self.columns
        {
            return Err(OperatorError::DimensionMismatch);
        }
        debug_assert_eq!(row_offsets.len(), self.rows + 1);
        debug_assert_eq!(column_indices.len(), self.num_nonzero);
        debug_assert_eq!(values.len(), self.num_nonzero);
        output.copy_from_slice(fixed_point_input);
        for (row, &incoming) in transpose_input.iter().enumerate() {
            let start = row_offsets[row] as usize;
            let end = row_offsets[row + 1] as usize;
            for index in start..end {
                output[column_indices[index] as usize] -= values[index] * incoming;
            }
        }
        Ok(())
    }

    pub fn apply_value_tangent(
        &self,
        row_offsets: &[i32],
        column_indices: &[i32],
        value_tangent: &[f64],
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_data(row_offsets, column_indices, value_tangent)?;
        self.apply(row_offsets, column_indices, value_tangent, input, output)
    }
}

fn validate_csr_structure(
    rows: usize,
    columns: usize,
    row_offsets: &[i32],
    column_indices: &[i32],
) -> Result<(), OperatorError> {
    if row_offsets.len() != rows + 1
        || row_offsets.first().copied() != Some(0)
        || row_offsets.last().copied() != i32::try_from(column_indices.len()).ok()
        || row_offsets
            .windows(2)
            .any(|window| window[0] < 0 || window[0] > window[1])
    {
        return Err(OperatorError::InvalidRowOffsets);
    }
    if column_indices
        .iter()
        .any(|&column| column < 0 || column as usize >= columns)
    {
        return Err(OperatorError::ColumnOutOfBounds);
    }
    Ok(())
}

/// Dense rectangular blocks arranged diagonally.
///
/// Each block maps one contiguous incoming-radiance range to one contiguous
/// outgoing-source range. Values are row-major within each block.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockDiagonalMatrix {
    output_size: usize,
    input_size: usize,
    output_offsets: Vec<usize>,
    input_offsets: Vec<usize>,
    value_offsets: Vec<usize>,
    values: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ScatteringOperator {
    Dense(BlockDiagonalMatrix),
    ScalarCoefficients(ScalarCoefficientScattering),
    VectorCoefficients(Box<VectorCoefficientScattering>),
}

impl ScatteringOperator {
    #[inline]
    pub fn output_size(&self) -> usize {
        match self {
            Self::Dense(matrix) => matrix.output_size(),
            Self::ScalarCoefficients(matrix) => matrix.output_size(),
            Self::VectorCoefficients(matrix) => matrix.output_size(),
        }
    }

    #[inline]
    pub fn input_size(&self) -> usize {
        match self {
            Self::Dense(matrix) => matrix.input_size(),
            Self::ScalarCoefficients(matrix) => matrix.input_size(),
            Self::VectorCoefficients(matrix) => matrix.input_size(),
        }
    }

    pub fn set_dense_values(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        match self {
            Self::Dense(matrix) => matrix.set_values(values),
            Self::ScalarCoefficients(matrix) => matrix.set_dense_values(values),
            Self::VectorCoefficients(matrix) => matrix.set_dense_values(values),
        }
    }

    pub fn set_coefficients(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        match self {
            Self::Dense(_) => Err(OperatorError::UnsupportedOperator),
            Self::ScalarCoefficients(matrix) => matrix.set_coefficients(values),
            Self::VectorCoefficients(matrix) => matrix.set_coefficients(values),
        }
    }

    #[inline]
    pub fn coefficient_value_size(&self) -> usize {
        match self {
            Self::Dense(_) => 0,
            Self::ScalarCoefficients(matrix) => matrix.coefficient_value_size(),
            Self::VectorCoefficients(matrix) => matrix.coefficient_value_size(),
        }
    }

    #[inline]
    pub fn dense_value_size(&self) -> usize {
        match self {
            Self::Dense(matrix) => matrix.value_size(),
            Self::ScalarCoefficients(matrix) => matrix.dense_value_size(),
            Self::VectorCoefficients(matrix) => matrix.dense_value_size(),
        }
    }

    #[inline]
    pub fn transpose_scratch_size(&self) -> usize {
        match self {
            Self::ScalarCoefficients(matrix) => matrix.transpose_scratch_size(),
            Self::Dense(_) | Self::VectorCoefficients(_) => 0,
        }
    }

    pub fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        match self {
            Self::Dense(matrix) => matrix.apply(input, output),
            Self::ScalarCoefficients(matrix) => matrix.apply(input, output),
            Self::VectorCoefficients(matrix) => matrix.apply(input, output),
        }
    }

    pub(crate) fn batch4_active_num_coefficients(
        &self,
        coefficients: &[f64],
    ) -> Result<usize, OperatorError> {
        match self {
            Self::ScalarCoefficients(matrix) => matrix.batch4_active_num_coefficients(coefficients),
            Self::Dense(_) | Self::VectorCoefficients(_) => Err(OperatorError::UnsupportedOperator),
        }
    }

    pub(crate) fn batch4_scratch_size(
        &self,
        active_num_coefficients: usize,
    ) -> Result<usize, OperatorError> {
        match self {
            Self::ScalarCoefficients(matrix) => {
                Ok(matrix.batch4_scratch_size(active_num_coefficients))
            }
            Self::Dense(_) | Self::VectorCoefficients(_) => Err(OperatorError::UnsupportedOperator),
        }
    }

    pub(crate) fn apply_batch4(
        &self,
        input: &[f64],
        coefficients: &[f64],
        dense_values: &[f64],
        active_num_coefficients: usize,
        output: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        match self {
            Self::ScalarCoefficients(matrix) => matrix.apply_batch4(
                input,
                coefficients,
                dense_values,
                active_num_coefficients,
                output,
                scratch,
            ),
            Self::Dense(_) | Self::VectorCoefficients(_) => Err(OperatorError::UnsupportedOperator),
        }
    }

    pub fn apply_transpose(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        let mut scratch = vec![0.0; self.transpose_scratch_size()];
        self.apply_transpose_with_scratch(input, output, &mut scratch)
    }

    pub fn apply_transpose_with_scratch(
        &self,
        input: &[f64],
        output: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if scratch.len() != self.transpose_scratch_size() {
            return Err(OperatorError::DimensionMismatch);
        }
        match self {
            Self::Dense(matrix) => matrix.apply_transpose(input, output),
            Self::ScalarCoefficients(matrix) => {
                matrix.apply_transpose_with_scratch(input, output, scratch)
            }
            Self::VectorCoefficients(matrix) => matrix.apply_transpose(input, output),
        }
    }

    pub(crate) fn apply_transpose_batch4(
        &self,
        input: &[f64],
        coefficients: &[f64],
        dense_values: &[f64],
        active_num_coefficients: usize,
        output: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        match self {
            Self::ScalarCoefficients(matrix) => matrix.apply_transpose_batch4(
                input,
                coefficients,
                dense_values,
                active_num_coefficients,
                output,
                scratch,
            ),
            Self::Dense(_) | Self::VectorCoefficients(_) => Err(OperatorError::UnsupportedOperator),
        }
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
        match self {
            Self::Dense(matrix) => {
                if !coefficient_tangent.is_empty() {
                    return Err(OperatorError::DimensionMismatch);
                }
                matrix.apply_jvp(input, input_tangent, dense_value_tangent, output_tangent)
            }
            Self::ScalarCoefficients(matrix) => matrix.apply_jvp(
                input,
                input_tangent,
                coefficient_tangent,
                dense_value_tangent,
                output_tangent,
            ),
            Self::VectorCoefficients(matrix) => matrix.apply_jvp(
                input,
                input_tangent,
                coefficient_tangent,
                dense_value_tangent,
                output_tangent,
            ),
        }
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
        match self {
            Self::Dense(matrix) => {
                if !coefficient_gradient.is_empty() {
                    return Err(OperatorError::DimensionMismatch);
                }
                matrix.apply_vjp(
                    input,
                    output_cotangent,
                    input_cotangent,
                    dense_value_gradient,
                )
            }
            Self::ScalarCoefficients(matrix) => matrix.apply_vjp(
                input,
                output_cotangent,
                input_cotangent,
                coefficient_gradient,
                dense_value_gradient,
            ),
            Self::VectorCoefficients(matrix) => matrix.apply_vjp(
                input,
                output_cotangent,
                input_cotangent,
                coefficient_gradient,
                dense_value_gradient,
            ),
        }
    }
}

impl From<BlockDiagonalMatrix> for ScatteringOperator {
    fn from(value: BlockDiagonalMatrix) -> Self {
        Self::Dense(value)
    }
}

impl From<ScalarCoefficientScattering> for ScatteringOperator {
    fn from(value: ScalarCoefficientScattering) -> Self {
        Self::ScalarCoefficients(value)
    }
}

impl From<VectorCoefficientScattering> for ScatteringOperator {
    fn from(value: VectorCoefficientScattering) -> Self {
        Self::VectorCoefficients(Box::new(value))
    }
}

impl BlockDiagonalMatrix {
    pub fn new(
        output_offsets: Vec<usize>,
        input_offsets: Vec<usize>,
        value_offsets: Vec<usize>,
        values: Vec<f64>,
    ) -> Result<Self, OperatorError> {
        let num_offsets = output_offsets.len();
        if num_offsets < 2
            || input_offsets.len() != num_offsets
            || value_offsets.len() != num_offsets
            || output_offsets[0] != 0
            || input_offsets[0] != 0
            || value_offsets[0] != 0
            || output_offsets
                .windows(2)
                .any(|window| window[0] > window[1])
            || input_offsets.windows(2).any(|window| window[0] > window[1])
            || value_offsets.windows(2).any(|window| window[0] > window[1])
            || value_offsets.last().copied() != Some(values.len())
        {
            return Err(OperatorError::InvalidBlockOffsets);
        }
        for block in 0..num_offsets - 1 {
            let rows = output_offsets[block + 1] - output_offsets[block];
            let columns = input_offsets[block + 1] - input_offsets[block];
            if value_offsets[block + 1] - value_offsets[block] != rows * columns {
                return Err(OperatorError::InvalidBlockOffsets);
            }
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        Ok(Self {
            output_size: *output_offsets.last().unwrap(),
            input_size: *input_offsets.last().unwrap(),
            output_offsets,
            input_offsets,
            value_offsets,
            values,
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
    pub fn value_size(&self) -> usize {
        self.values.len()
    }

    pub fn set_values(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        if values.len() != self.values.len() {
            return Err(OperatorError::DimensionMismatch);
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        self.values.copy_from_slice(values);
        Ok(())
    }

    pub fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        if input.len() != self.input_size || output.len() != self.output_size {
            return Err(OperatorError::DimensionMismatch);
        }
        for block in 0..self.output_offsets.len() - 1 {
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let value_start = self.value_offsets[block];
            let columns = input_end - input_start;
            for (local_row, result) in output[output_start..output_end].iter_mut().enumerate() {
                let row_start = value_start + local_row * columns;
                *result = self.values[row_start..row_start + columns]
                    .iter()
                    .zip(&input[input_start..input_end])
                    .map(|(&value, &incoming)| value * incoming)
                    .sum();
            }
        }
        Ok(())
    }

    pub fn apply_transpose(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        if input.len() != self.output_size || output.len() != self.input_size {
            return Err(OperatorError::DimensionMismatch);
        }
        output.fill(0.0);
        for block in 0..self.output_offsets.len() - 1 {
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let value_start = self.value_offsets[block];
            let columns = input_end - input_start;
            for (local_row, &input_value) in input[output_start..output_end].iter().enumerate() {
                let row_start = value_start + local_row * columns;
                for local_column in 0..columns {
                    output[input_start + local_column] +=
                        self.values[row_start + local_column] * input_value;
                }
            }
        }
        Ok(())
    }

    pub fn apply_jvp(
        &self,
        input: &[f64],
        input_tangent: &[f64],
        value_tangent: &[f64],
        output_tangent: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.input_size
            || input_tangent.len() != self.input_size
            || value_tangent.len() != self.values.len()
            || output_tangent.len() != self.output_size
        {
            return Err(OperatorError::DimensionMismatch);
        }
        for block in 0..self.output_offsets.len() - 1 {
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let value_start = self.value_offsets[block];
            let columns = input_end - input_start;
            for (local_row, result) in output_tangent[output_start..output_end]
                .iter_mut()
                .enumerate()
            {
                let row_start = value_start + local_row * columns;
                *result = self.values[row_start..row_start + columns]
                    .iter()
                    .zip(&input_tangent[input_start..input_end])
                    .zip(&value_tangent[row_start..row_start + columns])
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

    pub fn apply_vjp(
        &self,
        input: &[f64],
        output_cotangent: &[f64],
        input_cotangent: &mut [f64],
        value_gradient: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.input_size
            || output_cotangent.len() != self.output_size
            || input_cotangent.len() != self.input_size
            || value_gradient.len() != self.values.len()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        input_cotangent.fill(0.0);
        value_gradient.fill(0.0);
        for block in 0..self.output_offsets.len() - 1 {
            let output_start = self.output_offsets[block];
            let output_end = self.output_offsets[block + 1];
            let input_start = self.input_offsets[block];
            let input_end = self.input_offsets[block + 1];
            let value_start = self.value_offsets[block];
            let columns = input_end - input_start;
            for (local_row, &outgoing_cotangent) in output_cotangent[output_start..output_end]
                .iter()
                .enumerate()
            {
                let row_start = value_start + local_row * columns;
                for local_column in 0..columns {
                    input_cotangent[input_start + local_column] +=
                        self.values[row_start + local_column] * outgoing_cotangent;
                    value_gradient[row_start + local_column] +=
                        input[input_start + local_column] * outgoing_cotangent;
                }
            }
        }
        Ok(())
    }
}

/// The fixed-point map `G(x) = S (b + T x)`.
#[derive(Debug, Clone, PartialEq)]
pub struct FixedPointProblem {
    transport: CsrMatrix,
    scattering: ScatteringOperator,
}

impl FixedPointProblem {
    pub fn new<S: Into<ScatteringOperator>>(
        transport: CsrMatrix,
        scattering: S,
    ) -> Result<Self, OperatorError> {
        let scattering = scattering.into();
        if transport.rows() != scattering.input_size()
            || transport.columns() != scattering.output_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        Ok(Self {
            transport,
            scattering,
        })
    }

    #[inline]
    pub fn state_size(&self) -> usize {
        self.scattering.output_size()
    }

    #[inline]
    pub fn incoming_size(&self) -> usize {
        self.transport.rows()
    }

    pub fn set_scattering_values(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        self.scattering.set_dense_values(values)
    }

    pub fn set_scattering_coefficients(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        self.scattering.set_coefficients(values)
    }

    #[inline]
    pub fn transport_value_size(&self) -> usize {
        self.transport.num_nonzero()
    }

    #[inline]
    pub fn scattering_coefficient_size(&self) -> usize {
        self.scattering.coefficient_value_size()
    }

    #[inline]
    pub fn dense_scattering_value_size(&self) -> usize {
        self.scattering.dense_value_size()
    }

    #[inline]
    pub fn scattering_transpose_scratch_size(&self) -> usize {
        self.scattering.transpose_scratch_size()
    }

    pub fn validate_iteration_data(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
    ) -> Result<(), OperatorError> {
        self.transport.validate_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
        )?;
        if forcing.len() != self.transport.rows() {
            return Err(OperatorError::DimensionMismatch);
        }
        if forcing.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn validate_batch4_data(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        scattering_coefficients: &[f64],
        dense_scattering_values: &[f64],
        forcing: &[f64],
    ) -> Result<usize, OperatorError> {
        validate_csr_structure(
            self.transport.rows,
            self.transport.columns,
            transport_row_offsets,
            transport_column_indices,
        )?;
        if transport_column_indices.len() != self.transport.num_nonzero
            || transport_values.len() != self.transport.num_nonzero * LANES
            || dense_scattering_values.len() != self.scattering.dense_value_size() * LANES
            || forcing.len() != self.transport.rows() * LANES
            || transport_values.iter().any(|value| !value.is_finite())
            || dense_scattering_values
                .iter()
                .any(|value| !value.is_finite())
            || forcing.iter().any(|value| !value.is_finite())
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.scattering
            .batch4_active_num_coefficients(scattering_coefficients)
    }

    pub(crate) fn batch4_scattering_scratch_size(
        &self,
        active_num_coefficients: usize,
    ) -> Result<usize, OperatorError> {
        self.scattering.batch4_scratch_size(active_num_coefficients)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn apply_batch4(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        scattering_coefficients: &[f64],
        dense_scattering_values: &[f64],
        forcing: &[f64],
        active_num_coefficients: usize,
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
        scattering_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.transport.apply_batch4(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
        )?;
        for element in 0..self.incoming_size() {
            (Batch4::load(incoming_scratch, element) + Batch4::load(forcing, element))
                .store(incoming_scratch, element);
        }
        self.scattering.apply_batch4(
            incoming_scratch,
            scattering_coefficients,
            dense_scattering_values,
            active_num_coefficients,
            output,
            scattering_scratch,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn apply_linear_system_batch4(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        scattering_coefficients: &[f64],
        dense_scattering_values: &[f64],
        active_num_coefficients: usize,
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
        scattering_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.transport.apply_batch4(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
        )?;
        self.scattering.apply_batch4(
            incoming_scratch,
            scattering_coefficients,
            dense_scattering_values,
            active_num_coefficients,
            output,
            scattering_scratch,
        )?;
        for element in 0..self.state_size() {
            (Batch4::load(state, element) - Batch4::load(output, element)).store(output, element);
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn apply_linear_system_transpose_batch4(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        scattering_coefficients: &[f64],
        dense_scattering_values: &[f64],
        active_num_coefficients: usize,
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
        scattering_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if state.len() != self.state_size() * LANES
            || output.len() != self.state_size() * LANES
            || incoming_scratch.len() != self.incoming_size() * LANES
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.scattering.apply_transpose_batch4(
            state,
            scattering_coefficients,
            dense_scattering_values,
            active_num_coefficients,
            incoming_scratch,
            scattering_scratch,
        )?;
        self.transport.apply_transpose_system_batch4(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
            output,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn apply_linear_transpose_batch4(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        scattering_coefficients: &[f64],
        dense_scattering_values: &[f64],
        active_num_coefficients: usize,
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
        scattering_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if state.len() != self.state_size() * LANES
            || output.len() != self.state_size() * LANES
            || incoming_scratch.len() != self.incoming_size() * LANES
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.scattering.apply_transpose_batch4(
            state,
            scattering_coefficients,
            dense_scattering_values,
            active_num_coefficients,
            incoming_scratch,
            scattering_scratch,
        )?;
        self.transport.apply_transpose_batch4(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            incoming_scratch,
            output,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.transport.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
        )?;
        for (incoming, forcing) in incoming_scratch.iter_mut().zip(forcing) {
            *incoming += forcing;
        }
        self.scattering.apply(incoming_scratch, output)
    }

    pub fn apply_linear(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if state.len() != self.state_size()
            || output.len() != self.state_size()
            || incoming_scratch.len() != self.incoming_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.transport.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
        )?;
        self.scattering.apply(incoming_scratch, output)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply_linear_transpose(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
        scattering_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if state.len() != self.state_size()
            || output.len() != self.state_size()
            || incoming_scratch.len() != self.incoming_size()
            || scattering_scratch.len() != self.scattering_transpose_scratch_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.scattering.apply_transpose_with_scratch(
            state,
            incoming_scratch,
            scattering_scratch,
        )?;
        self.transport.apply_transpose(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            incoming_scratch,
            output,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply_linear_system_transpose(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
        scattering_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if state.len() != self.state_size()
            || output.len() != self.state_size()
            || incoming_scratch.len() != self.incoming_size()
            || scattering_scratch.len() != self.scattering_transpose_scratch_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.scattering.apply_transpose_with_scratch(
            state,
            incoming_scratch,
            scattering_scratch,
        )?;
        self.transport.apply_transpose_system(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
            output,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply_jvp(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        transport_value_tangent: &[f64],
        forcing: &[f64],
        forcing_tangent: &[f64],
        scattering_coefficient_tangent: &[f64],
        dense_scattering_value_tangent: &[f64],
        state: &[f64],
        state_tangent: &[f64],
        output: &mut [f64],
        output_tangent: &mut [f64],
        incoming_scratch: &mut [f64],
        incoming_tangent_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_iteration_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
        )?;
        if transport_value_tangent.len() != self.transport.num_nonzero()
            || forcing_tangent.len() != self.incoming_size()
            || scattering_coefficient_tangent.len() != self.scattering.coefficient_value_size()
            || dense_scattering_value_tangent.len() != self.scattering.dense_value_size()
            || state.len() != self.state_size()
            || state_tangent.len() != self.state_size()
            || output.len() != self.state_size()
            || output_tangent.len() != self.state_size()
            || incoming_scratch.len() != self.incoming_size()
            || incoming_tangent_scratch.len() != self.incoming_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.transport.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
        )?;
        self.transport.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state_tangent,
            incoming_tangent_scratch,
        )?;
        let mut transport_parameter_tangent = vec![0.0; self.incoming_size()];
        self.transport.apply_value_tangent(
            transport_row_offsets,
            transport_column_indices,
            transport_value_tangent,
            state,
            &mut transport_parameter_tangent,
        )?;
        for index in 0..self.incoming_size() {
            incoming_scratch[index] += forcing[index];
            incoming_tangent_scratch[index] +=
                transport_parameter_tangent[index] + forcing_tangent[index];
        }
        self.scattering.apply(incoming_scratch, output)?;
        self.scattering.apply_jvp(
            incoming_scratch,
            incoming_tangent_scratch,
            scattering_coefficient_tangent,
            dense_scattering_value_tangent,
            output_tangent,
        )
    }

    /// Parameter-only tangent of `G(x) = S(b + Tx)` at a fixed state.
    ///
    /// Unlike `apply_jvp`, this omits the zero state-tangent transport and the
    /// unused primal output when forming the right-hand side of an implicit
    /// JVP solve.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_parameter_jvp(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        transport_value_tangent: &[f64],
        forcing: &[f64],
        forcing_tangent: &[f64],
        scattering_coefficient_tangent: &[f64],
        dense_scattering_value_tangent: &[f64],
        state: &[f64],
        output_tangent: &mut [f64],
        incoming_scratch: &mut [f64],
        incoming_tangent_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_iteration_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
        )?;
        if transport_value_tangent.len() != self.transport.num_nonzero()
            || forcing_tangent.len() != self.incoming_size()
            || scattering_coefficient_tangent.len() != self.scattering.coefficient_value_size()
            || dense_scattering_value_tangent.len() != self.scattering.dense_value_size()
            || state.len() != self.state_size()
            || output_tangent.len() != self.state_size()
            || incoming_scratch.len() != self.incoming_size()
            || incoming_tangent_scratch.len() != self.incoming_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.transport.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
        )?;
        self.transport.apply_value_tangent(
            transport_row_offsets,
            transport_column_indices,
            transport_value_tangent,
            state,
            incoming_tangent_scratch,
        )?;
        for index in 0..self.incoming_size() {
            incoming_scratch[index] += forcing[index];
            incoming_tangent_scratch[index] += forcing_tangent[index];
        }
        self.scattering.apply_jvp(
            incoming_scratch,
            incoming_tangent_scratch,
            scattering_coefficient_tangent,
            dense_scattering_value_tangent,
            output_tangent,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply_vjp(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
        state: &[f64],
        output_cotangent: &[f64],
        state_cotangent: &mut [f64],
        transport_value_gradient: &mut [f64],
        scattering_coefficient_gradient: &mut [f64],
        dense_scattering_value_gradient: &mut [f64],
        forcing_gradient: &mut [f64],
        incoming_scratch: &mut [f64],
        incoming_cotangent_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_iteration_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
        )?;
        if state.len() != self.state_size()
            || output_cotangent.len() != self.state_size()
            || state_cotangent.len() != self.state_size()
            || (!transport_value_gradient.is_empty()
                && transport_value_gradient.len() != self.transport.num_nonzero())
            || scattering_coefficient_gradient.len() != self.scattering.coefficient_value_size()
            || dense_scattering_value_gradient.len() != self.scattering.dense_value_size()
            || forcing_gradient.len() != self.incoming_size()
            || incoming_scratch.len() != self.incoming_size()
            || incoming_cotangent_scratch.len() != self.incoming_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.transport.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
        )?;
        for (incoming, forcing) in incoming_scratch.iter_mut().zip(forcing) {
            *incoming += forcing;
        }
        let mut local_coefficient_gradient = vec![0.0; scattering_coefficient_gradient.len()];
        let mut local_dense_gradient = vec![0.0; dense_scattering_value_gradient.len()];
        self.scattering.apply_vjp(
            incoming_scratch,
            output_cotangent,
            incoming_cotangent_scratch,
            &mut local_coefficient_gradient,
            &mut local_dense_gradient,
        )?;
        for (gradient, local) in scattering_coefficient_gradient
            .iter_mut()
            .zip(local_coefficient_gradient)
        {
            *gradient += local;
        }
        for (gradient, local) in dense_scattering_value_gradient
            .iter_mut()
            .zip(local_dense_gradient)
        {
            *gradient += local;
        }
        for (gradient, &incoming_cotangent) in forcing_gradient
            .iter_mut()
            .zip(incoming_cotangent_scratch.iter())
        {
            *gradient += incoming_cotangent;
        }
        self.transport.apply_transpose(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            incoming_cotangent_scratch,
            state_cotangent,
        )?;
        if !transport_value_gradient.is_empty() {
            for (row, &incoming_cotangent) in incoming_cotangent_scratch.iter().enumerate() {
                let start = transport_row_offsets[row] as usize;
                let end = transport_row_offsets[row + 1] as usize;
                for index in start..end {
                    transport_value_gradient[index] +=
                        incoming_cotangent * state[transport_column_indices[index] as usize];
                }
            }
        }
        Ok(())
    }

    /// Parameter-only reverse derivative at a fixed state.
    ///
    /// The implicit adjoint has already propagated the state cotangent, so
    /// this omits the otherwise unused `T^T` state pullback.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_parameter_vjp(
        &self,
        transport_row_offsets: &[i32],
        transport_column_indices: &[i32],
        transport_values: &[f64],
        forcing: &[f64],
        state: &[f64],
        output_cotangent: &[f64],
        transport_value_gradient: &mut [f64],
        scattering_coefficient_gradient: &mut [f64],
        dense_scattering_value_gradient: &mut [f64],
        forcing_gradient: &mut [f64],
        incoming_scratch: &mut [f64],
        incoming_cotangent_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_iteration_data(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            forcing,
        )?;
        if state.len() != self.state_size()
            || output_cotangent.len() != self.state_size()
            || (!transport_value_gradient.is_empty()
                && transport_value_gradient.len() != self.transport.num_nonzero())
            || scattering_coefficient_gradient.len() != self.scattering.coefficient_value_size()
            || dense_scattering_value_gradient.len() != self.scattering.dense_value_size()
            || forcing_gradient.len() != self.incoming_size()
            || incoming_scratch.len() != self.incoming_size()
            || incoming_cotangent_scratch.len() != self.incoming_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.transport.apply(
            transport_row_offsets,
            transport_column_indices,
            transport_values,
            state,
            incoming_scratch,
        )?;
        for (incoming, forcing) in incoming_scratch.iter_mut().zip(forcing) {
            *incoming += forcing;
        }
        self.scattering.apply_vjp(
            incoming_scratch,
            output_cotangent,
            incoming_cotangent_scratch,
            scattering_coefficient_gradient,
            dense_scattering_value_gradient,
        )?;
        forcing_gradient.copy_from_slice(incoming_cotangent_scratch);
        if !transport_value_gradient.is_empty() {
            for (row, &incoming_cotangent) in incoming_cotangent_scratch.iter().enumerate() {
                let start = transport_row_offsets[row] as usize;
                let end = transport_row_offsets[row + 1] as usize;
                for index in start..end {
                    transport_value_gradient[index] =
                        incoming_cotangent * state[transport_column_indices[index] as usize];
                }
            }
        }
        Ok(())
    }
}
