use std::fmt::{Display, Formatter};

use super::ScalarCoefficientScattering;

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
/// retaining a second, widened copy of the transport operator.
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
}

impl ScatteringOperator {
    #[inline]
    pub fn output_size(&self) -> usize {
        match self {
            Self::Dense(matrix) => matrix.output_size(),
            Self::ScalarCoefficients(matrix) => matrix.output_size(),
        }
    }

    #[inline]
    pub fn input_size(&self) -> usize {
        match self {
            Self::Dense(matrix) => matrix.input_size(),
            Self::ScalarCoefficients(matrix) => matrix.input_size(),
        }
    }

    pub fn set_dense_values(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        match self {
            Self::Dense(matrix) => matrix.set_values(values),
            Self::ScalarCoefficients(matrix) => matrix.set_dense_values(values),
        }
    }

    pub fn set_coefficients(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        match self {
            Self::Dense(_) => Err(OperatorError::UnsupportedOperator),
            Self::ScalarCoefficients(matrix) => matrix.set_coefficients(values),
        }
    }

    pub fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        match self {
            Self::Dense(matrix) => matrix.apply(input, output),
            Self::ScalarCoefficients(matrix) => matrix.apply(input, output),
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
}
