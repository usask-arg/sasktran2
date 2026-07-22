use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperatorError {
    InvalidRowOffsets,
    InvalidBlockOffsets,
    ColumnOutOfBounds,
    DimensionMismatch,
    NonFiniteValue,
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
        }
    }
}

impl std::error::Error for OperatorError {}

/// Sparse transport operator with immutable structure and replaceable values.
#[derive(Debug, Clone, PartialEq)]
pub struct CsrMatrix {
    rows: usize,
    columns: usize,
    row_offsets: Vec<usize>,
    column_indices: Vec<usize>,
    values: Vec<f64>,
}

impl CsrMatrix {
    pub fn new(
        rows: usize,
        columns: usize,
        row_offsets: Vec<usize>,
        column_indices: Vec<usize>,
        values: Vec<f64>,
    ) -> Result<Self, OperatorError> {
        if row_offsets.len() != rows + 1
            || row_offsets.first().copied() != Some(0)
            || row_offsets.last().copied() != Some(values.len())
            || row_offsets.windows(2).any(|window| window[0] > window[1])
            || column_indices.len() != values.len()
        {
            return Err(OperatorError::InvalidRowOffsets);
        }
        if column_indices.iter().any(|&column| column >= columns) {
            return Err(OperatorError::ColumnOutOfBounds);
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        Ok(Self {
            rows,
            columns,
            row_offsets,
            column_indices,
            values,
        })
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
        if input.len() != self.columns || output.len() != self.rows {
            return Err(OperatorError::DimensionMismatch);
        }
        for (row, result) in output.iter_mut().enumerate() {
            let start = self.row_offsets[row];
            let end = self.row_offsets[row + 1];
            *result = self.column_indices[start..end]
                .iter()
                .zip(&self.values[start..end])
                .map(|(&column, &value)| value * input[column])
                .sum();
        }
        Ok(())
    }
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
    scattering: BlockDiagonalMatrix,
    forcing: Vec<f64>,
}

impl FixedPointProblem {
    pub fn new(
        transport: CsrMatrix,
        scattering: BlockDiagonalMatrix,
        forcing: Vec<f64>,
    ) -> Result<Self, OperatorError> {
        if transport.rows() != scattering.input_size()
            || transport.columns() != scattering.output_size()
            || forcing.len() != transport.rows()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        if forcing.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        Ok(Self {
            transport,
            scattering,
            forcing,
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

    pub fn set_transport_values(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        self.transport.set_values(values)
    }

    pub fn set_scattering_values(&mut self, values: &[f64]) -> Result<(), OperatorError> {
        self.scattering.set_values(values)
    }

    pub fn set_forcing(&mut self, forcing: &[f64]) -> Result<(), OperatorError> {
        if forcing.len() != self.forcing.len() {
            return Err(OperatorError::DimensionMismatch);
        }
        if forcing.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        self.forcing.copy_from_slice(forcing);
        Ok(())
    }

    pub fn apply(
        &self,
        state: &[f64],
        output: &mut [f64],
        incoming_scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.transport.apply(state, incoming_scratch)?;
        for (incoming, forcing) in incoming_scratch.iter_mut().zip(&self.forcing) {
            *incoming += forcing;
        }
        self.scattering.apply(incoming_scratch, output)
    }
}
