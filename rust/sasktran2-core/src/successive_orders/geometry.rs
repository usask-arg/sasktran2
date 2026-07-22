use super::OperatorError;

/// A source-grid point expressed entirely in global Cartesian coordinates.
///
/// `frame_x`, `frame_y`, and `frame_z` define the local angular basis. Keeping
/// the frame explicit avoids embedding 1-D solar-relative assumptions in the
/// transport and scattering operators.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SourceNode {
    pub position: [f64; 3],
    pub frame_x: [f64; 3],
    pub frame_y: [f64; 3],
    pub frame_z: [f64; 3],
    pub boundary: bool,
}

/// Variable-width interpolation stencils in compressed-row form.
///
/// A 1-D linear atmosphere normally contributes two entries per row, while a
/// structured 2-D cell contributes four. The representation does not assign
/// meaning to those widths.
#[derive(Debug, Clone, PartialEq)]
pub struct CompressedStencils {
    num_columns: usize,
    row_offsets: Vec<usize>,
    column_indices: Vec<usize>,
    weights: Vec<f64>,
}

impl CompressedStencils {
    pub fn new(
        num_columns: usize,
        row_offsets: Vec<usize>,
        column_indices: Vec<usize>,
        weights: Vec<f64>,
    ) -> Result<Self, OperatorError> {
        if row_offsets.is_empty() || row_offsets[0] != 0 {
            return Err(OperatorError::InvalidRowOffsets);
        }
        if column_indices.len() != weights.len()
            || row_offsets.last().copied() != Some(weights.len())
            || row_offsets.windows(2).any(|window| window[0] > window[1])
        {
            return Err(OperatorError::InvalidRowOffsets);
        }
        if column_indices.iter().any(|&index| index >= num_columns) {
            return Err(OperatorError::ColumnOutOfBounds);
        }
        if weights.iter().any(|weight| !weight.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        Ok(Self {
            num_columns,
            row_offsets,
            column_indices,
            weights,
        })
    }

    #[inline]
    pub fn num_rows(&self) -> usize {
        self.row_offsets.len() - 1
    }

    #[inline]
    pub fn num_columns(&self) -> usize {
        self.num_columns
    }

    pub fn interpolate(&self, values: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        if values.len() != self.num_columns || output.len() != self.num_rows() {
            return Err(OperatorError::DimensionMismatch);
        }
        for (row, result) in output.iter_mut().enumerate() {
            let start = self.row_offsets[row];
            let end = self.row_offsets[row + 1];
            *result = self.column_indices[start..end]
                .iter()
                .zip(&self.weights[start..end])
                .map(|(&column, &weight)| weight * values[column])
                .sum();
        }
        Ok(())
    }
}
