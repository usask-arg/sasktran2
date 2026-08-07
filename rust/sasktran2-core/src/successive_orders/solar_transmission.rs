use super::{
    OperatorError,
    simd::batch::{Batch4, LANES, interleave_wavelength_major},
};

#[derive(Debug, Clone, PartialEq)]
struct OwnedCsrMatrix {
    rows: usize,
    columns: usize,
    row_offsets: Vec<u32>,
    column_indices: Vec<u32>,
    values: Vec<f64>,
}

impl OwnedCsrMatrix {
    fn new(
        rows: usize,
        columns: usize,
        row_offsets: &[u32],
        column_indices: &[u32],
        values: &[f64],
    ) -> Result<Self, OperatorError> {
        if rows == 0
            || columns == 0
            || row_offsets.len() != rows + 1
            || row_offsets[0] != 0
            || row_offsets.windows(2).any(|window| window[0] > window[1])
            || row_offsets.last().copied().unwrap_or_default() as usize != values.len()
            || column_indices.len() != values.len()
        {
            return Err(OperatorError::InvalidRowOffsets);
        }
        if column_indices
            .iter()
            .any(|&column| column as usize >= columns)
        {
            return Err(OperatorError::ColumnOutOfBounds);
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        Ok(Self {
            rows,
            columns,
            row_offsets: row_offsets.to_vec(),
            column_indices: column_indices.to_vec(),
            values: values.to_vec(),
        })
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        if input.len() != self.columns || output.len() != self.rows {
            return Err(OperatorError::DimensionMismatch);
        }
        for (row, result) in output.iter_mut().enumerate() {
            let start = self.row_offsets[row] as usize;
            let end = self.row_offsets[row + 1] as usize;
            *result = self.column_indices[start..end]
                .iter()
                .zip(&self.values[start..end])
                .map(|(&column, &value)| value * input[column as usize])
                .sum();
        }
        Ok(())
    }

    fn apply_batch4_interleaved(
        &self,
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.columns * LANES || output.len() != self.rows * LANES {
            return Err(OperatorError::DimensionMismatch);
        }
        for row in 0..self.rows {
            let start = self.row_offsets[row] as usize;
            let end = self.row_offsets[row + 1] as usize;
            let mut result = Batch4::splat(0.0);
            for (&column, &value) in self.column_indices[start..end]
                .iter()
                .zip(&self.values[start..end])
            {
                result = result + Batch4::splat(value) * Batch4::load(input, column as usize);
            }
            result.store(output, row);
        }
        Ok(())
    }

    fn apply_transpose_add(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        if input.len() != self.rows || output.len() != self.columns {
            return Err(OperatorError::DimensionMismatch);
        }
        for (row, &row_value) in input.iter().enumerate() {
            let start = self.row_offsets[row] as usize;
            let end = self.row_offsets[row + 1] as usize;
            for (&column, &value) in self.column_indices[start..end]
                .iter()
                .zip(&self.values[start..end])
            {
                output[column as usize] += value * row_value;
            }
        }
        Ok(())
    }

    fn apply_transpose_add_batch4_interleaved(
        &self,
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        if input.len() != self.rows * LANES || output.len() != self.columns * LANES {
            return Err(OperatorError::DimensionMismatch);
        }
        for row in 0..self.rows {
            let row_value = Batch4::load(input, row);
            let start = self.row_offsets[row] as usize;
            let end = self.row_offsets[row + 1] as usize;
            for (&column, &value) in self.column_indices[start..end]
                .iter()
                .zip(&self.values[start..end])
            {
                let column = column as usize;
                (Batch4::load(output, column) + Batch4::splat(value) * row_value)
                    .store(output, column);
            }
        }
        Ok(())
    }

    fn storage_bytes(&self) -> usize {
        self.row_offsets.capacity() * size_of::<u32>()
            + self.column_indices.capacity() * size_of::<u32>()
            + self.values.capacity() * size_of::<f64>()
    }
}

#[derive(Debug, Clone, PartialEq)]
enum SolarPathMatrix {
    Dense {
        rows: usize,
        columns: usize,
        values: Vec<f64>,
    },
    Sparse(OwnedCsrMatrix),
}

impl SolarPathMatrix {
    fn rows(&self) -> usize {
        match self {
            Self::Dense { rows, .. } => *rows,
            Self::Sparse(matrix) => matrix.rows,
        }
    }

    fn columns(&self) -> usize {
        match self {
            Self::Dense { columns, .. } => *columns,
            Self::Sparse(matrix) => matrix.columns,
        }
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        match self {
            Self::Dense {
                rows,
                columns,
                values,
            } => {
                if input.len() != *columns || output.len() != *rows {
                    return Err(OperatorError::DimensionMismatch);
                }
                for (row, result) in output.iter_mut().enumerate() {
                    let row_start = row * columns;
                    *result = values[row_start..row_start + columns]
                        .iter()
                        .zip(input)
                        .map(|(&value, &input)| value * input)
                        .sum();
                }
                Ok(())
            }
            Self::Sparse(matrix) => matrix.apply(input, output),
        }
    }

    fn apply_batch4_interleaved(
        &self,
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        match self {
            Self::Dense {
                rows,
                columns,
                values,
            } => {
                if input.len() != columns * LANES || output.len() != rows * LANES {
                    return Err(OperatorError::DimensionMismatch);
                }
                for row in 0..*rows {
                    let row_start = row * columns;
                    let mut result = Batch4::splat(0.0);
                    for (column, &value) in
                        values[row_start..row_start + columns].iter().enumerate()
                    {
                        result = result + Batch4::splat(value) * Batch4::load(input, column);
                    }
                    result.store(output, row);
                }
                Ok(())
            }
            Self::Sparse(matrix) => matrix.apply_batch4_interleaved(input, output),
        }
    }

    fn apply_transpose_add(&self, input: &[f64], output: &mut [f64]) -> Result<(), OperatorError> {
        match self {
            Self::Dense {
                rows,
                columns,
                values,
            } => {
                if input.len() != *rows || output.len() != *columns {
                    return Err(OperatorError::DimensionMismatch);
                }
                for (row, &row_value) in input.iter().enumerate() {
                    let row_start = row * columns;
                    for (result, &value) in output
                        .iter_mut()
                        .zip(&values[row_start..row_start + columns])
                    {
                        *result += value * row_value;
                    }
                }
                Ok(())
            }
            Self::Sparse(matrix) => matrix.apply_transpose_add(input, output),
        }
    }

    fn apply_transpose_add_batch4_interleaved(
        &self,
        input: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        match self {
            Self::Dense {
                rows,
                columns,
                values,
            } => {
                if input.len() != rows * LANES || output.len() != columns * LANES {
                    return Err(OperatorError::DimensionMismatch);
                }
                for row in 0..*rows {
                    let row_value = Batch4::load(input, row);
                    let row_start = row * columns;
                    for (column, &value) in
                        values[row_start..row_start + columns].iter().enumerate()
                    {
                        (Batch4::load(output, column) + Batch4::splat(value) * row_value)
                            .store(output, column);
                    }
                }
                Ok(())
            }
            Self::Sparse(matrix) => matrix.apply_transpose_add_batch4_interleaved(input, output),
        }
    }

    fn storage_bytes(&self) -> usize {
        match self {
            Self::Dense { values, .. } => values.capacity() * size_of::<f64>(),
            Self::Sparse(matrix) => matrix.storage_bytes(),
        }
    }
}

/// Packed solar optical-depth operator for the delegated scalar
/// successive-orders source.
///
/// A table-based 1-D source uses a dense atmosphere-to-table path followed by
/// sparse interpolation to ray endpoints. Exact 2-D sources use one sparse
/// path directly. Primal, JVP, and VJP share the same immutable geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct SolarTransmissionOperator {
    path: SolarPathMatrix,
    interpolation: Option<OwnedCsrMatrix>,
    ground_hit: Vec<u8>,
}

impl SolarTransmissionOperator {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        path_rows: usize,
        path_columns: usize,
        dense_path_values: &[f64],
        path_row_offsets: &[u32],
        path_column_indices: &[u32],
        sparse_path_values: &[f64],
        interpolation_rows: usize,
        interpolation_columns: usize,
        interpolation_row_offsets: &[u32],
        interpolation_column_indices: &[u32],
        interpolation_values: &[f64],
        ground_hit: &[u8],
    ) -> Result<Self, OperatorError> {
        let path = if !dense_path_values.is_empty()
            && path_row_offsets.is_empty()
            && path_column_indices.is_empty()
            && sparse_path_values.is_empty()
        {
            if path_rows == 0
                || path_columns == 0
                || dense_path_values.len()
                    != path_rows
                        .checked_mul(path_columns)
                        .ok_or(OperatorError::DimensionMismatch)?
                || dense_path_values.iter().any(|value| !value.is_finite())
            {
                return Err(OperatorError::DimensionMismatch);
            }
            SolarPathMatrix::Dense {
                rows: path_rows,
                columns: path_columns,
                values: dense_path_values.to_vec(),
            }
        } else if dense_path_values.is_empty() {
            SolarPathMatrix::Sparse(OwnedCsrMatrix::new(
                path_rows,
                path_columns,
                path_row_offsets,
                path_column_indices,
                sparse_path_values,
            )?)
        } else {
            return Err(OperatorError::DimensionMismatch);
        };

        let interpolation = if interpolation_row_offsets.is_empty()
            && interpolation_column_indices.is_empty()
            && interpolation_values.is_empty()
        {
            if interpolation_rows != 0 || interpolation_columns != 0 {
                return Err(OperatorError::DimensionMismatch);
            }
            None
        } else {
            let matrix = OwnedCsrMatrix::new(
                interpolation_rows,
                interpolation_columns,
                interpolation_row_offsets,
                interpolation_column_indices,
                interpolation_values,
            )?;
            if matrix.columns != path.rows() {
                return Err(OperatorError::DimensionMismatch);
            }
            Some(matrix)
        };
        let output_size = interpolation
            .as_ref()
            .map_or_else(|| path.rows(), |matrix| matrix.rows);
        if ground_hit.len() != output_size || ground_hit.iter().any(|&value| value > 1) {
            return Err(OperatorError::DimensionMismatch);
        }
        Ok(Self {
            path,
            interpolation,
            ground_hit: ground_hit.to_vec(),
        })
    }

    #[inline]
    pub fn input_size(&self) -> usize {
        self.path.columns()
    }

    #[inline]
    pub fn output_size(&self) -> usize {
        self.ground_hit.len()
    }

    #[inline]
    pub fn scratch_size(&self) -> usize {
        self.output_size() + self.interpolation.as_ref().map_or(0, |_| self.path.rows())
    }

    #[inline]
    pub fn forward_scratch_size(&self) -> usize {
        self.interpolation.as_ref().map_or(0, |_| self.path.rows())
    }

    pub fn storage_bytes(&self) -> usize {
        self.path.storage_bytes()
            + self
                .interpolation
                .as_ref()
                .map_or(0, OwnedCsrMatrix::storage_bytes)
            + self.ground_hit.capacity() * size_of::<u8>()
    }

    pub fn calculate(
        &self,
        extinction: &[f64],
        solar_irradiance: f64,
        transmission: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_runtime(extinction, transmission, scratch)?;
        self.apply_optical_depth(extinction, transmission, scratch)?;
        for (value, &ground_hit) in transmission.iter_mut().zip(&self.ground_hit) {
            *value = if ground_hit == 0 {
                (-*value).exp() * solar_irradiance
            } else {
                0.0
            };
        }
        Ok(())
    }

    pub fn calculate_batch4(
        &self,
        extinction: &[f64],
        solar_irradiance: &[f64],
        transmission: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if extinction.len() != self.input_size() * LANES
            || solar_irradiance.len() != LANES
            || transmission.len() != self.output_size() * LANES
            || scratch.len() < self.forward_scratch_size() * LANES
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let mut interleaved_extinction = vec![0.0; extinction.len()];
        interleave_wavelength_major(extinction, self.input_size(), &mut interleaved_extinction);
        if let Some(interpolation) = &self.interpolation {
            let path_size = self.path.rows() * LANES;
            let path_optical_depth = &mut scratch[..path_size];
            self.path
                .apply_batch4_interleaved(&interleaved_extinction, path_optical_depth)?;
            interpolation.apply_batch4_interleaved(path_optical_depth, transmission)?;
        } else {
            self.path
                .apply_batch4_interleaved(&interleaved_extinction, transmission)?;
        }
        let irradiance = Batch4::from_array(solar_irradiance.try_into().unwrap());
        for (element, &ground_hit) in self.ground_hit.iter().enumerate() {
            let value = if ground_hit == 0 {
                (-Batch4::load(transmission, element)).exp() * irradiance
            } else {
                Batch4::splat(0.0)
            };
            value.store(transmission, element);
        }
        Ok(())
    }

    pub fn calculate_jvp(
        &self,
        extinction_tangent: &[f64],
        transmission: &[f64],
        transmission_tangent: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_runtime(extinction_tangent, transmission_tangent, scratch)?;
        if transmission.len() != self.output_size() {
            return Err(OperatorError::DimensionMismatch);
        }
        self.apply_optical_depth(extinction_tangent, transmission_tangent, scratch)?;
        for ((tangent, &value), &ground_hit) in transmission_tangent
            .iter_mut()
            .zip(transmission)
            .zip(&self.ground_hit)
        {
            *tangent = if ground_hit == 0 {
                -value * *tangent
            } else {
                0.0
            };
        }
        Ok(())
    }

    /// Calculates four wavelength JVPs together. All arrays use an
    /// element-major layout with four contiguous wavelength lanes.
    pub fn calculate_jvp_batch4(
        &self,
        extinction_tangent: &[f64],
        transmission: &[f64],
        transmission_tangent: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if extinction_tangent.len() != self.input_size() * LANES
            || transmission.len() != self.output_size() * LANES
            || transmission_tangent.len() != transmission.len()
            || scratch.len() < self.forward_scratch_size() * LANES
        {
            return Err(OperatorError::DimensionMismatch);
        }
        if let Some(interpolation) = &self.interpolation {
            let path_size = self.path.rows() * LANES;
            let path_optical_depth_tangent = &mut scratch[..path_size];
            self.path
                .apply_batch4_interleaved(extinction_tangent, path_optical_depth_tangent)?;
            interpolation
                .apply_batch4_interleaved(path_optical_depth_tangent, transmission_tangent)?;
        } else {
            self.path
                .apply_batch4_interleaved(extinction_tangent, transmission_tangent)?;
        }
        for (element, &ground_hit) in self.ground_hit.iter().enumerate() {
            let tangent = if ground_hit == 0 {
                -Batch4::load(transmission, element) * Batch4::load(transmission_tangent, element)
            } else {
                Batch4::splat(0.0)
            };
            tangent.store(transmission_tangent, element);
        }
        Ok(())
    }

    pub fn accumulate_vjp(
        &self,
        transmission: &[f64],
        transmission_cotangent: &[f64],
        extinction_gradient: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if transmission.len() != self.output_size()
            || transmission_cotangent.len() != self.output_size()
            || extinction_gradient.len() != self.input_size()
            || scratch.len() < self.scratch_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let (output_cotangent, remaining_scratch) = scratch.split_at_mut(self.output_size());
        for (((result, &transmission), &cotangent), &ground_hit) in output_cotangent
            .iter_mut()
            .zip(transmission)
            .zip(transmission_cotangent)
            .zip(&self.ground_hit)
        {
            *result = if ground_hit == 0 {
                -transmission * cotangent
            } else {
                0.0
            };
        }
        if let Some(interpolation) = &self.interpolation {
            let path_cotangent = &mut remaining_scratch[..self.path.rows()];
            path_cotangent.fill(0.0);
            interpolation.apply_transpose_add(output_cotangent, path_cotangent)?;
            self.path
                .apply_transpose_add(path_cotangent, extinction_gradient)?;
        } else {
            self.path
                .apply_transpose_add(output_cotangent, extinction_gradient)?;
        }
        Ok(())
    }

    /// Accumulates four wavelength VJPs together. The transmission
    /// cotangent is consumed in place so the caller can reuse its large JVP
    /// product buffer instead of retaining another four-lane allocation.
    pub fn accumulate_vjp_batch4(
        &self,
        transmission: &[f64],
        transmission_cotangent: &mut [f64],
        extinction_gradient: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if transmission.len() != self.output_size() * LANES
            || transmission_cotangent.len() != transmission.len()
            || extinction_gradient.len() != self.input_size() * LANES
            || scratch.len() < self.forward_scratch_size() * LANES
        {
            return Err(OperatorError::DimensionMismatch);
        }
        for (element, &ground_hit) in self.ground_hit.iter().enumerate() {
            let cotangent = if ground_hit == 0 {
                -Batch4::load(transmission, element) * Batch4::load(transmission_cotangent, element)
            } else {
                Batch4::splat(0.0)
            };
            cotangent.store(transmission_cotangent, element);
        }
        extinction_gradient.fill(0.0);
        if let Some(interpolation) = &self.interpolation {
            let path_cotangent = &mut scratch[..self.path.rows() * LANES];
            path_cotangent.fill(0.0);
            interpolation
                .apply_transpose_add_batch4_interleaved(transmission_cotangent, path_cotangent)?;
            self.path
                .apply_transpose_add_batch4_interleaved(path_cotangent, extinction_gradient)?;
        } else {
            self.path.apply_transpose_add_batch4_interleaved(
                transmission_cotangent,
                extinction_gradient,
            )?;
        }
        Ok(())
    }

    fn validate_runtime(
        &self,
        extinction: &[f64],
        output: &[f64],
        scratch: &[f64],
    ) -> Result<(), OperatorError> {
        if extinction.len() != self.input_size()
            || output.len() != self.output_size()
            || scratch.len() < self.forward_scratch_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        Ok(())
    }

    fn apply_optical_depth(
        &self,
        extinction: &[f64],
        optical_depth: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), OperatorError> {
        if let Some(interpolation) = &self.interpolation {
            let path_optical_depth = &mut scratch[..self.path.rows()];
            self.path.apply(extinction, path_optical_depth)?;
            interpolation.apply(path_optical_depth, optical_depth)
        } else {
            self.path.apply(extinction, optical_depth)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{LANES, SolarTransmissionOperator};

    fn chained_operator() -> SolarTransmissionOperator {
        SolarTransmissionOperator::new(
            2,
            2,
            &[1.0, 2.0, 3.0, 4.0],
            &[],
            &[],
            &[],
            3,
            2,
            &[0, 1, 3, 4],
            &[0, 0, 1, 1],
            &[1.0, 0.25, 0.75, 1.0],
            &[0, 0, 1],
        )
        .unwrap()
    }

    #[test]
    fn chained_primal_matches_dense_reference() {
        let operator = chained_operator();
        let mut transmission = [0.0; 3];
        let mut scratch = vec![0.0; operator.scratch_size()];
        operator
            .calculate(&[0.1, 0.2], 2.0, &mut transmission, &mut scratch)
            .unwrap();
        let table_od = [0.5_f64, 1.1_f64];
        let expected = [
            (-table_od[0]).exp() * 2.0,
            (-(0.25 * table_od[0] + 0.75 * table_od[1])).exp() * 2.0,
            0.0,
        ];
        for (&actual, expected) in transmission.iter().zip(expected) {
            assert!((actual - expected).abs() < 1.0e-14);
        }
    }

    #[test]
    fn batch4_primal_matches_scalar_wavelengths() {
        let operator = chained_operator();
        let extinction_lanes = [[0.1, 0.2], [0.2, 0.1], [0.3, 0.4], [0.05, 0.15]];
        let irradiance = [2.0, 1.5, 3.0, 0.75];
        let extinction: Vec<_> = extinction_lanes.into_iter().flatten().collect();
        let mut transmission = vec![0.0; operator.output_size() * LANES];
        let mut scratch = vec![0.0; operator.forward_scratch_size() * LANES];
        operator
            .calculate_batch4(&extinction, &irradiance, &mut transmission, &mut scratch)
            .unwrap();

        for lane in 0..LANES {
            let mut expected = vec![0.0; operator.output_size()];
            let mut scalar_scratch = vec![0.0; operator.forward_scratch_size()];
            operator
                .calculate(
                    &extinction_lanes[lane],
                    irradiance[lane],
                    &mut expected,
                    &mut scalar_scratch,
                )
                .unwrap();
            for (element, expected) in expected.into_iter().enumerate() {
                assert!(
                    (transmission[element * LANES + lane] - expected).abs()
                        <= 5.0e-14 * expected.abs().max(1.0)
                );
            }
        }
    }

    #[test]
    fn jvp_and_vjp_are_dual() {
        let operator = chained_operator();
        let mut transmission = [0.0; 3];
        let mut scratch = vec![0.0; operator.scratch_size()];
        operator
            .calculate(&[0.1, 0.2], 2.0, &mut transmission, &mut scratch)
            .unwrap();
        let extinction_tangent = [0.4, -0.3];
        let mut transmission_tangent = [0.0; 3];
        operator
            .calculate_jvp(
                &extinction_tangent,
                &transmission,
                &mut transmission_tangent,
                &mut scratch,
            )
            .unwrap();
        let transmission_cotangent = [0.7, -0.2, 3.0];
        let mut extinction_gradient = [0.0; 2];
        operator
            .accumulate_vjp(
                &transmission,
                &transmission_cotangent,
                &mut extinction_gradient,
                &mut scratch,
            )
            .unwrap();
        let forward_dot = transmission_tangent
            .iter()
            .zip(transmission_cotangent)
            .map(|(&left, right)| left * right)
            .sum::<f64>();
        let reverse_dot = extinction_tangent
            .iter()
            .zip(extinction_gradient)
            .map(|(&left, right)| left * right)
            .sum::<f64>();
        assert!((forward_dot - reverse_dot).abs() < 1.0e-12);
    }

    #[test]
    fn batch4_products_match_scalar_wavelengths() {
        let operator = chained_operator();
        let extinction_lanes = [[0.1, 0.2], [0.2, 0.1], [0.3, 0.4], [0.05, 0.15]];
        let tangent_lanes = [[0.4, -0.3], [-0.2, 0.1], [0.3, 0.2], [-0.1, -0.4]];
        let irradiance = [2.0, 1.5, 3.0, 0.75];
        let extinction: Vec<_> = extinction_lanes.into_iter().flatten().collect();
        let mut interleaved_tangent = vec![0.0; operator.input_size() * LANES];
        for element in 0..operator.input_size() {
            for lane in 0..LANES {
                interleaved_tangent[element * LANES + lane] = tangent_lanes[lane][element];
            }
        }
        let mut transmission = vec![0.0; operator.output_size() * LANES];
        let mut scratch = vec![0.0; operator.forward_scratch_size() * LANES];
        operator
            .calculate_batch4(&extinction, &irradiance, &mut transmission, &mut scratch)
            .unwrap();
        let mut transmission_tangent = vec![0.0; transmission.len()];
        operator
            .calculate_jvp_batch4(
                &interleaved_tangent,
                &transmission,
                &mut transmission_tangent,
                &mut scratch,
            )
            .unwrap();

        let cotangent_lanes = [
            [0.7, -0.2, 3.0],
            [-0.1, 0.5, 0.8],
            [0.3, 0.2, -0.7],
            [0.9, -0.4, 0.1],
        ];
        let mut transmission_cotangent = vec![0.0; transmission.len()];
        for element in 0..operator.output_size() {
            for lane in 0..LANES {
                transmission_cotangent[element * LANES + lane] = cotangent_lanes[lane][element];
            }
        }
        let mut extinction_gradient = vec![0.0; operator.input_size() * LANES];
        operator
            .accumulate_vjp_batch4(
                &transmission,
                &mut transmission_cotangent,
                &mut extinction_gradient,
                &mut scratch,
            )
            .unwrap();

        for lane in 0..LANES {
            let mut scalar_transmission = vec![0.0; operator.output_size()];
            let mut scalar_scratch = vec![0.0; operator.scratch_size()];
            operator
                .calculate(
                    &extinction_lanes[lane],
                    irradiance[lane],
                    &mut scalar_transmission,
                    &mut scalar_scratch,
                )
                .unwrap();
            let mut scalar_tangent = vec![0.0; operator.output_size()];
            operator
                .calculate_jvp(
                    &tangent_lanes[lane],
                    &scalar_transmission,
                    &mut scalar_tangent,
                    &mut scalar_scratch,
                )
                .unwrap();
            let mut scalar_gradient = vec![0.0; operator.input_size()];
            operator
                .accumulate_vjp(
                    &scalar_transmission,
                    &cotangent_lanes[lane],
                    &mut scalar_gradient,
                    &mut scalar_scratch,
                )
                .unwrap();
            for (element, &expected) in scalar_tangent.iter().enumerate() {
                assert!(
                    (transmission_tangent[element * LANES + lane] - expected).abs()
                        <= 5.0e-14 * expected.abs().max(1.0)
                );
            }
            for (element, &expected) in scalar_gradient.iter().enumerate() {
                assert!(
                    (extinction_gradient[element * LANES + lane] - expected).abs()
                        <= 5.0e-14 * expected.abs().max(1.0)
                );
            }
        }
    }

    #[test]
    fn sparse_direct_operator_is_supported() {
        let operator = SolarTransmissionOperator::new(
            2,
            2,
            &[],
            &[0, 1, 3],
            &[0, 0, 1],
            &[2.0, 1.0, 3.0],
            0,
            0,
            &[],
            &[],
            &[],
            &[0, 0],
        )
        .unwrap();
        let mut transmission = [0.0; 2];
        let mut scratch = vec![0.0; operator.scratch_size()];
        operator
            .calculate(&[0.5, 0.25], 1.0, &mut transmission, &mut scratch)
            .unwrap();
        assert!((transmission[0] - (-1.0_f64).exp()).abs() < 1.0e-14);
        assert!((transmission[1] - (-1.25_f64).exp()).abs() < 1.0e-14);
    }
}
