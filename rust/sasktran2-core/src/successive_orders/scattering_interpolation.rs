use super::OperatorError;

/// Contracts the atmospheric phase-coefficient derivative tensor with a
/// derivative-group tangent at one wavelength.
///
/// The input tensor follows Eigen's `(coefficient, atmosphere, wavelength,
/// group)` storage order.  The output is atmosphere-major and may retain only
/// a prefix of each coefficient row, which is useful when single scattering
/// consumes fewer moments than the atmosphere stores.
#[allow(clippy::too_many_arguments)]
pub fn atmospheric_coefficient_jvp(
    atmosphere_coefficient_derivatives: &[f64],
    atmosphere_coefficient_stride: usize,
    num_atmosphere_points: usize,
    num_wavelengths: usize,
    wavelength_index: usize,
    derivative_group_tangent: &[f64],
    output: &mut [f64],
) -> Result<(), OperatorError> {
    if atmosphere_coefficient_stride == 0
        || num_atmosphere_points == 0
        || num_wavelengths == 0
        || wavelength_index >= num_wavelengths
        || !derivative_group_tangent
            .len()
            .is_multiple_of(num_atmosphere_points)
        || !output.len().is_multiple_of(num_atmosphere_points)
    {
        return Err(OperatorError::DimensionMismatch);
    }
    let output_stride = output.len() / num_atmosphere_points;
    let num_derivative_groups = derivative_group_tangent.len() / num_atmosphere_points;
    if output_stride == 0 || output_stride > atmosphere_coefficient_stride {
        return Err(OperatorError::DimensionMismatch);
    }
    let required_derivatives = atmosphere_coefficient_stride
        .checked_mul(num_atmosphere_points)
        .and_then(|size| size.checked_mul(num_wavelengths))
        .and_then(|size| size.checked_mul(num_derivative_groups))
        .ok_or(OperatorError::DimensionMismatch)?;
    if atmosphere_coefficient_derivatives.len() != required_derivatives {
        return Err(OperatorError::DimensionMismatch);
    }

    output.fill(0.0);
    for group in 0..num_derivative_groups {
        let tangent_start = group * num_atmosphere_points;
        let tangent =
            &derivative_group_tangent[tangent_start..tangent_start + num_atmosphere_points];
        if tangent.iter().all(|&value| value == 0.0) {
            continue;
        }
        let wavelength_start = (group * num_wavelengths + wavelength_index)
            * num_atmosphere_points
            * atmosphere_coefficient_stride;
        for (atmosphere_index, &scale) in tangent.iter().enumerate() {
            if scale == 0.0 {
                continue;
            }
            let derivative_start =
                wavelength_start + atmosphere_index * atmosphere_coefficient_stride;
            let output_start = atmosphere_index * output_stride;
            for (result, &derivative) in output[output_start..output_start + output_stride]
                .iter_mut()
                .zip(
                    &atmosphere_coefficient_derivatives
                        [derivative_start..derivative_start + output_stride],
                )
            {
                *result += scale * derivative;
            }
        }
    }
    Ok(())
}

/// Pulls an atmosphere-major phase-coefficient cotangent back through the
/// Eigen atmospheric derivative tensor at one wavelength.
#[allow(clippy::too_many_arguments)]
pub fn accumulate_atmospheric_coefficient_vjp(
    coefficient_gradient: &[f64],
    atmosphere_coefficient_derivatives: &[f64],
    atmosphere_coefficient_stride: usize,
    num_atmosphere_points: usize,
    num_wavelengths: usize,
    wavelength_index: usize,
    active_derivative_groups: &[i32],
    atmosphere_gradient: &mut [f64],
) -> Result<(), OperatorError> {
    if atmosphere_coefficient_stride == 0
        || num_atmosphere_points == 0
        || num_wavelengths == 0
        || wavelength_index >= num_wavelengths
        || !coefficient_gradient
            .len()
            .is_multiple_of(num_atmosphere_points)
        || !atmosphere_gradient
            .len()
            .is_multiple_of(num_atmosphere_points)
    {
        return Err(OperatorError::DimensionMismatch);
    }
    let coefficient_count = coefficient_gradient.len() / num_atmosphere_points;
    let num_derivative_groups = atmosphere_gradient.len() / num_atmosphere_points;
    if coefficient_count == 0
        || coefficient_count > atmosphere_coefficient_stride
        || active_derivative_groups
            .iter()
            .any(|&group| group < 0 || group as usize >= num_derivative_groups)
        || active_derivative_groups
            .windows(2)
            .any(|window| window[0] >= window[1])
    {
        return Err(OperatorError::DimensionMismatch);
    }
    let required_derivatives = atmosphere_coefficient_stride
        .checked_mul(num_atmosphere_points)
        .and_then(|size| size.checked_mul(num_wavelengths))
        .and_then(|size| size.checked_mul(num_derivative_groups))
        .ok_or(OperatorError::DimensionMismatch)?;
    if atmosphere_coefficient_derivatives.len() != required_derivatives {
        return Err(OperatorError::DimensionMismatch);
    }

    for &group in active_derivative_groups {
        let group = group as usize;
        let wavelength_start = (group * num_wavelengths + wavelength_index)
            * num_atmosphere_points
            * atmosphere_coefficient_stride;
        let gradient_start = group * num_atmosphere_points;
        for atmosphere_index in 0..num_atmosphere_points {
            let coefficient_start = atmosphere_index * coefficient_count;
            let derivative_start =
                wavelength_start + atmosphere_index * atmosphere_coefficient_stride;
            atmosphere_gradient[gradient_start + atmosphere_index] += coefficient_gradient
                [coefficient_start..coefficient_start + coefficient_count]
                .iter()
                .zip(
                    &atmosphere_coefficient_derivatives
                        [derivative_start..derivative_start + coefficient_count],
                )
                .map(|(&gradient, &derivative)| gradient * derivative)
                .sum::<f64>();
        }
    }
    Ok(())
}

/// Geometry-neutral interpolation from atmospheric grid points to the
/// successive-orders scattering nodes.
///
/// Atmospheric coefficients use the Eigen tensor layout supplied by the C++
/// atmosphere: coefficient is the contiguous dimension, followed by geometry,
/// wavelength, and derivative group. Output coefficients are point-major.
#[derive(Debug, Clone, PartialEq)]
pub struct ScatteringCoefficientInterpolator {
    point_offsets: Vec<u32>,
    atmosphere_indices: Vec<u32>,
    weights: Vec<f64>,
    num_atmosphere_points: usize,
    num_output_coefficients: usize,
}

impl ScatteringCoefficientInterpolator {
    pub fn new(
        point_offsets: &[u32],
        atmosphere_indices: &[u32],
        weights: &[f64],
        num_atmosphere_points: usize,
        num_output_coefficients: usize,
    ) -> Result<Self, OperatorError> {
        if point_offsets.is_empty()
            || point_offsets[0] != 0
            || point_offsets.windows(2).any(|window| window[0] > window[1])
            || point_offsets.last().copied().unwrap_or_default() as usize
                != atmosphere_indices.len()
            || atmosphere_indices.len() != weights.len()
        {
            return Err(OperatorError::InvalidRowOffsets);
        }
        if num_atmosphere_points == 0 || num_output_coefficients == 0 {
            return Err(OperatorError::DimensionMismatch);
        }
        if atmosphere_indices
            .iter()
            .any(|&index| index as usize >= num_atmosphere_points)
        {
            return Err(OperatorError::ColumnOutOfBounds);
        }
        if weights.iter().any(|weight| !weight.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }
        Ok(Self {
            point_offsets: point_offsets.to_vec(),
            atmosphere_indices: atmosphere_indices.to_vec(),
            weights: weights.to_vec(),
            num_atmosphere_points,
            num_output_coefficients,
        })
    }

    #[inline]
    pub fn num_points(&self) -> usize {
        self.point_offsets.len() - 1
    }

    #[inline]
    pub fn num_weights(&self) -> usize {
        self.weights.len()
    }

    pub fn storage_bytes(&self) -> usize {
        self.point_offsets.capacity() * size_of::<u32>()
            + self.atmosphere_indices.capacity() * size_of::<u32>()
            + self.weights.capacity() * size_of::<f64>()
    }

    pub fn interpolate(
        &self,
        atmosphere_coefficients: &[f64],
        atmosphere_coefficient_stride: usize,
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_coefficient_stride(atmosphere_coefficient_stride)?;
        let required_input = atmosphere_coefficient_stride
            .checked_mul(self.num_atmosphere_points)
            .ok_or(OperatorError::DimensionMismatch)?;
        if atmosphere_coefficients.len() != required_input || output.len() != self.output_size() {
            return Err(OperatorError::DimensionMismatch);
        }

        output.fill(0.0);
        for point in 0..self.num_points() {
            let stencil_start = self.point_offsets[point] as usize;
            let stencil_end = self.point_offsets[point + 1] as usize;
            let output_start = point * self.num_output_coefficients;
            let point_output =
                &mut output[output_start..output_start + self.num_output_coefficients];
            for (&atmosphere_index, &weight) in self.atmosphere_indices[stencil_start..stencil_end]
                .iter()
                .zip(&self.weights[stencil_start..stencil_end])
            {
                let input_start = atmosphere_index as usize * atmosphere_coefficient_stride;
                let point_input = &atmosphere_coefficients
                    [input_start..input_start + self.num_output_coefficients];
                for (result, &value) in point_output.iter_mut().zip(point_input) {
                    *result += weight * value;
                }
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn interpolate_jvp(
        &self,
        atmosphere_coefficient_derivatives: &[f64],
        atmosphere_coefficient_stride: usize,
        num_wavelengths: usize,
        wavelength_index: usize,
        derivative_group_tangent: &[f64],
        output: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_coefficient_stride(atmosphere_coefficient_stride)?;
        if num_wavelengths == 0
            || wavelength_index >= num_wavelengths
            || !derivative_group_tangent
                .len()
                .is_multiple_of(self.num_atmosphere_points)
            || output.len() != self.output_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let num_derivative_groups = derivative_group_tangent.len() / self.num_atmosphere_points;
        let required_derivatives = self.derivative_tensor_size(
            atmosphere_coefficient_stride,
            num_wavelengths,
            num_derivative_groups,
        )?;
        if atmosphere_coefficient_derivatives.len() != required_derivatives {
            return Err(OperatorError::DimensionMismatch);
        }

        output.fill(0.0);
        for group in 0..num_derivative_groups {
            let tangent_start = group * self.num_atmosphere_points;
            let tangent = &derivative_group_tangent
                [tangent_start..tangent_start + self.num_atmosphere_points];
            if tangent.iter().all(|&value| value == 0.0) {
                continue;
            }
            let wavelength_start = (group * num_wavelengths + wavelength_index)
                * self.num_atmosphere_points
                * atmosphere_coefficient_stride;
            for point in 0..self.num_points() {
                let stencil_start = self.point_offsets[point] as usize;
                let stencil_end = self.point_offsets[point + 1] as usize;
                let output_start = point * self.num_output_coefficients;
                let point_output =
                    &mut output[output_start..output_start + self.num_output_coefficients];
                for (&atmosphere_index, &weight) in self.atmosphere_indices
                    [stencil_start..stencil_end]
                    .iter()
                    .zip(&self.weights[stencil_start..stencil_end])
                {
                    let atmosphere_index = atmosphere_index as usize;
                    let scale = weight * tangent[atmosphere_index];
                    if scale == 0.0 {
                        continue;
                    }
                    let derivative_start =
                        wavelength_start + atmosphere_index * atmosphere_coefficient_stride;
                    let point_derivative = &atmosphere_coefficient_derivatives
                        [derivative_start..derivative_start + self.num_output_coefficients];
                    for (result, &derivative) in point_output.iter_mut().zip(point_derivative) {
                        *result += scale * derivative;
                    }
                }
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_vjp(
        &self,
        coefficient_gradient: &[f64],
        atmosphere_coefficient_derivatives: &[f64],
        atmosphere_coefficient_stride: usize,
        num_wavelengths: usize,
        wavelength_index: usize,
        active_derivative_groups: &[i32],
        atmosphere_gradient: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_coefficient_stride(atmosphere_coefficient_stride)?;
        if num_wavelengths == 0
            || wavelength_index >= num_wavelengths
            || !atmosphere_gradient
                .len()
                .is_multiple_of(self.num_atmosphere_points)
            || coefficient_gradient.len() != self.output_size()
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let num_derivative_groups = atmosphere_gradient.len() / self.num_atmosphere_points;
        let required_derivatives = self.derivative_tensor_size(
            atmosphere_coefficient_stride,
            num_wavelengths,
            num_derivative_groups,
        )?;
        if atmosphere_coefficient_derivatives.len() != required_derivatives
            || active_derivative_groups
                .iter()
                .any(|&group| group < 0 || group as usize >= num_derivative_groups)
            || active_derivative_groups
                .windows(2)
                .any(|window| window[0] >= window[1])
        {
            return Err(OperatorError::DimensionMismatch);
        }

        for &group in active_derivative_groups {
            let group = group as usize;
            let wavelength_start = (group * num_wavelengths + wavelength_index)
                * self.num_atmosphere_points
                * atmosphere_coefficient_stride;
            let gradient_start = group * self.num_atmosphere_points;
            for point in 0..self.num_points() {
                let coefficient_start = point * self.num_output_coefficients;
                let point_gradient = &coefficient_gradient
                    [coefficient_start..coefficient_start + self.num_output_coefficients];
                let stencil_start = self.point_offsets[point] as usize;
                let stencil_end = self.point_offsets[point + 1] as usize;
                for (&atmosphere_index, &weight) in self.atmosphere_indices
                    [stencil_start..stencil_end]
                    .iter()
                    .zip(&self.weights[stencil_start..stencil_end])
                {
                    let atmosphere_index = atmosphere_index as usize;
                    let derivative_start =
                        wavelength_start + atmosphere_index * atmosphere_coefficient_stride;
                    let point_derivative = &atmosphere_coefficient_derivatives
                        [derivative_start..derivative_start + self.num_output_coefficients];
                    atmosphere_gradient[gradient_start + atmosphere_index] += weight
                        * point_gradient
                            .iter()
                            .zip(point_derivative)
                            .map(|(&gradient, &derivative)| gradient * derivative)
                            .sum::<f64>();
                }
            }
        }
        Ok(())
    }

    #[inline]
    fn output_size(&self) -> usize {
        self.num_points() * self.num_output_coefficients
    }

    fn validate_coefficient_stride(
        &self,
        atmosphere_coefficient_stride: usize,
    ) -> Result<(), OperatorError> {
        if atmosphere_coefficient_stride < self.num_output_coefficients {
            return Err(OperatorError::DimensionMismatch);
        }
        Ok(())
    }

    fn derivative_tensor_size(
        &self,
        atmosphere_coefficient_stride: usize,
        num_wavelengths: usize,
        num_derivative_groups: usize,
    ) -> Result<usize, OperatorError> {
        atmosphere_coefficient_stride
            .checked_mul(self.num_atmosphere_points)
            .and_then(|size| size.checked_mul(num_wavelengths))
            .and_then(|size| size.checked_mul(num_derivative_groups))
            .ok_or(OperatorError::DimensionMismatch)
    }
}

#[cfg(test)]
mod tests {
    use super::{
        ScatteringCoefficientInterpolator, accumulate_atmospheric_coefficient_vjp,
        atmospheric_coefficient_jvp,
    };

    fn interpolator() -> ScatteringCoefficientInterpolator {
        ScatteringCoefficientInterpolator::new(
            &[0, 2, 4],
            &[0, 1, 1, 2],
            &[0.25, 0.75, 0.5, 0.5],
            3,
            2,
        )
        .unwrap()
    }

    #[test]
    fn interpolates_strided_atmosphere_coefficients() {
        let interpolator = interpolator();
        let atmosphere = [1.0, 2.0, 99.0, 3.0, 4.0, 99.0, 5.0, 6.0, 99.0];
        let mut output = [0.0; 4];
        interpolator
            .interpolate(&atmosphere, 3, &mut output)
            .unwrap();
        assert_eq!(output, [2.5, 3.5, 4.0, 5.0]);
    }

    #[test]
    fn jvp_and_vjp_are_dual_for_tensor_layout() {
        let interpolator = interpolator();
        let num_wavelengths = 2;
        let num_groups = 2;
        let stride = 3;
        let mut derivatives = vec![0.0; stride * 3 * num_wavelengths * num_groups];
        for group in 0..num_groups {
            for wavelength in 0..num_wavelengths {
                for geometry in 0..3 {
                    for coefficient in 0..stride {
                        let index = coefficient
                            + stride * (geometry + 3 * (wavelength + num_wavelengths * group));
                        derivatives[index] = 0.1
                            * (1 + coefficient + 3 * geometry + 9 * wavelength + 18 * group) as f64;
                    }
                }
            }
        }
        let tangent = [1.0, -0.5, 0.25, 0.0, 2.0, -1.0];
        let mut output_tangent = [0.0; 4];
        interpolator
            .interpolate_jvp(
                &derivatives,
                stride,
                num_wavelengths,
                1,
                &tangent,
                &mut output_tangent,
            )
            .unwrap();

        let output_cotangent = [0.3, -0.7, 1.2, 0.4];
        let mut atmosphere_gradient = [0.0; 6];
        interpolator
            .accumulate_vjp(
                &output_cotangent,
                &derivatives,
                stride,
                num_wavelengths,
                1,
                &[0, 1],
                &mut atmosphere_gradient,
            )
            .unwrap();
        let forward_dot = output_tangent
            .iter()
            .zip(output_cotangent)
            .map(|(&left, right)| left * right)
            .sum::<f64>();
        let reverse_dot = tangent
            .iter()
            .zip(atmosphere_gradient)
            .map(|(&left, right)| left * right)
            .sum::<f64>();
        assert!((forward_dot - reverse_dot).abs() < 1.0e-12);
    }

    #[test]
    fn atmospheric_phase_products_are_dual_for_tensor_layout() {
        let num_atmosphere_points = 3;
        let num_wavelengths = 2;
        let num_groups = 2;
        let stride = 4;
        let retained_coefficients = 3;
        let mut derivatives =
            vec![0.0; stride * num_atmosphere_points * num_wavelengths * num_groups];
        for group in 0..num_groups {
            for wavelength in 0..num_wavelengths {
                for atmosphere in 0..num_atmosphere_points {
                    for coefficient in 0..stride {
                        let index = coefficient
                            + stride
                                * (atmosphere
                                    + num_atmosphere_points
                                        * (wavelength + num_wavelengths * group));
                        derivatives[index] = (index as f64 + 1.0) * 0.01;
                    }
                }
            }
        }
        let tangent = [0.4, -0.2, 0.7, 0.0, 0.3, -0.5];
        let mut coefficient_tangent = vec![0.0; num_atmosphere_points * retained_coefficients];
        atmospheric_coefficient_jvp(
            &derivatives,
            stride,
            num_atmosphere_points,
            num_wavelengths,
            1,
            &tangent,
            &mut coefficient_tangent,
        )
        .unwrap();

        let coefficient_cotangent = [0.3, -0.1, 0.8, 0.5, -0.4, 0.2, -0.7, 0.6, 0.9];
        let mut atmosphere_gradient = [0.0; 6];
        accumulate_atmospheric_coefficient_vjp(
            &coefficient_cotangent,
            &derivatives,
            stride,
            num_atmosphere_points,
            num_wavelengths,
            1,
            &[0, 1],
            &mut atmosphere_gradient,
        )
        .unwrap();
        let forward_dot = coefficient_tangent
            .iter()
            .zip(coefficient_cotangent)
            .map(|(&left, right)| left * right)
            .sum::<f64>();
        let reverse_dot = tangent
            .iter()
            .zip(atmosphere_gradient)
            .map(|(&left, right)| left * right)
            .sum::<f64>();
        assert!((forward_dot - reverse_dot).abs() < 1.0e-12);
    }

    #[test]
    fn rejects_invalid_geometry_and_runtime_shapes() {
        assert!(ScatteringCoefficientInterpolator::new(&[0, 2], &[0], &[1.0], 1, 1).is_err());
        let interpolator = interpolator();
        assert!(
            interpolator
                .interpolate(&[0.0; 9], 1, &mut [0.0; 4])
                .is_err()
        );
        assert!(
            interpolator
                .interpolate_jvp(&[], 3, 2, 2, &[], &mut [0.0; 4])
                .is_err()
        );
    }
}
