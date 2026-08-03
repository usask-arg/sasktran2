use super::OperatorError;

/// Geometry-only interpolation from the converged successive-orders state to
/// line-of-sight source samples.
///
/// The interpolation topology is immutable and shared by primal, JVP, and VJP
/// evaluations.  Wavelength-local radiance and atmosphere arrays remain owned
/// by the caller.  `solution_stride` permits direct reads from the
/// element-major wavelength-batch layout used by the engine without staging a
/// contiguous solution for each lane.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceInterpolator {
    ray_layer_offsets: Vec<u32>,
    layer_atmosphere_offsets: Vec<u32>,
    atmosphere_indices: Vec<u32>,
    atmosphere_weights: Vec<f64>,
    optical_depth_weights: Vec<f64>,
    layer_source_offsets: Vec<u32>,
    source_indices: Vec<u32>,
    source_weights: Vec<f64>,
    ray_ground_offsets: Vec<u32>,
    ground_source_indices: Vec<u32>,
    ground_source_weights: Vec<f64>,
    num_atmosphere_points: usize,
    num_solution_points: usize,
    num_stokes: usize,
}

impl SourceInterpolator {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ray_layer_offsets: &[u32],
        layer_atmosphere_offsets: &[u32],
        atmosphere_indices: &[u32],
        atmosphere_weights: &[f64],
        optical_depth_weights: &[f64],
        layer_source_offsets: &[u32],
        source_indices: &[u32],
        source_weights: &[f64],
        ray_ground_offsets: &[u32],
        ground_source_indices: &[u32],
        ground_source_weights: &[f64],
        num_atmosphere_points: usize,
        num_solution_points: usize,
        num_stokes: usize,
    ) -> Result<Self, OperatorError> {
        if num_atmosphere_points == 0
            || num_solution_points == 0
            || !matches!(num_stokes, 1 | 3)
            || !valid_offsets(ray_layer_offsets, None)
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let num_rays = ray_layer_offsets.len() - 1;
        let num_layers = ray_layer_offsets.last().copied().unwrap_or_default() as usize;
        if !valid_offsets(layer_atmosphere_offsets, Some(atmosphere_indices.len()))
            || layer_atmosphere_offsets.len() != num_layers + 1
            || atmosphere_indices.len() != atmosphere_weights.len()
            || atmosphere_indices.len() != optical_depth_weights.len()
            || !valid_offsets(layer_source_offsets, Some(source_indices.len()))
            || layer_source_offsets.len() != num_layers + 1
            || source_indices.len() != source_weights.len()
            || !valid_offsets(ray_ground_offsets, Some(ground_source_indices.len()))
            || ray_ground_offsets.len() != num_rays + 1
            || ground_source_indices.len() != ground_source_weights.len()
            || atmosphere_indices
                .iter()
                .any(|&index| index as usize >= num_atmosphere_points)
            || source_indices
                .iter()
                .chain(ground_source_indices)
                .any(|&index| index as usize >= num_solution_points)
            || atmosphere_weights
                .iter()
                .chain(optical_depth_weights)
                .chain(source_weights)
                .chain(ground_source_weights)
                .any(|value| !value.is_finite())
        {
            return Err(OperatorError::DimensionMismatch);
        }

        Ok(Self {
            ray_layer_offsets: ray_layer_offsets.to_vec(),
            layer_atmosphere_offsets: layer_atmosphere_offsets.to_vec(),
            atmosphere_indices: atmosphere_indices.to_vec(),
            atmosphere_weights: atmosphere_weights.to_vec(),
            optical_depth_weights: optical_depth_weights.to_vec(),
            layer_source_offsets: layer_source_offsets.to_vec(),
            source_indices: source_indices.to_vec(),
            source_weights: source_weights.to_vec(),
            ray_ground_offsets: ray_ground_offsets.to_vec(),
            ground_source_indices: ground_source_indices.to_vec(),
            ground_source_weights: ground_source_weights.to_vec(),
            num_atmosphere_points,
            num_solution_points,
            num_stokes,
        })
    }

    #[inline]
    pub fn num_rays(&self) -> usize {
        self.ray_layer_offsets.len() - 1
    }

    #[inline]
    pub fn num_layers(&self) -> usize {
        self.layer_atmosphere_offsets.len() - 1
    }

    pub fn storage_bytes(&self) -> usize {
        (self.ray_layer_offsets.capacity()
            + self.layer_atmosphere_offsets.capacity()
            + self.atmosphere_indices.capacity()
            + self.layer_source_offsets.capacity()
            + self.source_indices.capacity()
            + self.ray_ground_offsets.capacity()
            + self.ground_source_indices.capacity())
            * size_of::<u32>()
            + (self.atmosphere_weights.capacity()
                + self.optical_depth_weights.capacity()
                + self.source_weights.capacity()
                + self.ground_source_weights.capacity())
                * size_of::<f64>()
    }

    /// Integrates the complete successive-orders contribution along one LOS.
    ///
    /// Atmosphere values are wavelength-major, with atmosphere point varying
    /// fastest within each lane. The converged solution and output use the
    /// engine's element-major wavelength layout. Keeping both layouts explicit
    /// avoids staging either array at the C++/Rust boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_ray_value(
        &self,
        ray_index: usize,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        wavelength_count: usize,
        solution: &[f64],
        solution_stride: usize,
        source: &mut [f64],
        source_stride: usize,
    ) -> Result<(), OperatorError> {
        let layers = self.ray_layers(ray_index)?;
        self.validate_batched_runtime(
            extinction,
            single_scatter_albedo,
            wavelength_count,
            solution,
            solution_stride,
            source,
            source_stride,
        )?;

        for lane in 0..wavelength_count {
            let atmosphere_start = lane * self.num_atmosphere_points;
            let extinction =
                &extinction[atmosphere_start..atmosphere_start + self.num_atmosphere_points];
            let single_scatter_albedo = &single_scatter_albedo
                [atmosphere_start..atmosphere_start + self.num_atmosphere_points];
            let mut state = [0.0; 3];
            self.interpolate_ground(ray_index, solution, solution_stride, lane, &mut state);

            for layer in layers.clone() {
                let attenuation = (-self.interpolate_optical_depth(layer, extinction)).exp();
                if !attenuation.is_finite() {
                    return Err(OperatorError::NonFiniteValue);
                }
                let omega = self.interpolate_albedo(layer, single_scatter_albedo);
                let factor = 1.0 - attenuation;
                for (stokes, state_value) in state.iter_mut().enumerate().take(self.num_stokes) {
                    let outgoing =
                        self.interpolate_solution(layer, stokes, solution, solution_stride, lane);
                    *state_value = attenuation * *state_value + omega * factor * outgoing;
                }
            }

            for stokes in 0..self.num_stokes {
                source[stokes * source_stride + lane] += state[stokes];
            }
        }
        Ok(())
    }

    /// Integrates one LOS and adds its direct dense atmosphere Jacobian.
    ///
    /// The converged-state derivative is intentionally not included here; the
    /// Rust source exposes its complete Jacobian through the native JVP/VJP
    /// products.  This matches the previous per-layer dense callback while
    /// avoiding repeated scaling of the full derivative matrix.
    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_ray_jacobian(
        &self,
        ray_index: usize,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        wavelength_count: usize,
        solution: &[f64],
        solution_stride: usize,
        source: &mut [f64],
        source_stride: usize,
        source_derivative: &mut [f64],
        num_derivatives: usize,
        derivative_stride: usize,
        single_scatter_albedo_derivative_start: usize,
    ) -> Result<(), OperatorError> {
        let layers = self.ray_layers(ray_index)?;
        self.validate_batched_runtime(
            extinction,
            single_scatter_albedo,
            wavelength_count,
            solution,
            solution_stride,
            source,
            source_stride,
        )?;
        let required_derivatives = self
            .num_stokes
            .checked_mul(num_derivatives)
            .and_then(|size| size.checked_mul(derivative_stride))
            .ok_or(OperatorError::DimensionMismatch)?;
        if derivative_stride < wavelength_count
            || source_derivative.len() < required_derivatives
            || self.num_atmosphere_points > num_derivatives
            || single_scatter_albedo_derivative_start
                .checked_add(self.num_atmosphere_points)
                .is_none_or(|end| end > num_derivatives)
        {
            return Err(OperatorError::DimensionMismatch);
        }

        let num_ray_layers = layers.len();
        let vector_values = num_ray_layers * self.num_stokes;
        let outgoing_offset = vector_values;
        let attenuation_offset = outgoing_offset + vector_values;
        let omega_offset = attenuation_offset + num_ray_layers;
        let mut workspace = vec![0.0; omega_offset + num_ray_layers];

        for lane in 0..wavelength_count {
            let atmosphere_start = lane * self.num_atmosphere_points;
            let extinction =
                &extinction[atmosphere_start..atmosphere_start + self.num_atmosphere_points];
            let single_scatter_albedo = &single_scatter_albedo
                [atmosphere_start..atmosphere_start + self.num_atmosphere_points];
            let mut state = [0.0; 3];
            self.interpolate_ground(ray_index, solution, solution_stride, lane, &mut state);

            for (local_layer, layer) in layers.clone().enumerate() {
                let vector_start = local_layer * self.num_stokes;
                workspace[vector_start..vector_start + self.num_stokes]
                    .copy_from_slice(&state[..self.num_stokes]);
                let attenuation = (-self.interpolate_optical_depth(layer, extinction)).exp();
                if !attenuation.is_finite() {
                    return Err(OperatorError::NonFiniteValue);
                }
                let omega = self.interpolate_albedo(layer, single_scatter_albedo);
                workspace[attenuation_offset + local_layer] = attenuation;
                workspace[omega_offset + local_layer] = omega;
                let factor = 1.0 - attenuation;
                for (stokes, state_value) in state.iter_mut().enumerate().take(self.num_stokes) {
                    let outgoing =
                        self.interpolate_solution(layer, stokes, solution, solution_stride, lane);
                    workspace[outgoing_offset + vector_start + stokes] = outgoing;
                    *state_value = attenuation * *state_value + omega * factor * outgoing;
                }
            }

            for stokes in 0..self.num_stokes {
                source[stokes * source_stride + lane] += state[stokes];
            }

            let mut suffix_attenuation = 1.0;
            for (local_layer, layer) in layers.clone().enumerate().rev() {
                let vector_start = local_layer * self.num_stokes;
                let attenuation = workspace[attenuation_offset + local_layer];
                let omega = workspace[omega_offset + local_layer];
                let factor = 1.0 - attenuation;
                let atmosphere_start = self.layer_atmosphere_offsets[layer] as usize;
                let atmosphere_end = self.layer_atmosphere_offsets[layer + 1] as usize;
                for weight_index in atmosphere_start..atmosphere_end {
                    let atmosphere_index = self.atmosphere_indices[weight_index] as usize;
                    let extinction_weight = self.optical_depth_weights[weight_index];
                    let albedo_weight = self.atmosphere_weights[weight_index];
                    for stokes in 0..self.num_stokes {
                        let outgoing = workspace[outgoing_offset + vector_start + stokes];
                        let value_before = workspace[vector_start + stokes];
                        let extinction_output = ((atmosphere_index * self.num_stokes + stokes)
                            * derivative_stride)
                            + lane;
                        source_derivative[extinction_output] += suffix_attenuation
                            * extinction_weight
                            * attenuation
                            * (omega * outgoing - value_before);

                        let albedo_derivative =
                            single_scatter_albedo_derivative_start + atmosphere_index;
                        let albedo_output = ((albedo_derivative * self.num_stokes + stokes)
                            * derivative_stride)
                            + lane;
                        source_derivative[albedo_output] +=
                            suffix_attenuation * albedo_weight * factor * outgoing;
                    }
                }
                suffix_attenuation *= attenuation;
            }
        }
        Ok(())
    }

    /// Integrates the primal and one tangent together for one wavelength.
    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_ray_jvp(
        &self,
        ray_index: usize,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        extinction_tangent: &[f64],
        single_scatter_albedo_tangent: &[f64],
        solution: &[f64],
        solution_tangent: &[f64],
        source: &mut [f64],
        source_tangent: &mut [f64],
    ) -> Result<(), OperatorError> {
        let layers = self.ray_layers(ray_index)?;
        self.validate_runtime(single_scatter_albedo, solution, 1, 0)?;
        self.validate_runtime(single_scatter_albedo_tangent, solution_tangent, 1, 0)?;
        if extinction.len() != self.num_atmosphere_points
            || extinction_tangent.len() != self.num_atmosphere_points
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.validate_output(source, 1, 0)?;
        self.validate_output(source_tangent, 1, 0)?;

        let mut state = [0.0; 3];
        let mut state_tangent = [0.0; 3];
        self.interpolate_ground(ray_index, solution, 1, 0, &mut state);
        self.interpolate_ground(ray_index, solution_tangent, 1, 0, &mut state_tangent);

        for layer in layers {
            let optical_depth = self.interpolate_optical_depth(layer, extinction);
            let optical_depth_tangent =
                self.interpolate_optical_depth_tangent(layer, extinction_tangent);
            let attenuation = (-optical_depth).exp();
            if !attenuation.is_finite() {
                return Err(OperatorError::NonFiniteValue);
            }
            let omega = self.interpolate_albedo(layer, single_scatter_albedo);
            let omega_tangent = self.interpolate_albedo(layer, single_scatter_albedo_tangent);
            let factor = 1.0 - attenuation;
            let factor_tangent = attenuation * optical_depth_tangent;
            for stokes in 0..self.num_stokes {
                let value_before = state[stokes];
                let outgoing = self.interpolate_solution(layer, stokes, solution, 1, 0);
                let outgoing_tangent =
                    self.interpolate_solution(layer, stokes, solution_tangent, 1, 0);
                state[stokes] = attenuation * value_before + omega * factor * outgoing;
                state_tangent[stokes] = attenuation
                    * (state_tangent[stokes] - optical_depth_tangent * value_before)
                    + (omega_tangent * factor + omega * factor_tangent) * outgoing
                    + omega * factor * outgoing_tangent;
            }
        }

        for stokes in 0..self.num_stokes {
            source[stokes] += state[stokes];
            source_tangent[stokes] += state_tangent[stokes];
        }
        Ok(())
    }

    /// Pulls a LOS cotangent through the complete integration recurrence.
    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_ray_vjp(
        &self,
        ray_index: usize,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        solution: &[f64],
        cotangent: &[f64],
        solution_cotangent: &mut [f64],
        single_scatter_albedo_gradient: &mut [f64],
        extinction_gradient: &mut [f64],
    ) -> Result<(), OperatorError> {
        let layers = self.ray_layers(ray_index)?;
        self.validate_runtime(single_scatter_albedo, solution, 1, 0)?;
        if extinction.len() != self.num_atmosphere_points
            || cotangent.len() != self.num_stokes
            || solution_cotangent.len() != self.num_solution_points * self.num_stokes
            || single_scatter_albedo_gradient.len() != self.num_atmosphere_points
            || extinction_gradient.len() != self.num_atmosphere_points
        {
            return Err(OperatorError::DimensionMismatch);
        }
        if cotangent.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }

        let num_ray_layers = layers.len();
        let vector_values = num_ray_layers * self.num_stokes;
        let outgoing_offset = vector_values;
        let attenuation_offset = outgoing_offset + vector_values;
        let omega_offset = attenuation_offset + num_ray_layers;
        let mut workspace = vec![0.0; omega_offset + num_ray_layers];
        let mut state = [0.0; 3];
        self.interpolate_ground(ray_index, solution, 1, 0, &mut state);

        for (local_layer, layer) in layers.clone().enumerate() {
            let vector_start = local_layer * self.num_stokes;
            workspace[vector_start..vector_start + self.num_stokes]
                .copy_from_slice(&state[..self.num_stokes]);
            let attenuation = (-self.interpolate_optical_depth(layer, extinction)).exp();
            if !attenuation.is_finite() {
                return Err(OperatorError::NonFiniteValue);
            }
            let omega = self.interpolate_albedo(layer, single_scatter_albedo);
            workspace[attenuation_offset + local_layer] = attenuation;
            workspace[omega_offset + local_layer] = omega;
            let factor = 1.0 - attenuation;
            for (stokes, state_value) in state.iter_mut().enumerate().take(self.num_stokes) {
                let outgoing = self.interpolate_solution(layer, stokes, solution, 1, 0);
                workspace[outgoing_offset + vector_start + stokes] = outgoing;
                *state_value = attenuation * *state_value + omega * factor * outgoing;
            }
        }

        let mut state_cotangent = [0.0; 3];
        state_cotangent[..self.num_stokes].copy_from_slice(cotangent);
        for (local_layer, layer) in layers.enumerate().rev() {
            let vector_start = local_layer * self.num_stokes;
            let attenuation = workspace[attenuation_offset + local_layer];
            let omega = workspace[omega_offset + local_layer];
            let factor = 1.0 - attenuation;
            let mut albedo_cotangent = 0.0;
            let mut optical_depth_cotangent = 0.0;
            for stokes in 0..self.num_stokes {
                let outgoing = workspace[outgoing_offset + vector_start + stokes];
                let value_before = workspace[vector_start + stokes];
                albedo_cotangent += factor * state_cotangent[stokes] * outgoing;
                optical_depth_cotangent +=
                    attenuation * state_cotangent[stokes] * (omega * outgoing - value_before);
            }

            let source_start = self.layer_source_offsets[layer] as usize;
            let source_end = self.layer_source_offsets[layer + 1] as usize;
            for (&source_index, &weight) in self.source_indices[source_start..source_end]
                .iter()
                .zip(&self.source_weights[source_start..source_end])
            {
                let output_start = source_index as usize * self.num_stokes;
                for stokes in 0..self.num_stokes {
                    solution_cotangent[output_start + stokes] +=
                        weight * omega * factor * state_cotangent[stokes];
                }
            }

            let atmosphere_start = self.layer_atmosphere_offsets[layer] as usize;
            let atmosphere_end = self.layer_atmosphere_offsets[layer + 1] as usize;
            for index in atmosphere_start..atmosphere_end {
                let atmosphere_index = self.atmosphere_indices[index] as usize;
                single_scatter_albedo_gradient[atmosphere_index] +=
                    self.atmosphere_weights[index] * albedo_cotangent;
                extinction_gradient[atmosphere_index] +=
                    self.optical_depth_weights[index] * optical_depth_cotangent;
            }
            for value in &mut state_cotangent[..self.num_stokes] {
                *value *= attenuation;
            }
        }

        self.accumulate_ground_cotangent(ray_index, &state_cotangent, solution_cotangent);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_layer_value(
        &self,
        ray_index: usize,
        layer_index: usize,
        single_scatter_albedo: &[f64],
        attenuation: f64,
        solution: &[f64],
        solution_stride: usize,
        solution_lane: usize,
        source: &mut [f64],
        source_stride: usize,
        source_lane: usize,
    ) -> Result<(), OperatorError> {
        let layer = self.layer_index(ray_index, layer_index)?;
        self.validate_runtime(
            single_scatter_albedo,
            solution,
            solution_stride,
            solution_lane,
        )?;
        self.validate_output(source, source_stride, source_lane)?;
        if !attenuation.is_finite() {
            return Err(OperatorError::NonFiniteValue);
        }

        let omega = self.interpolate_albedo(layer, single_scatter_albedo);
        let amplitude = omega * (1.0 - attenuation);
        for stokes in 0..self.num_stokes {
            source[stokes * source_stride + source_lane] += amplitude
                * self.interpolate_solution(
                    layer,
                    stokes,
                    solution,
                    solution_stride,
                    solution_lane,
                );
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_layer_jvp(
        &self,
        ray_index: usize,
        layer_index: usize,
        single_scatter_albedo: &[f64],
        single_scatter_albedo_tangent: &[f64],
        extinction_tangent: &[f64],
        attenuation: f64,
        solution: &[f64],
        solution_tangent: &[f64],
        source: &mut [f64],
        source_tangent: &mut [f64],
    ) -> Result<(), OperatorError> {
        let layer = self.layer_index(ray_index, layer_index)?;
        self.validate_runtime(single_scatter_albedo, solution, 1, 0)?;
        self.validate_runtime(single_scatter_albedo_tangent, solution_tangent, 1, 0)?;
        if extinction_tangent.len() != self.num_atmosphere_points {
            return Err(OperatorError::DimensionMismatch);
        }
        self.validate_output(source, 1, 0)?;
        self.validate_output(source_tangent, 1, 0)?;
        if !attenuation.is_finite() {
            return Err(OperatorError::NonFiniteValue);
        }

        let omega = self.interpolate_albedo(layer, single_scatter_albedo);
        let omega_tangent = self.interpolate_albedo(layer, single_scatter_albedo_tangent);
        let optical_depth_tangent =
            self.interpolate_optical_depth_tangent(layer, extinction_tangent);
        let factor = 1.0 - attenuation;
        let factor_tangent = attenuation * optical_depth_tangent;
        for stokes in 0..self.num_stokes {
            let outgoing = self.interpolate_solution(layer, stokes, solution, 1, 0);
            let outgoing_tangent = self.interpolate_solution(layer, stokes, solution_tangent, 1, 0);
            source[stokes] += omega * factor * outgoing;
            source_tangent[stokes] += (omega_tangent * factor + omega * factor_tangent) * outgoing
                + omega * factor * outgoing_tangent;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_layer_vjp(
        &self,
        ray_index: usize,
        layer_index: usize,
        single_scatter_albedo: &[f64],
        attenuation: f64,
        solution: &[f64],
        cotangent: &[f64],
        source_value: &mut [f64],
        solution_cotangent: &mut [f64],
        single_scatter_albedo_gradient: &mut [f64],
        extinction_gradient: &mut [f64],
    ) -> Result<(), OperatorError> {
        let layer = self.layer_index(ray_index, layer_index)?;
        self.validate_runtime(single_scatter_albedo, solution, 1, 0)?;
        self.validate_output(source_value, 1, 0)?;
        if cotangent.len() != self.num_stokes
            || solution_cotangent.len() != self.num_solution_points * self.num_stokes
            || single_scatter_albedo_gradient.len() != self.num_atmosphere_points
            || extinction_gradient.len() != self.num_atmosphere_points
        {
            return Err(OperatorError::DimensionMismatch);
        }
        if !attenuation.is_finite() || cotangent.iter().any(|value| !value.is_finite()) {
            return Err(OperatorError::NonFiniteValue);
        }

        let omega = self.interpolate_albedo(layer, single_scatter_albedo);
        let factor = 1.0 - attenuation;
        let mut interpolated = [0.0; 3];
        let mut amplitude_cotangent = 0.0;
        for stokes in 0..self.num_stokes {
            interpolated[stokes] = self.interpolate_solution(layer, stokes, solution, 1, 0);
            source_value[stokes] += omega * factor * interpolated[stokes];
            amplitude_cotangent += cotangent[stokes] * interpolated[stokes];
        }

        let source_start = self.layer_source_offsets[layer] as usize;
        let source_end = self.layer_source_offsets[layer + 1] as usize;
        for (&source_index, &weight) in self.source_indices[source_start..source_end]
            .iter()
            .zip(&self.source_weights[source_start..source_end])
        {
            let output_start = source_index as usize * self.num_stokes;
            for stokes in 0..self.num_stokes {
                solution_cotangent[output_start + stokes] +=
                    weight * omega * factor * cotangent[stokes];
            }
        }

        let atmosphere_start = self.layer_atmosphere_offsets[layer] as usize;
        let atmosphere_end = self.layer_atmosphere_offsets[layer + 1] as usize;
        for (&atmosphere_index, &weight) in self.atmosphere_indices
            [atmosphere_start..atmosphere_end]
            .iter()
            .zip(&self.atmosphere_weights[atmosphere_start..atmosphere_end])
        {
            single_scatter_albedo_gradient[atmosphere_index as usize] +=
                weight * factor * amplitude_cotangent;
        }

        let optical_depth_cotangent = omega * attenuation * amplitude_cotangent;
        self.accumulate_optical_depth_gradient(layer, optical_depth_cotangent, extinction_gradient);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_layer_jacobian(
        &self,
        ray_index: usize,
        layer_index: usize,
        single_scatter_albedo: &[f64],
        attenuation: f64,
        solution: &[f64],
        solution_stride: usize,
        solution_lane: usize,
        source: &mut [f64],
        source_stride: usize,
        source_lane: usize,
        source_derivative: &mut [f64],
        num_derivatives: usize,
        derivative_stride: usize,
        single_scatter_albedo_derivative_start: usize,
    ) -> Result<(), OperatorError> {
        let layer = self.layer_index(ray_index, layer_index)?;
        self.validate_runtime(
            single_scatter_albedo,
            solution,
            solution_stride,
            solution_lane,
        )?;
        self.validate_output(source, source_stride, source_lane)?;
        let required_derivatives = self
            .num_stokes
            .checked_mul(num_derivatives)
            .and_then(|size| size.checked_mul(derivative_stride))
            .ok_or(OperatorError::DimensionMismatch)?;
        if derivative_stride == 0
            || source_lane >= derivative_stride
            || source_derivative.len() < required_derivatives
            || single_scatter_albedo_derivative_start
                .checked_add(self.num_atmosphere_points)
                .is_none_or(|end| end > num_derivatives)
            || self.num_atmosphere_points > num_derivatives
        {
            return Err(OperatorError::DimensionMismatch);
        }
        if !attenuation.is_finite() {
            return Err(OperatorError::NonFiniteValue);
        }

        let omega = self.interpolate_albedo(layer, single_scatter_albedo);
        let factor = 1.0 - attenuation;
        let atmosphere_start = self.layer_atmosphere_offsets[layer] as usize;
        let atmosphere_end = self.layer_atmosphere_offsets[layer + 1] as usize;
        for stokes in 0..self.num_stokes {
            let outgoing =
                self.interpolate_solution(layer, stokes, solution, solution_stride, solution_lane);
            source[stokes * source_stride + source_lane] += omega * factor * outgoing;

            for index in atmosphere_start..atmosphere_end {
                let atmosphere_index = self.atmosphere_indices[index] as usize;
                let extinction_output = ((atmosphere_index * self.num_stokes + stokes)
                    * derivative_stride)
                    + source_lane;
                source_derivative[extinction_output] +=
                    self.optical_depth_weights[index] * omega * attenuation * outgoing;

                let albedo_derivative = single_scatter_albedo_derivative_start + atmosphere_index;
                let albedo_output = ((albedo_derivative * self.num_stokes + stokes)
                    * derivative_stride)
                    + source_lane;
                source_derivative[albedo_output] +=
                    self.atmosphere_weights[index] * factor * outgoing;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn accumulate_ground_value(
        &self,
        ray_index: usize,
        solution: &[f64],
        solution_stride: usize,
        solution_lane: usize,
        source: &mut [f64],
        source_stride: usize,
        source_lane: usize,
    ) -> Result<(), OperatorError> {
        self.validate_ray(ray_index)?;
        self.validate_solution(solution, solution_stride, solution_lane)?;
        self.validate_output(source, source_stride, source_lane)?;
        let start = self.ray_ground_offsets[ray_index] as usize;
        let end = self.ray_ground_offsets[ray_index + 1] as usize;
        for (&source_index, &weight) in self.ground_source_indices[start..end]
            .iter()
            .zip(&self.ground_source_weights[start..end])
        {
            let output_start = source_index as usize * self.num_stokes;
            for stokes in 0..self.num_stokes {
                source[stokes * source_stride + source_lane] +=
                    weight * solution[(output_start + stokes) * solution_stride + solution_lane];
            }
        }
        Ok(())
    }

    pub fn accumulate_ground_jvp(
        &self,
        ray_index: usize,
        solution: &[f64],
        solution_tangent: &[f64],
        source: &mut [f64],
        source_tangent: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.accumulate_ground_value(ray_index, solution, 1, 0, source, 1, 0)?;
        self.accumulate_ground_value(ray_index, solution_tangent, 1, 0, source_tangent, 1, 0)
    }

    pub fn accumulate_ground_vjp(
        &self,
        ray_index: usize,
        cotangent: &[f64],
        solution_cotangent: &mut [f64],
    ) -> Result<(), OperatorError> {
        self.validate_ray(ray_index)?;
        if cotangent.len() != self.num_stokes
            || solution_cotangent.len() != self.num_solution_points * self.num_stokes
            || cotangent.iter().any(|value| !value.is_finite())
        {
            return Err(OperatorError::DimensionMismatch);
        }
        let start = self.ray_ground_offsets[ray_index] as usize;
        let end = self.ray_ground_offsets[ray_index + 1] as usize;
        for (&source_index, &weight) in self.ground_source_indices[start..end]
            .iter()
            .zip(&self.ground_source_weights[start..end])
        {
            let output_start = source_index as usize * self.num_stokes;
            for stokes in 0..self.num_stokes {
                solution_cotangent[output_start + stokes] += weight * cotangent[stokes];
            }
        }
        Ok(())
    }

    fn validate_runtime(
        &self,
        single_scatter_albedo: &[f64],
        solution: &[f64],
        solution_stride: usize,
        solution_lane: usize,
    ) -> Result<(), OperatorError> {
        if single_scatter_albedo.len() != self.num_atmosphere_points {
            return Err(OperatorError::DimensionMismatch);
        }
        self.validate_solution(solution, solution_stride, solution_lane)
    }

    #[allow(clippy::too_many_arguments)]
    fn validate_batched_runtime(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        wavelength_count: usize,
        solution: &[f64],
        solution_stride: usize,
        source: &[f64],
        source_stride: usize,
    ) -> Result<(), OperatorError> {
        let atmosphere_values = self
            .num_atmosphere_points
            .checked_mul(wavelength_count)
            .ok_or(OperatorError::DimensionMismatch)?;
        if wavelength_count == 0
            || extinction.len() < atmosphere_values
            || single_scatter_albedo.len() < atmosphere_values
            || solution_stride < wavelength_count
            || source_stride < wavelength_count
        {
            return Err(OperatorError::DimensionMismatch);
        }
        self.validate_solution(solution, solution_stride, wavelength_count - 1)?;
        self.validate_output(source, source_stride, wavelength_count - 1)
    }

    fn validate_output(
        &self,
        output: &[f64],
        output_stride: usize,
        output_lane: usize,
    ) -> Result<(), OperatorError> {
        let required = self
            .num_stokes
            .checked_mul(output_stride)
            .ok_or(OperatorError::DimensionMismatch)?;
        if output_stride == 0 || output_lane >= output_stride || output.len() < required {
            return Err(OperatorError::DimensionMismatch);
        }
        Ok(())
    }

    fn validate_solution(
        &self,
        solution: &[f64],
        solution_stride: usize,
        solution_lane: usize,
    ) -> Result<(), OperatorError> {
        let required = self
            .num_solution_points
            .checked_mul(self.num_stokes)
            .and_then(|size| size.checked_mul(solution_stride))
            .ok_or(OperatorError::DimensionMismatch)?;
        if solution_stride == 0 || solution_lane >= solution_stride || solution.len() < required {
            return Err(OperatorError::DimensionMismatch);
        }
        Ok(())
    }

    fn validate_ray(&self, ray_index: usize) -> Result<(), OperatorError> {
        if ray_index >= self.num_rays() {
            return Err(OperatorError::DimensionMismatch);
        }
        Ok(())
    }

    fn ray_layers(&self, ray_index: usize) -> Result<std::ops::Range<usize>, OperatorError> {
        self.validate_ray(ray_index)?;
        Ok(self.ray_layer_offsets[ray_index] as usize
            ..self.ray_layer_offsets[ray_index + 1] as usize)
    }

    fn layer_index(&self, ray_index: usize, layer_index: usize) -> Result<usize, OperatorError> {
        self.validate_ray(ray_index)?;
        let start = self.ray_layer_offsets[ray_index] as usize;
        let end = self.ray_layer_offsets[ray_index + 1] as usize;
        if layer_index >= end - start {
            return Err(OperatorError::DimensionMismatch);
        }
        Ok(start + layer_index)
    }

    fn interpolate_albedo(&self, layer: usize, single_scatter_albedo: &[f64]) -> f64 {
        let start = self.layer_atmosphere_offsets[layer] as usize;
        let end = self.layer_atmosphere_offsets[layer + 1] as usize;
        self.atmosphere_indices[start..end]
            .iter()
            .zip(&self.atmosphere_weights[start..end])
            .map(|(&index, &weight)| single_scatter_albedo[index as usize] * weight)
            .sum()
    }

    fn interpolate_optical_depth_tangent(&self, layer: usize, extinction_tangent: &[f64]) -> f64 {
        let start = self.layer_atmosphere_offsets[layer] as usize;
        let end = self.layer_atmosphere_offsets[layer + 1] as usize;
        self.atmosphere_indices[start..end]
            .iter()
            .zip(&self.optical_depth_weights[start..end])
            .map(|(&index, &weight)| extinction_tangent[index as usize] * weight)
            .sum()
    }

    fn interpolate_optical_depth(&self, layer: usize, extinction: &[f64]) -> f64 {
        let start = self.layer_atmosphere_offsets[layer] as usize;
        let end = self.layer_atmosphere_offsets[layer + 1] as usize;
        self.atmosphere_indices[start..end]
            .iter()
            .zip(&self.optical_depth_weights[start..end])
            .map(|(&index, &weight)| extinction[index as usize] * weight)
            .sum()
    }

    fn accumulate_optical_depth_gradient(
        &self,
        layer: usize,
        cotangent: f64,
        extinction_gradient: &mut [f64],
    ) {
        let start = self.layer_atmosphere_offsets[layer] as usize;
        let end = self.layer_atmosphere_offsets[layer + 1] as usize;
        for (&index, &weight) in self.atmosphere_indices[start..end]
            .iter()
            .zip(&self.optical_depth_weights[start..end])
        {
            extinction_gradient[index as usize] += weight * cotangent;
        }
    }

    fn interpolate_solution(
        &self,
        layer: usize,
        stokes: usize,
        solution: &[f64],
        solution_stride: usize,
        solution_lane: usize,
    ) -> f64 {
        let start = self.layer_source_offsets[layer] as usize;
        let end = self.layer_source_offsets[layer + 1] as usize;
        self.source_indices[start..end]
            .iter()
            .zip(&self.source_weights[start..end])
            .map(|(&index, &weight)| {
                let solution_index = index as usize * self.num_stokes + stokes;
                weight * solution[solution_index * solution_stride + solution_lane]
            })
            .sum()
    }

    fn interpolate_ground(
        &self,
        ray_index: usize,
        solution: &[f64],
        solution_stride: usize,
        solution_lane: usize,
        output: &mut [f64; 3],
    ) {
        let start = self.ray_ground_offsets[ray_index] as usize;
        let end = self.ray_ground_offsets[ray_index + 1] as usize;
        for (&source_index, &weight) in self.ground_source_indices[start..end]
            .iter()
            .zip(&self.ground_source_weights[start..end])
        {
            let output_start = source_index as usize * self.num_stokes;
            for stokes in 0..self.num_stokes {
                output[stokes] +=
                    weight * solution[(output_start + stokes) * solution_stride + solution_lane];
            }
        }
    }

    fn accumulate_ground_cotangent(
        &self,
        ray_index: usize,
        cotangent: &[f64; 3],
        solution_cotangent: &mut [f64],
    ) {
        let start = self.ray_ground_offsets[ray_index] as usize;
        let end = self.ray_ground_offsets[ray_index + 1] as usize;
        for (&source_index, &weight) in self.ground_source_indices[start..end]
            .iter()
            .zip(&self.ground_source_weights[start..end])
        {
            let output_start = source_index as usize * self.num_stokes;
            for stokes in 0..self.num_stokes {
                solution_cotangent[output_start + stokes] += weight * cotangent[stokes];
            }
        }
    }
}

fn valid_offsets(offsets: &[u32], expected_values: Option<usize>) -> bool {
    !offsets.is_empty()
        && offsets[0] == 0
        && offsets.windows(2).all(|window| window[0] <= window[1])
        && expected_values
            .is_none_or(|expected| offsets.last().copied().unwrap_or_default() as usize == expected)
}

#[cfg(test)]
mod tests {
    use super::SourceInterpolator;

    fn interpolator(num_stokes: usize) -> SourceInterpolator {
        SourceInterpolator::new(
            &[0, 1],
            &[0, 2],
            &[0, 1],
            &[0.25, 0.75],
            &[0.4, 0.6],
            &[0, 2],
            &[0, 1],
            &[0.4, 0.6],
            &[0, 1],
            &[1],
            &[0.8],
            2,
            2,
            num_stokes,
        )
        .unwrap()
    }

    fn two_layer_interpolator(num_stokes: usize) -> SourceInterpolator {
        SourceInterpolator::new(
            &[0, 2],
            &[0, 2, 4],
            &[0, 1, 1, 2],
            &[0.25, 0.75, 0.4, 0.6],
            &[0.4, 0.6, 0.2, 0.8],
            &[0, 2, 4],
            &[0, 1, 1, 2],
            &[0.4, 0.6, 0.3, 0.7],
            &[0, 1],
            &[2],
            &[0.8],
            3,
            3,
            num_stokes,
        )
        .unwrap()
    }

    #[test]
    fn scalar_primal_jvp_and_vjp_are_consistent() {
        let interpolation = interpolator(1);
        let albedo = [0.2, 0.6];
        let albedo_tangent = [0.1, -0.3];
        let solution = [2.0, 5.0];
        let solution_tangent = [-0.4, 0.7];
        let extinction_tangent = [0.5, 0.0];
        let attenuation = 0.3;
        let mut value = [0.0];
        let mut tangent = [0.0];
        interpolation
            .accumulate_layer_jvp(
                0,
                0,
                &albedo,
                &albedo_tangent,
                &extinction_tangent,
                attenuation,
                &solution,
                &solution_tangent,
                &mut value,
                &mut tangent,
            )
            .unwrap();

        let cotangent = [0.9];
        let mut reverse_value = [0.0];
        let mut solution_gradient = [0.0; 2];
        let mut albedo_gradient = [0.0; 2];
        let mut extinction_gradient = [0.0; 2];
        interpolation
            .accumulate_layer_vjp(
                0,
                0,
                &albedo,
                attenuation,
                &solution,
                &cotangent,
                &mut reverse_value,
                &mut solution_gradient,
                &mut albedo_gradient,
                &mut extinction_gradient,
            )
            .unwrap();
        let reverse_dot = albedo_tangent
            .iter()
            .zip(albedo_gradient)
            .map(|(&left, right)| left * right)
            .sum::<f64>()
            + solution_tangent
                .iter()
                .zip(solution_gradient)
                .map(|(&left, right)| left * right)
                .sum::<f64>()
            + extinction_tangent
                .iter()
                .zip(extinction_gradient)
                .map(|(&left, right)| left * right)
                .sum::<f64>();
        assert!((value[0] - reverse_value[0]).abs() < 1.0e-14);
        assert!((tangent[0] * cotangent[0] - reverse_dot).abs() < 1.0e-14);

        let mut jacobian_value = [0.0];
        let mut jacobian = [0.0; 4];
        interpolation
            .accumulate_layer_jacobian(
                0,
                0,
                &albedo,
                attenuation,
                &solution,
                1,
                0,
                &mut jacobian_value,
                1,
                0,
                &mut jacobian,
                4,
                1,
                2,
            )
            .unwrap();
        let direct_tangent = [
            extinction_tangent[0],
            extinction_tangent[1],
            albedo_tangent[0],
            albedo_tangent[1],
        ];
        let jacobian_jvp = jacobian
            .iter()
            .zip(direct_tangent)
            .map(|(&left, right)| left * right)
            .sum::<f64>();
        let mut direct_value = [0.0];
        let mut direct_jvp = [0.0];
        interpolation
            .accumulate_layer_jvp(
                0,
                0,
                &albedo,
                &albedo_tangent,
                &extinction_tangent,
                attenuation,
                &solution,
                &[0.0; 2],
                &mut direct_value,
                &mut direct_jvp,
            )
            .unwrap();
        assert!((jacobian_value[0] - direct_value[0]).abs() < 1.0e-14);
        assert!((jacobian_jvp - direct_jvp[0]).abs() < 1.0e-14);
    }

    #[test]
    fn vector_batch_and_ground_interpolation_use_requested_lane() {
        let interpolation = interpolator(3);
        let solution_batch = [
            1.0, 10.0, 2.0, 20.0, 3.0, 30.0, 4.0, 40.0, 5.0, 50.0, 6.0, 60.0,
        ];
        let mut layer = [0.0; 3];
        interpolation
            .accumulate_layer_value(
                0,
                0,
                &[0.2, 0.6],
                0.25,
                &solution_batch,
                2,
                1,
                &mut layer,
                1,
                0,
            )
            .unwrap();
        let omega = 0.25 * 0.2 + 0.75 * 0.6;
        for (actual, expected) in layer.into_iter().zip([28.0, 38.0, 48.0]) {
            assert!((actual - omega * 0.75 * expected).abs() < 1.0e-14);
        }

        let mut ground = [0.0; 3];
        interpolation
            .accumulate_ground_value(0, &solution_batch, 2, 1, &mut ground, 1, 0)
            .unwrap();
        for (actual, expected) in ground.into_iter().zip([40.0, 50.0, 60.0]) {
            assert!((actual - 0.8 * expected).abs() < 1.0e-14);
        }
    }

    #[test]
    fn whole_ray_primal_jacobian_jvp_and_vjp_are_consistent() {
        let interpolation = two_layer_interpolator(1);
        let extinction = [0.1, 0.2, 0.05];
        let albedo = [0.3, 0.6, 0.4];
        let solution = [2.0, 5.0, 7.0];
        let extinction_tangent = [0.2, -0.1, 0.4];
        let albedo_tangent = [-0.3, 0.25, 0.1];
        let solution_tangent = [0.5, -0.2, 0.3];

        let mut value = [0.0];
        interpolation
            .accumulate_ray_value(0, &extinction, &albedo, 1, &solution, 1, &mut value, 1)
            .unwrap();

        let attenuation_0 = (-(0.4 * extinction[0] + 0.6 * extinction[1])).exp();
        let attenuation_1 = (-(0.2 * extinction[1] + 0.8 * extinction[2])).exp();
        let omega_0 = 0.25 * albedo[0] + 0.75 * albedo[1];
        let omega_1 = 0.4 * albedo[1] + 0.6 * albedo[2];
        let outgoing_0 = 0.4 * solution[0] + 0.6 * solution[1];
        let outgoing_1 = 0.3 * solution[1] + 0.7 * solution[2];
        let expected = attenuation_1
            * (attenuation_0 * 0.8 * solution[2] + omega_0 * (1.0 - attenuation_0) * outgoing_0)
            + omega_1 * (1.0 - attenuation_1) * outgoing_1;
        assert!((value[0] - expected).abs() < 1.0e-14);

        let mut jvp_value = [0.0];
        let mut tangent = [0.0];
        interpolation
            .accumulate_ray_jvp(
                0,
                &extinction,
                &albedo,
                &extinction_tangent,
                &albedo_tangent,
                &solution,
                &solution_tangent,
                &mut jvp_value,
                &mut tangent,
            )
            .unwrap();
        assert!((jvp_value[0] - value[0]).abs() < 1.0e-14);

        let cotangent = [0.7];
        let mut solution_gradient = [0.0; 3];
        let mut albedo_gradient = [0.0; 3];
        let mut extinction_gradient = [0.0; 3];
        interpolation
            .accumulate_ray_vjp(
                0,
                &extinction,
                &albedo,
                &solution,
                &cotangent,
                &mut solution_gradient,
                &mut albedo_gradient,
                &mut extinction_gradient,
            )
            .unwrap();
        let reverse_dot = extinction_tangent
            .iter()
            .zip(extinction_gradient)
            .map(|(&left, right)| left * right)
            .sum::<f64>()
            + albedo_tangent
                .iter()
                .zip(albedo_gradient)
                .map(|(&left, right)| left * right)
                .sum::<f64>()
            + solution_tangent
                .iter()
                .zip(solution_gradient)
                .map(|(&left, right)| left * right)
                .sum::<f64>();
        assert!((tangent[0] * cotangent[0] - reverse_dot).abs() < 1.0e-13);

        let mut jacobian_value = [0.0];
        let mut jacobian = [0.0; 6];
        interpolation
            .accumulate_ray_jacobian(
                0,
                &extinction,
                &albedo,
                1,
                &solution,
                1,
                &mut jacobian_value,
                1,
                &mut jacobian,
                6,
                1,
                3,
            )
            .unwrap();
        let direct_tangent = [
            extinction_tangent[0],
            extinction_tangent[1],
            extinction_tangent[2],
            albedo_tangent[0],
            albedo_tangent[1],
            albedo_tangent[2],
        ];
        let jacobian_jvp = jacobian
            .iter()
            .zip(direct_tangent)
            .map(|(&left, right)| left * right)
            .sum::<f64>();
        let mut direct_value = [0.0];
        let mut direct_jvp = [0.0];
        interpolation
            .accumulate_ray_jvp(
                0,
                &extinction,
                &albedo,
                &extinction_tangent,
                &albedo_tangent,
                &solution,
                &[0.0; 3],
                &mut direct_value,
                &mut direct_jvp,
            )
            .unwrap();
        assert!((jacobian_value[0] - value[0]).abs() < 1.0e-14);
        assert!((jacobian_jvp - direct_jvp[0]).abs() < 1.0e-13);
    }
}
