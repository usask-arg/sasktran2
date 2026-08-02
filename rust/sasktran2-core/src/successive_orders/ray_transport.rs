//! Scalar transport assembly for traced successive-orders rays.
//!
//! The ray tracer remains responsible for reducing geometry to interpolation
//! stencils.  This module owns the wavelength-independent packed stencils and
//! applies wavelength-dependent extinction and single-scatter albedo directly
//! to the successive-orders CSR value array.

use anyhow::{Result, anyhow, bail};

/// Wavelength-independent scalar ray geometry compiled for transport assembly.
///
/// Layers for ray `r` occupy
/// `ray_layer_offsets[r]..ray_layer_offsets[r + 1]`.  Every layer has a shared
/// atmospheric index stencil with independent optical-depth and albedo weights.
/// Source and ground entries contain absolute indices into the already-compiled
/// CSR value array, so wavelength assembly does not need to search columns.
pub struct ScalarRayTransport {
    ray_layer_offsets: Vec<u32>,
    layer_atmosphere_offsets: Vec<u32>,
    atmosphere_indices: Vec<u32>,
    optical_depth_weights: Vec<f64>,
    albedo_weights: Vec<f64>,
    entrance_weights: Vec<f64>,
    exit_weights: Vec<f64>,
    layer_distance: Vec<f64>,
    layer_start_fraction: Vec<f64>,
    ray_phase_basis: Vec<f64>,
    num_phase_moments: usize,
    solar_transmission_on_atmosphere_grid: bool,
    layer_source_offsets: Vec<u32>,
    ray_transport_value_offsets: Vec<u32>,
    source_value_inner_indices: Vec<u16>,
    source_weights: Vec<f64>,
    ray_ground_offsets: Vec<u32>,
    ground_value_inner_indices: Vec<u16>,
    ground_weights: Vec<f64>,
    required_atmosphere_len: usize,
    required_transport_len: usize,
}

impl ScalarRayTransport {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ray_layer_offsets: &[u32],
        layer_atmosphere_offsets: &[u32],
        atmosphere_indices: &[u32],
        optical_depth_weights: &[f64],
        albedo_weights: &[f64],
        entrance_weights: &[f64],
        exit_weights: &[f64],
        layer_distance: &[f64],
        layer_start_fraction: &[f64],
        layer_end_fraction: &[f64],
        ray_scattering_cosine: &[f64],
        num_phase_moments: usize,
        solar_transmission_on_atmosphere_grid: bool,
        layer_source_offsets: &[u32],
        ray_transport_value_offsets: &[u32],
        source_value_inner_indices: &[u16],
        source_weights: &[f64],
        ray_ground_offsets: &[u32],
        ground_value_inner_indices: &[u16],
        ground_weights: &[f64],
    ) -> Result<Self> {
        validate_offsets("ray layer", ray_layer_offsets)?;
        validate_offsets("layer atmosphere", layer_atmosphere_offsets)?;
        validate_offsets("layer source", layer_source_offsets)?;
        validate_offsets("ray ground", ray_ground_offsets)?;

        let num_rays = ray_layer_offsets.len().saturating_sub(1);
        let num_layers = usize::try_from(*ray_layer_offsets.last().unwrap_or(&0))?;
        if layer_atmosphere_offsets.len() != num_layers + 1 {
            bail!(
                "layer-atmosphere offsets have length {}, expected {}",
                layer_atmosphere_offsets.len(),
                num_layers + 1
            );
        }
        if layer_source_offsets.len() != num_layers + 1 {
            bail!(
                "layer-source offsets have length {}, expected {}",
                layer_source_offsets.len(),
                num_layers + 1
            );
        }
        if ray_ground_offsets.len() != num_rays + 1 {
            bail!(
                "ray-ground offsets have length {}, expected {}",
                ray_ground_offsets.len(),
                num_rays + 1
            );
        }
        if ray_transport_value_offsets.len() != num_rays {
            bail!(
                "ray transport value offsets have length {}, expected {}",
                ray_transport_value_offsets.len(),
                num_rays
            );
        }
        if layer_distance.len() != num_layers
            || layer_start_fraction.len() != num_layers
            || layer_end_fraction.len() != num_layers
        {
            bail!(
                "layer geometry arrays have lengths ({}, {}, {}), expected {}",
                layer_distance.len(),
                layer_start_fraction.len(),
                layer_end_fraction.len(),
                num_layers
            );
        }
        if layer_start_fraction
            .iter()
            .zip(layer_end_fraction)
            .any(|(&start, &end)| {
                !start.is_finite() || !end.is_finite() || (start + end - 1.0).abs() > 1.0e-10
            })
        {
            bail!("layer endpoint quadrature fractions must sum to one");
        }
        if num_phase_moments == 0 {
            bail!("scalar first-order forcing requires at least one phase moment");
        }
        if ray_scattering_cosine.len() != num_rays {
            bail!(
                "ray scattering cosine array has length {}, expected {}",
                ray_scattering_cosine.len(),
                num_rays
            );
        }

        let num_atmosphere_weights = checked_last(layer_atmosphere_offsets)?;
        if atmosphere_indices.len() != num_atmosphere_weights
            || optical_depth_weights.len() != num_atmosphere_weights
            || albedo_weights.len() != num_atmosphere_weights
            || entrance_weights.len() != num_atmosphere_weights
            || exit_weights.len() != num_atmosphere_weights
        {
            bail!(
                "atmosphere stencil arrays have inconsistent lengths: indices {}, optical depth {}, albedo {}, entrance {}, exit {}, expected {}",
                atmosphere_indices.len(),
                optical_depth_weights.len(),
                albedo_weights.len(),
                entrance_weights.len(),
                exit_weights.len(),
                num_atmosphere_weights
            );
        }

        let num_source_weights = checked_last(layer_source_offsets)?;
        if source_value_inner_indices.len() != num_source_weights
            || source_weights.len() != num_source_weights
        {
            bail!(
                "source stencil arrays have inconsistent lengths: indices {}, weights {}, expected {}",
                source_value_inner_indices.len(),
                source_weights.len(),
                num_source_weights
            );
        }

        let num_ground_weights = checked_last(ray_ground_offsets)?;
        if ground_value_inner_indices.len() != num_ground_weights
            || ground_weights.len() != num_ground_weights
        {
            bail!(
                "ground stencil arrays have inconsistent lengths: indices {}, weights {}, expected {}",
                ground_value_inner_indices.len(),
                ground_weights.len(),
                num_ground_weights
            );
        }

        let required_atmosphere_len = atmosphere_indices
            .iter()
            .copied()
            .max()
            .map_or(0, |index| index as usize + 1);
        let mut required_transport_len = 0;
        for ray_index in 0..num_rays {
            let row_offset = ray_transport_value_offsets[ray_index] as usize;
            let layer_start = ray_layer_offsets[ray_index] as usize;
            let layer_end = ray_layer_offsets[ray_index + 1] as usize;
            let source_start = layer_source_offsets[layer_start] as usize;
            let source_end = layer_source_offsets[layer_end] as usize;
            for &inner_index in &source_value_inner_indices[source_start..source_end] {
                required_transport_len =
                    required_transport_len.max(row_offset + inner_index as usize + 1);
            }
            let ground_start = ray_ground_offsets[ray_index] as usize;
            let ground_end = ray_ground_offsets[ray_index + 1] as usize;
            for &inner_index in &ground_value_inner_indices[ground_start..ground_end] {
                required_transport_len =
                    required_transport_len.max(row_offset + inner_index as usize + 1);
            }
        }
        let mut ray_phase_basis = Vec::with_capacity(num_rays * num_phase_moments);
        for &cosine in ray_scattering_cosine {
            if !cosine.is_finite() {
                bail!("ray scattering cosine is not finite");
            }
            append_legendre_basis(
                cosine.clamp(-1.0, 1.0),
                num_phase_moments,
                &mut ray_phase_basis,
            );
        }

        Ok(Self {
            ray_layer_offsets: ray_layer_offsets.to_vec(),
            layer_atmosphere_offsets: layer_atmosphere_offsets.to_vec(),
            atmosphere_indices: atmosphere_indices.to_vec(),
            optical_depth_weights: optical_depth_weights.to_vec(),
            albedo_weights: albedo_weights.to_vec(),
            entrance_weights: entrance_weights.to_vec(),
            exit_weights: exit_weights.to_vec(),
            layer_distance: layer_distance.to_vec(),
            layer_start_fraction: layer_start_fraction.to_vec(),
            ray_phase_basis,
            num_phase_moments,
            solar_transmission_on_atmosphere_grid,
            layer_source_offsets: layer_source_offsets.to_vec(),
            ray_transport_value_offsets: ray_transport_value_offsets.to_vec(),
            source_value_inner_indices: source_value_inner_indices.to_vec(),
            source_weights: source_weights.to_vec(),
            ray_ground_offsets: ray_ground_offsets.to_vec(),
            ground_value_inner_indices: ground_value_inner_indices.to_vec(),
            ground_weights: ground_weights.to_vec(),
            required_atmosphere_len,
            required_transport_len,
        })
    }

    pub fn num_rays(&self) -> usize {
        self.ray_layer_offsets.len() - 1
    }

    pub fn num_layers(&self) -> usize {
        self.layer_atmosphere_offsets.len() - 1
    }

    pub fn num_atmosphere_weights(&self) -> usize {
        self.atmosphere_indices.len()
    }

    pub fn num_source_weights(&self) -> usize {
        self.source_value_inner_indices.len()
    }

    pub fn num_ground_weights(&self) -> usize {
        self.ground_value_inner_indices.len()
    }

    /// Persistent heap storage owned by this geometry, excluding `Vec` headers.
    pub fn storage_bytes(&self) -> usize {
        self.ray_layer_offsets.capacity() * size_of::<u32>()
            + self.layer_atmosphere_offsets.capacity() * size_of::<u32>()
            + self.atmosphere_indices.capacity() * size_of::<u32>()
            + self.optical_depth_weights.capacity() * size_of::<f64>()
            + self.albedo_weights.capacity() * size_of::<f64>()
            + self.entrance_weights.capacity() * size_of::<f64>()
            + self.exit_weights.capacity() * size_of::<f64>()
            + self.layer_distance.capacity() * size_of::<f64>()
            + self.layer_start_fraction.capacity() * size_of::<f64>()
            + self.ray_phase_basis.capacity() * size_of::<f64>()
            + self.layer_source_offsets.capacity() * size_of::<u32>()
            + self.ray_transport_value_offsets.capacity() * size_of::<u32>()
            + self.source_value_inner_indices.capacity() * size_of::<u16>()
            + self.source_weights.capacity() * size_of::<f64>()
            + self.ray_ground_offsets.capacity() * size_of::<u32>()
            + self.ground_value_inner_indices.capacity() * size_of::<u16>()
            + self.ground_weights.capacity() * size_of::<f64>()
    }

    /// Assemble one scalar wavelength without allocating.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        transport_values: &mut [f64],
        layer_optical_depth: &mut [f64],
        layer_attenuation: &mut [f64],
        layer_prefix_attenuation: &mut [f64],
        ray_end_attenuation: &mut [f64],
    ) -> Result<()> {
        self.validate_assembly_buffers(
            extinction,
            single_scatter_albedo,
            transport_values,
            layer_optical_depth,
            layer_attenuation,
            layer_prefix_attenuation,
            ray_end_attenuation,
        )?;
        self.assemble_impl(
            extinction,
            single_scatter_albedo,
            transport_values,
            layer_optical_depth,
            layer_attenuation,
            layer_prefix_attenuation,
            ray_end_attenuation,
            None,
        );
        Ok(())
    }

    /// Fuses transport assembly with scalar exact single-scatter layer
    /// forcing. Solar-path attenuation is supplied by the existing source for
    /// now; endpoint optical-property interpolation, phase evaluation, layer
    /// quadrature, and observer-path attenuation are all evaluated here.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_with_first_order(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        legendre_coefficients: &[f64],
        maximum_order: &[i32],
        solar_transmission: &[f64],
        transport_values: &mut [f64],
        first_order_radiance: &mut [f64],
        layer_optical_depth: &mut [f64],
        layer_attenuation: &mut [f64],
        layer_prefix_attenuation: &mut [f64],
        ray_end_attenuation: &mut [f64],
    ) -> Result<()> {
        self.validate_assembly_buffers(
            extinction,
            single_scatter_albedo,
            transport_values,
            layer_optical_depth,
            layer_attenuation,
            layer_prefix_attenuation,
            ray_end_attenuation,
        )?;
        if maximum_order.len() != extinction.len() {
            bail!(
                "maximum phase order length {} does not match {} atmosphere values",
                maximum_order.len(),
                extinction.len()
            );
        }
        let required_coefficients = extinction
            .len()
            .checked_mul(self.num_phase_moments)
            .ok_or_else(|| anyhow!("phase coefficient size overflow"))?;
        if legendre_coefficients.len() != required_coefficients {
            bail!(
                "phase coefficient array has length {}, expected {}",
                legendre_coefficients.len(),
                required_coefficients
            );
        }
        if maximum_order.iter().any(|&order| {
            order < 1 || usize::try_from(order).map_or(true, |value| value > self.num_phase_moments)
        }) {
            bail!(
                "maximum phase orders must lie in 1..={}",
                self.num_phase_moments
            );
        }
        let required_solar_transmission = if self.solar_transmission_on_atmosphere_grid {
            self.required_atmosphere_len
        } else {
            self.num_layers() + self.num_rays()
        };
        if solar_transmission.len() < required_solar_transmission {
            bail!(
                "solar transmission has length {}, but forcing requires {}",
                solar_transmission.len(),
                required_solar_transmission
            );
        }
        if first_order_radiance.len() != self.num_rays() {
            bail!(
                "first-order radiance length {} does not match {} rays",
                first_order_radiance.len(),
                self.num_rays()
            );
        }

        let forcing = FirstOrderInputs {
            legendre_coefficients,
            maximum_order,
            solar_transmission,
        };
        self.assemble_impl(
            extinction,
            single_scatter_albedo,
            transport_values,
            layer_optical_depth,
            layer_attenuation,
            layer_prefix_attenuation,
            ray_end_attenuation,
            Some((forcing, first_order_radiance)),
        );
        Ok(())
    }

    /// Directional derivative of the fused scalar transport and atmospheric
    /// first-order forcing. All output buffers are overwritten and no
    /// wavelength-sized allocation is performed.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_jvp_with_first_order(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        legendre_coefficients: &[f64],
        maximum_order: &[i32],
        solar_transmission: &[f64],
        extinction_tangent: &[f64],
        single_scatter_albedo_tangent: &[f64],
        legendre_coefficient_tangent: &[f64],
        solar_transmission_tangent: &[f64],
        transport_values: &mut [f64],
        transport_value_tangent: &mut [f64],
        first_order_radiance_tangent: &mut [f64],
        ray_end_attenuation: &mut [f64],
        ray_end_attenuation_tangent: &mut [f64],
    ) -> Result<()> {
        if extinction.len() < self.required_atmosphere_len
            || single_scatter_albedo.len() != extinction.len()
            || extinction_tangent.len() != extinction.len()
            || single_scatter_albedo_tangent.len() != extinction.len()
        {
            bail!("scalar JVP atmosphere arrays have inconsistent lengths");
        }
        let required_coefficients = extinction
            .len()
            .checked_mul(self.num_phase_moments)
            .ok_or_else(|| anyhow!("phase coefficient size overflow"))?;
        if legendre_coefficients.len() != required_coefficients
            || legendre_coefficient_tangent.len() != required_coefficients
            || maximum_order.len() != extinction.len()
        {
            bail!("scalar JVP phase arrays have inconsistent lengths");
        }
        if maximum_order.iter().any(|&order| {
            order < 1 || usize::try_from(order).map_or(true, |value| value > self.num_phase_moments)
        }) {
            bail!(
                "maximum phase orders must lie in 1..={}",
                self.num_phase_moments
            );
        }
        let required_solar_transmission = if self.solar_transmission_on_atmosphere_grid {
            self.required_atmosphere_len
        } else {
            self.num_layers() + self.num_rays()
        };
        if solar_transmission.len() < required_solar_transmission
            || solar_transmission_tangent.len() != solar_transmission.len()
        {
            bail!("scalar JVP solar-transmission arrays have inconsistent lengths");
        }
        if transport_values.len() < self.required_transport_len
            || transport_value_tangent.len() != transport_values.len()
        {
            bail!("scalar JVP transport arrays have inconsistent lengths");
        }
        if first_order_radiance_tangent.len() != self.num_rays()
            || ray_end_attenuation.len() != self.num_rays()
            || ray_end_attenuation_tangent.len() != self.num_rays()
        {
            bail!("scalar JVP ray output arrays have inconsistent lengths");
        }

        transport_values.fill(0.0);
        transport_value_tangent.fill(0.0);
        first_order_radiance_tangent.fill(0.0);
        let forcing = FirstOrderJvpInputs {
            legendre_coefficients,
            maximum_order,
            solar_transmission,
            legendre_coefficient_tangent,
            solar_transmission_tangent,
        };

        for ray_index in 0..self.num_rays() {
            let layer_start = self.ray_layer_offsets[ray_index] as usize;
            let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
            let transport_offset = self.ray_transport_value_offsets[ray_index] as usize;
            let mut current_attenuation = 1.0;
            let mut current_attenuation_tangent = 0.0;
            let mut ray_first_order_tangent = 0.0;
            let mut shared_endpoint_source = None;

            for layer_index in (layer_start..layer_end).rev() {
                let atmosphere_start = self.layer_atmosphere_offsets[layer_index] as usize;
                let atmosphere_end = self.layer_atmosphere_offsets[layer_index + 1] as usize;
                let mut optical_depth = 0.0;
                let mut optical_depth_tangent = 0.0;
                let mut albedo = 0.0;
                let mut albedo_tangent = 0.0;
                for stencil_index in atmosphere_start..atmosphere_end {
                    let atmosphere_index = self.atmosphere_indices[stencil_index] as usize;
                    let optical_depth_weight = self.optical_depth_weights[stencil_index];
                    let albedo_weight = self.albedo_weights[stencil_index];
                    optical_depth += optical_depth_weight * extinction[atmosphere_index];
                    optical_depth_tangent +=
                        optical_depth_weight * extinction_tangent[atmosphere_index];
                    albedo += albedo_weight * single_scatter_albedo[atmosphere_index];
                    albedo_tangent +=
                        albedo_weight * single_scatter_albedo_tangent[atmosphere_index];
                }
                let attenuation = (-optical_depth).exp();
                let attenuation_tangent = -attenuation * optical_depth_tangent;
                let source_factor = albedo * (1.0 - attenuation) * current_attenuation;
                let source_factor_tangent =
                    albedo_tangent * (1.0 - attenuation) * current_attenuation
                        - albedo * attenuation_tangent * current_attenuation
                        + albedo * (1.0 - attenuation) * current_attenuation_tangent;
                let source_start = self.layer_source_offsets[layer_index] as usize;
                let source_end = self.layer_source_offsets[layer_index + 1] as usize;
                for source_index in source_start..source_end {
                    let value_index =
                        transport_offset + self.source_value_inner_indices[source_index] as usize;
                    let weight = self.source_weights[source_index];
                    transport_values[value_index] += weight * source_factor;
                    transport_value_tangent[value_index] += weight * source_factor_tangent;
                }

                let start_source = shared_endpoint_source.unwrap_or_else(|| {
                    self.endpoint_source_jvp(
                        ray_index,
                        layer_index,
                        true,
                        &forcing,
                        extinction,
                        single_scatter_albedo,
                        extinction_tangent,
                        single_scatter_albedo_tangent,
                    )
                });
                let end_source = self.endpoint_source_jvp(
                    ray_index,
                    layer_index,
                    false,
                    &forcing,
                    extinction,
                    single_scatter_albedo,
                    extinction_tangent,
                    single_scatter_albedo_tangent,
                );
                if self.layer_distance[layer_index] >= 1.0e-4 {
                    let integration_factor = constant_source_factor(optical_depth, attenuation);
                    let integration_factor_tangent =
                        constant_source_factor_derivative(optical_depth, attenuation)
                            * optical_depth_tangent;
                    let endpoint_quadrature = start_source.value
                        * self.layer_start_fraction[layer_index]
                        + end_source.value * (1.0 - self.layer_start_fraction[layer_index]);
                    let endpoint_quadrature_tangent = start_source.tangent
                        * self.layer_start_fraction[layer_index]
                        + end_source.tangent * (1.0 - self.layer_start_fraction[layer_index]);
                    let distance = self.layer_distance[layer_index];
                    ray_first_order_tangent += distance
                        * (current_attenuation_tangent * integration_factor * endpoint_quadrature
                            + current_attenuation
                                * (integration_factor_tangent * endpoint_quadrature
                                    + integration_factor * endpoint_quadrature_tangent));
                }
                shared_endpoint_source = Some(end_source);

                current_attenuation_tangent = current_attenuation_tangent * attenuation
                    + current_attenuation * attenuation_tangent;
                current_attenuation *= attenuation;
            }

            let ground_start = self.ray_ground_offsets[ray_index] as usize;
            let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
            for ground_index in ground_start..ground_end {
                let value_index =
                    transport_offset + self.ground_value_inner_indices[ground_index] as usize;
                let weight = self.ground_weights[ground_index];
                transport_values[value_index] += weight * current_attenuation;
                transport_value_tangent[value_index] += weight * current_attenuation_tangent;
            }
            first_order_radiance_tangent[ray_index] = ray_first_order_tangent;
            ray_end_attenuation[ray_index] = current_attenuation;
            ray_end_attenuation_tangent[ray_index] = current_attenuation_tangent;
        }
        Ok(())
    }

    /// Reverse product for the fused scalar transport and atmospheric
    /// first-order forcing. The end-of-ray source is still evaluated in C++;
    /// its cotangent with respect to total observer attenuation is supplied in
    /// `ray_end_attenuation_cotangent`.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_vjp_with_first_order(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        legendre_coefficients: &[f64],
        maximum_order: &[i32],
        solar_transmission: &[f64],
        transport_value_gradient: &[f64],
        transport_column_indices: &[i32],
        solution: &[f64],
        first_order_radiance_gradient: &[f64],
        ray_end_attenuation_cotangent: &[f64],
        layer_optical_depth: &[f64],
        layer_attenuation: &[f64],
        layer_prefix_attenuation: &[f64],
        extinction_gradient: &mut [f64],
        single_scatter_albedo_gradient: &mut [f64],
        legendre_coefficient_gradient: &mut [f64],
        solar_transmission_gradient: &mut [f64],
    ) -> Result<()> {
        if extinction.len() < self.required_atmosphere_len
            || single_scatter_albedo.len() != extinction.len()
            || extinction_gradient.len() != extinction.len()
            || single_scatter_albedo_gradient.len() != extinction.len()
        {
            bail!("scalar VJP atmosphere arrays have inconsistent lengths");
        }
        let required_coefficients = extinction
            .len()
            .checked_mul(self.num_phase_moments)
            .ok_or_else(|| anyhow!("phase coefficient size overflow"))?;
        if legendre_coefficients.len() != required_coefficients
            || legendre_coefficient_gradient.len() != required_coefficients
            || maximum_order.len() != extinction.len()
        {
            bail!("scalar VJP phase arrays have inconsistent lengths");
        }
        if maximum_order.iter().any(|&order| {
            order < 1 || usize::try_from(order).map_or(true, |value| value > self.num_phase_moments)
        }) {
            bail!(
                "maximum phase orders must lie in 1..={}",
                self.num_phase_moments
            );
        }
        let required_solar_transmission = if self.solar_transmission_on_atmosphere_grid {
            self.required_atmosphere_len
        } else {
            self.num_layers() + self.num_rays()
        };
        if solar_transmission.len() < required_solar_transmission
            || solar_transmission_gradient.len() != solar_transmission.len()
        {
            bail!("scalar VJP solar-transmission arrays have inconsistent lengths");
        }
        let compact_transport_gradient = transport_value_gradient.is_empty();
        if !compact_transport_gradient
            && transport_value_gradient.len() < self.required_transport_len
        {
            bail!("scalar VJP transport gradient is too short");
        }
        if compact_transport_gradient {
            if transport_column_indices.len() < self.required_transport_len {
                bail!("scalar compact VJP transport columns are too short");
            }
            for &column in &transport_column_indices[..self.required_transport_len] {
                if column < 0 || column as usize >= solution.len() {
                    bail!("scalar compact VJP transport column is out of range");
                }
            }
        }
        if first_order_radiance_gradient.len() != self.num_rays()
            || ray_end_attenuation_cotangent.len() != self.num_rays()
        {
            bail!("scalar VJP ray gradient arrays have inconsistent lengths");
        }
        if layer_optical_depth.len() != self.num_layers()
            || layer_attenuation.len() != self.num_layers()
            || layer_prefix_attenuation.len() != self.num_layers()
        {
            bail!("scalar VJP layer scratch arrays have inconsistent lengths");
        }

        extinction_gradient.fill(0.0);
        single_scatter_albedo_gradient.fill(0.0);
        legendre_coefficient_gradient.fill(0.0);
        solar_transmission_gradient.fill(0.0);
        let forcing = FirstOrderInputs {
            legendre_coefficients,
            maximum_order,
            solar_transmission,
        };

        for ray_index in 0..self.num_rays() {
            let layer_start = self.ray_layer_offsets[ray_index] as usize;
            let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
            let transport_offset = self.ray_transport_value_offsets[ray_index] as usize;
            let forcing_cotangent = first_order_radiance_gradient[ray_index];
            let mut current_attenuation_cotangent = ray_end_attenuation_cotangent[ray_index];

            let ground_start = self.ray_ground_offsets[ray_index] as usize;
            let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
            for ground_index in ground_start..ground_end {
                let value_index =
                    transport_offset + self.ground_value_inner_indices[ground_index] as usize;
                let value_gradient = if compact_transport_gradient {
                    forcing_cotangent * solution[transport_column_indices[value_index] as usize]
                } else {
                    transport_value_gradient[value_index]
                };
                current_attenuation_cotangent += self.ground_weights[ground_index] * value_gradient;
            }

            // Reverse the observer-to-far-end forward traversal.
            for layer_index in layer_start..layer_end {
                let optical_depth = layer_optical_depth[layer_index];
                let attenuation = layer_attenuation[layer_index];
                let current_attenuation = layer_prefix_attenuation[layer_index];
                let atmosphere_start = self.layer_atmosphere_offsets[layer_index] as usize;
                let atmosphere_end = self.layer_atmosphere_offsets[layer_index + 1] as usize;

                let mut albedo = 0.0;
                for stencil_index in atmosphere_start..atmosphere_end {
                    let atmosphere_index = self.atmosphere_indices[stencil_index] as usize;
                    let albedo_weight = self.albedo_weights[stencil_index];
                    albedo += albedo_weight * single_scatter_albedo[atmosphere_index];
                }

                let mut source_factor_cotangent = 0.0;
                let source_start = self.layer_source_offsets[layer_index] as usize;
                let source_end = self.layer_source_offsets[layer_index + 1] as usize;
                for source_index in source_start..source_end {
                    let value_index =
                        transport_offset + self.source_value_inner_indices[source_index] as usize;
                    let value_gradient = if compact_transport_gradient {
                        forcing_cotangent * solution[transport_column_indices[value_index] as usize]
                    } else {
                        transport_value_gradient[value_index]
                    };
                    source_factor_cotangent += self.source_weights[source_index] * value_gradient;
                }
                let albedo_cotangent =
                    source_factor_cotangent * (1.0 - attenuation) * current_attenuation;
                let mut attenuation_cotangent =
                    -source_factor_cotangent * albedo * current_attenuation;
                let mut prefix_cotangent = source_factor_cotangent * albedo * (1.0 - attenuation);
                let mut optical_depth_cotangent = 0.0;

                if self.layer_distance[layer_index] >= 1.0e-4 && forcing_cotangent != 0.0 {
                    let start_endpoint = self.endpoint_primal(
                        ray_index,
                        layer_index,
                        true,
                        &forcing,
                        extinction,
                        single_scatter_albedo,
                    );
                    let end_endpoint = self.endpoint_primal(
                        ray_index,
                        layer_index,
                        false,
                        &forcing,
                        extinction,
                        single_scatter_albedo,
                    );
                    let endpoint_quadrature = start_endpoint.source
                        * self.layer_start_fraction[layer_index]
                        + end_endpoint.source * (1.0 - self.layer_start_fraction[layer_index]);
                    let integration_factor = constant_source_factor(optical_depth, attenuation);
                    let distance = self.layer_distance[layer_index];
                    prefix_cotangent +=
                        forcing_cotangent * integration_factor * endpoint_quadrature * distance;
                    optical_depth_cotangent += forcing_cotangent
                        * current_attenuation
                        * endpoint_quadrature
                        * distance
                        * constant_source_factor_derivative(optical_depth, attenuation);
                    let quadrature_cotangent =
                        forcing_cotangent * current_attenuation * integration_factor * distance;
                    self.accumulate_endpoint_vjp(
                        ray_index,
                        layer_index,
                        true,
                        &forcing,
                        start_endpoint,
                        quadrature_cotangent * self.layer_start_fraction[layer_index],
                        extinction_gradient,
                        single_scatter_albedo_gradient,
                        legendre_coefficient_gradient,
                        solar_transmission_gradient,
                    );
                    self.accumulate_endpoint_vjp(
                        ray_index,
                        layer_index,
                        false,
                        &forcing,
                        end_endpoint,
                        quadrature_cotangent * (1.0 - self.layer_start_fraction[layer_index]),
                        extinction_gradient,
                        single_scatter_albedo_gradient,
                        legendre_coefficient_gradient,
                        solar_transmission_gradient,
                    );
                }

                // A_next = A * attenuation.
                prefix_cotangent += current_attenuation_cotangent * attenuation;
                attenuation_cotangent += current_attenuation_cotangent * current_attenuation;
                // attenuation = exp(-optical_depth).
                optical_depth_cotangent -= attenuation * attenuation_cotangent;

                for stencil_index in atmosphere_start..atmosphere_end {
                    let atmosphere_index = self.atmosphere_indices[stencil_index] as usize;
                    let optical_depth_weight = self.optical_depth_weights[stencil_index];
                    let albedo_weight = self.albedo_weights[stencil_index];
                    extinction_gradient[atmosphere_index] +=
                        optical_depth_weight * optical_depth_cotangent;
                    single_scatter_albedo_gradient[atmosphere_index] +=
                        albedo_weight * albedo_cotangent;
                }
                current_attenuation_cotangent = prefix_cotangent;
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn validate_assembly_buffers(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        transport_values: &[f64],
        layer_optical_depth: &[f64],
        layer_attenuation: &[f64],
        layer_prefix_attenuation: &[f64],
        ray_end_attenuation: &[f64],
    ) -> Result<()> {
        if extinction.len() < self.required_atmosphere_len
            || single_scatter_albedo.len() < self.required_atmosphere_len
        {
            bail!(
                "atmosphere has {} extinction and {} albedo values, but ray stencils require {}",
                extinction.len(),
                single_scatter_albedo.len(),
                self.required_atmosphere_len
            );
        }
        if transport_values.len() < self.required_transport_len {
            bail!(
                "transport value array has length {}, but ray stencils require {}",
                transport_values.len(),
                self.required_transport_len
            );
        }
        let no_layer_scratch = layer_optical_depth.is_empty()
            && layer_attenuation.is_empty()
            && layer_prefix_attenuation.is_empty();
        let full_layer_scratch = layer_optical_depth.len() == self.num_layers()
            && layer_attenuation.len() == self.num_layers()
            && layer_prefix_attenuation.len() == self.num_layers();
        if !no_layer_scratch && !full_layer_scratch {
            bail!(
                "layer scratch lengths ({}, {}, {}) must all be zero or match {} packed layers",
                layer_optical_depth.len(),
                layer_attenuation.len(),
                layer_prefix_attenuation.len(),
                self.num_layers()
            );
        }
        if ray_end_attenuation.len() != self.num_rays() {
            bail!(
                "ray-end scratch length {} does not match {} packed rays",
                ray_end_attenuation.len(),
                self.num_rays()
            );
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble_impl(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        transport_values: &mut [f64],
        layer_optical_depth: &mut [f64],
        layer_attenuation: &mut [f64],
        layer_prefix_attenuation: &mut [f64],
        ray_end_attenuation: &mut [f64],
        mut first_order: Option<(FirstOrderInputs<'_>, &mut [f64])>,
    ) {
        transport_values.fill(0.0);
        if let Some((_, radiance)) = first_order.as_mut() {
            radiance.fill(0.0);
        }
        let retain_layer_scratch = !layer_optical_depth.is_empty();

        for ray_index in 0..self.num_rays() {
            let layer_start = self.ray_layer_offsets[ray_index] as usize;
            let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
            let transport_offset = self.ray_transport_value_offsets[ray_index] as usize;
            let mut current_attenuation = 1.0;
            let mut ray_first_order = 0.0;
            let mut shared_endpoint_source = None;

            for layer_index in (layer_start..layer_end).rev() {
                let atmosphere_start = self.layer_atmosphere_offsets[layer_index] as usize;
                let atmosphere_end = self.layer_atmosphere_offsets[layer_index + 1] as usize;

                let mut optical_depth = 0.0;
                let mut albedo = 0.0;
                for stencil_index in atmosphere_start..atmosphere_end {
                    let atmosphere_index = self.atmosphere_indices[stencil_index] as usize;
                    let optical_depth_weight = self.optical_depth_weights[stencil_index];
                    let albedo_weight = self.albedo_weights[stencil_index];
                    optical_depth += optical_depth_weight * extinction[atmosphere_index];
                    albedo += albedo_weight * single_scatter_albedo[atmosphere_index];
                }

                let attenuation = (-optical_depth).exp();
                if retain_layer_scratch {
                    layer_optical_depth[layer_index] = optical_depth;
                    layer_attenuation[layer_index] = attenuation;
                    layer_prefix_attenuation[layer_index] = current_attenuation;
                }

                let source_factor = albedo * (1.0 - attenuation) * current_attenuation;
                let source_start = self.layer_source_offsets[layer_index] as usize;
                let source_end = self.layer_source_offsets[layer_index + 1] as usize;
                for source_index in source_start..source_end {
                    let value_index =
                        transport_offset + self.source_value_inner_indices[source_index] as usize;
                    transport_values[value_index] +=
                        self.source_weights[source_index] * source_factor;
                }

                if let Some((forcing, _)) = first_order.as_ref() {
                    // Consecutive straight-ray layers share a boundary. In
                    // observer-to-far-end traversal, the previous exit is the
                    // next entrance, so only one new endpoint is evaluated per
                    // layer after the first.
                    let start_source = shared_endpoint_source.unwrap_or_else(|| {
                        self.endpoint_source(
                            ray_index,
                            layer_index,
                            true,
                            forcing,
                            extinction,
                            single_scatter_albedo,
                        )
                    });
                    let end_source = self.endpoint_source(
                        ray_index,
                        layer_index,
                        false,
                        forcing,
                        extinction,
                        single_scatter_albedo,
                    );
                    if self.layer_distance[layer_index] >= 1.0e-4 {
                        let integration_factor = constant_source_factor(optical_depth, attenuation);
                        ray_first_order += current_attenuation
                            * integration_factor
                            * (start_source * self.layer_start_fraction[layer_index]
                                + end_source * (1.0 - self.layer_start_fraction[layer_index]))
                            * self.layer_distance[layer_index];
                    }
                    shared_endpoint_source = Some(end_source);
                }

                current_attenuation *= attenuation;
            }

            let ground_start = self.ray_ground_offsets[ray_index] as usize;
            let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
            for ground_index in ground_start..ground_end {
                let value_index =
                    transport_offset + self.ground_value_inner_indices[ground_index] as usize;
                transport_values[value_index] +=
                    self.ground_weights[ground_index] * current_attenuation;
            }
            ray_end_attenuation[ray_index] = current_attenuation;
            if let Some((_, radiance)) = first_order.as_mut() {
                radiance[ray_index] = ray_first_order;
            }
        }
    }

    fn endpoint_source(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
        forcing: &FirstOrderInputs<'_>,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
    ) -> f64 {
        self.endpoint_primal(
            ray_index,
            layer_index,
            entrance,
            forcing,
            extinction,
            single_scatter_albedo,
        )
        .source
    }

    fn endpoint_primal(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
        forcing: &FirstOrderInputs<'_>,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
    ) -> EndpointPrimal {
        let atmosphere_start = self.layer_atmosphere_offsets[layer_index] as usize;
        let atmosphere_end = self.layer_atmosphere_offsets[layer_index + 1] as usize;
        let endpoint_weights = if entrance {
            &self.entrance_weights
        } else {
            &self.exit_weights
        };
        let phase_basis_start = ray_index * self.num_phase_moments;
        let phase_basis =
            &self.ray_phase_basis[phase_basis_start..phase_basis_start + self.num_phase_moments];

        let mut endpoint_extinction = 0.0;
        let mut endpoint_albedo = 0.0;
        let mut endpoint_phase = 0.0;
        let mut endpoint_solar_transmission = 0.0;
        for (&atmosphere_index, &weight) in self.atmosphere_indices
            [atmosphere_start..atmosphere_end]
            .iter()
            .zip(&endpoint_weights[atmosphere_start..atmosphere_end])
        {
            let atmosphere_index = atmosphere_index as usize;
            if weight == 0.0 {
                continue;
            }
            endpoint_extinction += weight * extinction[atmosphere_index];
            endpoint_albedo += weight * single_scatter_albedo[atmosphere_index];
            let coefficient_start = atmosphere_index * self.num_phase_moments;
            let order = forcing.maximum_order[atmosphere_index] as usize;
            endpoint_phase += weight
                * forcing.legendre_coefficients[coefficient_start..coefficient_start + order]
                    .iter()
                    .zip(&phase_basis[..order])
                    .map(|(coefficient, basis)| coefficient * basis)
                    .sum::<f64>();
            if self.solar_transmission_on_atmosphere_grid {
                endpoint_solar_transmission +=
                    weight * forcing.solar_transmission[atmosphere_index];
            }
        }
        if !self.solar_transmission_on_atmosphere_grid {
            // Every ray contributes one more endpoint than layer. Since flat
            // layer indices already include preceding rays' layers, adding the
            // ray index gives the matching endpoint offset.
            let exit_solar_index = layer_index + ray_index;
            endpoint_solar_transmission =
                forcing.solar_transmission[exit_solar_index + usize::from(entrance)];
        }

        let source =
            endpoint_extinction * endpoint_albedo * endpoint_solar_transmission * endpoint_phase
                / (4.0 * std::f64::consts::PI);
        EndpointPrimal {
            extinction: endpoint_extinction,
            albedo: endpoint_albedo,
            phase: endpoint_phase,
            solar_transmission: endpoint_solar_transmission,
            source,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn accumulate_endpoint_vjp(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
        forcing: &FirstOrderInputs<'_>,
        endpoint: EndpointPrimal,
        source_cotangent: f64,
        extinction_gradient: &mut [f64],
        single_scatter_albedo_gradient: &mut [f64],
        legendre_coefficient_gradient: &mut [f64],
        solar_transmission_gradient: &mut [f64],
    ) {
        if source_cotangent == 0.0 {
            return;
        }
        let scale = source_cotangent / (4.0 * std::f64::consts::PI);
        let extinction_cotangent =
            scale * endpoint.albedo * endpoint.solar_transmission * endpoint.phase;
        let albedo_cotangent =
            scale * endpoint.extinction * endpoint.solar_transmission * endpoint.phase;
        let solar_cotangent = scale * endpoint.extinction * endpoint.albedo * endpoint.phase;
        let phase_cotangent =
            scale * endpoint.extinction * endpoint.albedo * endpoint.solar_transmission;
        let atmosphere_start = self.layer_atmosphere_offsets[layer_index] as usize;
        let atmosphere_end = self.layer_atmosphere_offsets[layer_index + 1] as usize;
        let endpoint_weights = if entrance {
            &self.entrance_weights
        } else {
            &self.exit_weights
        };
        let phase_basis_start = ray_index * self.num_phase_moments;
        let phase_basis =
            &self.ray_phase_basis[phase_basis_start..phase_basis_start + self.num_phase_moments];
        for (&atmosphere_index, &weight) in self.atmosphere_indices
            [atmosphere_start..atmosphere_end]
            .iter()
            .zip(&endpoint_weights[atmosphere_start..atmosphere_end])
        {
            let atmosphere_index = atmosphere_index as usize;
            if weight == 0.0 {
                continue;
            }
            extinction_gradient[atmosphere_index] += weight * extinction_cotangent;
            single_scatter_albedo_gradient[atmosphere_index] += weight * albedo_cotangent;
            let coefficient_start = atmosphere_index * self.num_phase_moments;
            let order = forcing.maximum_order[atmosphere_index] as usize;
            for degree in 0..order {
                legendre_coefficient_gradient[coefficient_start + degree] +=
                    weight * phase_cotangent * phase_basis[degree];
            }
            if self.solar_transmission_on_atmosphere_grid {
                solar_transmission_gradient[atmosphere_index] += weight * solar_cotangent;
            }
        }
        if !self.solar_transmission_on_atmosphere_grid {
            let exit_solar_index = layer_index + ray_index;
            solar_transmission_gradient[exit_solar_index + usize::from(entrance)] +=
                solar_cotangent;
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn endpoint_source_jvp(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
        forcing: &FirstOrderJvpInputs<'_>,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        extinction_tangent: &[f64],
        single_scatter_albedo_tangent: &[f64],
    ) -> ValueTangent {
        let atmosphere_start = self.layer_atmosphere_offsets[layer_index] as usize;
        let atmosphere_end = self.layer_atmosphere_offsets[layer_index + 1] as usize;
        let endpoint_weights = if entrance {
            &self.entrance_weights
        } else {
            &self.exit_weights
        };
        let phase_basis_start = ray_index * self.num_phase_moments;
        let phase_basis =
            &self.ray_phase_basis[phase_basis_start..phase_basis_start + self.num_phase_moments];
        let mut endpoint_extinction = ValueTangent::default();
        let mut endpoint_albedo = ValueTangent::default();
        let mut endpoint_phase = ValueTangent::default();
        let mut endpoint_solar_transmission = ValueTangent::default();
        for (&atmosphere_index, &weight) in self.atmosphere_indices
            [atmosphere_start..atmosphere_end]
            .iter()
            .zip(&endpoint_weights[atmosphere_start..atmosphere_end])
        {
            let atmosphere_index = atmosphere_index as usize;
            if weight == 0.0 {
                continue;
            }
            endpoint_extinction.value += weight * extinction[atmosphere_index];
            endpoint_extinction.tangent += weight * extinction_tangent[atmosphere_index];
            endpoint_albedo.value += weight * single_scatter_albedo[atmosphere_index];
            endpoint_albedo.tangent += weight * single_scatter_albedo_tangent[atmosphere_index];
            let coefficient_start = atmosphere_index * self.num_phase_moments;
            let order = forcing.maximum_order[atmosphere_index] as usize;
            endpoint_phase.value += weight
                * forcing.legendre_coefficients[coefficient_start..coefficient_start + order]
                    .iter()
                    .zip(&phase_basis[..order])
                    .map(|(coefficient, basis)| coefficient * basis)
                    .sum::<f64>();
            endpoint_phase.tangent += weight
                * forcing.legendre_coefficient_tangent
                    [coefficient_start..coefficient_start + order]
                    .iter()
                    .zip(&phase_basis[..order])
                    .map(|(coefficient, basis)| coefficient * basis)
                    .sum::<f64>();
            if self.solar_transmission_on_atmosphere_grid {
                endpoint_solar_transmission.value +=
                    weight * forcing.solar_transmission[atmosphere_index];
                endpoint_solar_transmission.tangent +=
                    weight * forcing.solar_transmission_tangent[atmosphere_index];
            }
        }
        if !self.solar_transmission_on_atmosphere_grid {
            let exit_solar_index = layer_index + ray_index;
            let solar_index = exit_solar_index + usize::from(entrance);
            endpoint_solar_transmission.value = forcing.solar_transmission[solar_index];
            endpoint_solar_transmission.tangent = forcing.solar_transmission_tangent[solar_index];
        }
        let scale = 1.0 / (4.0 * std::f64::consts::PI);
        let value = endpoint_extinction.value
            * endpoint_albedo.value
            * endpoint_solar_transmission.value
            * endpoint_phase.value
            * scale;
        let tangent = (((endpoint_extinction.tangent * endpoint_albedo.value
            + endpoint_extinction.value * endpoint_albedo.tangent)
            * endpoint_solar_transmission.value
            + endpoint_extinction.value
                * endpoint_albedo.value
                * endpoint_solar_transmission.tangent)
            * endpoint_phase.value
            + endpoint_extinction.value
                * endpoint_albedo.value
                * endpoint_solar_transmission.value
                * endpoint_phase.tangent)
            * scale;
        ValueTangent { value, tangent }
    }
}

struct FirstOrderInputs<'a> {
    legendre_coefficients: &'a [f64],
    maximum_order: &'a [i32],
    solar_transmission: &'a [f64],
}

struct FirstOrderJvpInputs<'a> {
    legendre_coefficients: &'a [f64],
    maximum_order: &'a [i32],
    solar_transmission: &'a [f64],
    legendre_coefficient_tangent: &'a [f64],
    solar_transmission_tangent: &'a [f64],
}

#[derive(Clone, Copy, Default)]
struct EndpointPrimal {
    extinction: f64,
    albedo: f64,
    phase: f64,
    solar_transmission: f64,
    source: f64,
}

#[derive(Clone, Copy, Default)]
struct ValueTangent {
    value: f64,
    tangent: f64,
}

fn append_legendre_basis(cosine: f64, count: usize, output: &mut Vec<f64>) {
    output.push(1.0);
    if count == 1 {
        return;
    }
    output.push(cosine);
    for degree in 2..count {
        let degree_f64 = degree as f64;
        let previous = output[output.len() - 1];
        let previous_previous = output[output.len() - 2];
        output.push(
            ((2.0 * degree_f64 - 1.0) * cosine * previous - (degree_f64 - 1.0) * previous_previous)
                / degree_f64,
        );
    }
}

fn constant_source_factor(optical_depth: f64, attenuation: f64) -> f64 {
    if optical_depth.abs() < 1.0e-4 {
        let optical_depth_squared = optical_depth * optical_depth;
        let optical_depth_cubed = optical_depth_squared * optical_depth;
        let optical_depth_fourth = optical_depth_squared * optical_depth_squared;
        return 1.0 - optical_depth / 2.0 + optical_depth_squared / 6.0
            - optical_depth_cubed / 24.0
            + optical_depth_fourth / 120.0;
    }
    (1.0 - attenuation) / optical_depth
}

fn constant_source_factor_derivative(optical_depth: f64, attenuation: f64) -> f64 {
    if optical_depth.abs() < 1.0e-4 {
        let optical_depth_squared = optical_depth * optical_depth;
        let optical_depth_cubed = optical_depth_squared * optical_depth;
        return -0.5 + optical_depth / 3.0 - optical_depth_squared / 8.0
            + optical_depth_cubed / 30.0;
    }
    let value = (1.0 - attenuation) / optical_depth;
    1.0 / optical_depth - value * (1.0 + 1.0 / optical_depth)
}

fn validate_offsets(name: &str, offsets: &[u32]) -> Result<()> {
    if offsets.is_empty() {
        bail!("{name} offsets are empty");
    }
    if offsets[0] != 0 {
        bail!("{name} offsets must start at zero");
    }
    if offsets.windows(2).any(|pair| pair[0] > pair[1]) {
        bail!("{name} offsets are not monotonic");
    }
    Ok(())
}

fn checked_last(offsets: &[u32]) -> Result<usize> {
    offsets
        .last()
        .copied()
        .map(|value| value as usize)
        .ok_or_else(|| anyhow!("offset array is empty"))
}

#[cfg(test)]
mod tests {
    use super::ScalarRayTransport;

    fn geometry() -> ScalarRayTransport {
        // Ray zero has layers 0 and 1; ray one has layer 2.  Layers are stored
        // from the far end toward the observer and traversed in reverse.
        ScalarRayTransport::new(
            &[0, 2, 3],
            &[0, 1, 3, 4],
            &[0, 0, 1, 1],
            &[0.2, 0.1, 0.3, 0.4],
            &[1.0, 0.25, 0.75, 1.0],
            &[1.0, 0.25, 0.75, 1.0],
            &[1.0, 0.25, 0.75, 1.0],
            &[0.2, 0.4, 0.4],
            &[0.5, 0.5, 0.5],
            &[0.5, 0.5, 0.5],
            &[0.2, -0.4],
            3,
            false,
            &[0, 1, 3, 4],
            &[0, 2],
            &[0, 0, 1, 0],
            &[1.0, 0.4, 0.6, 1.0],
            &[0, 1, 1],
            &[1],
            &[0.5],
        )
        .unwrap()
    }

    #[test]
    fn assembles_transport_and_reusable_attenuation() {
        let geometry = geometry();
        let extinction = [2.0, 0.5];
        let albedo = [0.8, 0.4];
        let mut values = [7.0; 3];
        let mut optical_depth = [0.0; 3];
        let mut attenuation = [0.0; 3];
        let mut prefix = [0.0; 3];
        let mut end = [0.0; 2];

        geometry
            .assemble(
                &extinction,
                &albedo,
                &mut values,
                &mut optical_depth,
                &mut attenuation,
                &mut prefix,
                &mut end,
            )
            .unwrap();

        let od0: f64 = 0.2 * 2.0;
        let od1: f64 = 0.1 * 2.0 + 0.3 * 0.5;
        let od2: f64 = 0.4 * 0.5;
        let a0 = (-od0).exp();
        let a1 = (-od1).exp();
        let a2 = (-od2).exp();
        let w0 = 0.8;
        let w1 = 0.25 * 0.8 + 0.75 * 0.4;
        let w2 = 0.4;

        for (actual, expected) in optical_depth.iter().zip([od0, od1, od2]) {
            assert!((actual - expected).abs() < 1.0e-14);
        }
        for (actual, expected) in attenuation.iter().zip([a0, a1, a2]) {
            assert!((actual - expected).abs() < 1.0e-14);
        }
        for (actual, expected) in prefix.iter().zip([a1, 1.0, 1.0]) {
            assert!((actual - expected).abs() < 1.0e-14);
        }
        for (actual, expected) in end.iter().zip([a1 * a0, a2]) {
            assert!((actual - expected).abs() < 1.0e-14);
        }

        let expected0 = w0 * (1.0 - a0) * a1 + 0.4 * w1 * (1.0 - a1);
        let expected1 = 0.6 * w1 * (1.0 - a1) + 0.5 * a1 * a0;
        let expected2 = w2 * (1.0 - a2);
        assert!((values[0] - expected0).abs() < 1.0e-14);
        assert!((values[1] - expected1).abs() < 1.0e-14);
        assert!((values[2] - expected2).abs() < 1.0e-14);

        let expected_values = values;
        geometry
            .assemble(
                &extinction,
                &albedo,
                &mut values,
                &mut [],
                &mut [],
                &mut [],
                &mut end,
            )
            .unwrap();
        assert_eq!(values, expected_values);
    }

    #[test]
    fn validates_runtime_buffer_sizes() {
        let geometry = geometry();
        let mut values = [0.0; 3];
        let mut layer = [0.0; 3];
        let mut end = [0.0; 2];
        let error = geometry
            .assemble(
                &[1.0],
                &[1.0],
                &mut values,
                &mut layer.clone(),
                &mut layer.clone(),
                &mut layer,
                &mut end,
            )
            .unwrap_err();
        assert!(error.to_string().contains("ray stencils require 2"));
    }

    #[test]
    fn fuses_scalar_first_order_forcing() {
        let geometry = ScalarRayTransport::new(
            &[0, 1],
            &[0, 2],
            &[0, 1],
            &[0.2, 0.3],
            &[0.5, 0.5],
            &[0.25, 0.75],
            &[0.75, 0.25],
            &[0.5],
            &[0.4],
            &[0.6],
            &[0.3],
            3,
            false,
            &[0, 1],
            &[0],
            &[0],
            &[1.0],
            &[0, 0],
            &[],
            &[],
        )
        .unwrap();
        let extinction = [2.0, 1.0];
        let albedo = [0.8, 0.6];
        let coefficients = [1.0, 0.2, 0.1, 1.0, 0.4, 0.0];
        let maximum_order = [3, 2];
        let solar_transmission = [0.8, 0.9];
        let mut transport = [0.0];
        let mut first_order = [0.0];
        let mut optical_depth = [0.0];
        let mut attenuation = [0.0];
        let mut prefix = [0.0];
        let mut end = [0.0];

        geometry
            .assemble_with_first_order(
                &extinction,
                &albedo,
                &coefficients,
                &maximum_order,
                &solar_transmission,
                &mut transport,
                &mut first_order,
                &mut optical_depth,
                &mut attenuation,
                &mut prefix,
                &mut end,
            )
            .unwrap();

        let phase_0 = 1.0 + 0.2 * 0.3 + 0.1 * -0.365;
        let phase_1 = 1.0 + 0.4 * 0.3;
        let start = (0.25 * 2.0 + 0.75 * 1.0)
            * (0.25 * 0.8 + 0.75 * 0.6)
            * 0.9
            * (0.25 * phase_0 + 0.75 * phase_1)
            / (4.0 * std::f64::consts::PI);
        let finish = (0.75 * 2.0 + 0.25 * 1.0)
            * (0.75 * 0.8 + 0.25 * 0.6)
            * 0.8
            * (0.75 * phase_0 + 0.25 * phase_1)
            / (4.0 * std::f64::consts::PI);
        let od: f64 = 0.2 * 2.0 + 0.3 * 1.0;
        let layer_attenuation = (-od).exp();
        let integration_factor = (1.0 - layer_attenuation) / od;
        let expected = integration_factor * (0.4 * start + 0.6 * finish) * 0.5;
        assert!((first_order[0] - expected).abs() < 1.0e-14);

        let extinction_tangent = [0.1, -0.2];
        let albedo_tangent = [-0.03, 0.04];
        let coefficient_tangent = [0.02, -0.01, 0.03, -0.02, 0.01, 0.0];
        let solar_tangent = [-0.04, 0.05];
        let base_transport = transport;
        let base_first_order = first_order;
        let base_end = end;
        let mut transport_tangent = [0.0];
        let mut first_order_tangent = [0.0];
        let mut end_tangent = [0.0];
        geometry
            .assemble_jvp_with_first_order(
                &extinction,
                &albedo,
                &coefficients,
                &maximum_order,
                &solar_transmission,
                &extinction_tangent,
                &albedo_tangent,
                &coefficient_tangent,
                &solar_tangent,
                &mut transport,
                &mut transport_tangent,
                &mut first_order_tangent,
                &mut end,
                &mut end_tangent,
            )
            .unwrap();

        let epsilon = 1.0e-6;
        let perturbed_extinction: [f64; 2] =
            std::array::from_fn(|index| extinction[index] + epsilon * extinction_tangent[index]);
        let perturbed_albedo: [f64; 2] =
            std::array::from_fn(|index| albedo[index] + epsilon * albedo_tangent[index]);
        let perturbed_coefficients: [f64; 6] =
            std::array::from_fn(|index| coefficients[index] + epsilon * coefficient_tangent[index]);
        let perturbed_solar: [f64; 2] =
            std::array::from_fn(|index| solar_transmission[index] + epsilon * solar_tangent[index]);
        let mut perturbed_transport = [0.0];
        let mut perturbed_first_order = [0.0];
        let mut perturbed_optical_depth = [0.0];
        let mut perturbed_attenuation = [0.0];
        let mut perturbed_prefix = [0.0];
        let mut perturbed_end = [0.0];
        geometry
            .assemble_with_first_order(
                &perturbed_extinction,
                &perturbed_albedo,
                &perturbed_coefficients,
                &maximum_order,
                &perturbed_solar,
                &mut perturbed_transport,
                &mut perturbed_first_order,
                &mut perturbed_optical_depth,
                &mut perturbed_attenuation,
                &mut perturbed_prefix,
                &mut perturbed_end,
            )
            .unwrap();
        assert!(
            (transport_tangent[0] - (perturbed_transport[0] - base_transport[0]) / epsilon).abs()
                < 1.0e-7
        );
        assert!(
            (first_order_tangent[0] - (perturbed_first_order[0] - base_first_order[0]) / epsilon)
                .abs()
                < 1.0e-7
        );
        assert!((end_tangent[0] - (perturbed_end[0] - base_end[0]) / epsilon).abs() < 1.0e-7);

        let transport_cotangent = [0.75];
        let first_order_cotangent = [-0.5];
        let end_cotangent = [0.2];
        let mut extinction_gradient = [0.0; 2];
        let mut albedo_gradient = [0.0; 2];
        let mut coefficient_gradient = [0.0; 6];
        let mut solar_gradient = [0.0; 2];
        geometry
            .assemble_vjp_with_first_order(
                &extinction,
                &albedo,
                &coefficients,
                &maximum_order,
                &solar_transmission,
                &transport_cotangent,
                &[],
                &[],
                &first_order_cotangent,
                &end_cotangent,
                &optical_depth,
                &attenuation,
                &prefix,
                &mut extinction_gradient,
                &mut albedo_gradient,
                &mut coefficient_gradient,
                &mut solar_gradient,
            )
            .unwrap();
        let mut compact_extinction_gradient = [0.0; 2];
        let mut compact_albedo_gradient = [0.0; 2];
        let mut compact_coefficient_gradient = [0.0; 6];
        let mut compact_solar_gradient = [0.0; 2];
        geometry
            .assemble_vjp_with_first_order(
                &extinction,
                &albedo,
                &coefficients,
                &maximum_order,
                &solar_transmission,
                &[],
                &[0],
                &[transport_cotangent[0] / first_order_cotangent[0]],
                &first_order_cotangent,
                &end_cotangent,
                &optical_depth,
                &attenuation,
                &prefix,
                &mut compact_extinction_gradient,
                &mut compact_albedo_gradient,
                &mut compact_coefficient_gradient,
                &mut compact_solar_gradient,
            )
            .unwrap();
        assert_eq!(compact_extinction_gradient, extinction_gradient);
        assert_eq!(compact_albedo_gradient, albedo_gradient);
        assert_eq!(compact_coefficient_gradient, coefficient_gradient);
        assert_eq!(compact_solar_gradient, solar_gradient);
        let output_dot = transport_tangent[0] * transport_cotangent[0]
            + first_order_tangent[0] * first_order_cotangent[0]
            + end_tangent[0] * end_cotangent[0];
        let input_dot = extinction_gradient
            .iter()
            .zip(extinction_tangent)
            .map(|(gradient, tangent)| gradient * tangent)
            .sum::<f64>()
            + albedo_gradient
                .iter()
                .zip(albedo_tangent)
                .map(|(gradient, tangent)| gradient * tangent)
                .sum::<f64>()
            + coefficient_gradient
                .iter()
                .zip(coefficient_tangent)
                .map(|(gradient, tangent)| gradient * tangent)
                .sum::<f64>()
            + solar_gradient
                .iter()
                .zip(solar_tangent)
                .map(|(gradient, tangent)| gradient * tangent)
                .sum::<f64>();
        assert!((input_dot - output_dot).abs() < 1.0e-13);
    }
}
