//! Transport assembly for traced successive-orders rays.
//!
//! The ray tracer remains responsible for reducing geometry to interpolation
//! stencils.  This module owns the wavelength-independent packed stencils and
//! applies wavelength-dependent extinction and single-scatter albedo directly
//! to the successive-orders CSR value array.

use std::collections::HashMap;

use anyhow::{Result, anyhow, bail};

use crate::math::wigner::WignerDCalculator;

#[derive(Debug, Hash, PartialEq, Eq)]
enum EndpointStencilKey {
    Inline {
        len: u8,
        atmosphere_indices: [u32; 4],
        weight_bits: [u64; 4],
    },
    Heap(Box<[(u32, u64)]>),
}

impl EndpointStencilKey {
    fn new(atmosphere_indices: &[u32], weights: &[f64]) -> Self {
        debug_assert_eq!(atmosphere_indices.len(), weights.len());
        if atmosphere_indices.len() <= 4 {
            let mut inline_indices = [0; 4];
            let mut inline_weights = [0; 4];
            inline_indices[..atmosphere_indices.len()].copy_from_slice(atmosphere_indices);
            for (output, weight) in inline_weights.iter_mut().zip(weights) {
                *output = weight.to_bits();
            }
            Self::Inline {
                len: atmosphere_indices.len() as u8,
                atmosphere_indices: inline_indices,
                weight_bits: inline_weights,
            }
        } else {
            Self::Heap(
                atmosphere_indices
                    .iter()
                    .copied()
                    .zip(weights.iter().map(|weight| weight.to_bits()))
                    .collect(),
            )
        }
    }
}

/// Wavelength-independent ray geometry compiled for transport assembly.
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
    layer_exit_endpoint_indices: Vec<u32>,
    ray_terminal_endpoint_indices: Vec<u32>,
    layer_exit_phase_endpoint_slots: Vec<u16>,
    ray_terminal_phase_endpoint_slots: Vec<u16>,
    endpoint_atmosphere_offsets: Vec<u32>,
    endpoint_atmosphere_indices: Vec<u32>,
    endpoint_weights: Vec<f64>,
    layer_distance: Vec<f64>,
    layer_start_fraction: Vec<f64>,
    phase_geometry_ray_offsets: Vec<u32>,
    phase_geometry_ray_indices: Vec<u32>,
    phase_geometry_endpoint_offsets: Vec<u32>,
    phase_geometry_endpoint_indices: Vec<u32>,
    phase_endpoint_cache_size: usize,
    phase_geometry_endpoint_evaluations: Vec<usize>,
    ray_phase_basis: Vec<f64>,
    ray_polarized_phase_basis: Vec<f64>,
    phase_q_factor: Vec<f64>,
    phase_u_factor: Vec<f64>,
    num_phase_moments: usize,
    solar_transmission_on_atmosphere_grid: bool,
    layer_source_offsets: Vec<u32>,
    ray_transport_value_offsets: Vec<u32>,
    ray_transport_row_nnz: Vec<u32>,
    num_stokes: usize,
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
        ray_phase_q_factor: &[f64],
        ray_phase_u_factor: &[f64],
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
        if (!ray_phase_q_factor.is_empty() && ray_phase_q_factor.len() != num_rays)
            || (!ray_phase_u_factor.is_empty() && ray_phase_u_factor.len() != num_rays)
            || ray_phase_q_factor.len() != ray_phase_u_factor.len()
            || ray_phase_q_factor
                .iter()
                .chain(ray_phase_u_factor)
                .any(|factor| !factor.is_finite())
        {
            bail!("polarized ray phase factors have inconsistent lengths or values");
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
        let mut ray_transport_row_nnz = Vec::with_capacity(num_rays);
        for ray_index in 0..num_rays {
            let row_offset = ray_transport_value_offsets[ray_index] as usize;
            let layer_start = ray_layer_offsets[ray_index] as usize;
            let layer_end = ray_layer_offsets[ray_index + 1] as usize;
            let source_start = layer_source_offsets[layer_start] as usize;
            let source_end = layer_source_offsets[layer_end] as usize;
            let mut row_nnz = 0;
            for &inner_index in &source_value_inner_indices[source_start..source_end] {
                row_nnz = row_nnz.max(inner_index as usize + 1);
                required_transport_len =
                    required_transport_len.max(row_offset + inner_index as usize + 1);
            }
            let ground_start = ray_ground_offsets[ray_index] as usize;
            let ground_end = ray_ground_offsets[ray_index + 1] as usize;
            for &inner_index in &ground_value_inner_indices[ground_start..ground_end] {
                row_nnz = row_nnz.max(inner_index as usize + 1);
                required_transport_len =
                    required_transport_len.max(row_offset + inner_index as usize + 1);
            }
            ray_transport_row_nnz.push(u32::try_from(row_nnz)?);
        }
        let mut endpoint_map: HashMap<EndpointStencilKey, u32> = HashMap::new();
        let mut endpoint_atmosphere_offsets = vec![0];
        let mut endpoint_atmosphere_indices = Vec::new();
        let mut endpoint_weights = Vec::new();
        let mut layer_exit_endpoint_indices = Vec::with_capacity(num_layers);
        for layer_index in 0..num_layers {
            let atmosphere_start = layer_atmosphere_offsets[layer_index] as usize;
            let atmosphere_end = layer_atmosphere_offsets[layer_index + 1] as usize;
            layer_exit_endpoint_indices.push(intern_endpoint_stencil(
                &atmosphere_indices[atmosphere_start..atmosphere_end],
                &exit_weights[atmosphere_start..atmosphere_end],
                &mut endpoint_map,
                &mut endpoint_atmosphere_offsets,
                &mut endpoint_atmosphere_indices,
                &mut endpoint_weights,
            )?);
        }

        let mut ray_terminal_endpoint_indices = Vec::with_capacity(num_rays);
        for ray_index in 0..num_rays {
            let layer_end = ray_layer_offsets[ray_index + 1] as usize;
            if layer_end > ray_layer_offsets[ray_index] as usize {
                let terminal_layer = layer_end - 1;
                let atmosphere_start = layer_atmosphere_offsets[terminal_layer] as usize;
                let atmosphere_end = layer_atmosphere_offsets[terminal_layer + 1] as usize;
                ray_terminal_endpoint_indices.push(intern_endpoint_stencil(
                    &atmosphere_indices[atmosphere_start..atmosphere_end],
                    &entrance_weights[atmosphere_start..atmosphere_end],
                    &mut endpoint_map,
                    &mut endpoint_atmosphere_offsets,
                    &mut endpoint_atmosphere_indices,
                    &mut endpoint_weights,
                )?);
            } else {
                ray_terminal_endpoint_indices.push(u32::MAX);
            }
        }

        let polarized = !ray_phase_q_factor.is_empty();
        let mut phase_geometry_map = HashMap::new();
        let mut phase_geometry_rays: Vec<Vec<u32>> = Vec::new();
        let mut unique_scattering_cosines = Vec::new();
        let mut phase_q_factor = Vec::new();
        let mut phase_u_factor = Vec::new();
        for ray_index in 0..num_rays {
            let cosine = ray_scattering_cosine[ray_index].clamp(-1.0, 1.0);
            if !cosine.is_finite() {
                bail!("ray scattering cosine is not finite");
            }
            let q_factor = if polarized {
                ray_phase_q_factor[ray_index]
            } else {
                0.0
            };
            let u_factor = if polarized {
                ray_phase_u_factor[ray_index]
            } else {
                0.0
            };
            const PHASE_GEOMETRY_KEY_SCALE: f64 = 1.0e12;
            let key = [
                (cosine * PHASE_GEOMETRY_KEY_SCALE).round() as i64,
                (q_factor * PHASE_GEOMETRY_KEY_SCALE).round() as i64,
                (u_factor * PHASE_GEOMETRY_KEY_SCALE).round() as i64,
            ];
            let phase_geometry = match phase_geometry_map.get(&key).copied() {
                Some(index) => index,
                None => {
                    let index = unique_scattering_cosines.len();
                    phase_geometry_map.insert(key, index);
                    unique_scattering_cosines.push(cosine);
                    if polarized {
                        phase_q_factor.push(q_factor);
                        phase_u_factor.push(u_factor);
                    }
                    phase_geometry_rays.push(Vec::new());
                    index
                }
            };
            phase_geometry_rays[phase_geometry].push(u32::try_from(ray_index)?);
        }
        let mut phase_geometry_ray_offsets = Vec::with_capacity(phase_geometry_rays.len() + 1);
        let mut phase_geometry_ray_indices = Vec::with_capacity(num_rays);
        let mut phase_geometry_endpoint_offsets = Vec::with_capacity(phase_geometry_rays.len() + 1);
        let mut phase_geometry_endpoint_indices = Vec::new();
        let mut phase_geometry_endpoint_evaluations = Vec::with_capacity(phase_geometry_rays.len());
        let mut endpoint_seen = vec![false; endpoint_atmosphere_offsets.len() - 1];
        let mut endpoint_slot = vec![u32::MAX; endpoint_atmosphere_offsets.len() - 1];
        let mut layer_exit_phase_endpoint_slots_u32 = vec![u32::MAX; num_layers];
        let mut ray_terminal_phase_endpoint_slots_u32 = vec![u32::MAX; num_rays];
        let mut max_phase_geometry_endpoints = 0;
        phase_geometry_ray_offsets.push(0);
        phase_geometry_endpoint_offsets.push(0);
        for rays in phase_geometry_rays {
            let mut endpoint_evaluations = 0usize;
            let phase_endpoint_start = phase_geometry_endpoint_indices.len();
            endpoint_seen.fill(false);
            for &ray_index in &rays {
                let ray_index = ray_index as usize;
                let layer_start = ray_layer_offsets[ray_index] as usize;
                let layer_end = ray_layer_offsets[ray_index + 1] as usize;
                for layer_index in layer_start..layer_end {
                    let endpoint_index = layer_exit_endpoint_indices[layer_index] as usize;
                    if !endpoint_seen[endpoint_index] {
                        endpoint_seen[endpoint_index] = true;
                        endpoint_slot[endpoint_index] = u32::try_from(
                            phase_geometry_endpoint_indices.len() - phase_endpoint_start,
                        )?;
                        phase_geometry_endpoint_indices.push(endpoint_index as u32);
                    }
                    layer_exit_phase_endpoint_slots_u32[layer_index] =
                        endpoint_slot[endpoint_index];
                }
                endpoint_evaluations = endpoint_evaluations.saturating_add(
                    layer_atmosphere_offsets[layer_end] as usize
                        - layer_atmosphere_offsets[layer_start] as usize,
                );
                if layer_start < layer_end {
                    let endpoint_index = ray_terminal_endpoint_indices[ray_index] as usize;
                    if !endpoint_seen[endpoint_index] {
                        endpoint_seen[endpoint_index] = true;
                        endpoint_slot[endpoint_index] = u32::try_from(
                            phase_geometry_endpoint_indices.len() - phase_endpoint_start,
                        )?;
                        phase_geometry_endpoint_indices.push(endpoint_index as u32);
                    }
                    ray_terminal_phase_endpoint_slots_u32[ray_index] =
                        endpoint_slot[endpoint_index];
                    endpoint_evaluations = endpoint_evaluations.saturating_add(
                        layer_atmosphere_offsets[layer_end] as usize
                            - layer_atmosphere_offsets[layer_end - 1] as usize,
                    );
                }
            }
            phase_geometry_endpoint_evaluations.push(endpoint_evaluations);
            phase_geometry_ray_indices.extend(rays);
            phase_geometry_ray_offsets.push(u32::try_from(phase_geometry_ray_indices.len())?);
            phase_geometry_endpoint_offsets
                .push(u32::try_from(phase_geometry_endpoint_indices.len())?);
            max_phase_geometry_endpoints = max_phase_geometry_endpoints
                .max(phase_geometry_endpoint_indices.len() - phase_endpoint_start);
        }
        let compact_phase_endpoint_cache = endpoint_atmosphere_offsets.len() - 1
            > max_phase_geometry_endpoints.saturating_mul(4)
            && max_phase_geometry_endpoints < u16::MAX as usize;
        let phase_endpoint_cache_size = if compact_phase_endpoint_cache {
            max_phase_geometry_endpoints
        } else {
            endpoint_atmosphere_offsets.len() - 1
        };
        let layer_exit_phase_endpoint_slots = if compact_phase_endpoint_cache {
            layer_exit_phase_endpoint_slots_u32
                .into_iter()
                .map(|slot| slot as u16)
                .collect()
        } else {
            Vec::new()
        };
        let ray_terminal_phase_endpoint_slots = if compact_phase_endpoint_cache {
            ray_terminal_phase_endpoint_slots_u32
                .into_iter()
                .map(|slot| slot as u16)
                .collect()
        } else {
            Vec::new()
        };

        let num_phase_geometries = unique_scattering_cosines.len();
        let mut ray_phase_basis = Vec::with_capacity(num_phase_geometries * num_phase_moments);
        let mut ray_polarized_phase_basis = if polarized {
            Vec::with_capacity(num_phase_geometries * num_phase_moments)
        } else {
            Vec::new()
        };
        let polarized_calculator = polarized.then(|| WignerDCalculator::new(0, 2));
        for &cosine in &unique_scattering_cosines {
            append_legendre_basis(cosine, num_phase_moments, &mut ray_phase_basis);
            if let Some(calculator) = &polarized_calculator {
                let basis_start = ray_polarized_phase_basis.len();
                ray_polarized_phase_basis.resize(basis_start + num_phase_moments, 0.0);
                calculator.vector_d(
                    cosine.acos(),
                    &mut ray_polarized_phase_basis[basis_start..basis_start + num_phase_moments],
                );
            }
        }

        Ok(Self {
            ray_layer_offsets: ray_layer_offsets.to_vec(),
            layer_atmosphere_offsets: layer_atmosphere_offsets.to_vec(),
            atmosphere_indices: atmosphere_indices.to_vec(),
            optical_depth_weights: optical_depth_weights.to_vec(),
            albedo_weights: albedo_weights.to_vec(),
            layer_exit_endpoint_indices,
            ray_terminal_endpoint_indices,
            layer_exit_phase_endpoint_slots,
            ray_terminal_phase_endpoint_slots,
            endpoint_atmosphere_offsets,
            endpoint_atmosphere_indices,
            endpoint_weights,
            layer_distance: layer_distance.to_vec(),
            layer_start_fraction: layer_start_fraction.to_vec(),
            phase_geometry_ray_offsets,
            phase_geometry_ray_indices,
            phase_geometry_endpoint_offsets,
            phase_geometry_endpoint_indices,
            phase_endpoint_cache_size,
            phase_geometry_endpoint_evaluations,
            ray_phase_basis,
            ray_polarized_phase_basis,
            phase_q_factor,
            phase_u_factor,
            num_phase_moments,
            solar_transmission_on_atmosphere_grid,
            layer_source_offsets: layer_source_offsets.to_vec(),
            ray_transport_value_offsets: ray_transport_value_offsets.to_vec(),
            ray_transport_row_nnz,
            num_stokes: 1,
            source_value_inner_indices: source_value_inner_indices.to_vec(),
            source_weights: source_weights.to_vec(),
            ray_ground_offsets: ray_ground_offsets.to_vec(),
            ground_value_inner_indices: ground_value_inner_indices.to_vec(),
            ground_weights: ground_weights.to_vec(),
            required_atmosphere_len,
            required_transport_len,
        })
    }

    /// Configures the packed CSR row layout for a Stokes-vector transport.
    /// The transport attenuation is scalar, but each Stokes row receives the
    /// same interpolation weights at its own CSR value offset.
    pub fn with_stokes_layout(
        mut self,
        num_stokes: usize,
        ray_transport_row_nnz: &[u32],
    ) -> Result<Self> {
        if num_stokes == 0 {
            bail!("ray transport requires at least one Stokes component");
        }
        if ray_transport_row_nnz.len() != self.num_rays() {
            bail!(
                "ray transport row sizes have length {}, expected {}",
                ray_transport_row_nnz.len(),
                self.num_rays()
            );
        }

        let mut required_transport_len = 0;
        for (ray_index, &ray_row_nnz) in ray_transport_row_nnz.iter().enumerate() {
            let row_nnz = ray_row_nnz as usize;
            let layer_start = self.ray_layer_offsets[ray_index] as usize;
            let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
            let source_start = self.layer_source_offsets[layer_start] as usize;
            let source_end = self.layer_source_offsets[layer_end] as usize;
            let ground_start = self.ray_ground_offsets[ray_index] as usize;
            let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
            if self.source_value_inner_indices[source_start..source_end]
                .iter()
                .chain(&self.ground_value_inner_indices[ground_start..ground_end])
                .any(|&inner| inner as usize >= row_nnz)
            {
                bail!("ray transport inner index exceeds its CSR row size");
            }
            let ray_end = (self.ray_transport_value_offsets[ray_index] as usize)
                .checked_add(
                    row_nnz
                        .checked_mul(num_stokes)
                        .ok_or_else(|| anyhow!("ray transport size overflow"))?,
                )
                .ok_or_else(|| anyhow!("ray transport size overflow"))?;
            required_transport_len = required_transport_len.max(ray_end);
        }
        self.num_stokes = num_stokes;
        self.ray_transport_row_nnz = ray_transport_row_nnz.to_vec();
        self.required_transport_len = required_transport_len;
        if num_stokes == 3 && self.phase_q_factor.len() != self.num_phase_geometries() {
            bail!("polarized first-order geometry is missing phase rotation factors");
        }
        if num_stokes != 1 && num_stokes != 3 {
            bail!("first-order transport supports one or three Stokes components");
        }
        if num_stokes == 1 {
            self.ray_polarized_phase_basis.clear();
            self.ray_polarized_phase_basis.shrink_to_fit();
            self.phase_q_factor.clear();
            self.phase_q_factor.shrink_to_fit();
            self.phase_u_factor.clear();
            self.phase_u_factor.shrink_to_fit();
        }
        Ok(self)
    }

    pub fn num_rays(&self) -> usize {
        self.ray_layer_offsets.len() - 1
    }

    pub fn num_layers(&self) -> usize {
        self.layer_atmosphere_offsets.len() - 1
    }

    /// Number of CSR values required by this packed transport geometry.
    pub fn transport_value_size(&self) -> usize {
        self.required_transport_len
    }

    /// Number of Stokes rows assembled for each traced ray.
    pub fn num_stokes(&self) -> usize {
        self.num_stokes
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

    pub fn num_phase_geometries(&self) -> usize {
        self.phase_geometry_ray_offsets.len() - 1
    }

    /// Returns whether at least one phase geometry benefits from evaluating
    /// phase functions on the complete atmosphere grid instead of only at the
    /// endpoint interpolation stencils that use that geometry.
    ///
    /// The grouped path is especially effective for 1D atmospheres, where
    /// many rays share a small set of scattering geometries.  In 2D the
    /// geometry count can grow much faster than the number of endpoint
    /// stencils; the direct path avoids forming that mostly-unused Cartesian
    /// product and its wavelength-local scratch arrays.
    pub fn uses_grouped_phase_evaluation(&self, num_atmosphere_points: usize) -> bool {
        self.num_grouped_phase_geometries(num_atmosphere_points) != 0
    }

    pub fn num_grouped_phase_geometries(&self, num_atmosphere_points: usize) -> usize {
        (0..self.num_phase_geometries())
            .filter(|&phase_geometry| {
                self.uses_grouped_phase_geometry(phase_geometry, num_atmosphere_points)
            })
            .count()
    }

    pub fn num_unique_endpoint_stencils(&self) -> usize {
        self.endpoint_atmosphere_offsets.len() - 1
    }

    /// Persistent heap storage owned by this geometry, excluding `Vec` headers.
    pub fn storage_bytes(&self) -> usize {
        self.ray_layer_offsets.capacity() * size_of::<u32>()
            + self.layer_atmosphere_offsets.capacity() * size_of::<u32>()
            + self.atmosphere_indices.capacity() * size_of::<u32>()
            + self.optical_depth_weights.capacity() * size_of::<f64>()
            + self.albedo_weights.capacity() * size_of::<f64>()
            + self.layer_exit_endpoint_indices.capacity() * size_of::<u32>()
            + self.ray_terminal_endpoint_indices.capacity() * size_of::<u32>()
            + self.layer_exit_phase_endpoint_slots.capacity() * size_of::<u16>()
            + self.ray_terminal_phase_endpoint_slots.capacity() * size_of::<u16>()
            + self.endpoint_atmosphere_offsets.capacity() * size_of::<u32>()
            + self.endpoint_atmosphere_indices.capacity() * size_of::<u32>()
            + self.endpoint_weights.capacity() * size_of::<f64>()
            + self.layer_distance.capacity() * size_of::<f64>()
            + self.layer_start_fraction.capacity() * size_of::<f64>()
            + self.phase_geometry_ray_offsets.capacity() * size_of::<u32>()
            + self.phase_geometry_ray_indices.capacity() * size_of::<u32>()
            + self.phase_geometry_endpoint_offsets.capacity() * size_of::<u32>()
            + self.phase_geometry_endpoint_indices.capacity() * size_of::<u32>()
            + self.phase_geometry_endpoint_evaluations.capacity() * size_of::<usize>()
            + self.ray_phase_basis.capacity() * size_of::<f64>()
            + self.ray_polarized_phase_basis.capacity() * size_of::<f64>()
            + self.phase_q_factor.capacity() * size_of::<f64>()
            + self.phase_u_factor.capacity() * size_of::<f64>()
            + self.layer_source_offsets.capacity() * size_of::<u32>()
            + self.ray_transport_value_offsets.capacity() * size_of::<u32>()
            + self.ray_transport_row_nnz.capacity() * size_of::<u32>()
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

    /// Assembles the transport operator and its directional derivative for
    /// any Stokes dimension. Layer attenuation scratch is retained for the
    /// C++ polarized first-order source during the staged migration.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_transport_jvp(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        extinction_tangent: &[f64],
        single_scatter_albedo_tangent: &[f64],
        transport_values: &mut [f64],
        transport_value_tangent: &mut [f64],
        layer_optical_depth: &mut [f64],
        layer_attenuation: &mut [f64],
        layer_prefix_attenuation: &mut [f64],
        layer_prefix_attenuation_tangent: &mut [f64],
        ray_end_attenuation: &mut [f64],
        ray_end_attenuation_tangent: &mut [f64],
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
        if extinction_tangent.len() != extinction.len()
            || single_scatter_albedo_tangent.len() != single_scatter_albedo.len()
            || transport_value_tangent.len() != transport_values.len()
            || layer_optical_depth.len() != self.num_layers()
            || layer_prefix_attenuation_tangent.len() != self.num_layers()
            || ray_end_attenuation_tangent.len() != self.num_rays()
        {
            bail!("ray transport JVP arrays have inconsistent lengths");
        }

        transport_values.fill(0.0);
        transport_value_tangent.fill(0.0);
        for ray_index in 0..self.num_rays() {
            let layer_start = self.ray_layer_offsets[ray_index] as usize;
            let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
            let transport_offset = self.ray_transport_value_offsets[ray_index] as usize;
            let mut current_attenuation = 1.0;
            let mut current_attenuation_tangent = 0.0;

            for layer_index in (layer_start..layer_end).rev() {
                let atmosphere_start = self.layer_atmosphere_offsets[layer_index] as usize;
                let atmosphere_end = self.layer_atmosphere_offsets[layer_index + 1] as usize;
                let mut optical_depth = 0.0;
                let mut optical_depth_tangent = 0.0;
                let mut albedo = 0.0;
                let mut albedo_tangent = 0.0;
                for stencil_index in atmosphere_start..atmosphere_end {
                    let atmosphere_index = self.atmosphere_indices[stencil_index] as usize;
                    optical_depth +=
                        self.optical_depth_weights[stencil_index] * extinction[atmosphere_index];
                    optical_depth_tangent += self.optical_depth_weights[stencil_index]
                        * extinction_tangent[atmosphere_index];
                    albedo += self.albedo_weights[stencil_index]
                        * single_scatter_albedo[atmosphere_index];
                    albedo_tangent += self.albedo_weights[stencil_index]
                        * single_scatter_albedo_tangent[atmosphere_index];
                }
                let attenuation = (-optical_depth).exp();
                let attenuation_tangent = -attenuation * optical_depth_tangent;
                layer_optical_depth[layer_index] = optical_depth;
                layer_attenuation[layer_index] = attenuation;
                layer_prefix_attenuation[layer_index] = current_attenuation;
                layer_prefix_attenuation_tangent[layer_index] = current_attenuation_tangent;

                let source_factor = albedo * (1.0 - attenuation) * current_attenuation;
                let source_factor_tangent =
                    albedo_tangent * (1.0 - attenuation) * current_attenuation
                        - albedo * attenuation_tangent * current_attenuation
                        + albedo * (1.0 - attenuation) * current_attenuation_tangent;
                let source_start = self.layer_source_offsets[layer_index] as usize;
                let source_end = self.layer_source_offsets[layer_index + 1] as usize;
                if self.num_stokes == 1 {
                    for source_index in source_start..source_end {
                        let value_index = transport_offset
                            + self.source_value_inner_indices[source_index] as usize;
                        let weight = self.source_weights[source_index];
                        add_packed_value(transport_values, value_index, weight * source_factor);
                        add_packed_value(
                            transport_value_tangent,
                            value_index,
                            weight * source_factor_tangent,
                        );
                    }
                } else {
                    for source_index in source_start..source_end {
                        let weight = self.source_weights[source_index];
                        for stokes in 0..self.num_stokes {
                            let value_index = self.transport_value_index(
                                ray_index,
                                transport_offset,
                                stokes,
                                self.source_value_inner_indices[source_index],
                            );
                            transport_values[value_index] += weight * source_factor;
                            transport_value_tangent[value_index] += weight * source_factor_tangent;
                        }
                    }
                }

                current_attenuation_tangent = current_attenuation_tangent * attenuation
                    + current_attenuation * attenuation_tangent;
                current_attenuation *= attenuation;
            }

            let ground_start = self.ray_ground_offsets[ray_index] as usize;
            let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
            if self.num_stokes == 1 {
                for ground_index in ground_start..ground_end {
                    let value_index =
                        transport_offset + self.ground_value_inner_indices[ground_index] as usize;
                    let weight = self.ground_weights[ground_index];
                    add_packed_value(transport_values, value_index, weight * current_attenuation);
                    add_packed_value(
                        transport_value_tangent,
                        value_index,
                        weight * current_attenuation_tangent,
                    );
                }
            } else {
                for ground_index in ground_start..ground_end {
                    let weight = self.ground_weights[ground_index];
                    for stokes in 0..self.num_stokes {
                        let value_index = self.transport_value_index(
                            ray_index,
                            transport_offset,
                            stokes,
                            self.ground_value_inner_indices[ground_index],
                        );
                        transport_values[value_index] += weight * current_attenuation;
                        transport_value_tangent[value_index] +=
                            weight * current_attenuation_tangent;
                    }
                }
            }
            ray_end_attenuation[ray_index] = current_attenuation;
            ray_end_attenuation_tangent[ray_index] = current_attenuation_tangent;
        }
        Ok(())
    }

    /// Pulls a cotangent on the assembled CSR transport values back to
    /// extinction and single-scatter albedo for any Stokes dimension.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_transport_vjp(
        &self,
        single_scatter_albedo: &[f64],
        transport_value_gradient: &[f64],
        transport_column_indices: &[i32],
        solution: &[f64],
        first_order_radiance_gradient: &[f64],
        layer_optical_depth: &[f64],
        layer_attenuation: &[f64],
        layer_prefix_attenuation: &[f64],
        extinction_gradient: &mut [f64],
        single_scatter_albedo_gradient: &mut [f64],
    ) -> Result<()> {
        let compact_transport_gradient = transport_value_gradient.is_empty();
        if single_scatter_albedo.len() < self.required_atmosphere_len
            || (!compact_transport_gradient
                && transport_value_gradient.len() < self.required_transport_len)
            || layer_optical_depth.len() != self.num_layers()
            || layer_attenuation.len() != self.num_layers()
            || layer_prefix_attenuation.len() != self.num_layers()
            || extinction_gradient.len() != single_scatter_albedo.len()
            || single_scatter_albedo_gradient.len() != single_scatter_albedo.len()
        {
            bail!("ray transport VJP arrays have inconsistent lengths");
        }
        if compact_transport_gradient {
            if transport_column_indices.len() < self.required_transport_len
                || first_order_radiance_gradient.len() != self.num_rays() * self.num_stokes
            {
                bail!("compact ray transport VJP arrays have inconsistent lengths");
            }
            for &column in &transport_column_indices[..self.required_transport_len] {
                if column < 0 || column as usize >= solution.len() {
                    bail!("compact ray transport VJP column is out of range");
                }
            }
        }

        extinction_gradient.fill(0.0);
        single_scatter_albedo_gradient.fill(0.0);
        for ray_index in 0..self.num_rays() {
            let layer_start = self.ray_layer_offsets[ray_index] as usize;
            let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
            let transport_offset = self.ray_transport_value_offsets[ray_index] as usize;
            let mut current_attenuation_cotangent = 0.0;

            let ground_start = self.ray_ground_offsets[ray_index] as usize;
            let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
            if self.num_stokes == 1 {
                let forcing_cotangent = first_order_radiance_gradient[ray_index];
                for ground_index in ground_start..ground_end {
                    let value_index =
                        transport_offset + self.ground_value_inner_indices[ground_index] as usize;
                    let value_cotangent = if compact_transport_gradient {
                        forcing_cotangent
                            * packed_value(
                                solution,
                                packed_column(transport_column_indices, value_index),
                            )
                    } else {
                        packed_value(transport_value_gradient, value_index)
                    };
                    current_attenuation_cotangent +=
                        self.ground_weights[ground_index] * value_cotangent;
                }
            } else {
                for ground_index in ground_start..ground_end {
                    let mut value_cotangent = 0.0;
                    for stokes in 0..self.num_stokes {
                        let value_index = self.transport_value_index(
                            ray_index,
                            transport_offset,
                            stokes,
                            self.ground_value_inner_indices[ground_index],
                        );
                        value_cotangent += if compact_transport_gradient {
                            first_order_radiance_gradient[ray_index * self.num_stokes + stokes]
                                * solution[transport_column_indices[value_index] as usize]
                        } else {
                            transport_value_gradient[value_index]
                        };
                    }
                    current_attenuation_cotangent +=
                        self.ground_weights[ground_index] * value_cotangent;
                }
            }

            for layer_index in layer_start..layer_end {
                let attenuation = layer_attenuation[layer_index];
                let current_attenuation = layer_prefix_attenuation[layer_index];
                let atmosphere_start = self.layer_atmosphere_offsets[layer_index] as usize;
                let atmosphere_end = self.layer_atmosphere_offsets[layer_index + 1] as usize;
                let mut albedo = 0.0;
                for stencil_index in atmosphere_start..atmosphere_end {
                    let atmosphere_index = self.atmosphere_indices[stencil_index] as usize;
                    albedo += self.albedo_weights[stencil_index]
                        * single_scatter_albedo[atmosphere_index];
                }

                let mut source_factor_cotangent = 0.0;
                let source_start = self.layer_source_offsets[layer_index] as usize;
                let source_end = self.layer_source_offsets[layer_index + 1] as usize;
                if self.num_stokes == 1 {
                    let forcing_cotangent = first_order_radiance_gradient[ray_index];
                    for source_index in source_start..source_end {
                        let value_index = transport_offset
                            + self.source_value_inner_indices[source_index] as usize;
                        let value_cotangent = if compact_transport_gradient {
                            forcing_cotangent
                                * packed_value(
                                    solution,
                                    packed_column(transport_column_indices, value_index),
                                )
                        } else {
                            packed_value(transport_value_gradient, value_index)
                        };
                        source_factor_cotangent +=
                            self.source_weights[source_index] * value_cotangent;
                    }
                } else {
                    for source_index in source_start..source_end {
                        let mut value_cotangent = 0.0;
                        for stokes in 0..self.num_stokes {
                            let value_index = self.transport_value_index(
                                ray_index,
                                transport_offset,
                                stokes,
                                self.source_value_inner_indices[source_index],
                            );
                            value_cotangent += if compact_transport_gradient {
                                first_order_radiance_gradient[ray_index * self.num_stokes + stokes]
                                    * solution[transport_column_indices[value_index] as usize]
                            } else {
                                transport_value_gradient[value_index]
                            };
                        }
                        source_factor_cotangent +=
                            self.source_weights[source_index] * value_cotangent;
                    }
                }

                let albedo_cotangent =
                    source_factor_cotangent * (1.0 - attenuation) * current_attenuation;
                let mut attenuation_cotangent =
                    -source_factor_cotangent * albedo * current_attenuation;
                let mut prefix_cotangent = source_factor_cotangent * albedo * (1.0 - attenuation);
                prefix_cotangent += current_attenuation_cotangent * attenuation;
                attenuation_cotangent += current_attenuation_cotangent * current_attenuation;
                let optical_depth_cotangent = -attenuation * attenuation_cotangent;

                for stencil_index in atmosphere_start..atmosphere_end {
                    let atmosphere_index = self.atmosphere_indices[stencil_index] as usize;
                    extinction_gradient[atmosphere_index] +=
                        self.optical_depth_weights[stencil_index] * optical_depth_cotangent;
                    single_scatter_albedo_gradient[atmosphere_index] +=
                        self.albedo_weights[stencil_index] * albedo_cotangent;
                }
                current_attenuation_cotangent = prefix_cotangent;
            }
        }
        Ok(())
    }

    /// Fuses transport assembly with exact single-scatter layer
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
        end_of_ray_source: &[f64],
        transport_values: &mut [f64],
        first_order_radiance: &mut [f64],
        layer_optical_depth: &mut [f64],
        layer_attenuation: &mut [f64],
        layer_prefix_attenuation: &mut [f64],
        ray_end_attenuation: &mut [f64],
    ) -> Result<()> {
        self.validate_first_order()?;
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
            .checked_mul(self.num_phase_moments * self.phase_coefficient_families())
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
        if first_order_radiance.len() != self.num_rays() * self.num_stokes {
            bail!(
                "first-order radiance length {} does not match {} ray-Stokes values",
                first_order_radiance.len(),
                self.num_rays() * self.num_stokes
            );
        }
        if end_of_ray_source.len() != first_order_radiance.len() {
            bail!("end-of-ray source has an inconsistent length");
        }

        let forcing = FirstOrderInputs {
            legendre_coefficients,
            maximum_order,
            solar_transmission,
            end_of_ray_source,
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

    /// Directional derivative of the fused transport and atmospheric
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
        end_of_ray_source: &[f64],
        end_of_ray_source_tangent: &[f64],
        transport_values: &mut [f64],
        transport_value_tangent: &mut [f64],
        first_order_radiance_tangent: &mut [f64],
        ray_end_attenuation: &mut [f64],
        ray_end_attenuation_tangent: &mut [f64],
    ) -> Result<()> {
        self.validate_first_order()?;
        if extinction.len() < self.required_atmosphere_len
            || single_scatter_albedo.len() != extinction.len()
            || extinction_tangent.len() != extinction.len()
            || single_scatter_albedo_tangent.len() != extinction.len()
        {
            bail!("first-order JVP atmosphere arrays have inconsistent lengths");
        }
        let required_coefficients = extinction
            .len()
            .checked_mul(self.num_phase_moments * self.phase_coefficient_families())
            .ok_or_else(|| anyhow!("phase coefficient size overflow"))?;
        if legendre_coefficients.len() != required_coefficients
            || legendre_coefficient_tangent.len() != required_coefficients
            || maximum_order.len() != extinction.len()
        {
            bail!("first-order JVP phase arrays have inconsistent lengths");
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
            bail!("first-order JVP solar-transmission arrays have inconsistent lengths");
        }
        if transport_values.len() < self.required_transport_len
            || transport_value_tangent.len() != transport_values.len()
        {
            bail!("first-order JVP transport arrays have inconsistent lengths");
        }
        if first_order_radiance_tangent.len() != self.num_rays() * self.num_stokes
            || end_of_ray_source.len() != first_order_radiance_tangent.len()
            || end_of_ray_source_tangent.len() != first_order_radiance_tangent.len()
            || ray_end_attenuation.len() != self.num_rays()
            || ray_end_attenuation_tangent.len() != self.num_rays()
        {
            bail!("first-order JVP ray output arrays have inconsistent lengths");
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
            end_of_ray_source,
            end_of_ray_source_tangent,
        };
        let endpoint_medium = self.interpolate_endpoint_medium_jvp(
            extinction,
            single_scatter_albedo,
            forcing.solar_transmission,
            extinction_tangent,
            single_scatter_albedo_tangent,
            forcing.solar_transmission_tangent,
        );
        let grouped_phase = self.uses_grouped_phase_evaluation(maximum_order.len());
        let phase_value_size = if grouped_phase {
            maximum_order.len() * self.num_stokes
        } else {
            0
        };
        let mut phase_values = vec![0.0; phase_value_size];
        let mut phase_value_tangent = vec![0.0; phase_value_size];
        let endpoint_value_size = self.phase_endpoint_cache_size * self.num_stokes;
        let mut endpoint_scattering_values = vec![0.0; endpoint_value_size];
        let mut endpoint_scattering_value_tangent = vec![0.0; endpoint_value_size];
        for phase_geometry in 0..self.num_phase_geometries() {
            let grouped_phase_geometry =
                self.uses_grouped_phase_geometry(phase_geometry, maximum_order.len());
            if grouped_phase_geometry {
                self.fill_phase_values_jvp(
                    phase_geometry,
                    forcing.legendre_coefficients,
                    forcing.legendre_coefficient_tangent,
                    forcing.maximum_order,
                    &mut phase_values,
                    &mut phase_value_tangent,
                );
            }
            let phase_values_for_geometry = if grouped_phase_geometry {
                phase_values.as_slice()
            } else {
                &[]
            };
            let phase_tangent_for_geometry = if grouped_phase_geometry {
                phase_value_tangent.as_slice()
            } else {
                &[]
            };
            self.fill_endpoint_scattering_values_jvp(
                phase_geometry,
                &forcing,
                phase_values_for_geometry,
                phase_tangent_for_geometry,
                &endpoint_medium,
                &mut endpoint_scattering_values,
                &mut endpoint_scattering_value_tangent,
            );
            for &ray_index in self.phase_geometry_rays(phase_geometry) {
                let ray_index = ray_index as usize;
                let layer_start = self.ray_layer_offsets[ray_index] as usize;
                let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
                let transport_offset = self.ray_transport_value_offsets[ray_index] as usize;
                let mut current_attenuation = 1.0;
                let mut current_attenuation_tangent = 0.0;
                let mut ray_first_order_tangent = [0.0; 3];
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
                    if self.num_stokes == 1 {
                        for source_index in source_start..source_end {
                            let value_index = transport_offset
                                + self.source_value_inner_indices[source_index] as usize;
                            let weight = self.source_weights[source_index];
                            add_packed_value(transport_values, value_index, weight * source_factor);
                            add_packed_value(
                                transport_value_tangent,
                                value_index,
                                weight * source_factor_tangent,
                            );
                        }
                    } else {
                        for source_index in source_start..source_end {
                            let weight = self.source_weights[source_index];
                            for stokes in 0..self.num_stokes {
                                let value_index = self.transport_value_index(
                                    ray_index,
                                    transport_offset,
                                    stokes,
                                    self.source_value_inner_indices[source_index],
                                );
                                transport_values[value_index] += weight * source_factor;
                                transport_value_tangent[value_index] +=
                                    weight * source_factor_tangent;
                            }
                        }
                    }

                    let start_source = shared_endpoint_source.unwrap_or_else(|| {
                        self.endpoint_source_jvp(
                            ray_index,
                            layer_index,
                            true,
                            &forcing,
                            &endpoint_medium,
                            &endpoint_scattering_values,
                            &endpoint_scattering_value_tangent,
                        )
                    });
                    let end_source = self.endpoint_source_jvp(
                        ray_index,
                        layer_index,
                        false,
                        &forcing,
                        &endpoint_medium,
                        &endpoint_scattering_values,
                        &endpoint_scattering_value_tangent,
                    );
                    if self.layer_distance[layer_index] >= 1.0e-4 {
                        let (integration_factor, integration_factor_derivative) =
                            constant_source_factor_and_derivative(optical_depth, attenuation);
                        let integration_factor_tangent =
                            integration_factor_derivative * optical_depth_tangent;
                        let distance = self.layer_distance[layer_index];
                        if self.num_stokes == 1 {
                            let endpoint_quadrature = start_source.value[0]
                                * self.layer_start_fraction[layer_index]
                                + end_source.value[0]
                                    * (1.0 - self.layer_start_fraction[layer_index]);
                            let endpoint_quadrature_tangent = start_source.tangent[0]
                                * self.layer_start_fraction[layer_index]
                                + end_source.tangent[0]
                                    * (1.0 - self.layer_start_fraction[layer_index]);
                            ray_first_order_tangent[0] += distance
                                * (current_attenuation_tangent
                                    * integration_factor
                                    * endpoint_quadrature
                                    + current_attenuation
                                        * (integration_factor_tangent * endpoint_quadrature
                                            + integration_factor * endpoint_quadrature_tangent));
                        } else {
                            for (stokes, ray_tangent) in ray_first_order_tangent[..self.num_stokes]
                                .iter_mut()
                                .enumerate()
                            {
                                let endpoint_quadrature = start_source.value[stokes]
                                    * self.layer_start_fraction[layer_index]
                                    + end_source.value[stokes]
                                        * (1.0 - self.layer_start_fraction[layer_index]);
                                let endpoint_quadrature_tangent = start_source.tangent[stokes]
                                    * self.layer_start_fraction[layer_index]
                                    + end_source.tangent[stokes]
                                        * (1.0 - self.layer_start_fraction[layer_index]);
                                *ray_tangent += distance
                                    * (current_attenuation_tangent
                                        * integration_factor
                                        * endpoint_quadrature
                                        + current_attenuation
                                            * (integration_factor_tangent * endpoint_quadrature
                                                + integration_factor
                                                    * endpoint_quadrature_tangent));
                            }
                        }
                    }
                    shared_endpoint_source = Some(end_source);

                    current_attenuation_tangent = current_attenuation_tangent * attenuation
                        + current_attenuation * attenuation_tangent;
                    current_attenuation *= attenuation;
                }

                let ground_start = self.ray_ground_offsets[ray_index] as usize;
                let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
                if self.num_stokes == 1 {
                    for ground_index in ground_start..ground_end {
                        let value_index = transport_offset
                            + self.ground_value_inner_indices[ground_index] as usize;
                        let weight = self.ground_weights[ground_index];
                        add_packed_value(
                            transport_values,
                            value_index,
                            weight * current_attenuation,
                        );
                        add_packed_value(
                            transport_value_tangent,
                            value_index,
                            weight * current_attenuation_tangent,
                        );
                    }
                } else {
                    for ground_index in ground_start..ground_end {
                        let weight = self.ground_weights[ground_index];
                        for stokes in 0..self.num_stokes {
                            let value_index = self.transport_value_index(
                                ray_index,
                                transport_offset,
                                stokes,
                                self.ground_value_inner_indices[ground_index],
                            );
                            transport_values[value_index] += weight * current_attenuation;
                            transport_value_tangent[value_index] +=
                                weight * current_attenuation_tangent;
                        }
                    }
                }
                let forcing_start = ray_index * self.num_stokes;
                if self.num_stokes == 1 {
                    ray_first_order_tangent[0] += current_attenuation_tangent
                        * forcing.end_of_ray_source[forcing_start]
                        + current_attenuation * forcing.end_of_ray_source_tangent[forcing_start];
                } else {
                    for (stokes, ray_tangent) in ray_first_order_tangent[..self.num_stokes]
                        .iter_mut()
                        .enumerate()
                    {
                        *ray_tangent += current_attenuation_tangent
                            * forcing.end_of_ray_source[forcing_start + stokes]
                            + current_attenuation
                                * forcing.end_of_ray_source_tangent[forcing_start + stokes];
                    }
                }
                first_order_radiance_tangent[forcing_start..forcing_start + self.num_stokes]
                    .copy_from_slice(&ray_first_order_tangent[..self.num_stokes]);
                ray_end_attenuation[ray_index] = current_attenuation;
                ray_end_attenuation_tangent[ray_index] = current_attenuation_tangent;
            }
        }
        Ok(())
    }

    /// Reverse product for the fused transport and atmospheric
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
        ray_end_attenuation: &[f64],
        ray_end_attenuation_cotangent: &[f64],
        end_of_ray_source: &[f64],
        layer_optical_depth: &[f64],
        layer_attenuation: &[f64],
        layer_prefix_attenuation: &[f64],
        extinction_gradient: &mut [f64],
        single_scatter_albedo_gradient: &mut [f64],
        legendre_coefficient_gradient: &mut [f64],
        solar_transmission_gradient: &mut [f64],
        end_of_ray_source_gradient: &mut [f64],
    ) -> Result<()> {
        self.validate_first_order()?;
        if extinction.len() < self.required_atmosphere_len
            || single_scatter_albedo.len() != extinction.len()
            || extinction_gradient.len() != extinction.len()
            || single_scatter_albedo_gradient.len() != extinction.len()
        {
            bail!("first-order VJP atmosphere arrays have inconsistent lengths");
        }
        let required_coefficients = extinction
            .len()
            .checked_mul(self.num_phase_moments * self.phase_coefficient_families())
            .ok_or_else(|| anyhow!("phase coefficient size overflow"))?;
        if legendre_coefficients.len() != required_coefficients
            || legendre_coefficient_gradient.len() != required_coefficients
            || maximum_order.len() != extinction.len()
        {
            bail!("first-order VJP phase arrays have inconsistent lengths");
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
            bail!("first-order VJP solar-transmission arrays have inconsistent lengths");
        }
        let compact_transport_gradient = transport_value_gradient.is_empty();
        if !compact_transport_gradient
            && transport_value_gradient.len() < self.required_transport_len
        {
            bail!("first-order VJP transport gradient is too short");
        }
        if compact_transport_gradient {
            if transport_column_indices.len() < self.required_transport_len {
                bail!("first-order compact VJP transport columns are too short");
            }
            for &column in &transport_column_indices[..self.required_transport_len] {
                if column < 0 || column as usize >= solution.len() {
                    bail!("first-order compact VJP transport column is out of range");
                }
            }
        }
        if first_order_radiance_gradient.len() != self.num_rays() * self.num_stokes
            || ray_end_attenuation.len() != self.num_rays()
            || ray_end_attenuation_cotangent.len() != self.num_rays()
            || end_of_ray_source.len() != first_order_radiance_gradient.len()
            || end_of_ray_source_gradient.len() != first_order_radiance_gradient.len()
        {
            bail!("first-order VJP ray gradient arrays have inconsistent lengths");
        }
        if layer_optical_depth.len() != self.num_layers()
            || layer_attenuation.len() != self.num_layers()
            || layer_prefix_attenuation.len() != self.num_layers()
        {
            bail!("first-order VJP layer scratch arrays have inconsistent lengths");
        }
        extinction_gradient.fill(0.0);
        single_scatter_albedo_gradient.fill(0.0);
        legendre_coefficient_gradient.fill(0.0);
        solar_transmission_gradient.fill(0.0);
        end_of_ray_source_gradient.fill(0.0);
        let forcing = FirstOrderInputs {
            legendre_coefficients,
            maximum_order,
            solar_transmission,
            end_of_ray_source,
        };
        let endpoint_medium = self.interpolate_endpoint_medium(
            extinction,
            single_scatter_albedo,
            forcing.solar_transmission,
        );
        let num_endpoints = self.num_unique_endpoint_stencils();
        let mut endpoint_extinction_gradient = vec![0.0; num_endpoints];
        let mut endpoint_albedo_gradient = vec![0.0; num_endpoints];
        let mut endpoint_solar_gradient = if self.solar_transmission_on_atmosphere_grid {
            vec![0.0; num_endpoints]
        } else {
            Vec::new()
        };
        let phase_endpoint_value_size = self.phase_endpoint_cache_size * self.num_stokes;
        let mut endpoint_phase_gradient = vec![0.0; phase_endpoint_value_size];
        let mut endpoint_weighted_source_cotangent = vec![0.0; phase_endpoint_value_size];
        let mut endpoint_source_cotangent = if self.solar_transmission_on_atmosphere_grid {
            vec![0.0; phase_endpoint_value_size]
        } else {
            Vec::new()
        };
        let grouped_phase = self.uses_grouped_phase_evaluation(maximum_order.len());
        let phase_value_size = if grouped_phase {
            maximum_order.len() * self.num_stokes
        } else {
            0
        };
        let mut phase_values = vec![0.0; phase_value_size];
        let mut phase_value_gradient = vec![0.0; phase_value_size];
        let mut endpoint_phase_values = vec![0.0; phase_endpoint_value_size];
        for phase_geometry in 0..self.num_phase_geometries() {
            let grouped_phase_geometry =
                self.uses_grouped_phase_geometry(phase_geometry, maximum_order.len());
            if grouped_phase_geometry {
                self.fill_phase_values(
                    phase_geometry,
                    forcing.legendre_coefficients,
                    forcing.maximum_order,
                    &mut phase_values,
                );
            }
            phase_value_gradient.fill(0.0);
            endpoint_phase_gradient.fill(0.0);
            endpoint_weighted_source_cotangent.fill(0.0);
            endpoint_source_cotangent.fill(0.0);
            let phase_values_for_geometry = if grouped_phase_geometry {
                phase_values.as_slice()
            } else {
                &[]
            };
            self.fill_endpoint_phase_cache(
                phase_geometry,
                &forcing,
                phase_values_for_geometry,
                &mut endpoint_phase_values,
            );
            for &ray_index in self.phase_geometry_rays(phase_geometry) {
                let ray_index = ray_index as usize;
                let layer_start = self.ray_layer_offsets[ray_index] as usize;
                let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
                let transport_offset = self.ray_transport_value_offsets[ray_index] as usize;
                let forcing_start = ray_index * self.num_stokes;
                let mut forcing_cotangent = [0.0; 3];
                let mut current_attenuation_cotangent = ray_end_attenuation_cotangent[ray_index];
                let first_order_active = if self.num_stokes == 1 {
                    let forcing_value_cotangent = first_order_radiance_gradient[forcing_start];
                    forcing_cotangent[0] = forcing_value_cotangent;
                    current_attenuation_cotangent +=
                        forcing_value_cotangent * forcing.end_of_ray_source[forcing_start];
                    end_of_ray_source_gradient[forcing_start] =
                        ray_end_attenuation[ray_index] * forcing_value_cotangent;
                    forcing_value_cotangent != 0.0
                } else {
                    forcing_cotangent[..self.num_stokes].copy_from_slice(
                        &first_order_radiance_gradient
                            [forcing_start..forcing_start + self.num_stokes],
                    );
                    for (stokes, &forcing_value_cotangent) in
                        forcing_cotangent[..self.num_stokes].iter().enumerate()
                    {
                        current_attenuation_cotangent += forcing_value_cotangent
                            * forcing.end_of_ray_source[forcing_start + stokes];
                        end_of_ray_source_gradient[forcing_start + stokes] =
                            ray_end_attenuation[ray_index] * forcing_value_cotangent;
                    }
                    forcing_cotangent[..self.num_stokes]
                        .iter()
                        .any(|value| *value != 0.0)
                };
                let mut shared_endpoint = None;
                let mut shared_endpoint_cotangent = [0.0; 3];

                let ground_start = self.ray_ground_offsets[ray_index] as usize;
                let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
                if self.num_stokes == 1 {
                    for ground_index in ground_start..ground_end {
                        let value_index = transport_offset
                            + self.ground_value_inner_indices[ground_index] as usize;
                        let value_gradient = if compact_transport_gradient {
                            forcing_cotangent[0]
                                * packed_value(
                                    solution,
                                    packed_column(transport_column_indices, value_index),
                                )
                        } else {
                            packed_value(transport_value_gradient, value_index)
                        };
                        current_attenuation_cotangent +=
                            self.ground_weights[ground_index] * value_gradient;
                    }
                } else {
                    for ground_index in ground_start..ground_end {
                        for (stokes, &forcing_value_cotangent) in
                            forcing_cotangent[..self.num_stokes].iter().enumerate()
                        {
                            let value_index = self.transport_value_index(
                                ray_index,
                                transport_offset,
                                stokes,
                                self.ground_value_inner_indices[ground_index],
                            );
                            let value_gradient = if compact_transport_gradient {
                                forcing_value_cotangent
                                    * solution[transport_column_indices[value_index] as usize]
                            } else {
                                transport_value_gradient[value_index]
                            };
                            current_attenuation_cotangent +=
                                self.ground_weights[ground_index] * value_gradient;
                        }
                    }
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
                    if self.num_stokes == 1 {
                        for source_index in source_start..source_end {
                            let value_index = transport_offset
                                + self.source_value_inner_indices[source_index] as usize;
                            let value_gradient = if compact_transport_gradient {
                                forcing_cotangent[0]
                                    * packed_value(
                                        solution,
                                        packed_column(transport_column_indices, value_index),
                                    )
                            } else {
                                packed_value(transport_value_gradient, value_index)
                            };
                            source_factor_cotangent +=
                                self.source_weights[source_index] * value_gradient;
                        }
                    } else {
                        for source_index in source_start..source_end {
                            for (stokes, &forcing_value_cotangent) in
                                forcing_cotangent[..self.num_stokes].iter().enumerate()
                            {
                                let value_index = self.transport_value_index(
                                    ray_index,
                                    transport_offset,
                                    stokes,
                                    self.source_value_inner_indices[source_index],
                                );
                                let value_gradient = if compact_transport_gradient {
                                    forcing_value_cotangent
                                        * solution[transport_column_indices[value_index] as usize]
                                } else {
                                    transport_value_gradient[value_index]
                                };
                                source_factor_cotangent +=
                                    self.source_weights[source_index] * value_gradient;
                            }
                        }
                    }
                    let albedo_cotangent =
                        source_factor_cotangent * (1.0 - attenuation) * current_attenuation;
                    let mut attenuation_cotangent =
                        -source_factor_cotangent * albedo * current_attenuation;
                    let mut prefix_cotangent =
                        source_factor_cotangent * albedo * (1.0 - attenuation);
                    let mut optical_depth_cotangent = 0.0;

                    if first_order_active {
                        // The forward pass evaluates every shared straight-ray
                        // endpoint once: start(i) is end(i + 1). Reverse the same
                        // graph, accumulating the carried start cotangent into the
                        // next layer's end before visiting that endpoint once.
                        let end_endpoint = shared_endpoint.unwrap_or_else(|| {
                            self.endpoint_primal(
                                ray_index,
                                layer_index,
                                false,
                                &forcing,
                                &endpoint_medium,
                                &endpoint_phase_values,
                            )
                        });
                        let start_endpoint = if layer_index + 1 < layer_end {
                            self.endpoint_primal(
                                ray_index,
                                layer_index + 1,
                                false,
                                &forcing,
                                &endpoint_medium,
                                &endpoint_phase_values,
                            )
                        } else {
                            self.endpoint_primal(
                                ray_index,
                                layer_index,
                                true,
                                &forcing,
                                &endpoint_medium,
                                &endpoint_phase_values,
                            )
                        };
                        let mut start_cotangent = [0.0; 3];
                        let mut end_cotangent = shared_endpoint_cotangent;

                        if self.layer_distance[layer_index] >= 1.0e-4 {
                            let (integration_factor, integration_factor_derivative) =
                                constant_source_factor_and_derivative(optical_depth, attenuation);
                            let distance = self.layer_distance[layer_index];
                            let forcing_endpoint_dot = if self.num_stokes == 1 {
                                forcing_cotangent[0]
                                    * (start_endpoint[0] * self.layer_start_fraction[layer_index]
                                        + end_endpoint[0]
                                            * (1.0 - self.layer_start_fraction[layer_index]))
                            } else {
                                let mut endpoint_quadrature = [0.0; 3];
                                for (stokes, endpoint_value) in endpoint_quadrature
                                    [..self.num_stokes]
                                    .iter_mut()
                                    .enumerate()
                                {
                                    *endpoint_value = start_endpoint[stokes]
                                        * self.layer_start_fraction[layer_index]
                                        + end_endpoint[stokes]
                                            * (1.0 - self.layer_start_fraction[layer_index]);
                                }
                                dot_stokes(
                                    &forcing_cotangent,
                                    &endpoint_quadrature,
                                    self.num_stokes,
                                )
                            };
                            prefix_cotangent +=
                                forcing_endpoint_dot * integration_factor * distance;
                            optical_depth_cotangent += forcing_endpoint_dot
                                * current_attenuation
                                * distance
                                * integration_factor_derivative;
                            let quadrature_scale =
                                current_attenuation * integration_factor * distance;
                            if self.num_stokes == 1 {
                                let cotangent = forcing_cotangent[0] * quadrature_scale;
                                start_cotangent[0] =
                                    cotangent * self.layer_start_fraction[layer_index];
                                end_cotangent[0] +=
                                    cotangent * (1.0 - self.layer_start_fraction[layer_index]);
                            } else {
                                for stokes in 0..self.num_stokes {
                                    let cotangent = forcing_cotangent[stokes] * quadrature_scale;
                                    start_cotangent[stokes] =
                                        cotangent * self.layer_start_fraction[layer_index];
                                    end_cotangent[stokes] +=
                                        cotangent * (1.0 - self.layer_start_fraction[layer_index]);
                                }
                            }
                        }
                        self.accumulate_endpoint_source_cotangent(
                            ray_index,
                            layer_index,
                            false,
                            end_cotangent,
                            forcing.solar_transmission,
                            &endpoint_medium,
                            &endpoint_phase_values,
                            &mut endpoint_weighted_source_cotangent,
                            &mut endpoint_source_cotangent,
                            solar_transmission_gradient,
                        );
                        shared_endpoint = Some(start_endpoint);
                        shared_endpoint_cotangent = start_cotangent;
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
                if first_order_active && layer_start < layer_end {
                    self.accumulate_endpoint_source_cotangent(
                        ray_index,
                        layer_end - 1,
                        true,
                        shared_endpoint_cotangent,
                        forcing.solar_transmission,
                        &endpoint_medium,
                        &endpoint_phase_values,
                        &mut endpoint_weighted_source_cotangent,
                        &mut endpoint_source_cotangent,
                        solar_transmission_gradient,
                    );
                }
            }
            let endpoint_phase_gradient_for_geometry = if grouped_phase_geometry {
                endpoint_phase_gradient.as_mut_slice()
            } else {
                &mut []
            };
            self.finish_endpoint_source_vjp(
                phase_geometry,
                forcing.maximum_order,
                &endpoint_medium,
                &endpoint_phase_values,
                &endpoint_weighted_source_cotangent,
                &endpoint_source_cotangent,
                &mut endpoint_extinction_gradient,
                &mut endpoint_albedo_gradient,
                endpoint_phase_gradient_for_geometry,
                &mut endpoint_solar_gradient,
                legendre_coefficient_gradient,
            );
            if grouped_phase_geometry {
                for (endpoint_slot, &endpoint_index) in self
                    .phase_geometry_endpoints(phase_geometry)
                    .iter()
                    .enumerate()
                {
                    let endpoint_index = endpoint_index as usize;
                    let (atmosphere_indices, weights) =
                        self.endpoint_stencil_by_index(endpoint_index);
                    let endpoint_start = self
                        .phase_endpoint_cache_slot(endpoint_slot, endpoint_index)
                        * self.num_stokes;
                    for (&atmosphere_index, &weight) in atmosphere_indices.iter().zip(weights) {
                        let phase_start = atmosphere_index as usize * self.num_stokes;
                        if self.num_stokes == 1 {
                            phase_value_gradient[phase_start] +=
                                weight * endpoint_phase_gradient[endpoint_start];
                        } else {
                            for stokes in 0..self.num_stokes {
                                phase_value_gradient[phase_start + stokes] +=
                                    weight * endpoint_phase_gradient[endpoint_start + stokes];
                            }
                        }
                    }
                }
                self.accumulate_phase_values_vjp(
                    phase_geometry,
                    forcing.maximum_order,
                    &phase_value_gradient,
                    legendre_coefficient_gradient,
                );
            }
        }
        for endpoint_index in 0..num_endpoints {
            let (atmosphere_indices, weights) = self.endpoint_stencil_by_index(endpoint_index);
            for (&atmosphere_index, &weight) in atmosphere_indices.iter().zip(weights) {
                let atmosphere_index = atmosphere_index as usize;
                extinction_gradient[atmosphere_index] +=
                    weight * endpoint_extinction_gradient[endpoint_index];
                single_scatter_albedo_gradient[atmosphere_index] +=
                    weight * endpoint_albedo_gradient[endpoint_index];
                if self.solar_transmission_on_atmosphere_grid {
                    solar_transmission_gradient[atmosphere_index] +=
                        weight * endpoint_solar_gradient[endpoint_index];
                }
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
        first_order: Option<(FirstOrderInputs<'_>, &mut [f64])>,
    ) {
        transport_values.fill(0.0);
        if let Some((forcing, radiance)) = first_order {
            radiance.fill(0.0);
            let endpoint_medium = self.interpolate_endpoint_medium(
                extinction,
                single_scatter_albedo,
                forcing.solar_transmission,
            );
            let grouped_phase = self.uses_grouped_phase_evaluation(forcing.maximum_order.len());
            let phase_value_size = if grouped_phase {
                forcing.maximum_order.len() * self.num_stokes
            } else {
                0
            };
            let mut phase_values = vec![0.0; phase_value_size];
            let mut endpoint_scattering_values =
                vec![0.0; self.phase_endpoint_cache_size * self.num_stokes];
            for phase_geometry in 0..self.num_phase_geometries() {
                let grouped_phase_geometry =
                    self.uses_grouped_phase_geometry(phase_geometry, forcing.maximum_order.len());
                if grouped_phase_geometry {
                    self.fill_phase_values(
                        phase_geometry,
                        forcing.legendre_coefficients,
                        forcing.maximum_order,
                        &mut phase_values,
                    );
                }
                let phase_values_for_geometry = if grouped_phase_geometry {
                    phase_values.as_slice()
                } else {
                    &[]
                };
                self.fill_endpoint_scattering_values(
                    phase_geometry,
                    forcing.legendre_coefficients,
                    forcing.maximum_order,
                    phase_values_for_geometry,
                    &endpoint_medium,
                    &mut endpoint_scattering_values,
                );
                for &ray_index in self.phase_geometry_rays(phase_geometry) {
                    self.assemble_ray(
                        ray_index as usize,
                        extinction,
                        single_scatter_albedo,
                        transport_values,
                        layer_optical_depth,
                        layer_attenuation,
                        layer_prefix_attenuation,
                        ray_end_attenuation,
                        Some((&forcing, &endpoint_scattering_values)),
                        &endpoint_medium,
                        radiance,
                    );
                }
            }
        } else {
            let mut no_first_order_radiance = [];
            for ray_index in 0..self.num_rays() {
                self.assemble_ray(
                    ray_index,
                    extinction,
                    single_scatter_albedo,
                    transport_values,
                    layer_optical_depth,
                    layer_attenuation,
                    layer_prefix_attenuation,
                    ray_end_attenuation,
                    None,
                    &[],
                    &mut no_first_order_radiance,
                );
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble_ray(
        &self,
        ray_index: usize,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        transport_values: &mut [f64],
        layer_optical_depth: &mut [f64],
        layer_attenuation: &mut [f64],
        layer_prefix_attenuation: &mut [f64],
        ray_end_attenuation: &mut [f64],
        first_order: Option<(&FirstOrderInputs<'_>, &[f64])>,
        endpoint_medium: &[EndpointMedium],
        first_order_radiance: &mut [f64],
    ) {
        let retain_layer_scratch = !layer_optical_depth.is_empty();
        let layer_start = self.ray_layer_offsets[ray_index] as usize;
        let layer_end = self.ray_layer_offsets[ray_index + 1] as usize;
        let transport_offset = self.ray_transport_value_offsets[ray_index] as usize;
        let mut current_attenuation = 1.0;
        let mut ray_first_order = [0.0; 3];
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
            if self.num_stokes == 1 {
                for source_index in source_start..source_end {
                    let value_index =
                        transport_offset + self.source_value_inner_indices[source_index] as usize;
                    add_packed_value(
                        transport_values,
                        value_index,
                        self.source_weights[source_index] * source_factor,
                    );
                }
            } else {
                for source_index in source_start..source_end {
                    for stokes in 0..self.num_stokes {
                        let value_index = self.transport_value_index(
                            ray_index,
                            transport_offset,
                            stokes,
                            self.source_value_inner_indices[source_index],
                        );
                        transport_values[value_index] +=
                            self.source_weights[source_index] * source_factor;
                    }
                }
            }

            if let Some((forcing, endpoint_scattering_values)) = first_order {
                // Consecutive straight-ray layers share a boundary. In
                // observer-to-far-end traversal, the previous exit is the next
                // entrance, so only one new endpoint is evaluated per layer.
                let start_source = shared_endpoint_source.unwrap_or_else(|| {
                    self.endpoint_source(
                        ray_index,
                        layer_index,
                        true,
                        forcing,
                        endpoint_medium,
                        endpoint_scattering_values,
                    )
                });
                let end_source = self.endpoint_source(
                    ray_index,
                    layer_index,
                    false,
                    forcing,
                    endpoint_medium,
                    endpoint_scattering_values,
                );
                if self.layer_distance[layer_index] >= 1.0e-4 {
                    let integration_factor = constant_source_factor(optical_depth, attenuation);
                    let scale =
                        current_attenuation * integration_factor * self.layer_distance[layer_index];
                    if self.num_stokes == 1 {
                        ray_first_order[0] += scale
                            * (start_source[0] * self.layer_start_fraction[layer_index]
                                + end_source[0] * (1.0 - self.layer_start_fraction[layer_index]));
                    } else {
                        for stokes in 0..self.num_stokes {
                            ray_first_order[stokes] += scale
                                * (start_source[stokes] * self.layer_start_fraction[layer_index]
                                    + end_source[stokes]
                                        * (1.0 - self.layer_start_fraction[layer_index]));
                        }
                    }
                }
                shared_endpoint_source = Some(end_source);
            }

            current_attenuation *= attenuation;
        }

        let ground_start = self.ray_ground_offsets[ray_index] as usize;
        let ground_end = self.ray_ground_offsets[ray_index + 1] as usize;
        if self.num_stokes == 1 {
            for ground_index in ground_start..ground_end {
                let value_index =
                    transport_offset + self.ground_value_inner_indices[ground_index] as usize;
                add_packed_value(
                    transport_values,
                    value_index,
                    self.ground_weights[ground_index] * current_attenuation,
                );
            }
        } else {
            for ground_index in ground_start..ground_end {
                for stokes in 0..self.num_stokes {
                    let value_index = self.transport_value_index(
                        ray_index,
                        transport_offset,
                        stokes,
                        self.ground_value_inner_indices[ground_index],
                    );
                    transport_values[value_index] +=
                        self.ground_weights[ground_index] * current_attenuation;
                }
            }
        }
        ray_end_attenuation[ray_index] = current_attenuation;
        if let Some((forcing, _)) = first_order {
            let forcing_start = ray_index * self.num_stokes;
            if self.num_stokes == 1 {
                ray_first_order[0] +=
                    current_attenuation * forcing.end_of_ray_source[forcing_start];
            } else {
                for (stokes, ray_value) in ray_first_order[..self.num_stokes].iter_mut().enumerate()
                {
                    *ray_value +=
                        current_attenuation * forcing.end_of_ray_source[forcing_start + stokes];
                }
            }
            first_order_radiance[forcing_start..forcing_start + self.num_stokes]
                .copy_from_slice(&ray_first_order[..self.num_stokes]);
        }
    }

    fn validate_first_order(&self) -> Result<()> {
        if self.num_stokes != 1 && self.num_stokes != 3 {
            bail!("first-order forcing requires one or three Stokes components");
        }
        if self.ray_terminal_endpoint_indices.len() != self.num_rays()
            || self.layer_exit_endpoint_indices.len() != self.num_layers()
            || (!self.layer_exit_phase_endpoint_slots.is_empty()
                && (self.layer_exit_phase_endpoint_slots.len() != self.num_layers()
                    || self
                        .layer_exit_phase_endpoint_slots
                        .iter()
                        .any(|&slot| slot as usize >= self.phase_endpoint_cache_size)
                    || self.ray_terminal_phase_endpoint_slots.len() != self.num_rays()))
            || (self.layer_exit_phase_endpoint_slots.is_empty()
                && !self.ray_terminal_phase_endpoint_slots.is_empty())
            || self.endpoint_atmosphere_offsets.len() != self.num_unique_endpoint_stencils() + 1
            || self.endpoint_atmosphere_indices.len() != self.endpoint_weights.len()
            || self.layer_distance.len() != self.num_layers()
            || self.layer_start_fraction.len() != self.num_layers()
            || self.phase_geometry_ray_indices.len() != self.num_rays()
            || self.phase_geometry_endpoint_offsets.len() != self.num_phase_geometries() + 1
            || self
                .phase_geometry_endpoint_indices
                .iter()
                .any(|&endpoint| endpoint as usize >= self.num_unique_endpoint_stencils())
            || self.phase_geometry_endpoint_evaluations.len() != self.num_phase_geometries()
            || self.ray_phase_basis.len() != self.num_phase_geometries() * self.num_phase_moments
        {
            bail!("first-order ray geometry is incomplete");
        }
        if self.num_stokes == 3
            && (self.ray_polarized_phase_basis.len()
                != self.num_phase_geometries() * self.num_phase_moments
                || self.phase_q_factor.len() != self.num_phase_geometries()
                || self.phase_u_factor.len() != self.num_phase_geometries())
        {
            bail!("polarized first-order ray geometry is incomplete");
        }
        Ok(())
    }

    #[inline]
    fn phase_coefficient_families(&self) -> usize {
        if self.num_stokes == 1 { 1 } else { 4 }
    }

    #[inline]
    fn phase_geometry_rays(&self, phase_geometry: usize) -> &[u32] {
        let start = self.phase_geometry_ray_offsets[phase_geometry] as usize;
        let end = self.phase_geometry_ray_offsets[phase_geometry + 1] as usize;
        &self.phase_geometry_ray_indices[start..end]
    }

    #[inline]
    fn phase_geometry_endpoints(&self, phase_geometry: usize) -> &[u32] {
        let start = self.phase_geometry_endpoint_offsets[phase_geometry] as usize;
        let end = self.phase_geometry_endpoint_offsets[phase_geometry + 1] as usize;
        &self.phase_geometry_endpoint_indices[start..end]
    }

    #[inline]
    fn uses_grouped_phase_geometry(
        &self,
        phase_geometry: usize,
        num_atmosphere_points: usize,
    ) -> bool {
        num_atmosphere_points <= self.phase_geometry_endpoint_evaluations[phase_geometry]
    }

    fn endpoint_phase_direct(
        &self,
        phase_geometry: usize,
        atmosphere_indices: &[u32],
        endpoint_weights: &[f64],
        legendre_coefficients: &[f64],
        maximum_order: &[i32],
    ) -> [f64; 3] {
        let basis_start = phase_geometry * self.num_phase_moments;
        let intensity_basis =
            &self.ray_phase_basis[basis_start..basis_start + self.num_phase_moments];
        let polarized_basis = if self.num_stokes == 3 {
            &self.ray_polarized_phase_basis[basis_start..basis_start + self.num_phase_moments]
        } else {
            &[]
        };
        let coefficient_families = self.phase_coefficient_families();
        let mut phase = [0.0; 3];
        for (&atmosphere_index, &weight) in atmosphere_indices.iter().zip(endpoint_weights) {
            if weight == 0.0 {
                continue;
            }
            let atmosphere_index = atmosphere_index as usize;
            let order = maximum_order[atmosphere_index] as usize;
            let coefficient_start =
                atmosphere_index * self.num_phase_moments * coefficient_families;
            if self.num_stokes == 1 {
                let mut intensity = 0.0;
                for degree in 0..order {
                    intensity +=
                        legendre_coefficients[coefficient_start + degree] * intensity_basis[degree];
                }
                phase[0] += weight * intensity;
            } else {
                let mut intensity = 0.0;
                let mut polarized = 0.0;
                for degree in 0..order {
                    intensity += legendre_coefficients[coefficient_start + degree * 4]
                        * intensity_basis[degree];
                    polarized += legendre_coefficients[coefficient_start + degree * 4 + 3]
                        * polarized_basis[degree];
                }
                phase[0] += weight * intensity;
                phase[1] -= weight * self.phase_q_factor[phase_geometry] * polarized;
                phase[2] -= weight * self.phase_u_factor[phase_geometry] * polarized;
            }
        }
        phase
    }

    #[allow(clippy::too_many_arguments)]
    fn accumulate_endpoint_phase_vjp_direct(
        &self,
        phase_geometry: usize,
        atmosphere_index: usize,
        weight: f64,
        maximum_order: &[i32],
        phase_cotangent: [f64; 3],
        legendre_coefficient_gradient: &mut [f64],
    ) {
        let basis_start = phase_geometry * self.num_phase_moments;
        let intensity_basis =
            &self.ray_phase_basis[basis_start..basis_start + self.num_phase_moments];
        let order = maximum_order[atmosphere_index] as usize;
        let coefficient_start =
            atmosphere_index * self.num_phase_moments * self.phase_coefficient_families();
        if self.num_stokes == 1 {
            let cotangent = weight * phase_cotangent[0];
            for degree in 0..order {
                legendre_coefficient_gradient[coefficient_start + degree] +=
                    cotangent * intensity_basis[degree];
            }
        } else {
            let polarized_basis =
                &self.ray_polarized_phase_basis[basis_start..basis_start + self.num_phase_moments];
            let intensity_cotangent = weight * phase_cotangent[0];
            let polarized_cotangent = -weight
                * (self.phase_q_factor[phase_geometry] * phase_cotangent[1]
                    + self.phase_u_factor[phase_geometry] * phase_cotangent[2]);
            for degree in 0..order {
                legendre_coefficient_gradient[coefficient_start + degree * 4] +=
                    intensity_cotangent * intensity_basis[degree];
                legendre_coefficient_gradient[coefficient_start + degree * 4 + 3] +=
                    polarized_cotangent * polarized_basis[degree];
            }
        }
    }

    fn fill_phase_values(
        &self,
        phase_geometry: usize,
        legendre_coefficients: &[f64],
        maximum_order: &[i32],
        phase_values: &mut [f64],
    ) {
        debug_assert_eq!(phase_values.len(), maximum_order.len() * self.num_stokes);
        let basis_start = phase_geometry * self.num_phase_moments;
        let intensity_basis =
            &self.ray_phase_basis[basis_start..basis_start + self.num_phase_moments];
        let polarized_basis = if self.num_stokes == 3 {
            &self.ray_polarized_phase_basis[basis_start..basis_start + self.num_phase_moments]
        } else {
            &[]
        };
        for (atmosphere_index, &order) in maximum_order.iter().enumerate() {
            let order = order as usize;
            let coefficient_start =
                atmosphere_index * self.num_phase_moments * self.phase_coefficient_families();
            let output_start = atmosphere_index * self.num_stokes;
            if self.num_stokes == 1 {
                phase_values[output_start] = legendre_coefficients
                    [coefficient_start..coefficient_start + order]
                    .iter()
                    .zip(&intensity_basis[..order])
                    .map(|(coefficient, basis)| coefficient * basis)
                    .sum();
            } else {
                let mut intensity = 0.0;
                let mut polarized = 0.0;
                for degree in 0..order {
                    intensity += legendre_coefficients[coefficient_start + degree * 4]
                        * intensity_basis[degree];
                    polarized += legendre_coefficients[coefficient_start + degree * 4 + 3]
                        * polarized_basis[degree];
                }
                phase_values[output_start] = intensity;
                phase_values[output_start + 1] = -self.phase_q_factor[phase_geometry] * polarized;
                phase_values[output_start + 2] = -self.phase_u_factor[phase_geometry] * polarized;
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn fill_phase_values_jvp(
        &self,
        phase_geometry: usize,
        legendre_coefficients: &[f64],
        legendre_coefficient_tangent: &[f64],
        maximum_order: &[i32],
        phase_values: &mut [f64],
        phase_value_tangent: &mut [f64],
    ) {
        self.fill_phase_values(
            phase_geometry,
            legendre_coefficients,
            maximum_order,
            phase_values,
        );
        self.fill_phase_values(
            phase_geometry,
            legendre_coefficient_tangent,
            maximum_order,
            phase_value_tangent,
        );
    }

    fn accumulate_phase_values_vjp(
        &self,
        phase_geometry: usize,
        maximum_order: &[i32],
        phase_value_gradient: &[f64],
        legendre_coefficient_gradient: &mut [f64],
    ) {
        let basis_start = phase_geometry * self.num_phase_moments;
        let intensity_basis =
            &self.ray_phase_basis[basis_start..basis_start + self.num_phase_moments];
        let polarized_basis = if self.num_stokes == 3 {
            &self.ray_polarized_phase_basis[basis_start..basis_start + self.num_phase_moments]
        } else {
            &[]
        };
        for (atmosphere_index, &order) in maximum_order.iter().enumerate() {
            let order = order as usize;
            let coefficient_start =
                atmosphere_index * self.num_phase_moments * self.phase_coefficient_families();
            let gradient_start = atmosphere_index * self.num_stokes;
            if self.num_stokes == 1 {
                let phase_cotangent = phase_value_gradient[gradient_start];
                for degree in 0..order {
                    legendre_coefficient_gradient[coefficient_start + degree] +=
                        phase_cotangent * intensity_basis[degree];
                }
            } else {
                let polarized_cotangent = -self.phase_q_factor[phase_geometry]
                    * phase_value_gradient[gradient_start + 1]
                    - self.phase_u_factor[phase_geometry]
                        * phase_value_gradient[gradient_start + 2];
                for degree in 0..order {
                    legendre_coefficient_gradient[coefficient_start + degree * 4] +=
                        phase_value_gradient[gradient_start] * intensity_basis[degree];
                    legendre_coefficient_gradient[coefficient_start + degree * 4 + 3] +=
                        polarized_cotangent * polarized_basis[degree];
                }
            }
        }
    }

    #[inline]
    fn transport_value_index(
        &self,
        ray_index: usize,
        transport_offset: usize,
        stokes: usize,
        inner_index: u16,
    ) -> usize {
        transport_offset
            + stokes * self.ray_transport_row_nnz[ray_index] as usize
            + inner_index as usize
    }

    #[cfg(test)]
    #[inline(always)]
    fn endpoint_index(&self, ray_index: usize, layer_index: usize, entrance: bool) -> usize {
        if entrance {
            debug_assert_eq!(
                layer_index + 1,
                self.ray_layer_offsets[ray_index + 1] as usize
            );
            self.ray_terminal_endpoint_indices[ray_index] as usize
        } else {
            self.layer_exit_endpoint_indices[layer_index] as usize
        }
    }

    #[inline(always)]
    fn endpoint_location(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
    ) -> (usize, usize) {
        if entrance {
            debug_assert_eq!(
                layer_index + 1,
                self.ray_layer_offsets[ray_index + 1] as usize
            );
            let endpoint_index = self.ray_terminal_endpoint_indices[ray_index] as usize;
            let endpoint_slot = if self.layer_exit_phase_endpoint_slots.is_empty() {
                endpoint_index
            } else {
                self.ray_terminal_phase_endpoint_slots[ray_index] as usize
            };
            (endpoint_index, endpoint_slot)
        } else {
            let endpoint_index = self.layer_exit_endpoint_indices[layer_index] as usize;
            let endpoint_slot = if self.layer_exit_phase_endpoint_slots.is_empty() {
                endpoint_index
            } else {
                self.layer_exit_phase_endpoint_slots[layer_index] as usize
            };
            (endpoint_index, endpoint_slot)
        }
    }

    #[inline]
    fn phase_endpoint_cache_slot(&self, endpoint_slot: usize, endpoint_index: usize) -> usize {
        if self.layer_exit_phase_endpoint_slots.is_empty() {
            endpoint_index
        } else {
            endpoint_slot
        }
    }

    #[inline]
    fn endpoint_stencil_by_index(&self, endpoint_index: usize) -> (&[u32], &[f64]) {
        let start = self.endpoint_atmosphere_offsets[endpoint_index] as usize;
        let end = self.endpoint_atmosphere_offsets[endpoint_index + 1] as usize;
        (
            &self.endpoint_atmosphere_indices[start..end],
            &self.endpoint_weights[start..end],
        )
    }

    #[cfg(test)]
    fn endpoint_stencil(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
    ) -> (&[u32], &[f64]) {
        self.endpoint_stencil_by_index(self.endpoint_index(ray_index, layer_index, entrance))
    }

    fn interpolate_endpoint_medium(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        solar_transmission: &[f64],
    ) -> Vec<EndpointMedium> {
        let mut result = vec![EndpointMedium::default(); self.num_unique_endpoint_stencils()];
        for (endpoint_index, medium) in result.iter_mut().enumerate() {
            let (atmosphere_indices, weights) = self.endpoint_stencil_by_index(endpoint_index);
            for (&atmosphere_index, &weight) in atmosphere_indices.iter().zip(weights) {
                let atmosphere_index = atmosphere_index as usize;
                medium.extinction += weight * extinction[atmosphere_index];
                medium.albedo += weight * single_scatter_albedo[atmosphere_index];
                if self.solar_transmission_on_atmosphere_grid {
                    medium.solar_transmission += weight * solar_transmission[atmosphere_index];
                }
            }
        }
        result
    }

    fn endpoint_phase_by_index(
        &self,
        phase_geometry: usize,
        endpoint_index: usize,
        legendre_coefficients: &[f64],
        maximum_order: &[i32],
        phase_values: &[f64],
    ) -> [f64; 3] {
        let (atmosphere_indices, weights) = self.endpoint_stencil_by_index(endpoint_index);
        if phase_values.is_empty() {
            return self.endpoint_phase_direct(
                phase_geometry,
                atmosphere_indices,
                weights,
                legendre_coefficients,
                maximum_order,
            );
        }

        let mut phase = [0.0; 3];
        if self.num_stokes == 1 {
            for (&atmosphere_index, &weight) in atmosphere_indices.iter().zip(weights) {
                phase[0] += weight * phase_values[atmosphere_index as usize];
            }
        } else {
            for (&atmosphere_index, &weight) in atmosphere_indices.iter().zip(weights) {
                let phase_start = atmosphere_index as usize * self.num_stokes;
                for stokes in 0..self.num_stokes {
                    phase[stokes] += weight * phase_values[phase_start + stokes];
                }
            }
        }
        phase
    }

    #[allow(clippy::too_many_arguments)]
    fn fill_endpoint_scattering_values(
        &self,
        phase_geometry: usize,
        legendre_coefficients: &[f64],
        maximum_order: &[i32],
        phase_values: &[f64],
        endpoint_medium: &[EndpointMedium],
        endpoint_scattering_values: &mut [f64],
    ) {
        let inverse_sphere = 1.0 / (4.0 * std::f64::consts::PI);
        for (endpoint_slot, &endpoint_index) in self
            .phase_geometry_endpoints(phase_geometry)
            .iter()
            .enumerate()
        {
            let endpoint_index = endpoint_index as usize;
            let phase = self.endpoint_phase_by_index(
                phase_geometry,
                endpoint_index,
                legendre_coefficients,
                maximum_order,
                phase_values,
            );
            let medium = endpoint_medium[endpoint_index];
            let scale = medium.extinction * medium.albedo * inverse_sphere;
            let output_start =
                self.phase_endpoint_cache_slot(endpoint_slot, endpoint_index) * self.num_stokes;
            if self.num_stokes == 1 {
                endpoint_scattering_values[output_start] = scale * phase[0];
            } else {
                for stokes in 0..self.num_stokes {
                    endpoint_scattering_values[output_start + stokes] = scale * phase[stokes];
                }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn fill_endpoint_scattering_values_jvp(
        &self,
        phase_geometry: usize,
        forcing: &FirstOrderJvpInputs<'_>,
        phase_values: &[f64],
        phase_value_tangent: &[f64],
        endpoint_medium: &[EndpointMediumJvp],
        endpoint_scattering_values: &mut [f64],
        endpoint_scattering_value_tangent: &mut [f64],
    ) {
        let inverse_sphere = 1.0 / (4.0 * std::f64::consts::PI);
        for (endpoint_slot, &endpoint_index) in self
            .phase_geometry_endpoints(phase_geometry)
            .iter()
            .enumerate()
        {
            let endpoint_index = endpoint_index as usize;
            let phase = self.endpoint_phase_by_index(
                phase_geometry,
                endpoint_index,
                forcing.legendre_coefficients,
                forcing.maximum_order,
                phase_values,
            );
            let phase_tangent = self.endpoint_phase_by_index(
                phase_geometry,
                endpoint_index,
                forcing.legendre_coefficient_tangent,
                forcing.maximum_order,
                phase_value_tangent,
            );
            let medium = endpoint_medium[endpoint_index];
            let scale = medium.extinction.value * medium.albedo.value * inverse_sphere;
            let scale_tangent = (medium.extinction.tangent * medium.albedo.value
                + medium.extinction.value * medium.albedo.tangent)
                * inverse_sphere;
            let output_start =
                self.phase_endpoint_cache_slot(endpoint_slot, endpoint_index) * self.num_stokes;
            if self.num_stokes == 1 {
                endpoint_scattering_values[output_start] = scale * phase[0];
                endpoint_scattering_value_tangent[output_start] =
                    scale_tangent * phase[0] + scale * phase_tangent[0];
            } else {
                for stokes in 0..self.num_stokes {
                    endpoint_scattering_values[output_start + stokes] = scale * phase[stokes];
                    endpoint_scattering_value_tangent[output_start + stokes] =
                        scale_tangent * phase[stokes] + scale * phase_tangent[stokes];
                }
            }
        }
    }

    fn fill_endpoint_phase_cache(
        &self,
        phase_geometry: usize,
        forcing: &FirstOrderInputs<'_>,
        phase_values: &[f64],
        endpoint_phase_values: &mut [f64],
    ) {
        for (endpoint_slot, &endpoint_index) in self
            .phase_geometry_endpoints(phase_geometry)
            .iter()
            .enumerate()
        {
            let endpoint_index = endpoint_index as usize;
            let phase = self.endpoint_phase_by_index(
                phase_geometry,
                endpoint_index,
                forcing.legendre_coefficients,
                forcing.maximum_order,
                phase_values,
            );
            let output_start =
                self.phase_endpoint_cache_slot(endpoint_slot, endpoint_index) * self.num_stokes;
            endpoint_phase_values[output_start..output_start + self.num_stokes]
                .copy_from_slice(&phase[..self.num_stokes]);
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn interpolate_endpoint_medium_jvp(
        &self,
        extinction: &[f64],
        single_scatter_albedo: &[f64],
        solar_transmission: &[f64],
        extinction_tangent: &[f64],
        single_scatter_albedo_tangent: &[f64],
        solar_transmission_tangent: &[f64],
    ) -> Vec<EndpointMediumJvp> {
        let mut result = vec![EndpointMediumJvp::default(); self.num_unique_endpoint_stencils()];
        for (endpoint_index, medium) in result.iter_mut().enumerate() {
            let (atmosphere_indices, weights) = self.endpoint_stencil_by_index(endpoint_index);
            for (&atmosphere_index, &weight) in atmosphere_indices.iter().zip(weights) {
                let atmosphere_index = atmosphere_index as usize;
                medium.extinction.value += weight * extinction[atmosphere_index];
                medium.extinction.tangent += weight * extinction_tangent[atmosphere_index];
                medium.albedo.value += weight * single_scatter_albedo[atmosphere_index];
                medium.albedo.tangent += weight * single_scatter_albedo_tangent[atmosphere_index];
                if self.solar_transmission_on_atmosphere_grid {
                    medium.solar_transmission.value +=
                        weight * solar_transmission[atmosphere_index];
                    medium.solar_transmission.tangent +=
                        weight * solar_transmission_tangent[atmosphere_index];
                }
            }
        }
        result
    }

    #[allow(clippy::too_many_arguments)]
    fn endpoint_source(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
        forcing: &FirstOrderInputs<'_>,
        endpoint_medium: &[EndpointMedium],
        endpoint_scattering_values: &[f64],
    ) -> [f64; 3] {
        let (endpoint_index, endpoint_slot) =
            self.endpoint_location(ray_index, layer_index, entrance);
        let solar_transmission = if self.solar_transmission_on_atmosphere_grid {
            endpoint_medium[endpoint_index].solar_transmission
        } else {
            let exit_solar_index = layer_index + ray_index;
            forcing.solar_transmission[exit_solar_index + usize::from(entrance)]
        };
        let source_start = endpoint_slot * self.num_stokes;
        let mut source = [0.0; 3];
        if self.num_stokes == 1 {
            source[0] = endpoint_scattering_values[source_start] * solar_transmission;
        } else {
            for stokes in 0..self.num_stokes {
                source[stokes] =
                    endpoint_scattering_values[source_start + stokes] * solar_transmission;
            }
        }
        source
    }

    #[allow(clippy::too_many_arguments)]
    fn endpoint_primal(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
        forcing: &FirstOrderInputs<'_>,
        endpoint_medium: &[EndpointMedium],
        endpoint_phase_values: &[f64],
    ) -> [f64; 3] {
        let (endpoint_index, endpoint_slot) =
            self.endpoint_location(ray_index, layer_index, entrance);
        let medium = endpoint_medium[endpoint_index];
        let phase_start = endpoint_slot * self.num_stokes;
        let mut endpoint_phase = [0.0; 3];
        endpoint_phase[..self.num_stokes]
            .copy_from_slice(&endpoint_phase_values[phase_start..phase_start + self.num_stokes]);
        let endpoint_solar_transmission = if self.solar_transmission_on_atmosphere_grid {
            medium.solar_transmission
        } else {
            // Every ray contributes one more endpoint than layer. Since flat
            // layer indices already include preceding rays' layers, adding the
            // ray index gives the matching endpoint offset.
            let exit_solar_index = layer_index + ray_index;
            forcing.solar_transmission[exit_solar_index + usize::from(entrance)]
        };

        let amplitude = medium.extinction * medium.albedo * endpoint_solar_transmission
            / (4.0 * std::f64::consts::PI);
        let mut source = [0.0; 3];
        if self.num_stokes == 1 {
            source[0] = amplitude * endpoint_phase[0];
        } else {
            for stokes in 0..self.num_stokes {
                source[stokes] = amplitude * endpoint_phase[stokes];
            }
        }
        source
    }

    #[allow(clippy::too_many_arguments)]
    fn accumulate_endpoint_source_cotangent(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
        source_cotangent: [f64; 3],
        solar_transmission: &[f64],
        endpoint_medium: &[EndpointMedium],
        endpoint_phase_values: &[f64],
        endpoint_weighted_source_cotangent: &mut [f64],
        endpoint_source_cotangent: &mut [f64],
        solar_transmission_gradient: &mut [f64],
    ) {
        if source_cotangent[..self.num_stokes]
            .iter()
            .all(|value| *value == 0.0)
        {
            return;
        }
        let (endpoint_index, endpoint_slot) =
            self.endpoint_location(ray_index, layer_index, entrance);
        let endpoint_start = endpoint_slot * self.num_stokes;
        let endpoint_solar_transmission = if self.solar_transmission_on_atmosphere_grid {
            endpoint_medium[endpoint_index].solar_transmission
        } else {
            let exit_solar_index = layer_index + ray_index;
            solar_transmission[exit_solar_index + usize::from(entrance)]
        };
        for stokes in 0..self.num_stokes {
            endpoint_weighted_source_cotangent[endpoint_start + stokes] +=
                source_cotangent[stokes] * endpoint_solar_transmission;
            if self.solar_transmission_on_atmosphere_grid {
                endpoint_source_cotangent[endpoint_start + stokes] += source_cotangent[stokes];
            }
        }
        if !self.solar_transmission_on_atmosphere_grid {
            let medium = endpoint_medium[endpoint_index];
            let amplitude_cotangent = if self.num_stokes == 1 {
                source_cotangent[0] * endpoint_phase_values[endpoint_start]
            } else {
                source_cotangent[..self.num_stokes]
                    .iter()
                    .zip(&endpoint_phase_values[endpoint_start..endpoint_start + self.num_stokes])
                    .map(|(cotangent, phase)| cotangent * phase)
                    .sum()
            };
            let solar_cotangent = amplitude_cotangent * medium.extinction * medium.albedo
                / (4.0 * std::f64::consts::PI);
            let exit_solar_index = layer_index + ray_index;
            solar_transmission_gradient[exit_solar_index + usize::from(entrance)] +=
                solar_cotangent;
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_endpoint_source_vjp(
        &self,
        phase_geometry: usize,
        maximum_order: &[i32],
        endpoint_medium: &[EndpointMedium],
        endpoint_phase_values: &[f64],
        endpoint_weighted_source_cotangent: &[f64],
        endpoint_source_cotangent: &[f64],
        endpoint_extinction_gradient: &mut [f64],
        endpoint_albedo_gradient: &mut [f64],
        endpoint_phase_gradient: &mut [f64],
        endpoint_solar_gradient: &mut [f64],
        legendre_coefficient_gradient: &mut [f64],
    ) {
        let inverse_sphere = 1.0 / (4.0 * std::f64::consts::PI);
        for (endpoint_slot, &endpoint_index) in self
            .phase_geometry_endpoints(phase_geometry)
            .iter()
            .enumerate()
        {
            let endpoint_index = endpoint_index as usize;
            let endpoint_start =
                self.phase_endpoint_cache_slot(endpoint_slot, endpoint_index) * self.num_stokes;
            let weighted_cotangent = &endpoint_weighted_source_cotangent
                [endpoint_start..endpoint_start + self.num_stokes];
            let source_cotangent = if self.solar_transmission_on_atmosphere_grid {
                &endpoint_source_cotangent[endpoint_start..endpoint_start + self.num_stokes]
            } else {
                &[]
            };
            if weighted_cotangent.iter().all(|value| *value == 0.0)
                && source_cotangent.iter().all(|value| *value == 0.0)
            {
                continue;
            }
            let phase = &endpoint_phase_values[endpoint_start..endpoint_start + self.num_stokes];
            let amplitude_cotangent: f64 = weighted_cotangent
                .iter()
                .zip(phase)
                .map(|(cotangent, phase)| cotangent * phase)
                .sum();
            let medium = endpoint_medium[endpoint_index];
            endpoint_extinction_gradient[endpoint_index] +=
                amplitude_cotangent * medium.albedo * inverse_sphere;
            endpoint_albedo_gradient[endpoint_index] +=
                amplitude_cotangent * medium.extinction * inverse_sphere;
            let phase_scale = medium.extinction * medium.albedo * inverse_sphere;
            let mut phase_cotangent = [0.0; 3];
            for stokes in 0..self.num_stokes {
                phase_cotangent[stokes] = weighted_cotangent[stokes] * phase_scale;
            }
            if endpoint_phase_gradient.is_empty() {
                let (atmosphere_indices, weights) = self.endpoint_stencil_by_index(endpoint_index);
                for (&atmosphere_index, &weight) in atmosphere_indices.iter().zip(weights) {
                    self.accumulate_endpoint_phase_vjp_direct(
                        phase_geometry,
                        atmosphere_index as usize,
                        weight,
                        maximum_order,
                        phase_cotangent,
                        legendre_coefficient_gradient,
                    );
                }
            } else {
                for stokes in 0..self.num_stokes {
                    endpoint_phase_gradient[endpoint_start + stokes] += phase_cotangent[stokes];
                }
            }
            if self.solar_transmission_on_atmosphere_grid {
                let solar_amplitude_cotangent: f64 = source_cotangent
                    .iter()
                    .zip(phase)
                    .map(|(cotangent, phase)| cotangent * phase)
                    .sum();
                endpoint_solar_gradient[endpoint_index] += solar_amplitude_cotangent * phase_scale;
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn endpoint_source_jvp(
        &self,
        ray_index: usize,
        layer_index: usize,
        entrance: bool,
        forcing: &FirstOrderJvpInputs<'_>,
        endpoint_medium: &[EndpointMediumJvp],
        endpoint_scattering_values: &[f64],
        endpoint_scattering_value_tangent: &[f64],
    ) -> EndpointValueTangent {
        let (endpoint_index, endpoint_slot) =
            self.endpoint_location(ray_index, layer_index, entrance);
        let medium = endpoint_medium[endpoint_index];
        let endpoint_solar_transmission = if self.solar_transmission_on_atmosphere_grid {
            medium.solar_transmission
        } else {
            let exit_solar_index = layer_index + ray_index;
            let solar_index = exit_solar_index + usize::from(entrance);
            ValueTangent {
                value: forcing.solar_transmission[solar_index],
                tangent: forcing.solar_transmission_tangent[solar_index],
            }
        };
        let source_start = endpoint_slot * self.num_stokes;
        let mut result = EndpointValueTangent::default();
        if self.num_stokes == 1 {
            let scattering_value = endpoint_scattering_values[source_start];
            result.value[0] = scattering_value * endpoint_solar_transmission.value;
            result.tangent[0] = endpoint_scattering_value_tangent[source_start]
                * endpoint_solar_transmission.value
                + scattering_value * endpoint_solar_transmission.tangent;
        } else {
            for stokes in 0..self.num_stokes {
                let scattering_value = endpoint_scattering_values[source_start + stokes];
                result.value[stokes] = scattering_value * endpoint_solar_transmission.value;
                result.tangent[stokes] = endpoint_scattering_value_tangent[source_start + stokes]
                    * endpoint_solar_transmission.value
                    + scattering_value * endpoint_solar_transmission.tangent;
            }
        }
        result
    }
}

struct FirstOrderInputs<'a> {
    legendre_coefficients: &'a [f64],
    maximum_order: &'a [i32],
    solar_transmission: &'a [f64],
    end_of_ray_source: &'a [f64],
}

struct FirstOrderJvpInputs<'a> {
    legendre_coefficients: &'a [f64],
    maximum_order: &'a [i32],
    solar_transmission: &'a [f64],
    legendre_coefficient_tangent: &'a [f64],
    solar_transmission_tangent: &'a [f64],
    end_of_ray_source: &'a [f64],
    end_of_ray_source_tangent: &'a [f64],
}

#[derive(Clone, Copy, Default)]
struct EndpointMedium {
    extinction: f64,
    albedo: f64,
    solar_transmission: f64,
}

#[derive(Clone, Copy, Default)]
struct ValueTangent {
    value: f64,
    tangent: f64,
}

#[derive(Clone, Copy, Default)]
struct EndpointMediumJvp {
    extinction: ValueTangent,
    albedo: ValueTangent,
    solar_transmission: ValueTangent,
}

#[derive(Clone, Copy, Default)]
struct EndpointValueTangent {
    value: [f64; 3],
    tangent: [f64; 3],
}

#[inline]
fn dot_stokes(left: &[f64; 3], right: &[f64; 3], num_stokes: usize) -> f64 {
    left[..num_stokes]
        .iter()
        .zip(&right[..num_stokes])
        .map(|(left, right)| left * right)
        .sum()
}

#[inline(always)]
fn add_packed_value(values: &mut [f64], index: usize, increment: f64) {
    debug_assert!(index < values.len());
    // SAFETY: every public assembly entry point checks the buffer against
    // `required_transport_len`, which is the maximum packed index compiled
    // into this immutable geometry.
    unsafe {
        *values.get_unchecked_mut(index) += increment;
    }
}

#[inline(always)]
fn packed_value(values: &[f64], index: usize) -> f64 {
    debug_assert!(index < values.len());
    // SAFETY: callers validate the corresponding packed buffer against
    // `required_transport_len` before entering the hot ray traversal.
    unsafe { *values.get_unchecked(index) }
}

#[inline(always)]
fn packed_column(columns: &[i32], index: usize) -> usize {
    debug_assert!(index < columns.len());
    // SAFETY: compact VJP validation checks the full required column range
    // before any ray traversal and rejects negative column indices.
    unsafe { *columns.get_unchecked(index) as usize }
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

fn constant_source_factor_and_derivative(optical_depth: f64, attenuation: f64) -> (f64, f64) {
    if optical_depth.abs() < 1.0e-4 {
        let optical_depth_squared = optical_depth * optical_depth;
        let optical_depth_cubed = optical_depth_squared * optical_depth;
        let value = 1.0 - optical_depth / 2.0 + optical_depth_squared / 6.0
            - optical_depth_cubed / 24.0
            + optical_depth_squared * optical_depth_squared / 120.0;
        let derivative =
            -0.5 + optical_depth / 3.0 - optical_depth_squared / 8.0 + optical_depth_cubed / 30.0;
        return (value, derivative);
    }
    let inverse_optical_depth = 1.0 / optical_depth;
    let value = (1.0 - attenuation) * inverse_optical_depth;
    (
        value,
        inverse_optical_depth - value * (1.0 + inverse_optical_depth),
    )
}

fn intern_endpoint_stencil(
    atmosphere_indices: &[u32],
    weights: &[f64],
    endpoint_map: &mut HashMap<EndpointStencilKey, u32>,
    endpoint_offsets: &mut Vec<u32>,
    endpoint_atmosphere_indices: &mut Vec<u32>,
    endpoint_weights: &mut Vec<f64>,
) -> Result<u32> {
    debug_assert_eq!(atmosphere_indices.len(), weights.len());
    let key = EndpointStencilKey::new(atmosphere_indices, weights);
    if let Some(&endpoint_index) = endpoint_map.get(&key) {
        return Ok(endpoint_index);
    }

    let endpoint_index = u32::try_from(endpoint_offsets.len() - 1)?;
    endpoint_atmosphere_indices.extend_from_slice(atmosphere_indices);
    endpoint_weights.extend_from_slice(weights);
    endpoint_offsets.push(u32::try_from(endpoint_weights.len())?);
    endpoint_map.insert(key, endpoint_index);
    Ok(endpoint_index)
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
            &[],
            &[],
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

    fn compact_phase_cache_geometry() -> ScalarRayTransport {
        ScalarRayTransport::new(
            &[0, 1, 2, 3, 4, 5],
            &[0, 1, 2, 3, 4, 5],
            &[0, 1, 2, 3, 4],
            &[0.2, 0.25, 0.3, 0.35, 0.4],
            &[1.0; 5],
            &[1.0; 5],
            &[1.0; 5],
            &[0.5; 5],
            &[0.5; 5],
            &[0.5; 5],
            &[-0.8, -0.4, 0.0, 0.4, 0.8],
            &[],
            &[],
            2,
            true,
            &[0, 1, 2, 3, 4, 5],
            &[0, 1, 2, 3, 4],
            &[0; 5],
            &[1.0; 5],
            &[0; 6],
            &[],
            &[],
        )
        .unwrap()
    }

    #[test]
    fn compact_phase_cache_products_are_dual() {
        let geometry = compact_phase_cache_geometry();
        assert_eq!(geometry.phase_endpoint_cache_size, 1);
        assert!(!geometry.layer_exit_phase_endpoint_slots.is_empty());

        let extinction = [0.8, 1.0, 1.2, 1.4, 1.6];
        let albedo = [0.3, 0.4, 0.5, 0.6, 0.7];
        let coefficients = [1.0, 0.1, 1.0, 0.2, 1.0, 0.3, 1.0, 0.4, 1.0, 0.5];
        let maximum_order = [2; 5];
        let solar = [0.9, 0.85, 0.8, 0.75, 0.7];
        let end_source = [0.01, 0.02, 0.03, 0.04, 0.05];
        let mut transport = [0.0; 5];
        let mut first_order = [0.0; 5];
        let mut optical_depth = [0.0; 5];
        let mut attenuation = [0.0; 5];
        let mut prefix = [0.0; 5];
        let mut ray_end = [0.0; 5];
        geometry
            .assemble_with_first_order(
                &extinction,
                &albedo,
                &coefficients,
                &maximum_order,
                &solar,
                &end_source,
                &mut transport,
                &mut first_order,
                &mut optical_depth,
                &mut attenuation,
                &mut prefix,
                &mut ray_end,
            )
            .unwrap();

        let extinction_tangent = [0.01, -0.02, 0.03, -0.04, 0.05];
        let albedo_tangent = [-0.02, 0.01, -0.01, 0.02, -0.03];
        let coefficient_tangent = [
            0.01, -0.01, 0.02, -0.02, 0.03, -0.03, 0.04, -0.04, 0.05, -0.05,
        ];
        let solar_tangent = [0.02, -0.01, 0.03, -0.02, 0.01];
        let end_source_tangent = [-0.01, 0.02, -0.03, 0.04, -0.05];
        let mut transport_tangent = [0.0; 5];
        let mut first_order_tangent = [0.0; 5];
        let mut ray_end_tangent = [0.0; 5];
        geometry
            .assemble_jvp_with_first_order(
                &extinction,
                &albedo,
                &coefficients,
                &maximum_order,
                &solar,
                &extinction_tangent,
                &albedo_tangent,
                &coefficient_tangent,
                &solar_tangent,
                &end_source,
                &end_source_tangent,
                &mut transport,
                &mut transport_tangent,
                &mut first_order_tangent,
                &mut ray_end,
                &mut ray_end_tangent,
            )
            .unwrap();

        let transport_cotangent = [0.2, -0.3, 0.4, -0.5, 0.6];
        let first_order_cotangent = [-0.4, 0.3, -0.2, 0.1, 0.5];
        let ray_end_cotangent = [0.1, -0.2, 0.3, -0.4, 0.5];
        let mut extinction_gradient = [0.0; 5];
        let mut albedo_gradient = [0.0; 5];
        let mut coefficient_gradient = [0.0; 10];
        let mut solar_gradient = [0.0; 5];
        let mut end_source_gradient = [0.0; 5];
        geometry
            .assemble_vjp_with_first_order(
                &extinction,
                &albedo,
                &coefficients,
                &maximum_order,
                &solar,
                &transport_cotangent,
                &[],
                &[],
                &first_order_cotangent,
                &ray_end,
                &ray_end_cotangent,
                &end_source,
                &optical_depth,
                &attenuation,
                &prefix,
                &mut extinction_gradient,
                &mut albedo_gradient,
                &mut coefficient_gradient,
                &mut solar_gradient,
                &mut end_source_gradient,
            )
            .unwrap();

        let output_dot = transport_tangent
            .iter()
            .zip(transport_cotangent)
            .map(|(tangent, cotangent)| tangent * cotangent)
            .sum::<f64>()
            + first_order_tangent
                .iter()
                .zip(first_order_cotangent)
                .map(|(tangent, cotangent)| tangent * cotangent)
                .sum::<f64>()
            + ray_end_tangent
                .iter()
                .zip(ray_end_cotangent)
                .map(|(tangent, cotangent)| tangent * cotangent)
                .sum::<f64>();
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
                .sum::<f64>()
            + end_source_gradient
                .iter()
                .zip(end_source_tangent)
                .map(|(gradient, tangent)| gradient * tangent)
                .sum::<f64>();
        assert!((input_dot - output_dot).abs() < 1.0e-13);
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
    fn vector_transport_jvp_and_vjp_are_dual() {
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
            &[0.7],
            &[-0.2],
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
        .unwrap()
        .with_stokes_layout(3, &[1])
        .unwrap();
        let extinction = [2.0, 1.0];
        let albedo = [0.8, 0.6];
        let extinction_tangent = [0.1, -0.2];
        let albedo_tangent = [-0.03, 0.04];
        let mut transport = [0.0; 3];
        let mut transport_tangent = [0.0; 3];
        let mut optical_depth = [0.0];
        let mut attenuation = [0.0];
        let mut prefix = [0.0];
        let mut prefix_tangent = [0.0];
        let mut end = [0.0];
        let mut end_tangent = [0.0];
        geometry
            .assemble_transport_jvp(
                &extinction,
                &albedo,
                &extinction_tangent,
                &albedo_tangent,
                &mut transport,
                &mut transport_tangent,
                &mut optical_depth,
                &mut attenuation,
                &mut prefix,
                &mut prefix_tangent,
                &mut end,
                &mut end_tangent,
            )
            .unwrap();
        assert_eq!(transport[0], transport[1]);
        assert_eq!(transport[1], transport[2]);

        let cotangent = [0.3, -0.2, 0.7];
        let mut extinction_gradient = [0.0; 2];
        let mut albedo_gradient = [0.0; 2];
        geometry
            .assemble_transport_vjp(
                &albedo,
                &cotangent,
                &[],
                &[],
                &[],
                &optical_depth,
                &attenuation,
                &prefix,
                &mut extinction_gradient,
                &mut albedo_gradient,
            )
            .unwrap();
        let lhs: f64 = transport_tangent
            .iter()
            .zip(cotangent)
            .map(|(tangent, cotangent)| tangent * cotangent)
            .sum();
        let rhs: f64 = extinction_tangent
            .iter()
            .zip(extinction_gradient)
            .map(|(tangent, gradient)| tangent * gradient)
            .sum::<f64>()
            + albedo_tangent
                .iter()
                .zip(albedo_gradient)
                .map(|(tangent, gradient)| tangent * gradient)
                .sum::<f64>();
        assert!((lhs - rhs).abs() < 1.0e-14);

        let forcing_gradient = [0.5, -0.2, 1.0];
        let solution = [0.6, 0.7, -0.4];
        let dense_gradient = [
            forcing_gradient[0] * solution[0],
            forcing_gradient[1] * solution[1],
            forcing_gradient[2] * solution[2],
        ];
        let mut dense_extinction_gradient = [0.0; 2];
        let mut dense_albedo_gradient = [0.0; 2];
        geometry
            .assemble_transport_vjp(
                &albedo,
                &dense_gradient,
                &[],
                &[],
                &[],
                &optical_depth,
                &attenuation,
                &prefix,
                &mut dense_extinction_gradient,
                &mut dense_albedo_gradient,
            )
            .unwrap();
        let mut compact_extinction_gradient = [0.0; 2];
        let mut compact_albedo_gradient = [0.0; 2];
        geometry
            .assemble_transport_vjp(
                &albedo,
                &[],
                &[0, 1, 2],
                &solution,
                &forcing_gradient,
                &optical_depth,
                &attenuation,
                &prefix,
                &mut compact_extinction_gradient,
                &mut compact_albedo_gradient,
            )
            .unwrap();
        assert_eq!(compact_extinction_gradient, dense_extinction_gradient);
        assert_eq!(compact_albedo_gradient, dense_albedo_gradient);
    }

    #[test]
    fn adaptive_direct_phase_matches_grouped_phase_and_vjp() {
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
            &[0.7],
            &[-0.2],
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
        .unwrap()
        .with_stokes_layout(3, &[1])
        .unwrap();
        assert!(geometry.uses_grouped_phase_evaluation(2));
        assert!(!geometry.uses_grouped_phase_evaluation(5));

        let maximum_order = [3, 2, 1, 1, 1];
        let mut coefficients = vec![0.0; maximum_order.len() * 3 * 4];
        coefficients[..24].copy_from_slice(&[
            1.0, 0.0, 0.0, 0.2, 0.3, 0.0, 0.0, -0.1, 0.05, 0.0, 0.0, 0.08, 0.9, 0.0, 0.0, -0.15,
            0.2, 0.0, 0.0, 0.07, 0.0, 0.0, 0.0, 0.0,
        ]);
        let (endpoint_indices, endpoint_weights) = geometry.endpoint_stencil(0, 0, false);
        let direct = geometry.endpoint_phase_direct(
            0,
            endpoint_indices,
            endpoint_weights,
            &coefficients,
            &maximum_order,
        );
        let mut phase_values = vec![0.0; maximum_order.len() * 3];
        geometry.fill_phase_values(0, &coefficients, &maximum_order, &mut phase_values);
        let mut grouped = [0.0; 3];
        for (&atmosphere_index, &weight) in endpoint_indices.iter().zip(endpoint_weights) {
            for stokes in 0..3 {
                grouped[stokes] += weight * phase_values[atmosphere_index as usize * 3 + stokes];
            }
        }
        for (direct, grouped) in direct.into_iter().zip(grouped) {
            assert!((direct - grouped).abs() < 1.0e-15);
        }

        let phase_cotangent = [0.4, -0.3, 0.2];
        let mut direct_gradient = vec![0.0; coefficients.len()];
        for (atmosphere_index, &weight) in endpoint_weights.iter().enumerate() {
            geometry.accumulate_endpoint_phase_vjp_direct(
                0,
                atmosphere_index,
                weight,
                &maximum_order,
                phase_cotangent,
                &mut direct_gradient,
            );
        }
        let mut phase_value_gradient = vec![0.0; maximum_order.len() * 3];
        for (atmosphere_index, &weight) in endpoint_weights.iter().enumerate() {
            for stokes in 0..3 {
                phase_value_gradient[atmosphere_index * 3 + stokes] =
                    weight * phase_cotangent[stokes];
            }
        }
        let mut grouped_gradient = vec![0.0; coefficients.len()];
        geometry.accumulate_phase_values_vjp(
            0,
            &maximum_order,
            &phase_value_gradient,
            &mut grouped_gradient,
        );
        for (direct, grouped) in direct_gradient.into_iter().zip(grouped_gradient) {
            assert!((direct - grouped).abs() < 1.0e-15);
        }
    }

    #[test]
    fn fuses_vector_first_order_and_products_are_dual() {
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
            &[0.7],
            &[-0.2],
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
        .unwrap()
        .with_stokes_layout(3, &[1])
        .unwrap();
        let extinction = [2.0, 1.0];
        let albedo = [0.8, 0.6];
        let mut coefficients = [0.0; 24];
        for location in 0..2 {
            let start = location * 12;
            coefficients[start] = 1.0;
            coefficients[start + 3] = 0.15 + 0.05 * location as f64;
            coefficients[start + 4] = 0.2 - 0.03 * location as f64;
            coefficients[start + 7] = -0.1 + 0.02 * location as f64;
            coefficients[start + 8] = 0.04;
            coefficients[start + 11] = 0.08;
        }
        let maximum_order = [3, 3];
        let solar_transmission = [0.8, 0.9];
        let end_of_ray_source = [0.1, -0.2, 0.3];
        let mut transport = [0.0; 3];
        let mut first_order = [0.0; 3];
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
                &end_of_ray_source,
                &mut transport,
                &mut first_order,
                &mut optical_depth,
                &mut attenuation,
                &mut prefix,
                &mut end,
            )
            .unwrap();
        assert!(first_order.iter().all(|value| value.is_finite()));
        assert_ne!(first_order[1], 0.0);
        assert_ne!(first_order[2], 0.0);

        let extinction_tangent = [0.1, -0.2];
        let albedo_tangent = [-0.03, 0.04];
        let coefficient_tangent: [f64; 24] =
            std::array::from_fn(|index| (index as f64 - 9.0) * 0.002);
        let solar_tangent = [-0.04, 0.05];
        let end_of_ray_source_tangent = [0.04, 0.01, -0.03];
        let mut transport_tangent = [0.0; 3];
        let mut first_order_tangent = [0.0; 3];
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
                &end_of_ray_source,
                &end_of_ray_source_tangent,
                &mut transport,
                &mut transport_tangent,
                &mut first_order_tangent,
                &mut end,
                &mut end_tangent,
            )
            .unwrap();

        let transport_cotangent = [0.75, -0.2, 0.4];
        let first_order_cotangent = [-0.5, 0.3, 0.8];
        let end_cotangent = [0.2];
        let mut extinction_gradient = [0.0; 2];
        let mut albedo_gradient = [0.0; 2];
        let mut coefficient_gradient = [0.0; 24];
        let mut solar_gradient = [0.0; 2];
        let mut end_of_ray_source_gradient = [0.0; 3];
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
                &end,
                &end_cotangent,
                &end_of_ray_source,
                &optical_depth,
                &attenuation,
                &prefix,
                &mut extinction_gradient,
                &mut albedo_gradient,
                &mut coefficient_gradient,
                &mut solar_gradient,
                &mut end_of_ray_source_gradient,
            )
            .unwrap();
        let output_dot = transport_tangent
            .iter()
            .zip(transport_cotangent)
            .map(|(tangent, cotangent)| tangent * cotangent)
            .sum::<f64>()
            + first_order_tangent
                .iter()
                .zip(first_order_cotangent)
                .map(|(tangent, cotangent)| tangent * cotangent)
                .sum::<f64>()
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
                .sum::<f64>()
            + end_of_ray_source_gradient
                .iter()
                .zip(end_of_ray_source_tangent)
                .map(|(gradient, tangent)| gradient * tangent)
                .sum::<f64>();
        assert!((input_dot - output_dot).abs() < 1.0e-13);
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
            &[],
            &[],
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
        let end_of_ray_source = [0.12];
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
                &end_of_ray_source,
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
        let expected = integration_factor * (0.4 * start + 0.6 * finish) * 0.5
            + layer_attenuation * end_of_ray_source[0];
        assert!((first_order[0] - expected).abs() < 1.0e-14);

        let extinction_tangent = [0.1, -0.2];
        let albedo_tangent = [-0.03, 0.04];
        let coefficient_tangent = [0.02, -0.01, 0.03, -0.02, 0.01, 0.0];
        let solar_tangent = [-0.04, 0.05];
        let end_of_ray_source_tangent = [-0.07];
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
                &end_of_ray_source,
                &end_of_ray_source_tangent,
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
        let perturbed_end_of_ray_source =
            [end_of_ray_source[0] + epsilon * end_of_ray_source_tangent[0]];
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
                &perturbed_end_of_ray_source,
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
        let mut end_of_ray_source_gradient = [0.0];
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
                &end,
                &end_cotangent,
                &end_of_ray_source,
                &optical_depth,
                &attenuation,
                &prefix,
                &mut extinction_gradient,
                &mut albedo_gradient,
                &mut coefficient_gradient,
                &mut solar_gradient,
                &mut end_of_ray_source_gradient,
            )
            .unwrap();
        let mut compact_extinction_gradient = [0.0; 2];
        let mut compact_albedo_gradient = [0.0; 2];
        let mut compact_coefficient_gradient = [0.0; 6];
        let mut compact_solar_gradient = [0.0; 2];
        let mut compact_end_of_ray_source_gradient = [0.0];
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
                &end,
                &end_cotangent,
                &end_of_ray_source,
                &optical_depth,
                &attenuation,
                &prefix,
                &mut compact_extinction_gradient,
                &mut compact_albedo_gradient,
                &mut compact_coefficient_gradient,
                &mut compact_solar_gradient,
                &mut compact_end_of_ray_source_gradient,
            )
            .unwrap();
        assert_eq!(compact_extinction_gradient, extinction_gradient);
        assert_eq!(compact_albedo_gradient, albedo_gradient);
        assert_eq!(compact_coefficient_gradient, coefficient_gradient);
        assert_eq!(compact_solar_gradient, solar_gradient);
        assert_eq!(
            compact_end_of_ray_source_gradient,
            end_of_ray_source_gradient
        );
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
                .sum::<f64>()
            + end_of_ray_source_gradient[0] * end_of_ray_source_tangent[0];
        assert!((input_dot - output_dot).abs() < 1.0e-13);
    }
}
