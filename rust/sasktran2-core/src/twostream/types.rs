use std::fmt;

use crate::raytracer::{CellId, GeometryKind, TracedRay};

/// Source terms supported by the production C++ two-stream implementation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SourceMode {
    Solar,
    Thermal,
}

/// A plane-parallel upwelling line of sight.
#[derive(Clone, Copy, Debug)]
pub struct View {
    pub cosine: f64,
    pub relative_azimuth: f64,
}

/// Wavelength-independent two-stream geometry.
#[derive(Clone, Debug)]
pub struct Geometry {
    /// Layer thicknesses ordered from top of atmosphere to the surface.
    pub layer_thickness: Vec<f64>,
    /// Row-major matrix mapping layer optical depth to slant optical depth at
    /// each lower layer boundary. Shape `[num_layers, num_layers]`.
    pub chapman_factors: Vec<f64>,
    pub solar_cosine: f64,
    pub quadrature_cosine: f64,
}

/// Compiled wavelength-independent geometry for a batch of spherical lines of
/// sight. Segment arrays are ordered from the far end of each ray toward the
/// observer, matching [`TracedRay::layers`]. Atmospheric grid indices use the
/// two-stream top-down level order.
#[derive(Clone, Debug)]
pub struct SphericalRayGeometry {
    pub(crate) num_levels: usize,
    pub(crate) ray_offsets: Vec<usize>,
    pub(crate) ground_hit: Vec<bool>,
    pub(crate) ground_cos_sza: Vec<f64>,
    pub(crate) segment_layers: Vec<usize>,
    pub(crate) segment_fractions: Vec<f64>,
    pub(crate) segment_cosines: Vec<f64>,
    pub(crate) segment_relative_azimuths: Vec<f64>,
    pub(crate) segment_cos_sza: Vec<f64>,
    pub(crate) od_offsets: Vec<usize>,
    pub(crate) od_indices: Vec<usize>,
    pub(crate) od_weights: Vec<f64>,
}

impl SphericalRayGeometry {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        num_levels: usize,
        ray_offsets: Vec<usize>,
        ground_hit: Vec<bool>,
        ground_cos_sza: Vec<f64>,
        segment_layers: Vec<usize>,
        segment_fractions: Vec<f64>,
        segment_cosines: Vec<f64>,
        segment_relative_azimuths: Vec<f64>,
        segment_cos_sza: Vec<f64>,
        od_offsets: Vec<usize>,
        od_indices: Vec<usize>,
        od_weights: Vec<f64>,
    ) -> Result<Self, TwoStreamError> {
        let result = Self {
            num_levels,
            ray_offsets,
            ground_hit,
            ground_cos_sza,
            segment_layers,
            segment_fractions,
            segment_cosines,
            segment_relative_azimuths,
            segment_cos_sza,
            od_offsets,
            od_indices,
            od_weights,
        };
        result.validate()?;
        Ok(result)
    }

    /// Compile the existing pure-Rust one-dimensional traced-ray primitives.
    /// `altitudes_bottom_up` must be the altitude grid used by the ray tracer.
    pub fn from_traced_rays(
        altitudes_bottom_up: &[f64],
        rays: &[TracedRay],
    ) -> Result<Self, TwoStreamError> {
        use std::collections::BTreeMap;

        if altitudes_bottom_up.len() < 2
            || altitudes_bottom_up
                .windows(2)
                .any(|pair| !pair[0].is_finite() || pair[1] <= pair[0])
        {
            return Err(TwoStreamError::invalid(
                "spherical ray compilation requires increasing altitudes",
            ));
        }
        let num_levels = altitudes_bottom_up.len();
        let num_layers = num_levels - 1;
        let mut ray_offsets = Vec::with_capacity(rays.len() + 1);
        let mut ground_hit = Vec::with_capacity(rays.len());
        let mut ground_cos_sza = Vec::with_capacity(rays.len());
        let mut segment_layers = Vec::new();
        let mut segment_fractions = Vec::new();
        let mut segment_cosines = Vec::new();
        let mut segment_relative_azimuths = Vec::new();
        let mut segment_cos_sza = Vec::new();
        let mut od_offsets = Vec::new();
        let mut od_indices = Vec::new();
        let mut od_weights = Vec::new();
        ray_offsets.push(0);
        od_offsets.push(0);

        for ray in rays {
            if ray.geometry != GeometryKind::Spherical {
                return Err(TwoStreamError::invalid(
                    "two-stream spherical paths require spherical traced rays",
                ));
            }
            ground_hit.push(ray.ground_is_hit);
            ground_cos_sza.push(ray.layers.first().map_or(0.0, |layer| layer.cos_sza_exit));
            for layer in &ray.layers {
                let (start_weight, end_weight) = normalized_source_weights(
                    layer.od_quad_start_fraction,
                    layer.od_quad_end_fraction,
                );
                let altitude =
                    start_weight * layer.entrance.altitude + end_weight * layer.exit.altitude;
                let bottom_layer = match layer.cell {
                    Some(CellId::AltitudeLayer(index)) => index.min(num_layers - 1),
                    _ => altitudes_bottom_up
                        .partition_point(|candidate| *candidate <= altitude)
                        .saturating_sub(1)
                        .min(num_layers - 1),
                };
                let top_down_layer = num_layers - 1 - bottom_layer;
                let bottom = altitudes_bottom_up[bottom_layer];
                let top = altitudes_bottom_up[bottom_layer + 1];
                let fraction_from_top = ((top - altitude) / (top - bottom)).clamp(0.0, 1.0);
                let local_up = (layer.entrance.position * start_weight
                    + layer.exit.position * end_weight)
                    .normalized();
                let cosine = -layer.average_look_away.dot(local_up);
                let sin_azimuth =
                    start_weight * layer.saz_entrance.sin() + end_weight * layer.saz_exit.sin();
                let cos_azimuth =
                    start_weight * layer.saz_entrance.cos() + end_weight * layer.saz_exit.cos();

                segment_layers.push(top_down_layer);
                segment_fractions.push(fraction_from_top);
                segment_cosines.push(cosine.clamp(-1.0, 1.0));
                segment_relative_azimuths.push(sin_azimuth.atan2(cos_azimuth));
                segment_cos_sza.push(
                    (start_weight * layer.cos_sza_entrance + end_weight * layer.cos_sza_exit)
                        .clamp(-1.0, 1.0),
                );

                let mut weights = BTreeMap::<usize, f64>::new();
                for interpolation in layer.entrance.interpolation.iter() {
                    *weights.entry(interpolation.index).or_default() +=
                        layer.od_quad_start * interpolation.weight;
                }
                for interpolation in layer.exit.interpolation.iter() {
                    *weights.entry(interpolation.index).or_default() +=
                        layer.od_quad_end * interpolation.weight;
                }
                for (bottom_up_index, weight) in weights {
                    if weight != 0.0 {
                        od_indices.push(num_levels - 1 - bottom_up_index);
                        od_weights.push(weight);
                    }
                }
                od_offsets.push(od_indices.len());
            }
            ray_offsets.push(segment_layers.len());
        }

        Self::new(
            num_levels,
            ray_offsets,
            ground_hit,
            ground_cos_sza,
            segment_layers,
            segment_fractions,
            segment_cosines,
            segment_relative_azimuths,
            segment_cos_sza,
            od_offsets,
            od_indices,
            od_weights,
        )
    }

    pub fn num_views(&self) -> usize {
        self.ground_hit.len()
    }

    pub fn num_segments(&self) -> usize {
        self.segment_layers.len()
    }

    fn validate(&self) -> Result<(), TwoStreamError> {
        let num_views = self.ground_hit.len();
        let num_segments = self.segment_layers.len();
        if self.num_levels < 2
            || self.ray_offsets.len() != num_views + 1
            || self.ground_cos_sza.len() != num_views
            || self.ray_offsets.first() != Some(&0)
            || self.ray_offsets.last() != Some(&num_segments)
            || self.ray_offsets.windows(2).any(|pair| pair[0] > pair[1])
        {
            return Err(TwoStreamError::invalid("invalid spherical ray offsets"));
        }
        if [
            self.segment_fractions.len(),
            self.segment_cosines.len(),
            self.segment_relative_azimuths.len(),
            self.segment_cos_sza.len(),
        ]
        .into_iter()
        .any(|len| len != num_segments)
            || self.od_offsets.len() != num_segments + 1
            || self.od_offsets.first() != Some(&0)
            || self.od_offsets.last() != Some(&self.od_indices.len())
            || self.od_indices.len() != self.od_weights.len()
            || self.od_offsets.windows(2).any(|pair| pair[0] > pair[1])
        {
            return Err(TwoStreamError::invalid(
                "spherical segment arrays have inconsistent shapes",
            ));
        }
        if self
            .segment_layers
            .iter()
            .any(|&layer| layer + 1 >= self.num_levels)
            || self
                .od_indices
                .iter()
                .any(|&index| index >= self.num_levels)
            || self
                .segment_fractions
                .iter()
                .any(|&fraction| !fraction.is_finite() || !(0.0..=1.0).contains(&fraction))
            || self
                .segment_cosines
                .iter()
                .any(|&cosine| !cosine.is_finite() || !(-1.0..=1.0).contains(&cosine))
            || self
                .segment_cos_sza
                .iter()
                .any(|&cosine| !cosine.is_finite() || !(-1.0..=1.0).contains(&cosine))
            || self
                .segment_relative_azimuths
                .iter()
                .chain(&self.ground_cos_sza)
                .chain(&self.od_weights)
                .any(|value| !value.is_finite())
        {
            return Err(TwoStreamError::invalid(
                "spherical segment geometry contains invalid values",
            ));
        }
        Ok(())
    }
}

fn normalized_source_weights(start: f64, end: f64) -> (f64, f64) {
    let sum = start + end;
    if start.is_finite() && end.is_finite() && sum > 0.0 {
        (start / sum, end / sum)
    } else {
        (0.5, 0.5)
    }
}

/// Two-stream column geometries and the spherical lines of sight sampling
/// their local source fields.
#[derive(Clone, Debug)]
pub struct SphericalGeometry {
    pub(crate) sza_grid: Vec<f64>,
    pub(crate) columns: Vec<Geometry>,
    pub(crate) rays: SphericalRayGeometry,
}

impl SphericalGeometry {
    pub fn new(
        sza_grid: Vec<f64>,
        columns: Vec<Geometry>,
        rays: SphericalRayGeometry,
    ) -> Result<Self, TwoStreamError> {
        if columns.is_empty()
            || columns.len() != sza_grid.len()
            || sza_grid
                .iter()
                .any(|value| !value.is_finite() || !(-1.0..=1.0).contains(value))
            || sza_grid.windows(2).any(|pair| pair[1] <= pair[0])
        {
            return Err(TwoStreamError::invalid("invalid spherical SZA grid"));
        }
        for column in &columns {
            if column.num_layers() + 1 != rays.num_levels {
                return Err(TwoStreamError::invalid(
                    "spherical column and ray level counts differ",
                ));
            }
            if column.layer_thickness != columns[0].layer_thickness {
                return Err(TwoStreamError::invalid(
                    "spherical columns must use the same vertical grid",
                ));
            }
        }
        Ok(Self {
            sza_grid,
            columns,
            rays,
        })
    }

    pub fn rays(&self) -> &SphericalRayGeometry {
        &self.rays
    }
}

impl Geometry {
    pub fn new(layer_thickness: Vec<f64>, chapman_factors: Vec<f64>, solar_cosine: f64) -> Self {
        Self {
            layer_thickness,
            chapman_factors,
            solar_cosine,
            quadrature_cosine: 0.5,
        }
    }

    pub fn num_layers(&self) -> usize {
        self.layer_thickness.len()
    }
}

/// Atmospheric quantities on layer boundaries, ordered TOA to surface.
///
/// Level arrays have shape `[num_layers + 1, num_wavelengths]`. Surface and
/// irradiance arrays have length `num_wavelengths`.
#[derive(Clone, Debug)]
pub struct AtmosphereBatch {
    pub num_wavelengths: usize,
    pub extinction: Vec<f64>,
    pub single_scatter_albedo: Vec<f64>,
    pub first_legendre: Vec<f64>,
    pub emission: Option<Vec<f64>>,
    /// Lambertian albedo (not BRDF divided by pi).
    pub surface_albedo: Vec<f64>,
    pub surface_emission: Option<Vec<f64>>,
    pub solar_irradiance: Option<Vec<f64>>,
}

/// Prepared layer quantities. All layer arrays use `[layer, wavelength]`.
#[derive(Clone, Debug)]
pub struct LayerInputs {
    pub num_layers: usize,
    pub num_wavelengths: usize,
    pub optical_depth: Vec<f64>,
    pub single_scatter_albedo: Vec<f64>,
    pub first_legendre: Vec<f64>,
    pub transmission: Option<Vec<f64>>,
    pub average_secant: Option<Vec<f64>>,
    pub thermal_b0: Option<Vec<f64>>,
    pub thermal_b1: Option<Vec<f64>>,
    pub surface_albedo: Vec<f64>,
    pub surface_emission: Option<Vec<f64>>,
}

#[derive(Clone, Debug)]
pub struct RadianceBatch {
    pub num_views: usize,
    pub num_wavelengths: usize,
    /// Shape `[view, wavelength]`.
    pub values: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct LayerAdjoints {
    pub optical_depth: Vec<f64>,
    pub single_scatter_albedo: Vec<f64>,
    pub first_legendre: Vec<f64>,
    /// Shape `[layer boundary, wavelength]` for solar mode.
    pub transmission: Option<Vec<f64>>,
    pub average_secant: Option<Vec<f64>>,
    pub thermal_b0: Option<Vec<f64>>,
    pub thermal_b1: Option<Vec<f64>>,
    pub surface_albedo: Vec<f64>,
    pub surface_emission: Option<Vec<f64>>,
}

#[derive(Clone, Debug)]
pub struct AtmosphereAdjoints {
    pub extinction: Vec<f64>,
    pub single_scatter_albedo: Vec<f64>,
    pub first_legendre: Vec<f64>,
    pub emission: Option<Vec<f64>>,
    pub surface_albedo: Vec<f64>,
    pub surface_emission: Option<Vec<f64>>,
}

/// Per-view atmospheric Jacobians for a radiance batch.
///
/// Level quantities have shape `[view, level, wavelength]`; surface
/// quantities have shape `[view, wavelength]`.  The wavelength dimension is
/// contiguous so the engine can copy one wavelength/view result without
/// transposing the SIMD-friendly storage.
#[derive(Clone, Debug)]
pub struct AtmosphereJacobians {
    pub num_views: usize,
    pub num_levels: usize,
    pub num_wavelengths: usize,
    pub extinction: Vec<f64>,
    pub single_scatter_albedo: Vec<f64>,
    pub first_legendre: Vec<f64>,
    pub emission: Option<Vec<f64>>,
    pub surface_albedo: Vec<f64>,
    pub surface_emission: Option<Vec<f64>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TwoStreamError(pub(crate) String);

impl fmt::Display for TwoStreamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TwoStreamError {}

impl TwoStreamError {
    pub(crate) fn invalid(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}
