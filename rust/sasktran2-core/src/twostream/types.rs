use std::fmt;

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
