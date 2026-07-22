//! Experimental wavelength-batched two-stream radiative-transfer solver.
//!
//! The implementation intentionally has no dependency on the C++ engine.  Its
//! arrays are layer-major and wavelength-contiguous: element `(layer, wavel)`
//! is stored at `layer * num_wavelengths + wavel`.  That layout is part of the
//! API because the hot kernels vectorize over wavelength.

const THERMAL_MIN_OPTICAL_DEPTH: f64 = 1.0e-10;
const THERMAL_MIN_EMISSION: f64 = 1.0e-30;
const THERMAL_RELATIVE_DIFFERENCE: f64 = 1.0e-15;

#[inline]
fn thermal_profile_is_active(top: f64, bottom: f64, optical_depth: f64) -> bool {
    optical_depth > THERMAL_MIN_OPTICAL_DEPTH
        && top > THERMAL_MIN_EMISSION
        && bottom > THERMAL_MIN_EMISSION
        && (top - bottom).abs() > THERMAL_RELATIVE_DIFFERENCE * top.max(bottom)
}

/// Exponential thermal-source slope in optical-depth coordinates.
///
/// The thresholds match the standard discrete-ordinate layer preparation and
/// make transparent, vanishing-emission, and isothermal limits constant.
#[inline]
fn thermal_profile_slope(top: f64, bottom: f64, optical_depth: f64) -> f64 {
    if thermal_profile_is_active(top, bottom, optical_depth) {
        (top / bottom).ln() / optical_depth
    } else {
        0.0
    }
}

/// Reverse the active branch of [`thermal_profile_slope`].
#[inline]
fn thermal_profile_slope_adjoint(
    top: f64,
    bottom: f64,
    optical_depth: f64,
    slope: f64,
    adjoint: f64,
) -> (f64, f64, f64) {
    if thermal_profile_is_active(top, bottom, optical_depth) {
        (
            adjoint / (top * optical_depth),
            -adjoint / (bottom * optical_depth),
            -adjoint * slope / optical_depth,
        )
    } else {
        (0.0, 0.0, 0.0)
    }
}

#[cfg(not(test))]
mod cxx;
mod explicit;
#[cfg(test)]
mod reverse;
mod solver;
mod types;

#[cfg(feature = "simd")]
mod simd;

pub use solver::{ExecutionPolicy, TwoStreamSolver, Workspace};
pub use types::{
    AtmosphereAdjoints, AtmosphereBatch, AtmosphereJacobians, Geometry, LayerAdjoints, LayerInputs,
    RadianceBatch, SourceMode, SphericalGeometry, SphericalRayGeometry, TwoStreamError, View,
};

#[cfg(test)]
mod tests;
