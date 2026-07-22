//! Experimental wavelength-batched two-stream radiative-transfer solver.
//!
//! The implementation intentionally has no dependency on the C++ engine.  Its
//! arrays are layer-major and wavelength-contiguous: element `(layer, wavel)`
//! is stored at `layer * num_wavelengths + wavel`.  That layout is part of the
//! API because the hot kernels vectorize over wavelength.

const THERMAL_MIN_OPTICAL_DEPTH: f64 = 1.0e-10;
const THERMAL_MIN_EMISSION: f64 = 1.0e-30;

#[inline]
fn thermal_profile_is_active(top: f64, bottom: f64, optical_depth: f64) -> bool {
    optical_depth > THERMAL_MIN_OPTICAL_DEPTH
        && top > THERMAL_MIN_EMISSION
        && bottom > THERMAL_MIN_EMISSION
}

/// Exponential thermal-source slope in optical-depth coordinates.
///
/// Transparent and vanishing-emission layers use a constant profile. Positive
/// profiles use `ln1p` so the isothermal limit remains continuous and retains
/// its well-defined derivatives with respect to the endpoint emissions.
#[inline]
fn thermal_profile_slope(top: f64, bottom: f64, optical_depth: f64) -> f64 {
    if thermal_profile_is_active(top, bottom, optical_depth) {
        ((top - bottom) / bottom).ln_1p() / optical_depth
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

#[cfg(test)]
mod thermal_profile_tests {
    use super::{thermal_profile_slope, thermal_profile_slope_adjoint};

    #[test]
    fn isothermal_slope_has_continuous_endpoint_derivatives() {
        let emission = 2.0;
        let optical_depth = 0.5;
        let slope = thermal_profile_slope(emission, emission, optical_depth);
        assert_eq!(slope, 0.0);

        let (d_top, d_bottom, d_od) =
            thermal_profile_slope_adjoint(emission, emission, optical_depth, slope, 1.0);
        assert_eq!(d_top, 1.0);
        assert_eq!(d_bottom, -1.0);
        assert_eq!(d_od, 0.0);

        let epsilon = 1.0e-6;
        let numeric_top = (thermal_profile_slope(
            emission + epsilon,
            emission,
            optical_depth,
        ) - thermal_profile_slope(emission - epsilon, emission, optical_depth))
            / (2.0 * epsilon);
        let numeric_bottom = (thermal_profile_slope(
            emission,
            emission + epsilon,
            optical_depth,
        ) - thermal_profile_slope(emission, emission - epsilon, optical_depth))
            / (2.0 * epsilon);

        assert!((d_top - numeric_top).abs() < 1.0e-10);
        assert!((d_bottom - numeric_bottom).abs() < 1.0e-10);
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
