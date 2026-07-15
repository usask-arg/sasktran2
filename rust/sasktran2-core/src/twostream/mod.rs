//! Experimental wavelength-batched two-stream radiative-transfer solver.
//!
//! The implementation intentionally has no dependency on the C++ engine.  Its
//! arrays are layer-major and wavelength-contiguous: element `(layer, wavel)`
//! is stored at `layer * num_wavelengths + wavel`.  That layout is part of the
//! API because the hot kernels vectorize over wavelength.

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
