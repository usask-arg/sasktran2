//! Geometry-neutral operators and iterative solvers for successive orders of
//! scattering.
//!
//! Geometry builders are deliberately kept outside the numerical core. They
//! reduce one- or multi-dimensional traced geometry to source nodes and
//! compressed interpolation stencils. The transport, scattering, and solver
//! implementations therefore do not depend on altitude, solar-zenith angle,
//! or structured-cell conventions.

mod coefficient;
#[cfg(not(test))]
mod cxx;
mod geometry;
mod operator;
mod ray_transport;
mod solver;

pub use coefficient::{
    ScalarCoefficientBasis, ScalarCoefficientScattering, VectorCoefficientBasis,
    VectorCoefficientScattering,
};
pub use geometry::{CompressedStencils, SourceNode};
pub use operator::{
    BlockDiagonalMatrix, CsrMatrix, FixedPointProblem, OperatorError, ScatteringOperator,
};
pub use ray_transport::ScalarRayTransport;
pub use solver::{
    ConvergenceReason, SolverConfig, SolverDiagnostics, SolverError, SuccessiveOrdersSolver,
};

#[cfg(test)]
mod tests;
