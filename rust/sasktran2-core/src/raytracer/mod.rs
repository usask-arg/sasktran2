//! Pure Rust ray tracing primitives for 1D vertical and structured 2D grids.
//!
//! This module intentionally has no C++ or Python interface dependencies.  The
//! data model separates geometric intersections from medium topology so future
//! 2D/3D grids can map the same primitive intersections onto richer cells.

pub mod grid;
pub mod grid2d;
pub mod layer;
pub mod od_quadrature;
pub mod primitive;
pub mod ray;
pub mod refraction;
pub mod solar;
pub mod trace;
pub mod trace2d;
pub mod vec3;

pub mod cxx;

#[cfg(test)]
mod parity_tests;

pub use grid::{
    BoundarySet, BoundaryTag, CellId, GeometryKind, GridError, InterpolationMethod,
    InterpolationStencil, InterpolationWeight, MediumGrid, VerticalGrid1D,
};
pub use grid2d::{AngularBasis, StructuredGrid2D};
pub use layer::{Layer, LayerType, TraceEvent, TraceEventKind, TracePoint, TracedRay};
pub use primitive::{Crossing, Intersection, Primitive, PrimitiveId, PrimitiveKind};
pub use ray::{Ray, StraightPath};
pub use refraction::{PathIntegral, RefractiveProfile};
pub use solar::SolarContext;
pub use trace::{TraceOptions, TraceScratch, VerticalRayTracer};
pub use trace2d::{StructuredRayTracer2D, TraceOptions2D, TraceScratch2D};
pub use vec3::Vec3;
