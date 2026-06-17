//! Pure Rust ray tracing primitives and 1D vertical-grid tracing.
//!
//! This module intentionally has no C++ or Python interface dependencies.  The
//! data model separates geometric intersections from medium topology so future
//! 2D/3D grids can map the same primitive intersections onto richer cells.

pub mod grid;
pub mod layer;
pub mod od_quadrature;
pub mod primitive;
pub mod ray;
pub mod refraction;
pub mod solar;
pub mod trace;
pub mod vec3;

#[cfg(test)]
mod parity_tests;

pub use grid::{
    BoundarySet, BoundaryTag, CellId, GeometryKind, GridError, InterpolationMethod,
    InterpolationStencil, InterpolationWeight, MediumGrid, VerticalGrid1D,
};
pub use layer::{Layer, LayerType, TraceEvent, TraceEventKind, TracePoint, TracedRay};
pub use primitive::{Crossing, Intersection, Primitive, PrimitiveId, PrimitiveKind};
pub use ray::{Ray, StraightPath};
pub use refraction::{PathIntegral, RefractiveProfile};
pub use solar::SolarContext;
pub use trace::{TraceOptions, TraceScratch, VerticalRayTracer};
pub use vec3::Vec3;
