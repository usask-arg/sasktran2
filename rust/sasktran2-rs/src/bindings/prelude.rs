#[derive(Clone, Copy, Debug)]
pub enum Stokes {
    Stokes1,
    Stokes3,
}

impl Stokes {
    pub fn num_legendre(&self) -> usize {
        match self {
            Self::Stokes1 => 1,
            Self::Stokes3 => 4,
        }
    }

    pub fn num_stokes(&self) -> usize {
        match self {
            Self::Stokes1 => 1,
            Self::Stokes3 => 3,
        }
    }
}

pub use super::atmosphere::Atmosphere;
pub use super::config::{Config, MultipleScatterSource, SingleScatterSource};
pub use super::engine::Engine;
pub use super::geometry::Geometry1D;
pub use super::geometry::{GeometryType, InterpolationMethod};
pub use super::output::Output;
pub use super::viewing_geometry::ViewingGeometry;
pub use anyhow::{Result, anyhow};
