#[derive(Clone, Copy, Debug)]
pub enum Stokes {
    Stokes1,
    Stokes3
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


pub use anyhow::{Result, anyhow};
pub use crate::atmosphere::Atmosphere;
pub use crate::config::{Config, MultipleScatterSource, SingleScatterSource};
pub use crate::engine::Engine;
pub use crate::geometry::Geometry1D;
pub use crate::output::Output;
pub use crate::geometry::{InterpolationMethod, GeometryType};
pub use crate::viewing_geometry::ViewingGeometry;
