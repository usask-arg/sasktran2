use std::error::Error;
use std::fmt;

use super::primitive::{Intersection, Primitive, PrimitiveId};
use super::vec3::Vec3;

pub const MAX_BOUNDARY_TAGS: usize = 4;
pub const MAX_STENCIL_WEIGHTS: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GeometryKind {
    PlaneParallel,
    Spherical,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InterpolationMethod {
    Shell,
    Linear,
    Lower,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BoundaryTag {
    Surface { altitude_index: usize },
    TopOfAtmosphere { altitude_index: usize },
    Altitude { index: usize },
    Cross { index: usize },
    Custom { id: u32 },
}

impl BoundaryTag {
    #[inline(always)]
    pub fn altitude_index(self) -> Option<usize> {
        match self {
            Self::Surface { altitude_index } | Self::TopOfAtmosphere { altitude_index } => {
                Some(altitude_index)
            }
            Self::Altitude { index } => Some(index),
            Self::Cross { .. } | Self::Custom { .. } => None,
        }
    }

    #[inline(always)]
    pub fn is_surface(self) -> bool {
        matches!(self, Self::Surface { .. })
    }

    #[inline(always)]
    pub fn is_top_of_atmosphere(self) -> bool {
        matches!(self, Self::TopOfAtmosphere { .. })
    }

    #[inline(always)]
    pub fn is_vertical(self) -> bool {
        self.altitude_index().is_some()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundarySet {
    tags: [Option<BoundaryTag>; MAX_BOUNDARY_TAGS],
    len: u8,
}

impl Default for BoundarySet {
    fn default() -> Self {
        Self {
            tags: [None; MAX_BOUNDARY_TAGS],
            len: 0,
        }
    }
}

impl BoundarySet {
    #[inline(always)]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len as usize
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    pub fn push(&mut self, tag: BoundaryTag) {
        if self.contains(tag) {
            return;
        }
        let len = self.len();
        if len < MAX_BOUNDARY_TAGS {
            self.tags[len] = Some(tag);
            self.len += 1;
        }
    }

    #[inline]
    pub fn contains(&self, tag: BoundaryTag) -> bool {
        self.iter().any(|existing| existing == tag)
    }

    #[inline]
    pub fn contains_surface(&self) -> bool {
        self.iter().any(BoundaryTag::is_surface)
    }

    #[inline]
    pub fn contains_top_of_atmosphere(&self) -> bool {
        self.iter().any(BoundaryTag::is_top_of_atmosphere)
    }

    #[inline]
    pub fn contains_vertical(&self) -> bool {
        self.iter().any(BoundaryTag::is_vertical)
    }

    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = BoundaryTag> + '_ {
        self.tags[..self.len()].iter().filter_map(|tag| *tag)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CellId {
    AltitudeLayer(usize),
    Structured2D {
        altitude_index: usize,
        cross_index: usize,
    },
    Unstructured(u64),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterpolationWeight {
    pub index: usize,
    pub weight: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterpolationStencil {
    weights: [InterpolationWeight; MAX_STENCIL_WEIGHTS],
    len: u8,
}

impl Default for InterpolationStencil {
    fn default() -> Self {
        Self {
            weights: [InterpolationWeight {
                index: 0,
                weight: 0.0,
            }; MAX_STENCIL_WEIGHTS],
            len: 0,
        }
    }
}

impl InterpolationStencil {
    #[inline(always)]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len as usize
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    pub fn push(&mut self, index: usize, weight: f64) {
        let len = self.len();
        if len < MAX_STENCIL_WEIGHTS {
            self.weights[len] = InterpolationWeight { index, weight };
            self.len += 1;
        }
    }

    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = InterpolationWeight> + '_ {
        self.weights[..self.len()].iter().copied()
    }
}

pub trait MediumGrid {
    fn locate_cell(&self, point: Vec3) -> Option<CellId>;
    fn interpolation_weights(&self, point: Vec3) -> InterpolationStencil;
    fn classify_intersection(&self, intersection: &Intersection) -> BoundaryTag;
}

#[derive(Debug, Clone, PartialEq)]
pub enum GridError {
    TooFewAltitudes,
    NonIncreasingAltitude { index: usize },
}

impl fmt::Display for GridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooFewAltitudes => write!(f, "a vertical grid needs at least two altitudes"),
            Self::NonIncreasingAltitude { index } => {
                write!(
                    f,
                    "altitudes must be strictly increasing near index {index}"
                )
            }
        }
    }
}

impl Error for GridError {}

#[derive(Debug, Clone, PartialEq)]
pub struct VerticalGrid1D {
    earth_radius: f64,
    altitudes: Vec<f64>,
    interpolation_method: InterpolationMethod,
    geometry: GeometryKind,
}

impl VerticalGrid1D {
    pub fn new(
        earth_radius: f64,
        altitudes: Vec<f64>,
        interpolation_method: InterpolationMethod,
        geometry: GeometryKind,
    ) -> Result<Self, GridError> {
        if altitudes.len() < 2 {
            return Err(GridError::TooFewAltitudes);
        }
        for i in 1..altitudes.len() {
            if altitudes[i] <= altitudes[i - 1] {
                return Err(GridError::NonIncreasingAltitude { index: i });
            }
        }
        Ok(Self {
            earth_radius,
            altitudes,
            interpolation_method,
            geometry,
        })
    }

    #[inline(always)]
    pub fn earth_radius(&self) -> f64 {
        self.earth_radius
    }

    #[inline(always)]
    pub fn altitudes(&self) -> &[f64] {
        &self.altitudes
    }

    #[inline(always)]
    pub fn interpolation_method(&self) -> InterpolationMethod {
        self.interpolation_method
    }

    #[inline(always)]
    pub fn geometry(&self) -> GeometryKind {
        self.geometry
    }

    #[inline(always)]
    pub fn ground_altitude(&self) -> f64 {
        self.altitudes[0]
    }

    #[inline(always)]
    pub fn top_altitude(&self) -> f64 {
        self.altitudes[self.altitudes.len() - 1]
    }

    #[inline(always)]
    pub fn boundary_tag(&self, altitude_index: usize) -> BoundaryTag {
        if altitude_index == 0 {
            BoundaryTag::Surface { altitude_index }
        } else if altitude_index == self.altitudes.len() - 1 {
            BoundaryTag::TopOfAtmosphere { altitude_index }
        } else {
            BoundaryTag::Altitude {
                index: altitude_index,
            }
        }
    }

    pub fn vertical_primitives(&self) -> Vec<Primitive> {
        self.altitudes
            .iter()
            .enumerate()
            .map(|(index, altitude)| {
                let tag = self.boundary_tag(index);
                let id = PrimitiveId(index);
                match self.geometry {
                    GeometryKind::Spherical => {
                        Primitive::sphere(id, Vec3::ZERO, self.earth_radius + altitude, tag)
                    }
                    GeometryKind::PlaneParallel => Primitive::plane(
                        id,
                        Vec3::new(0.0, 0.0, self.earth_radius + altitude),
                        Vec3::UNIT_Z,
                        tag,
                    ),
                }
            })
            .collect()
    }

    #[inline]
    pub fn altitude_at(&self, point: Vec3) -> f64 {
        match self.geometry {
            GeometryKind::Spherical => point.norm() - self.earth_radius,
            GeometryKind::PlaneParallel => point.z - self.earth_radius,
        }
    }

    #[inline]
    pub fn radius_for_altitude(&self, altitude: f64) -> f64 {
        self.earth_radius + altitude
    }

    #[inline]
    pub fn is_inside_altitude(&self, altitude: f64, epsilon: f64) -> bool {
        altitude >= self.ground_altitude() - epsilon && altitude <= self.top_altitude() + epsilon
    }

    #[inline]
    pub fn on_exact_vertical_boundary(&self, altitude: f64, epsilon: f64) -> Option<usize> {
        self.altitudes
            .iter()
            .position(|boundary| (altitude - boundary).abs() <= epsilon)
    }

    fn upper_bound(&self, altitude: f64) -> usize {
        self.altitudes
            .partition_point(|boundary| *boundary <= altitude)
    }

    pub fn locate_altitude_layer(&self, altitude: f64, epsilon: f64) -> Option<usize> {
        if !self.is_inside_altitude(altitude, epsilon) {
            return None;
        }
        if let Some(index) = self.on_exact_vertical_boundary(altitude, epsilon) {
            return Some(index.min(self.altitudes.len() - 2));
        }
        let upper = self.upper_bound(altitude);
        if upper == 0 || upper >= self.altitudes.len() {
            None
        } else {
            Some(upper - 1)
        }
    }

    pub fn interpolation_weights_at_altitude(&self, altitude: f64) -> InterpolationStencil {
        let mut stencil = InterpolationStencil::new();

        if let Some(index) = self.on_exact_vertical_boundary(altitude, 1e-9) {
            stencil.push(index, 1.0);
            return stencil;
        }

        if altitude <= self.altitudes[0] {
            stencil.push(0, 1.0);
            return stencil;
        }

        let last = self.altitudes.len() - 1;
        if altitude >= self.altitudes[last] {
            stencil.push(last, 1.0);
            return stencil;
        }

        let upper = self.upper_bound(altitude);
        let lower = upper - 1;

        match self.interpolation_method {
            InterpolationMethod::Linear => {
                let width = self.altitudes[upper] - self.altitudes[lower];
                let upper_weight = (altitude - self.altitudes[lower]) / width;
                stencil.push(lower, 1.0 - upper_weight);
                stencil.push(upper, upper_weight);
            }
            InterpolationMethod::Shell | InterpolationMethod::Lower => {
                stencil.push(lower, 1.0);
            }
        }

        stencil
    }
}

impl MediumGrid for VerticalGrid1D {
    fn locate_cell(&self, point: Vec3) -> Option<CellId> {
        self.locate_altitude_layer(self.altitude_at(point), 1e-8)
            .map(CellId::AltitudeLayer)
    }

    fn interpolation_weights(&self, point: Vec3) -> InterpolationStencil {
        self.interpolation_weights_at_altitude(self.altitude_at(point))
    }

    fn classify_intersection(&self, intersection: &Intersection) -> BoundaryTag {
        intersection.boundary
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid() -> VerticalGrid1D {
        VerticalGrid1D::new(
            10.0,
            vec![0.0, 10.0, 20.0, 30.0],
            InterpolationMethod::Linear,
            GeometryKind::Spherical,
        )
        .unwrap()
    }

    #[test]
    fn vertical_grid_builds_sphere_primitives() {
        let grid = grid();
        let primitives = grid.vertical_primitives();
        assert_eq!(primitives.len(), 4);
        assert!(primitives[0].boundary().is_surface());
        assert!(primitives[3].boundary().is_top_of_atmosphere());
    }

    #[test]
    fn interpolation_weights_are_linear_in_altitude() {
        let stencil = grid().interpolation_weights_at_altitude(15.0);
        let weights: Vec<_> = stencil.iter().collect();
        assert_eq!(weights.len(), 2);
        assert_eq!(weights[0].index, 1);
        assert_eq!(weights[1].index, 2);
        assert!((weights[0].weight - 0.5).abs() < 1e-12);
        assert!((weights[1].weight - 0.5).abs() < 1e-12);
    }

    #[test]
    fn future_cell_ids_are_not_limited_to_altitude_layers() {
        let cell = CellId::Structured2D {
            altitude_index: 2,
            cross_index: 4,
        };
        assert_eq!(
            cell,
            CellId::Structured2D {
                altitude_index: 2,
                cross_index: 4
            }
        );
    }
}
