use super::grid::{BoundarySet, CellId, GeometryKind, InterpolationMethod, InterpolationStencil};
use super::vec3::Vec3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TraceEventKind {
    Observer,
    Boundary,
    Tangent,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TraceEvent {
    pub kind: TraceEventKind,
    pub boundaries: BoundarySet,
}

impl TraceEvent {
    #[inline(always)]
    pub fn observer() -> Self {
        Self {
            kind: TraceEventKind::Observer,
            boundaries: BoundarySet::new(),
        }
    }

    #[inline(always)]
    pub fn boundary(boundaries: BoundarySet) -> Self {
        Self {
            kind: TraceEventKind::Boundary,
            boundaries,
        }
    }

    #[inline(always)]
    pub fn tangent(boundaries: BoundarySet) -> Self {
        Self {
            kind: TraceEventKind::Tangent,
            boundaries,
        }
    }

    #[inline(always)]
    pub fn is_tangent(self) -> bool {
        self.kind == TraceEventKind::Tangent
    }

    #[inline(always)]
    pub fn is_observer(self) -> bool {
        self.kind == TraceEventKind::Observer
    }

    #[inline(always)]
    pub fn is_boundary(self) -> bool {
        self.kind == TraceEventKind::Boundary
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TracePoint {
    pub distance: f64,
    pub position: Vec3,
    pub altitude: f64,
    pub event: TraceEvent,
    pub interpolation: InterpolationStencil,
    pub cell: Option<CellId>,
    pub on_exact_vertical_boundary: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LayerType {
    Complete,
    Partial,
    Tangent,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Layer {
    pub entrance: TracePoint,
    pub exit: TracePoint,
    pub average_look_away: Vec3,
    pub geometric_distance: f64,
    pub curvature_factor: f64,
    pub od_quad_start: f64,
    pub od_quad_end: f64,
    pub od_quad_start_fraction: f64,
    pub od_quad_end_fraction: f64,
    pub saz_entrance: f64,
    pub saz_exit: f64,
    pub cos_sza_entrance: f64,
    pub cos_sza_exit: f64,
    pub layer_type: LayerType,
    pub cell: Option<CellId>,
}

impl Layer {
    pub fn new(entrance: TracePoint, exit: TracePoint, cell: Option<CellId>) -> Self {
        let delta = exit.position - entrance.position;
        let geometric_distance = delta.norm();
        let average_look_away = if geometric_distance == 0.0 {
            Vec3::ZERO
        } else {
            delta / geometric_distance
        };

        let layer_type = if entrance.event.is_tangent() || exit.event.is_tangent() {
            LayerType::Tangent
        } else if entrance.event.is_boundary()
            && exit.event.is_boundary()
            && entrance.event.boundaries.contains_vertical()
            && exit.event.boundaries.contains_vertical()
        {
            LayerType::Complete
        } else {
            LayerType::Partial
        };

        Self {
            entrance,
            exit,
            average_look_away,
            geometric_distance,
            curvature_factor: 1.0,
            od_quad_start: 0.0,
            od_quad_end: 0.0,
            od_quad_start_fraction: 0.0,
            od_quad_end_fraction: 0.0,
            saz_entrance: 0.0,
            saz_exit: 0.0,
            cos_sza_entrance: 0.0,
            cos_sza_exit: 0.0,
            layer_type,
            cell,
        }
    }

    #[inline(always)]
    pub fn effective_distance(&self) -> f64 {
        self.geometric_distance * self.curvature_factor
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TracedRay {
    pub observer: super::ray::Ray,
    pub is_straight: bool,
    pub ground_is_hit: bool,
    pub layers: Vec<Layer>,
    pub tangent_radius: f64,
    pub geometry: GeometryKind,
    pub interpolation_method: InterpolationMethod,
}

impl TracedRay {
    pub fn new(
        observer: super::ray::Ray,
        geometry: GeometryKind,
        interpolation_method: InterpolationMethod,
    ) -> Self {
        Self {
            observer,
            is_straight: true,
            ground_is_hit: false,
            layers: Vec::new(),
            tangent_radius: f64::NAN,
            geometry,
            interpolation_method,
        }
    }

    pub fn reset(
        &mut self,
        observer: super::ray::Ray,
        geometry: GeometryKind,
        interpolation_method: InterpolationMethod,
    ) {
        self.observer = observer;
        self.is_straight = true;
        self.ground_is_hit = false;
        self.layers.clear();
        self.tangent_radius = f64::NAN;
        self.geometry = geometry;
        self.interpolation_method = interpolation_method;
    }
}
