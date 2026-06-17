use super::grid::{BoundarySet, BoundaryTag, GeometryKind, MediumGrid, VerticalGrid1D};
use super::layer::{Layer, TraceEvent, TraceEventKind, TracePoint, TracedRay};
use super::od_quadrature::add_od_quadrature;
use super::primitive::{Intersection, Primitive};
use super::ray::Ray;
use super::solar::{SolarContext, add_solar_parameters};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TraceOptions {
    pub solar: Option<SolarContext>,
}

#[derive(Debug, Clone, Copy)]
struct EventSeed {
    distance: f64,
    kind: TraceEventKind,
    boundaries: BoundarySet,
}

#[derive(Debug, Default)]
pub struct TraceScratch {
    intersections: Vec<Intersection>,
    events: Vec<EventSeed>,
}

impl TraceScratch {
    pub fn with_capacity(num_primitives: usize) -> Self {
        Self {
            intersections: Vec::with_capacity(num_primitives * 2),
            events: Vec::with_capacity(num_primitives * 2 + 2),
        }
    }

    #[inline]
    pub fn clear(&mut self) {
        self.intersections.clear();
        self.events.clear();
    }
}

#[derive(Debug, Clone)]
pub struct VerticalRayTracer {
    grid: VerticalGrid1D,
    primitives: Vec<Primitive>,
    epsilon: f64,
}

impl VerticalRayTracer {
    pub fn new(grid: VerticalGrid1D) -> Self {
        let primitives = grid.vertical_primitives();
        Self {
            grid,
            primitives,
            epsilon: 1e-9,
        }
    }

    pub fn with_epsilon(mut self, epsilon: f64) -> Self {
        self.epsilon = epsilon;
        self
    }

    #[inline(always)]
    pub fn grid(&self) -> &VerticalGrid1D {
        &self.grid
    }

    #[inline(always)]
    pub fn primitives(&self) -> &[Primitive] {
        &self.primitives
    }

    pub fn trace(&self, ray: Ray, options: TraceOptions) -> TracedRay {
        let mut result =
            TracedRay::new(ray, self.grid.geometry(), self.grid.interpolation_method());
        let mut scratch = TraceScratch::with_capacity(self.primitives.len());
        self.trace_into(ray, &mut result, &mut scratch, options);
        result
    }

    pub fn trace_into(
        &self,
        ray: Ray,
        result: &mut TracedRay,
        scratch: &mut TraceScratch,
        options: TraceOptions,
    ) {
        scratch.clear();
        result.reset(ray, self.grid.geometry(), self.grid.interpolation_method());
        result.tangent_radius = match self.grid.geometry() {
            GeometryKind::Spherical => ray.spherical_tangent_radius(),
            GeometryKind::PlaneParallel => f64::NAN,
        };

        for primitive in &self.primitives {
            primitive.intersections(ray, self.epsilon, &mut scratch.intersections);
        }
        scratch
            .intersections
            .sort_unstable_by(|lhs, rhs| lhs.distance.total_cmp(&rhs.distance));

        let observer_altitude = self.grid.altitude_at(ray.origin);
        let observer_inside = self
            .grid
            .is_inside_altitude(observer_altitude, self.epsilon);

        let surface_distance =
            first_distance_with_tag(&scratch.intersections, BoundaryTag::is_surface);
        let top_distances =
            distances_with_tag(&scratch.intersections, BoundaryTag::is_top_of_atmosphere);

        let Some((start_distance, end_distance, ground_is_hit)) =
            self.atmospheric_interval(observer_inside, surface_distance, &top_distances)
        else {
            return;
        };
        result.ground_is_hit = ground_is_hit;

        if observer_inside {
            add_event_seed(
                &mut scratch.events,
                0.0,
                TraceEventKind::Observer,
                None,
                self.epsilon,
            );
        }

        for intersection in &scratch.intersections {
            if intersection.distance < start_distance - self.epsilon
                || intersection.distance > end_distance + self.epsilon
            {
                continue;
            }
            add_event_seed(
                &mut scratch.events,
                intersection.distance,
                TraceEventKind::Boundary,
                Some(self.grid.classify_intersection(intersection)),
                self.epsilon,
            );
        }

        self.add_tangent_event_if_needed(ray, start_distance, end_distance, ground_is_hit, scratch);

        scratch
            .events
            .sort_unstable_by(|lhs, rhs| lhs.distance.total_cmp(&rhs.distance));

        let mut points = Vec::with_capacity(scratch.events.len());
        for event in &scratch.events {
            if event.distance < start_distance - self.epsilon
                || event.distance > end_distance + self.epsilon
            {
                continue;
            }
            points.push(self.trace_point(ray, *event));
        }

        if points.len() < 2 {
            return;
        }

        result.layers.reserve(points.len() - 1);
        for window in points.windows(2).rev() {
            let entrance = window[0];
            let exit = window[1];
            if (exit.distance - entrance.distance).abs() <= self.epsilon {
                continue;
            }
            let midpoint = ray.point_at((entrance.distance + exit.distance) / 2.0);
            let cell = self.grid.locate_cell(midpoint);
            let mut layer = Layer::new(entrance, exit, cell);
            add_od_quadrature(
                &mut layer,
                self.grid.geometry(),
                self.grid.interpolation_method(),
            );
            if let Some(solar) = options.solar {
                add_solar_parameters(&mut layer, solar);
            }
            result.layers.push(layer);
        }
    }

    fn atmospheric_interval(
        &self,
        observer_inside: bool,
        surface_distance: Option<f64>,
        top_distances: &[f64],
    ) -> Option<(f64, f64, bool)> {
        if let Some(ground_distance) = surface_distance {
            let start_distance = if observer_inside {
                0.0
            } else {
                *top_distances
                    .iter()
                    .find(|distance| **distance <= ground_distance + self.epsilon)?
            };
            return Some((start_distance, ground_distance, true));
        }

        if observer_inside {
            let end_distance = *top_distances
                .iter()
                .find(|distance| **distance > self.epsilon)?;
            return Some((0.0, end_distance, false));
        }

        if top_distances.len() >= 2 {
            Some((top_distances[0], top_distances[1], false))
        } else {
            None
        }
    }

    fn add_tangent_event_if_needed(
        &self,
        ray: Ray,
        start_distance: f64,
        end_distance: f64,
        ground_is_hit: bool,
        scratch: &mut TraceScratch,
    ) {
        if self.grid.geometry() != GeometryKind::Spherical || ground_is_hit {
            return;
        }

        let tangent_distance = ray.spherical_tangent_distance();
        if tangent_distance <= start_distance + self.epsilon
            || tangent_distance >= end_distance - self.epsilon
        {
            return;
        }

        let tangent_altitude = ray.spherical_tangent_radius() - self.grid.earth_radius();
        if tangent_altitude <= self.grid.ground_altitude() + self.epsilon
            || tangent_altitude >= self.grid.top_altitude() - self.epsilon
        {
            return;
        }

        add_event_seed(
            &mut scratch.events,
            tangent_distance,
            TraceEventKind::Tangent,
            None,
            self.epsilon,
        );
    }

    fn trace_point(&self, ray: Ray, event: EventSeed) -> TracePoint {
        let distance = event.distance;
        let position = ray.point_at(distance);
        let altitude = self.grid.altitude_at(position);
        let event = match event.kind {
            TraceEventKind::Observer => TraceEvent::observer(),
            TraceEventKind::Boundary => TraceEvent::boundary(event.boundaries),
            TraceEventKind::Tangent => TraceEvent::tangent(event.boundaries),
        };

        TracePoint {
            distance,
            position,
            altitude,
            event,
            interpolation: self.grid.interpolation_weights(position),
            cell: self.grid.locate_cell(position),
            on_exact_vertical_boundary: self
                .grid
                .on_exact_vertical_boundary(altitude, 1e-8)
                .is_some(),
        }
    }
}

fn first_distance_with_tag(
    intersections: &[Intersection],
    predicate: impl Fn(BoundaryTag) -> bool,
) -> Option<f64> {
    intersections
        .iter()
        .find(|intersection| predicate(intersection.boundary))
        .map(|intersection| intersection.distance)
}

fn distances_with_tag(
    intersections: &[Intersection],
    predicate: impl Fn(BoundaryTag) -> bool,
) -> Vec<f64> {
    intersections
        .iter()
        .filter(|intersection| predicate(intersection.boundary))
        .map(|intersection| intersection.distance)
        .collect()
}

fn add_event_seed(
    events: &mut Vec<EventSeed>,
    distance: f64,
    kind: TraceEventKind,
    boundary: Option<BoundaryTag>,
    epsilon: f64,
) {
    if let Some(existing) = events
        .iter_mut()
        .find(|event| (event.distance - distance).abs() <= epsilon)
    {
        if kind == TraceEventKind::Tangent {
            existing.kind = TraceEventKind::Tangent;
        } else if existing.kind != TraceEventKind::Tangent {
            existing.kind = kind;
        }
        if let Some(boundary) = boundary {
            existing.boundaries.push(boundary);
        }
        return;
    }

    let mut boundaries = BoundarySet::new();
    if let Some(boundary) = boundary {
        boundaries.push(boundary);
    }
    events.push(EventSeed {
        distance,
        kind,
        boundaries,
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raytracer::grid::{GeometryKind, InterpolationMethod, VerticalGrid1D};
    use crate::raytracer::layer::LayerType;
    use crate::raytracer::vec3::Vec3;

    fn spherical_tracer() -> VerticalRayTracer {
        let grid = VerticalGrid1D::new(
            10.0,
            vec![0.0, 10.0, 20.0, 30.0],
            InterpolationMethod::Linear,
            GeometryKind::Spherical,
        )
        .unwrap();
        VerticalRayTracer::new(grid)
    }

    fn plane_tracer() -> VerticalRayTracer {
        let grid = VerticalGrid1D::new(
            10.0,
            vec![0.0, 10.0, 20.0, 30.0],
            InterpolationMethod::Linear,
            GeometryKind::PlaneParallel,
        )
        .unwrap();
        VerticalRayTracer::new(grid)
    }

    fn min_endpoint_altitude(ray: &TracedRay) -> f64 {
        ray.layers.iter().fold(f64::INFINITY, |min_alt, layer| {
            min_alt
                .min(layer.entrance.altitude)
                .min(layer.exit.altitude)
        })
    }

    #[test]
    fn spherical_observer_outside_ground_viewing_traces_to_surface() {
        let tracer = spherical_tracer();
        let ray = Ray::new(Vec3::new(0.0, 0.0, 50.0), Vec3::new(0.0, 0.0, -1.0));

        let traced = tracer.trace(ray, TraceOptions::default());

        assert!(traced.ground_is_hit);
        assert_eq!(traced.layers.len(), 3);
        assert!((min_endpoint_altitude(&traced) - 0.0).abs() < 1e-8);
        assert_eq!(traced.layers[0].layer_type, LayerType::Complete);
    }

    #[test]
    fn spherical_observer_outside_limb_viewing_inserts_tangent_event() {
        let tracer = spherical_tracer();
        let tangent_radius: f64 = 25.0;
        let origin_radius: f64 = 50.0;
        let ratio = tangent_radius / origin_radius;
        let dz = -(1.0_f64 - ratio * ratio).sqrt();
        let dx = ratio;
        let ray = Ray::new(Vec3::new(0.0, 0.0, origin_radius), Vec3::new(dx, 0.0, dz));

        let traced = tracer.trace(ray, TraceOptions::default());

        assert!(!traced.ground_is_hit);
        assert_eq!(traced.layers.len(), 4);
        assert!((min_endpoint_altitude(&traced) - 15.0).abs() < 1e-8);
        assert!(
            traced
                .layers
                .iter()
                .any(|layer| layer.layer_type == LayerType::Tangent)
        );
    }

    #[test]
    fn spherical_observer_inside_looking_up_starts_with_partial_layer() {
        let tracer = spherical_tracer();
        let ray = Ray::new(Vec3::new(0.0, 0.0, 25.0), Vec3::new(0.0, 0.0, 1.0));

        let traced = tracer.trace(ray, TraceOptions::default());

        assert!(!traced.ground_is_hit);
        assert_eq!(traced.layers.len(), 2);
        assert!(
            traced
                .layers
                .iter()
                .any(|layer| layer.layer_type == LayerType::Partial)
        );
    }

    #[test]
    fn plane_parallel_observer_outside_ground_viewing_traces_to_surface() {
        let tracer = plane_tracer();
        let ray = Ray::new(Vec3::new(0.0, 0.0, 50.0), Vec3::new(0.0, 0.0, -1.0));

        let traced = tracer.trace(ray, TraceOptions::default());

        assert!(traced.ground_is_hit);
        assert_eq!(traced.layers.len(), 3);
        assert!((min_endpoint_altitude(&traced) - 0.0).abs() < 1e-8);
    }
}
