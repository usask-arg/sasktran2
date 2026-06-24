use super::grid::{
    BoundarySet, BoundaryTag, CellId, GeometryKind, InterpolationStencil, MediumGrid,
    VerticalGrid1D,
};
use super::layer::{Layer, TraceEvent, TraceEventKind, TracePoint, TracedRay};
use super::od_quadrature::add_od_quadrature;
use super::primitive::{Intersection, Primitive};
use super::ray::Ray;
use super::refraction::RefractiveProfile;
use super::solar::{SolarContext, add_solar_parameters};
use super::vec3::Vec3;

pub const NADIR_VIEWING_CUTOFF: f64 = 0.999999;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TraceOptions<'a> {
    pub solar: Option<SolarContext>,
    pub refraction: Option<&'a RefractiveProfile>,
}

#[derive(Debug, Clone, Copy)]
struct EventSeed {
    distance: f64,
    altitude: f64,
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

    pub fn trace(&self, ray: Ray, options: TraceOptions<'_>) -> TracedRay {
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
        options: TraceOptions<'_>,
    ) {
        scratch.clear();
        result.reset(ray, self.grid.geometry(), self.grid.interpolation_method());

        let straight_tangent_radius = match self.grid.geometry() {
            GeometryKind::Spherical => ray.spherical_tangent_radius(),
            GeometryKind::PlaneParallel => f64::NAN,
        };
        let refraction_profile = self.active_refraction(ray, options.refraction);
        result.is_straight = refraction_profile.is_none();
        result.tangent_radius = match self.grid.geometry() {
            GeometryKind::Spherical => refraction_profile
                .map(|profile| profile.tangent_radius(straight_tangent_radius))
                .unwrap_or(straight_tangent_radius),
            GeometryKind::PlaneParallel => f64::NAN,
        };

        match self.grid.geometry() {
            GeometryKind::Spherical => {
                if !self.build_spherical_events(
                    ray,
                    result.tangent_radius,
                    straight_tangent_radius,
                    scratch,
                    result,
                ) {
                    return;
                }
                self.build_spherical_layers_from_events(scratch, result);
                self.finalize_spherical_layers(ray, result, refraction_profile);
            }
            GeometryKind::PlaneParallel => {
                if !self.build_plane_parallel_events(ray, scratch, result) {
                    return;
                }
                self.build_plane_parallel_layers_from_events(ray, scratch, result);
            }
        }

        for layer in &mut result.layers {
            add_od_quadrature(
                layer,
                self.grid.geometry(),
                self.grid.interpolation_method(),
            );
            if let Some(solar) = options.solar {
                add_solar_parameters(layer, solar);
            }
        }
    }

    fn active_refraction<'a>(
        &self,
        ray: Ray,
        profile: Option<&'a RefractiveProfile>,
    ) -> Option<&'a RefractiveProfile> {
        if self.grid.geometry() != GeometryKind::Spherical {
            return None;
        }
        if ray.cos_zenith_at_origin().abs() > NADIR_VIEWING_CUTOFF {
            return None;
        }
        profile
    }

    fn build_spherical_events(
        &self,
        ray: Ray,
        tangent_radius: f64,
        straight_tangent_radius: f64,
        scratch: &mut TraceScratch,
        result: &mut TracedRay,
    ) -> bool {
        let observer_altitude = self.grid.altitude_at(ray.origin);
        let observer_outside = observer_altitude >= self.grid.top_altitude();
        let cos_viewing = ray.cos_zenith_at_origin();
        let tangent_altitude =
            self.snap_to_vertical_boundary(tangent_radius - self.grid.earth_radius());
        let straight_tangent_altitude =
            self.snap_to_vertical_boundary(straight_tangent_radius - self.grid.earth_radius());

        if tangent_altitude >= self.grid.top_altitude() || (observer_outside && cos_viewing > 0.0) {
            result.ground_is_hit = false;
            return false;
        }

        let mut sequence = 0.0;
        let altitudes = self.grid.altitudes();

        if observer_outside {
            if tangent_altitude > self.grid.ground_altitude() {
                result.ground_is_hit = false;
                let above_tangent_idx = upper_bound(altitudes, tangent_altitude);
                if above_tangent_idx >= altitudes.len() {
                    return false;
                }
                for index in (above_tangent_idx..altitudes.len()).rev() {
                    self.push_boundary_event(scratch, &mut sequence, index);
                }
                self.push_tangent_event(scratch, &mut sequence, tangent_altitude);
                for index in above_tangent_idx..altitudes.len() {
                    self.push_boundary_event(scratch, &mut sequence, index);
                }
            } else {
                result.ground_is_hit = true;
                for index in (0..altitudes.len()).rev() {
                    self.push_boundary_event(scratch, &mut sequence, index);
                }
            }
            return scratch.events.len() >= 2;
        }

        self.push_observer_event(scratch, &mut sequence, observer_altitude);

        if cos_viewing > 0.0 {
            result.ground_is_hit = false;
            let start_index = upper_bound(altitudes, observer_altitude);
            for index in start_index..altitudes.len() {
                self.push_boundary_event(scratch, &mut sequence, index);
            }
        } else if tangent_altitude > self.grid.ground_altitude() {
            result.ground_is_hit = false;
            let observer_index = upper_bound(altitudes, observer_altitude);
            // C++ chooses the limb branch from the refracted tangent radius,
            // then sequences inside-limb layer boundaries from the straight
            // tangent altitude.
            let above_tangent_idx = upper_bound(altitudes, straight_tangent_altitude);
            for index in (above_tangent_idx..observer_index).rev() {
                self.push_boundary_event(scratch, &mut sequence, index);
            }
            self.push_tangent_event(scratch, &mut sequence, straight_tangent_altitude);
            for index in above_tangent_idx..altitudes.len() {
                self.push_boundary_event(scratch, &mut sequence, index);
            }
        } else {
            result.ground_is_hit = true;
            let start_index = upper_bound(altitudes, observer_altitude);
            if start_index == 0 {
                return false;
            }
            for index in (0..start_index).rev() {
                self.push_boundary_event(scratch, &mut sequence, index);
            }
        }

        scratch.events.len() >= 2
    }

    fn build_spherical_layers_from_events(&self, scratch: &TraceScratch, result: &mut TracedRay) {
        if scratch.events.len() < 2 {
            return;
        }

        result.layers.reserve(scratch.events.len() - 1);
        for window in scratch.events.windows(2).rev() {
            let entrance = self.trace_spherical_event_point(window[0]);
            let exit = self.trace_spherical_event_point(window[1]);
            if (exit.distance - entrance.distance).abs() <= self.epsilon {
                continue;
            }
            result.layers.push(Layer::unfinalized(entrance, exit, None));
        }
    }

    fn build_plane_parallel_layers_from_events(
        &self,
        ray: Ray,
        scratch: &mut TraceScratch,
        result: &mut TracedRay,
    ) {
        scratch
            .events
            .sort_unstable_by(|lhs, rhs| lhs.distance.total_cmp(&rhs.distance));

        if scratch.events.len() < 2 {
            return;
        }

        result.layers.reserve(scratch.events.len() - 1);
        for window in scratch.events.windows(2).rev() {
            let entrance = self.trace_point(ray, window[0]);
            let exit = self.trace_point(ray, window[1]);
            if (exit.distance - entrance.distance).abs() <= self.epsilon {
                continue;
            }
            let midpoint = ray.point_at((entrance.distance + exit.distance) / 2.0);
            let cell = self.grid.locate_cell(midpoint);
            result.layers.push(Layer::new(entrance, exit, cell));
        }
    }

    fn build_plane_parallel_events(
        &self,
        ray: Ray,
        scratch: &mut TraceScratch,
        result: &mut TracedRay,
    ) -> bool {
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
            return false;
        };
        result.ground_is_hit = ground_is_hit;

        if observer_inside {
            add_event_seed(
                &mut scratch.events,
                0.0,
                observer_altitude,
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
                self.grid.altitude_at(intersection.point),
                TraceEventKind::Boundary,
                Some(self.grid.classify_intersection(intersection)),
                self.epsilon,
            );
        }

        scratch.events.len() >= 2
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

    fn push_boundary_event(&self, scratch: &mut TraceScratch, sequence: &mut f64, index: usize) {
        scratch.events.push(EventSeed {
            distance: *sequence,
            altitude: self.grid.altitudes()[index],
            kind: TraceEventKind::Boundary,
            boundaries: boundary_set(self.grid.boundary_tag(index)),
        });
        *sequence += 1.0;
    }

    fn push_observer_event(
        &self,
        scratch: &mut TraceScratch,
        sequence: &mut f64,
        observer_altitude: f64,
    ) {
        scratch.events.push(EventSeed {
            distance: *sequence,
            altitude: observer_altitude,
            kind: TraceEventKind::Observer,
            boundaries: BoundarySet::new(),
        });
        *sequence += 1.0;
    }

    fn push_tangent_event(
        &self,
        scratch: &mut TraceScratch,
        sequence: &mut f64,
        tangent_altitude: f64,
    ) {
        let boundaries = self
            .grid
            .on_exact_vertical_boundary(tangent_altitude, self.epsilon)
            .map(|index| boundary_set(self.grid.boundary_tag(index)))
            .unwrap_or_default();
        scratch.events.push(EventSeed {
            distance: *sequence,
            altitude: tangent_altitude,
            kind: TraceEventKind::Tangent,
            boundaries,
        });
        *sequence += 1.0;
    }

    fn finalize_spherical_layers(
        &self,
        ray: Ray,
        result: &mut TracedRay,
        refraction_profile: Option<&RefractiveProfile>,
    ) {
        let mut total_deflection = 0.0;
        let tangent_radius = result.tangent_radius;
        let refraction = refraction_profile
            .map(|profile| (profile, profile.refractive_index_at_radius(tangent_radius)));
        let mut x_basis = Vec3::ZERO;
        let mut y_basis = Vec3::ZERO;
        let mut previous_exit = None;
        let observer_inside = self.grid.altitude_at(ray.origin) < self.grid.top_altitude();

        for layer_index in (0..result.layers.len()).rev() {
            let layer = &mut result.layers[layer_index];

            if let Some(exit) = previous_exit {
                layer.entrance = exit;
            } else {
                let entrance_position = if observer_inside {
                    ray.origin
                } else {
                    ray.point_at(self.distance_to_spherical_altitude(
                        ray,
                        self.grid.top_altitude(),
                        1.0,
                        1.0,
                    ))
                };
                self.update_point_position(&mut layer.entrance, entrance_position);
                x_basis = layer.entrance.position.normalized();
                y_basis = x_basis
                    .cross(ray.direction)
                    .normalized()
                    .cross(x_basis)
                    .normalized();
            }

            let r_entrance = self.grid.radius_for_altitude(layer.entrance.altitude);
            let r_exit = self.grid.radius_for_altitude(layer.exit.altitude);

            if let Some((profile, tangent_refractive_index)) = refraction {
                let refracted_path = profile.integrate_path_with_tangent_index(
                    tangent_radius,
                    tangent_refractive_index,
                    r_entrance,
                    r_exit,
                );
                total_deflection += refracted_path.deflection_angle;
                let exit_position =
                    r_exit * (x_basis * total_deflection.cos() + y_basis * total_deflection.sin());
                self.update_point_position(&mut layer.exit, exit_position);
                let delta = layer.exit.position - layer.entrance.position;
                layer.geometric_distance = delta.norm();
                layer.average_look_away = if layer.geometric_distance == 0.0 {
                    Vec3::ZERO
                } else {
                    delta / layer.geometric_distance
                };
                layer.curvature_factor = if layer.geometric_distance == 0.0 {
                    1.0
                } else {
                    refracted_path.path_length / layer.geometric_distance
                };
            } else {
                layer.curvature_factor = 1.0;
                layer.geometric_distance =
                    spherical_shell_distance(tangent_radius, r_entrance, r_exit);
                layer.average_look_away = ray.direction;
                let exit_position =
                    layer.entrance.position + layer.average_look_away * layer.geometric_distance;
                self.update_point_position(&mut layer.exit, exit_position);
            }

            layer.cell = self
                .grid
                .locate_altitude_layer((layer.entrance.altitude + layer.exit.altitude) / 2.0, 1e-8)
                .map(super::grid::CellId::AltitudeLayer);

            previous_exit = Some(layer.exit);
        }
    }

    fn update_point_position(&self, point: &mut TracePoint, position: super::vec3::Vec3) {
        point.position = position;
        if point.event.is_observer() || !point.event.boundaries.contains_vertical() {
            point.altitude = self.grid.altitude_at(position);
        }

        if let Some(index) = vertical_boundary_index(point.event.boundaries) {
            point.interpolation = exact_boundary_stencil(index);
            point.cell = Some(CellId::AltitudeLayer(
                index.min(self.grid.altitudes().len() - 2),
            ));
            point.on_exact_vertical_boundary = true;
        } else {
            point.interpolation = self.grid.interpolation_weights_at_altitude(point.altitude);
            point.cell = self
                .grid
                .locate_altitude_layer(point.altitude, 1e-8)
                .map(CellId::AltitudeLayer);
            point.on_exact_vertical_boundary = self
                .grid
                .on_exact_vertical_boundary(point.altitude, 1e-8)
                .is_some();
        }
    }

    fn distance_to_spherical_altitude(
        &self,
        ray: Ray,
        altitude: f64,
        direction: f64,
        side: f64,
    ) -> f64 {
        let cos_zenith = ray.cos_zenith_at_origin().abs();
        let ro = ray.origin.norm();
        let re = self.grid.earth_radius() + altitude;
        let rtsq = ro * ro * (1.0 - cos_zenith * cos_zenith);
        let tangent_distance = side * direction * ro * cos_zenith;
        let dist_from_tangent = if rtsq > re * re && (rtsq - re * re).abs() < 100.0 {
            0.0
        } else {
            side * direction * (re * re - rtsq).abs().sqrt()
        };

        if side > 0.0 {
            tangent_distance - dist_from_tangent
        } else {
            tangent_distance + dist_from_tangent
        }
    }

    fn snap_to_vertical_boundary(&self, altitude: f64) -> f64 {
        self.grid
            .on_exact_vertical_boundary(altitude, self.epsilon.max(1e-6))
            .map(|index| self.grid.altitudes()[index])
            .unwrap_or(altitude)
    }

    fn trace_spherical_event_point(&self, event: EventSeed) -> TracePoint {
        let event_kind = match event.kind {
            TraceEventKind::Observer => TraceEvent::observer(),
            TraceEventKind::Boundary => TraceEvent::boundary(event.boundaries),
            TraceEventKind::Tangent => TraceEvent::tangent(event.boundaries),
        };

        TracePoint {
            distance: event.distance,
            position: Vec3::ZERO,
            altitude: event.altitude,
            event: event_kind,
            interpolation: Default::default(),
            cell: None,
            on_exact_vertical_boundary: event.boundaries.contains_vertical()
                || self
                    .grid
                    .on_exact_vertical_boundary(event.altitude, 1e-8)
                    .is_some(),
        }
    }

    fn trace_point(&self, ray: Ray, event: EventSeed) -> TracePoint {
        let distance = event.distance;
        let position = ray.point_at(distance);
        let altitude = event.altitude;
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
            interpolation: self.grid.interpolation_weights_at_altitude(altitude),
            cell: self
                .grid
                .locate_altitude_layer(altitude, 1e-8)
                .map(super::grid::CellId::AltitudeLayer),
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
    altitude: f64,
    kind: TraceEventKind,
    boundary: Option<BoundaryTag>,
    epsilon: f64,
) {
    if let Some(existing) = events
        .iter_mut()
        .find(|event| (event.distance - distance).abs() <= epsilon)
    {
        existing.altitude = altitude;
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
        altitude,
        kind,
        boundaries,
    });
}

fn boundary_set(boundary: BoundaryTag) -> BoundarySet {
    let mut boundaries = BoundarySet::new();
    boundaries.push(boundary);
    boundaries
}

fn vertical_boundary_index(boundaries: BoundarySet) -> Option<usize> {
    boundaries.iter().find_map(BoundaryTag::altitude_index)
}

fn exact_boundary_stencil(index: usize) -> InterpolationStencil {
    let mut stencil = InterpolationStencil::new();
    stencil.push(index, 1.0);
    stencil
}

fn upper_bound(values: &[f64], value: f64) -> usize {
    values.partition_point(|candidate| *candidate <= value)
}

fn spherical_shell_distance(tangent_radius: f64, r_entrance: f64, r_exit: f64) -> f64 {
    ((r_entrance * r_entrance - tangent_radius * tangent_radius)
        .max(0.0)
        .sqrt()
        - (r_exit * r_exit - tangent_radius * tangent_radius)
            .max(0.0)
            .sqrt())
    .abs()
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
