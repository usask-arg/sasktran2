//! Straight and altitude-only refracted tracing through [`StructuredGrid2D`].
//!
//! This is currently a Rust-only API. It is intentionally not connected to the
//! CXX bridge, Python bindings, or the engine's default raytracer selection.

use super::grid::{BoundarySet, BoundaryTag, GeometryKind, MediumGrid, VerticalGrid1D};
use super::grid2d::StructuredGrid2D;
use super::layer::{Layer, TraceEvent, TraceEventKind, TracePoint, TracedRay};
use super::od_quadrature::add_od_quadrature;
use super::primitive::{Crossing, Intersection, Primitive};
use super::ray::Ray;
use super::refraction::{PathIntegral, RefractiveProfile};
use super::solar::{SolarContext, add_solar_parameters_reusing_exit};
use super::trace::{NADIR_VIEWING_CUTOFF, TraceOptions, TraceScratch, VerticalRayTracer};
use super::vec3::Vec3;

/// Optional per-ray inputs for structured 2D tracing.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TraceOptions2D<'a> {
    pub solar: Option<SolarContext>,
    /// Altitude-only refractive profile used for this ray.
    ///
    /// The profile is borrowed for the duration of the trace, allowing every
    /// ray to use a different column without rebuilding the 2D tracer. Its
    /// earth radius must match the structured grid's earth radius.
    pub refraction: Option<&'a RefractiveProfile>,
}

#[derive(Debug, Clone, Copy)]
struct Event2D {
    distance: f64,
    boundaries: BoundarySet,
    is_observer: bool,
    is_tangent: bool,
}

#[derive(Debug, Clone, Copy)]
struct CurvedEvent2D {
    angle: f64,
    path_distance: f64,
    point: TracePoint,
}

/// Reusable allocations for [`StructuredRayTracer2D::trace_into`].
#[derive(Debug, Default)]
pub struct TraceScratch2D {
    intersections: Vec<Intersection>,
    events: Vec<Event2D>,
    curved_events: Vec<CurvedEvent2D>,
    radial_result: Option<TracedRay>,
    radial_scratch: TraceScratch,
}

impl TraceScratch2D {
    pub fn with_capacity(num_primitives: usize) -> Self {
        Self {
            intersections: Vec::with_capacity(num_primitives * 2),
            events: Vec::with_capacity(num_primitives * 2 + 2),
            curved_events: Vec::with_capacity(num_primitives + 2),
            radial_result: None,
            radial_scratch: TraceScratch::with_capacity(num_primitives),
        }
    }

    #[inline]
    pub fn clear(&mut self) {
        self.intersections.clear();
        self.events.clear();
        self.curved_events.clear();
        self.radial_scratch.clear();
    }
}

/// Ray tracer for a finite spherical altitude-by-angle grid.
///
/// Refracted rays assume spherical symmetry for each individual trace: the
/// supplied refractive profile may differ between rays, but varies only with
/// altitude along that ray.
#[derive(Debug, Clone)]
pub struct StructuredRayTracer2D {
    grid: StructuredGrid2D,
    primitives: Vec<Primitive>,
    radial_tracer: VerticalRayTracer,
    epsilon: f64,
}

impl StructuredRayTracer2D {
    pub fn new(grid: StructuredGrid2D) -> Self {
        let primitives = grid.primitives();
        let radial_grid = VerticalGrid1D::new(
            grid.earth_radius(),
            grid.altitudes().to_vec(),
            grid.altitude_interpolation(),
            GeometryKind::Spherical,
        )
        .expect("a validated 2D grid always defines a valid radial grid");
        Self {
            grid,
            primitives,
            radial_tracer: VerticalRayTracer::new(radial_grid),
            epsilon: 1e-9,
        }
    }

    pub fn with_epsilon(mut self, epsilon: f64) -> Self {
        self.epsilon = epsilon;
        self.radial_tracer = self.radial_tracer.with_epsilon(epsilon);
        self
    }

    #[inline(always)]
    pub fn grid(&self) -> &StructuredGrid2D {
        &self.grid
    }

    #[inline(always)]
    pub fn primitives(&self) -> &[Primitive] {
        &self.primitives
    }

    pub fn trace(&self, ray: Ray, options: TraceOptions2D<'_>) -> TracedRay {
        let mut result = TracedRay::new(
            ray,
            GeometryKind::Spherical,
            self.grid.altitude_interpolation(),
        );
        let mut scratch = TraceScratch2D::with_capacity(self.primitives.len());
        self.trace_into(ray, &mut result, &mut scratch, options);
        result
    }

    pub fn trace_into(
        &self,
        ray: Ray,
        result: &mut TracedRay,
        scratch: &mut TraceScratch2D,
        options: TraceOptions2D<'_>,
    ) {
        scratch.clear();
        result.reset(
            ray,
            GeometryKind::Spherical,
            self.grid.altitude_interpolation(),
        );

        let active_refraction = options
            .refraction
            .filter(|_| ray.cos_zenith_at_origin().abs() <= NADIR_VIEWING_CUTOFF);
        if let Some(profile) = active_refraction {
            if profile
                .refractive_index()
                .iter()
                .all(|refractive_index| *refractive_index == 1.0)
            {
                // Preserve the 1D/C++ convention that a supplied profile marks
                // the result as refracted, but avoid quadrature dither for the
                // exactly straight unity-index path.
                self.trace_straight_geometry(ray, result, scratch);
                result.is_straight = false;
            } else {
                self.trace_refracted_geometry(ray, result, scratch, profile);
            }
        } else {
            self.trace_straight_geometry(ray, result, scratch);
        }

        let reuse_shared_boundary_solar = result.is_straight;
        let mut shared_boundary_solar = None;
        for layer in &mut result.layers {
            add_od_quadrature(
                layer,
                GeometryKind::Spherical,
                self.grid.altitude_interpolation(),
            );
            if let Some(solar) = options.solar {
                let exit_solar = if reuse_shared_boundary_solar {
                    shared_boundary_solar
                } else {
                    None
                };
                let entrance_solar = add_solar_parameters_reusing_exit(layer, solar, exit_solar);
                if reuse_shared_boundary_solar {
                    shared_boundary_solar = Some(entrance_solar);
                }
            }
        }
    }

    fn trace_straight_geometry(
        &self,
        ray: Ray,
        result: &mut TracedRay,
        scratch: &mut TraceScratch2D,
    ) {
        result.is_straight = true;
        result.tangent_radius = ray.spherical_tangent_radius();

        for primitive in &self.primitives {
            primitive.intersections(ray, self.epsilon, &mut scratch.intersections);
        }
        scratch.events.push(Event2D {
            distance: 0.0,
            boundaries: BoundarySet::new(),
            is_observer: true,
            is_tangent: false,
        });
        for intersection in &scratch.intersections {
            let mut boundaries = BoundarySet::new();
            boundaries.push(intersection.boundary);
            scratch.events.push(Event2D {
                distance: intersection.distance,
                boundaries,
                is_observer: false,
                is_tangent: intersection.crossing == Crossing::Tangent,
            });
        }

        let tangent_distance = ray.spherical_tangent_distance();
        if tangent_distance > self.epsilon {
            let tangent_position = ray.point_at(tangent_distance);
            if self.grid.locate_cell(tangent_position).is_some() {
                scratch.events.push(Event2D {
                    distance: tangent_distance,
                    boundaries: BoundarySet::new(),
                    is_observer: false,
                    is_tangent: true,
                });
            }
        }
        sort_and_merge_events(&mut scratch.events, self.epsilon);

        let surface_limit = scratch
            .intersections
            .iter()
            .filter(|intersection| {
                intersection.boundary.is_surface() && intersection.distance > self.epsilon
            })
            .map(|intersection| intersection.distance)
            .min_by(f64::total_cmp);

        let mut atmospheric_segment_started = false;
        let mut far_boundaries = BoundarySet::new();
        for window in scratch.events.windows(2) {
            let near = window[0];
            let far = window[1];
            if surface_limit.is_some_and(|limit| near.distance >= limit - self.epsilon) {
                break;
            }
            if far.distance - near.distance <= self.epsilon {
                continue;
            }

            let midpoint = ray.point_at((near.distance + far.distance) / 2.0);
            let Some(cell) = self.grid.locate_cell(midpoint) else {
                if atmospheric_segment_started {
                    break;
                }
                continue;
            };

            atmospheric_segment_started = true;
            let entrance = self.trace_point(ray, near);
            let exit = self.trace_point(ray, far);
            result.layers.push(Layer::new(entrance, exit, Some(cell)));
            far_boundaries = far.boundaries;

            if surface_limit.is_some_and(|limit| far.distance >= limit - self.epsilon) {
                break;
            }
        }

        result.ground_is_hit = !result.layers.is_empty() && far_boundaries.contains_surface();
        result.layers.reverse();
    }

    fn trace_refracted_geometry(
        &self,
        ray: Ray,
        result: &mut TracedRay,
        scratch: &mut TraceScratch2D,
        profile: &RefractiveProfile,
    ) {
        const ANGULAR_EVENT_EPSILON: f64 = 1e-8;

        let radial_result = scratch.radial_result.get_or_insert_with(|| {
            TracedRay::new(
                ray,
                GeometryKind::Spherical,
                self.grid.altitude_interpolation(),
            )
        });
        self.radial_tracer.trace_into(
            ray,
            radial_result,
            &mut scratch.radial_scratch,
            TraceOptions {
                solar: None,
                refraction: Some(profile),
            },
        );

        result.is_straight = radial_result.is_straight;
        result.tangent_radius = radial_result.tangent_radius;
        if radial_result.layers.is_empty() {
            return;
        }

        let ray_plane_normal = ray.origin.cross(ray.direction).normalized();
        if ray_plane_normal.norm_squared() == 0.0 {
            // Near-radial rays are normally routed through the straight path by
            // the viewing cutoff. Keep a defensive fallback for degenerate input.
            self.trace_straight_geometry(ray, result, scratch);
            return;
        }

        let tangent_index = profile.refractive_index_at_radius(result.tangent_radius);
        let mut cumulative_path = 0.0;
        let mut atmospheric_segment_started = false;
        let mut far_boundaries = BoundarySet::new();

        'radial_layers: for radial_layer in radial_result.layers.iter().rev() {
            scratch.curved_events.clear();

            let start = self.structured_point(radial_layer.entrance, cumulative_path);
            let layer_path = radial_layer.effective_distance();
            let end = self.structured_point(radial_layer.exit, cumulative_path + layer_path);
            let start_direction = start.position.normalized();
            let end_direction = end.position.normalized();
            let total_angle = forward_angle(
                start_direction,
                end_direction,
                ray_plane_normal,
                self.epsilon,
            );

            scratch.curved_events.push(CurvedEvent2D {
                angle: 0.0,
                path_distance: 0.0,
                point: start,
            });

            if total_angle > ANGULAR_EVENT_EPSILON {
                for (horizontal_index, horizontal_angle) in
                    self.grid.horizontal_angles().iter().enumerate()
                {
                    let Some(crossing_direction) = angular_crossing_direction(
                        ray_plane_normal,
                        self.grid.basis().angular_normal(*horizontal_angle),
                        self.grid.basis().radial_direction(*horizontal_angle),
                        self.epsilon,
                    ) else {
                        continue;
                    };
                    let crossing_angle = forward_angle(
                        start_direction,
                        crossing_direction,
                        ray_plane_normal,
                        self.epsilon,
                    );
                    if crossing_angle <= ANGULAR_EVENT_EPSILON
                        || crossing_angle >= total_angle - ANGULAR_EVENT_EPSILON
                    {
                        continue;
                    }

                    let integral = refracted_integral_at_angle(
                        profile,
                        result.tangent_radius,
                        tangent_index,
                        start.position.norm(),
                        end.position.norm(),
                        crossing_angle,
                    );
                    let radius = integral.radius;
                    let position = radius * crossing_direction;
                    let mut boundaries = BoundarySet::new();
                    boundaries.push(BoundaryTag::Horizontal {
                        index: horizontal_index,
                    });
                    scratch.curved_events.push(CurvedEvent2D {
                        angle: crossing_angle,
                        path_distance: integral.integral.path_length,
                        point: TracePoint {
                            distance: cumulative_path + integral.integral.path_length,
                            position,
                            altitude: self.grid.altitude_at(position),
                            event: TraceEvent::boundary(boundaries),
                            interpolation: self.grid.interpolation_weights_at(position),
                            cell: self.grid.locate_cell(position),
                            on_exact_vertical_boundary: false,
                        },
                    });
                }
            }

            scratch.curved_events.push(CurvedEvent2D {
                angle: total_angle,
                path_distance: layer_path,
                point: end,
            });
            scratch
                .curved_events
                .sort_unstable_by(|lhs, rhs| lhs.angle.total_cmp(&rhs.angle));

            for window in scratch.curved_events.windows(2) {
                let near = window[0];
                let far = window[1];
                let effective_distance = far.path_distance - near.path_distance;
                if effective_distance <= self.epsilon {
                    continue;
                }

                let midpoint_angle = (near.angle + far.angle) / 2.0;
                let midpoint_direction =
                    rotate_in_plane(start_direction, ray_plane_normal, midpoint_angle);
                let midpoint_radius =
                    (near.point.position.norm() + far.point.position.norm()) / 2.0;
                let midpoint = midpoint_radius * midpoint_direction;
                let Some(cell) = self.grid.locate_cell(midpoint) else {
                    if atmospheric_segment_started {
                        break 'radial_layers;
                    }
                    continue;
                };

                atmospheric_segment_started = true;
                let mut layer = Layer::new(near.point, far.point, Some(cell));
                layer.curvature_factor = if layer.geometric_distance == 0.0 {
                    1.0
                } else {
                    effective_distance / layer.geometric_distance
                };
                result.layers.push(layer);
                far_boundaries = far.point.event.boundaries;
            }

            cumulative_path += layer_path;
        }

        result.ground_is_hit = !result.layers.is_empty()
            && radial_result.ground_is_hit
            && far_boundaries.contains_surface();
        result.layers.reverse();
    }

    fn structured_point(&self, point: TracePoint, distance: f64) -> TracePoint {
        let mut boundaries = point.event.boundaries;
        let angle = self.grid.horizontal_angle_at(point.position);
        for (index, boundary) in self.grid.horizontal_angles().iter().enumerate() {
            if (angle - boundary).abs() <= 1e-8 {
                boundaries.push(BoundaryTag::Horizontal { index });
            }
        }
        let event = TraceEvent {
            kind: point.event.kind,
            boundaries,
        };

        TracePoint {
            distance,
            position: point.position,
            altitude: point.altitude,
            event,
            interpolation: self.grid.interpolation_weights_at(point.position),
            cell: self.grid.locate_cell(point.position),
            on_exact_vertical_boundary: point.on_exact_vertical_boundary,
        }
    }

    fn trace_point(&self, ray: Ray, event: Event2D) -> TracePoint {
        let position = ray.point_at(event.distance);
        let trace_event = TraceEvent {
            kind: if event.is_tangent {
                TraceEventKind::Tangent
            } else if event.is_observer {
                TraceEventKind::Observer
            } else {
                TraceEventKind::Boundary
            },
            boundaries: event.boundaries,
        };

        TracePoint {
            distance: event.distance,
            position,
            altitude: self.grid.altitude_at(position),
            event: trace_event,
            interpolation: self.grid.interpolation_weights_at(position),
            cell: self.grid.locate_cell(position),
            on_exact_vertical_boundary: event.boundaries.contains_vertical(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct RefractedIntegralAtRadius {
    radius: f64,
    integral: PathIntegral,
}

fn refracted_integral_at_angle(
    profile: &RefractiveProfile,
    tangent_radius: f64,
    tangent_index: f64,
    start_radius: f64,
    end_radius: f64,
    target_angle: f64,
) -> RefractedIntegralAtRadius {
    const MAX_ITERATIONS: usize = 40;
    const ANGLE_TOLERANCE: f64 = 1e-12;
    const RADIUS_TOLERANCE: f64 = 1e-7;

    let mut lower_fraction = 0.0;
    let mut upper_fraction = 1.0;
    let mut best = RefractedIntegralAtRadius {
        radius: start_radius,
        integral: PathIntegral {
            path_length: 0.0,
            deflection_angle: 0.0,
        },
    };

    for _ in 0..MAX_ITERATIONS {
        let fraction = (lower_fraction + upper_fraction) / 2.0;
        let radius = start_radius + fraction * (end_radius - start_radius);
        let integral = profile.integrate_path_with_tangent_index(
            tangent_radius,
            tangent_index,
            start_radius,
            radius,
        );
        best = RefractedIntegralAtRadius { radius, integral };

        if (integral.deflection_angle - target_angle).abs() <= ANGLE_TOLERANCE
            || (upper_fraction - lower_fraction) * (end_radius - start_radius).abs()
                <= RADIUS_TOLERANCE
        {
            break;
        }
        if integral.deflection_angle < target_angle {
            lower_fraction = fraction;
        } else {
            upper_fraction = fraction;
        }
    }
    best
}

fn angular_crossing_direction(
    ray_plane_normal: Vec3,
    boundary_normal: Vec3,
    selected_radial_direction: Vec3,
    epsilon: f64,
) -> Option<Vec3> {
    let mut direction = ray_plane_normal.cross(boundary_normal);
    if direction.norm() <= epsilon {
        return None;
    }
    direction = direction.normalized();
    let selected_half = direction.dot(selected_radial_direction);
    if selected_half.abs() <= epsilon {
        return None;
    }
    if selected_half < 0.0 {
        direction = -direction;
    }
    Some(direction)
}

fn forward_angle(from: Vec3, to: Vec3, plane_normal: Vec3, epsilon: f64) -> f64 {
    let mut angle = plane_normal.dot(from.cross(to)).atan2(from.dot(to));
    if angle < -epsilon {
        angle += 2.0 * std::f64::consts::PI;
    }
    angle.max(0.0)
}

fn rotate_in_plane(direction: Vec3, plane_normal: Vec3, angle: f64) -> Vec3 {
    direction * angle.cos() + plane_normal.cross(direction) * angle.sin()
}

fn sort_and_merge_events(events: &mut Vec<Event2D>, epsilon: f64) {
    events.sort_unstable_by(|lhs, rhs| lhs.distance.total_cmp(&rhs.distance));

    let mut write = 0;
    for read in 0..events.len() {
        let event = events[read];
        if write > 0 && (events[write - 1].distance - event.distance).abs() <= epsilon {
            for boundary in event.boundaries.iter() {
                events[write - 1].boundaries.push(boundary);
            }
            events[write - 1].is_observer |= event.is_observer;
            events[write - 1].is_tangent |= event.is_tangent;
        } else {
            events[write] = event;
            write += 1;
        }
    }
    events.truncate(write);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raytracer::grid::{BoundaryTag, CellId, InterpolationMethod};
    use crate::raytracer::grid2d::AngularBasis;
    use crate::raytracer::layer::LayerType;
    use crate::raytracer::solar::calculate_csz_saz;
    use crate::raytracer::trace::{TraceOptions, VerticalRayTracer};
    use crate::raytracer::vec3::Vec3;

    const EPSILON: f64 = 1e-9;

    fn grid() -> StructuredGrid2D {
        StructuredGrid2D::new(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![-0.5, 0.0, 0.5],
            InterpolationMethod::Linear,
        )
        .unwrap()
    }

    fn radial_direction(angle: f64) -> Vec3 {
        Vec3::new(angle.sin(), 0.0, angle.cos())
    }

    fn radial_ray(radius: f64, angle: f64, inward: bool) -> Ray {
        let direction = radial_direction(angle);
        Ray::new(
            radius * direction,
            if inward { -direction } else { direction },
        )
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= EPSILON,
            "actual {actual} != expected {expected}"
        );
    }

    fn unity_profile() -> RefractiveProfile {
        RefractiveProfile::new(10.0, vec![0.0, 10.0, 20.0], vec![1.0, 1.0, 1.0]).unwrap()
    }

    fn refractive_profile(scale: f64) -> RefractiveProfile {
        RefractiveProfile::new(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![1.0 + 0.02 * scale, 1.0 + 0.01 * scale, 1.0],
        )
        .unwrap()
    }

    #[test]
    fn outside_radial_ray_traces_to_surface() {
        let traced = StructuredRayTracer2D::new(grid())
            .trace(radial_ray(40.0, 0.25, true), TraceOptions2D::default());

        assert!(traced.ground_is_hit);
        assert_eq!(traced.layers.len(), 2);
        assert_eq!(
            traced.layers[0].cell,
            Some(CellId::Structured2D {
                altitude_index: 0,
                horizontal_index: 1,
            })
        );
        assert_eq!(
            traced.layers[1].cell,
            Some(CellId::Structured2D {
                altitude_index: 1,
                horizontal_index: 1,
            })
        );
        assert_close(traced.layers[0].geometric_distance, 10.0);
        assert_close(traced.layers[1].geometric_distance, 10.0);
    }

    #[test]
    fn out_of_plane_ray_is_traced_on_same_horizontal_coordinate() {
        let angle = 0.25;
        let traced = StructuredRayTracer2D::new(grid()).trace(
            Ray::new(
                40.0 * radial_direction(angle) + 2.0 * Vec3::UNIT_Y,
                -radial_direction(angle),
            ),
            TraceOptions2D::default(),
        );

        assert!(traced.ground_is_hit);
        assert!(!traced.layers.is_empty());
        for layer in &traced.layers {
            assert!(matches!(
                layer.cell,
                Some(CellId::Structured2D {
                    horizontal_index: 1,
                    ..
                })
            ));
        }
    }

    #[test]
    fn outside_observer_looking_up_is_empty() {
        let traced = StructuredRayTracer2D::new(grid())
            .trace(radial_ray(40.0, 0.25, false), TraceOptions2D::default());

        assert!(!traced.ground_is_hit);
        assert!(traced.layers.is_empty());
    }

    #[test]
    fn inside_observer_looking_up_starts_with_partial_layer() {
        let traced = StructuredRayTracer2D::new(grid())
            .trace(radial_ray(25.0, -0.25, false), TraceOptions2D::default());

        assert!(!traced.ground_is_hit);
        assert_eq!(traced.layers.len(), 1);
        assert_close(traced.layers[0].geometric_distance, 5.0);
        assert_eq!(
            traced.layers[0].entrance.event.kind,
            TraceEventKind::Observer
        );
    }

    #[test]
    fn ray_crosses_internal_horizontal_boundary() {
        let origin = 15.0 * radial_direction(-0.25);
        let traced = StructuredRayTracer2D::new(grid())
            .trace(Ray::new(origin, Vec3::UNIT_X), TraceOptions2D::default());

        assert_eq!(traced.layers.len(), 2);
        assert_eq!(
            traced.layers[0].cell,
            Some(CellId::Structured2D {
                altitude_index: 0,
                horizontal_index: 1,
            })
        );
        assert_eq!(
            traced.layers[1].cell,
            Some(CellId::Structured2D {
                altitude_index: 0,
                horizontal_index: 0,
            })
        );
        assert_close(
            traced.layers[0].entrance.position.x,
            traced.layers[1].exit.position.x,
        );
        assert_close(traced.layers[0].entrance.position.x, 0.0);
        assert_eq!(traced.layers[0].layer_type, LayerType::Tangent);
        assert_eq!(traced.layers[1].layer_type, LayerType::Tangent);
    }

    #[test]
    fn ray_can_enter_through_horizontal_domain_boundary() {
        let origin = 15.0 * radial_direction(-0.75);
        let traced = StructuredRayTracer2D::new(grid())
            .trace(Ray::new(origin, Vec3::UNIT_X), TraceOptions2D::default());

        assert!(!traced.layers.is_empty());
        let near_layer = traced.layers.last().unwrap();
        assert!(
            near_layer
                .entrance
                .event
                .boundaries
                .contains(BoundaryTag::Horizontal { index: 0 })
        );
    }

    #[test]
    fn altitude_and_horizontal_corner_hits_are_merged() {
        let boundary_angle = -0.5;
        let corner = 30.0 * radial_direction(boundary_angle);
        let inward_and_across = (-radial_direction(boundary_angle)
            + 0.3 * AngularBasis::CANONICAL.angular_normal(boundary_angle))
        .normalized();
        let ray = Ray::new(corner - 5.0 * inward_and_across, inward_and_across);
        let traced = StructuredRayTracer2D::new(grid()).trace(ray, TraceOptions2D::default());

        let near_layer = traced.layers.last().expect("ray enters the 2D domain");
        assert!(
            near_layer
                .entrance
                .event
                .boundaries
                .contains(BoundaryTag::TopOfAtmosphere { altitude_index: 2 })
        );
        assert!(
            near_layer
                .entrance
                .event
                .boundaries
                .contains(BoundaryTag::Horizontal { index: 0 })
        );
        assert!(near_layer.entrance.event.boundaries.contains_vertical());
        assert!(near_layer.entrance.event.boundaries.contains_horizontal());
        assert_eq!(near_layer.entrance.event.boundaries.len(), 2);
        assert_eq!(near_layer.layer_type, LayerType::Complete);
    }

    #[test]
    fn ray_missing_horizontal_domain_is_empty() {
        let traced = StructuredRayTracer2D::new(grid())
            .trace(radial_ray(40.0, 0.75, true), TraceOptions2D::default());

        assert!(traced.layers.is_empty());
        assert!(!traced.ground_is_hit);
    }

    #[test]
    fn surface_blocks_far_side_atmosphere() {
        let far_side_grid = StructuredGrid2D::new(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![2.8, 3.0, 3.2],
            InterpolationMethod::Linear,
        )
        .unwrap();
        let traced = StructuredRayTracer2D::new(far_side_grid)
            .trace(radial_ray(40.0, 0.0, true), TraceOptions2D::default());

        assert!(traced.layers.is_empty());
        assert!(!traced.ground_is_hit);
    }

    #[test]
    fn tangent_sphere_event_is_preserved() {
        let tracer = StructuredRayTracer2D::new(grid());
        let traced = tracer.trace(
            Ray::new(Vec3::new(-40.0, 0.0, 20.0), Vec3::UNIT_X),
            TraceOptions2D::default(),
        );

        assert!(traced.layers.iter().any(|layer| {
            layer.entrance.event.kind == TraceEventKind::Tangent
                || layer.exit.event.kind == TraceEventKind::Tangent
        }));
    }

    #[test]
    fn tangent_surface_contact_stops_the_ray() {
        let traced = StructuredRayTracer2D::new(grid()).trace(
            Ray::new(Vec3::new(-40.0, 0.0, 10.0), Vec3::UNIT_X),
            TraceOptions2D::default(),
        );

        assert!(traced.ground_is_hit);
        assert!(!traced.layers.is_empty());
        let farthest_layer = &traced.layers[0];
        assert_eq!(farthest_layer.exit.event.kind, TraceEventKind::Tangent);
        assert!(farthest_layer.exit.event.boundaries.contains_surface());
        assert_close(farthest_layer.exit.distance, 40.0);
    }

    #[test]
    fn horizontally_uniform_2d_trace_matches_1d_limb_distances() {
        let altitudes = vec![0.0, 10.0, 20.0];
        let vertical_grid = super::super::grid::VerticalGrid1D::new(
            10.0,
            altitudes.clone(),
            InterpolationMethod::Linear,
            GeometryKind::Spherical,
        )
        .unwrap();
        let structured_grid = StructuredGrid2D::new(
            10.0,
            altitudes,
            vec![-0.1, 2.5],
            InterpolationMethod::Linear,
        )
        .unwrap();
        let tangent_radius = 15.0;
        let observer_radius = 40.0;
        let ray = Ray::new(
            Vec3::new(0.0, 0.0, observer_radius),
            Vec3::new(
                tangent_radius / observer_radius,
                0.0,
                -(1.0 - (tangent_radius / observer_radius).powi(2)).sqrt(),
            ),
        );

        let vertical = VerticalRayTracer::new(vertical_grid).trace(ray, TraceOptions::default());
        let structured =
            StructuredRayTracer2D::new(structured_grid).trace(ray, TraceOptions2D::default());

        assert_eq!(structured.ground_is_hit, vertical.ground_is_hit);
        assert_eq!(structured.layers.len(), vertical.layers.len());
        for (structured, vertical) in structured.layers.iter().zip(&vertical.layers) {
            assert_close(structured.geometric_distance, vertical.geometric_distance);
            let Some(CellId::Structured2D { altitude_index, .. }) = structured.cell else {
                panic!("expected a structured 2D cell");
            };
            assert_eq!(vertical.cell, Some(CellId::AltitudeLayer(altitude_index)));
        }
    }

    #[test]
    fn constant_field_optical_depth_equals_total_path_length() {
        let traced = StructuredRayTracer2D::new(grid())
            .trace(radial_ray(40.0, 0.25, true), TraceOptions2D::default());

        let optical_depth: f64 = traced
            .layers
            .iter()
            .map(|layer| layer.od_quad_start + layer.od_quad_end)
            .sum();
        assert_close(optical_depth, 20.0);
    }

    #[test]
    fn adjacent_layers_share_exact_boundary_positions() {
        let origin = 15.0 * radial_direction(-0.25);
        let traced = StructuredRayTracer2D::new(grid())
            .trace(Ray::new(origin, Vec3::UNIT_X), TraceOptions2D::default());

        for adjacent in traced.layers.windows(2) {
            assert_eq!(adjacent[0].entrance.position, adjacent[1].exit.position);
        }
    }

    #[test]
    fn straight_solar_metadata_matches_and_reuses_endpoint_geometry() {
        let origin = 15.0 * radial_direction(-0.25);
        let traced = StructuredRayTracer2D::new(grid()).trace(
            Ray::new(origin, Vec3::UNIT_X),
            TraceOptions2D {
                solar: Some(SolarContext::new(Vec3::UNIT_Z, GeometryKind::Spherical)),
                ..TraceOptions2D::default()
            },
        );

        for layer in &traced.layers {
            let (entrance_csz, entrance_saz) = calculate_csz_saz(
                Vec3::UNIT_Z,
                layer.entrance.position,
                layer.average_look_away,
                GeometryKind::Spherical,
            );
            let (exit_csz, exit_saz) = calculate_csz_saz(
                Vec3::UNIT_Z,
                layer.exit.position,
                layer.average_look_away,
                GeometryKind::Spherical,
            );
            assert_close(layer.cos_sza_entrance, entrance_csz);
            assert_close(layer.saz_entrance, entrance_saz);
            assert_close(layer.cos_sza_exit, exit_csz);
            assert_close(layer.saz_exit, exit_saz);
        }
        for adjacent in traced.layers.windows(2) {
            assert_eq!(adjacent[0].cos_sza_entrance, adjacent[1].cos_sza_exit);
            assert_eq!(adjacent[0].saz_entrance, adjacent[1].saz_exit);
        }
    }

    #[test]
    fn trace_into_does_not_leak_previous_layers() {
        let tracer = StructuredRayTracer2D::new(grid());
        let mut result = tracer.trace(radial_ray(40.0, 0.25, true), TraceOptions2D::default());
        let mut scratch = TraceScratch2D::with_capacity(tracer.primitives().len());

        tracer.trace_into(
            radial_ray(40.0, 0.75, true),
            &mut result,
            &mut scratch,
            TraceOptions2D::default(),
        );

        assert!(result.layers.is_empty());
        assert!(!result.ground_is_hit);
    }

    #[test]
    fn rotated_basis_preserves_distances_and_cells() {
        let canonical_grid = grid();
        let rotated_basis = AngularBasis::new(Vec3::UNIT_X, -Vec3::UNIT_Z).unwrap();
        let rotated_grid = StructuredGrid2D::new_with_basis(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![-0.5, 0.0, 0.5],
            InterpolationMethod::Linear,
            rotated_basis,
        )
        .unwrap();

        let canonical = StructuredRayTracer2D::new(canonical_grid)
            .trace(radial_ray(40.0, 0.25, true), TraceOptions2D::default());
        let rotated_direction = rotated_basis.radial_direction(0.25);
        let rotated = StructuredRayTracer2D::new(rotated_grid).trace(
            Ray::new(40.0 * rotated_direction, -rotated_direction),
            TraceOptions2D::default(),
        );

        assert_eq!(canonical.layers.len(), rotated.layers.len());
        for (canonical, rotated) in canonical.layers.iter().zip(&rotated.layers) {
            assert_eq!(canonical.cell, rotated.cell);
            assert_close(canonical.geometric_distance, rotated.geometric_distance);
        }
    }

    #[test]
    fn near_radial_ray_ignores_supplied_refraction_consistently() {
        let tracer = StructuredRayTracer2D::new(grid());
        let profile = refractive_profile(1.5);
        let ray = radial_ray(40.0, 0.25, true);

        let straight = tracer.trace(ray, TraceOptions2D::default());
        let with_profile = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert!(with_profile.is_straight);
        assert_eq!(with_profile, straight);
    }

    #[test]
    fn refracted_ray_can_enter_through_horizontal_domain_boundary() {
        let tracer = StructuredRayTracer2D::new(grid());
        let profile = refractive_profile(1.0);
        let ray = Ray::new(15.0 * radial_direction(-0.75), Vec3::UNIT_X);
        let traced = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert!(!traced.is_straight);
        let near_layer = traced.layers.last().expect("ray enters the 2D domain");
        assert!(
            near_layer
                .entrance
                .event
                .boundaries
                .contains(BoundaryTag::Horizontal { index: 0 })
        );
    }

    #[test]
    fn unity_profile_altitude_and_horizontal_corner_is_single_event() {
        let boundary_angle = -0.5;
        let corner = 30.0 * radial_direction(boundary_angle);
        let inward_and_across = (-radial_direction(boundary_angle)
            + 0.3 * AngularBasis::CANONICAL.angular_normal(boundary_angle))
        .normalized();
        let ray = Ray::new(corner - 5.0 * inward_and_across, inward_and_across);
        let profile = unity_profile();
        let traced = StructuredRayTracer2D::new(grid()).trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        let corner_event = traced
            .layers
            .iter()
            .flat_map(|layer| [layer.entrance, layer.exit])
            .find(|point| {
                point
                    .event
                    .boundaries
                    .contains(BoundaryTag::TopOfAtmosphere { altitude_index: 2 })
                    && point
                        .event
                        .boundaries
                        .contains(BoundaryTag::Horizontal { index: 0 })
            })
            .expect("unity-profile trace preserves the corner event");
        assert_eq!(corner_event.event.boundaries.len(), 2);
    }

    #[test]
    fn refracted_inside_observer_looking_up_starts_at_observer() {
        let tracer = StructuredRayTracer2D::new(grid());
        let profile = refractive_profile(1.0);
        let ray = Ray::new(25.0 * radial_direction(0.0), Vec3::new(0.2, 0.0, 1.0));
        let traced = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert!(!traced.is_straight);
        assert!(!traced.ground_is_hit);
        assert!(!traced.layers.is_empty());
        assert!(traced.layers.last().unwrap().entrance.event.is_observer());
    }

    #[test]
    fn refracted_surface_blocks_far_side_horizontal_domain() {
        let far_side_grid = StructuredGrid2D::new(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![2.8, 3.0, 3.2],
            InterpolationMethod::Linear,
        )
        .unwrap();
        let profile = refractive_profile(1.0);
        let traced = StructuredRayTracer2D::new(far_side_grid).trace(
            Ray::new(Vec3::new(0.0, 0.0, 40.0), Vec3::new(0.1, 0.0, -1.0)),
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert!(!traced.is_straight);
        assert!(traced.layers.is_empty());
        assert!(!traced.ground_is_hit);
    }

    #[test]
    fn unity_refraction_matches_straight_2d_geometry() {
        let tracer = StructuredRayTracer2D::new(grid());
        let ray = Ray::new(15.0 * radial_direction(-0.25), Vec3::UNIT_X);
        let profile = unity_profile();

        let straight = tracer.trace(ray, TraceOptions2D::default());
        let refracted = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert!(!refracted.is_straight);
        assert_eq!(refracted.ground_is_hit, straight.ground_is_hit);
        assert_eq!(refracted.tangent_radius, straight.tangent_radius);
        assert_eq!(refracted.layers, straight.layers);
    }

    #[test]
    fn horizontally_uniform_2d_refraction_matches_1d() {
        let altitudes = vec![0.0, 10.0, 20.0];
        let vertical_grid = super::super::grid::VerticalGrid1D::new(
            10.0,
            altitudes.clone(),
            InterpolationMethod::Linear,
            GeometryKind::Spherical,
        )
        .unwrap();
        let structured_grid = StructuredGrid2D::new(
            10.0,
            altitudes,
            vec![-0.1, 2.5],
            InterpolationMethod::Linear,
        )
        .unwrap();
        let profile = refractive_profile(1.0);
        let ray = Ray::new(
            Vec3::new(0.0, 0.0, 40.0),
            Vec3::new(
                15.0 / 40.0,
                0.0,
                -(1.0_f64 - (15.0_f64 / 40.0).powi(2)).sqrt(),
            ),
        );

        let vertical = VerticalRayTracer::new(vertical_grid).trace(
            ray,
            TraceOptions {
                solar: None,
                refraction: Some(&profile),
            },
        );
        let structured = StructuredRayTracer2D::new(structured_grid).trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert_close(structured.tangent_radius, vertical.tangent_radius);
        assert_eq!(structured.ground_is_hit, vertical.ground_is_hit);
        assert_eq!(structured.layers.len(), vertical.layers.len());
        for (structured, vertical) in structured.layers.iter().zip(&vertical.layers) {
            let Some(CellId::Structured2D { altitude_index, .. }) = structured.cell else {
                panic!("expected a structured 2D cell");
            };
            assert_eq!(vertical.cell, Some(CellId::AltitudeLayer(altitude_index)));
            assert_close(
                structured.effective_distance(),
                vertical.effective_distance(),
            );
            assert!((structured.entrance.position - vertical.entrance.position).norm() < 1e-8);
            assert!((structured.exit.position - vertical.exit.position).norm() < 1e-8);
        }
    }

    #[test]
    fn each_ray_can_borrow_a_different_refractive_profile() {
        let tracer = StructuredRayTracer2D::new(grid());
        let ray = Ray::new(Vec3::new(0.0, 0.0, 40.0), Vec3::new(0.375, 0.0, -1.0));
        let weak = refractive_profile(0.25);
        let strong = refractive_profile(1.5);

        let weak_trace = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&weak),
            },
        );
        let strong_trace = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&strong),
            },
        );

        assert!(strong_trace.tangent_radius < weak_trace.tangent_radius);
        let weak_path: f64 = weak_trace
            .layers
            .iter()
            .map(Layer::effective_distance)
            .sum();
        let strong_path: f64 = strong_trace
            .layers
            .iter()
            .map(Layer::effective_distance)
            .sum();
        assert!((strong_path - weak_path).abs() > 1e-6);
    }

    #[test]
    fn reusable_storage_can_alternate_refractive_profiles() {
        let tracer = StructuredRayTracer2D::new(grid());
        let ray = Ray::new(Vec3::new(0.0, 0.0, 40.0), Vec3::new(0.375, 0.0, -1.0));
        let weak = refractive_profile(0.25);
        let strong = refractive_profile(1.5);
        let mut result = tracer.trace(ray, TraceOptions2D::default());
        let mut scratch = TraceScratch2D::with_capacity(tracer.primitives().len());

        tracer.trace_into(
            ray,
            &mut result,
            &mut scratch,
            TraceOptions2D {
                solar: None,
                refraction: Some(&strong),
            },
        );
        let strong_first = result.clone();
        tracer.trace_into(
            ray,
            &mut result,
            &mut scratch,
            TraceOptions2D {
                solar: None,
                refraction: Some(&weak),
            },
        );
        assert_ne!(result.tangent_radius, strong_first.tangent_radius);
        tracer.trace_into(
            ray,
            &mut result,
            &mut scratch,
            TraceOptions2D {
                solar: None,
                refraction: Some(&strong),
            },
        );

        assert_eq!(result, strong_first);
    }

    #[test]
    fn refracted_ray_is_split_at_horizontal_boundaries() {
        let tracer = StructuredRayTracer2D::new(grid());
        let profile = refractive_profile(1.0);
        let ray = Ray::new(15.0 * radial_direction(-0.25), Vec3::UNIT_X);
        let traced = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert!(traced.layers.len() >= 2);
        assert!(traced.layers.iter().any(|layer| {
            layer
                .entrance
                .event
                .boundaries
                .contains(BoundaryTag::Horizontal { index: 1 })
                || layer
                    .exit
                    .event
                    .boundaries
                    .contains(BoundaryTag::Horizontal { index: 1 })
        }));
        assert!(traced.layers.iter().any(|layer| {
            matches!(
                layer.cell,
                Some(CellId::Structured2D {
                    horizontal_index: 0,
                    ..
                })
            )
        }));
        assert!(traced.layers.iter().any(|layer| {
            matches!(
                layer.cell,
                Some(CellId::Structured2D {
                    horizontal_index: 1,
                    ..
                })
            )
        }));
        for adjacent in traced.layers.windows(2) {
            assert_eq!(adjacent[0].entrance.position, adjacent[1].exit.position);
        }
    }

    #[test]
    fn refracted_layers_preserve_path_and_interpolation_invariants() {
        let tracer = StructuredRayTracer2D::new(grid());
        let profile = refractive_profile(1.0);
        let ray = Ray::new(Vec3::new(0.0, 0.0, 40.0), Vec3::new(0.375, 0.0, -1.0));
        let traced = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert!(!traced.is_straight);
        for layer in &traced.layers {
            assert!(layer.geometric_distance.is_finite());
            assert!(layer.geometric_distance > 0.0);
            assert!(layer.curvature_factor.is_finite());
            assert!(layer.curvature_factor > 0.0);
            assert!(layer.effective_distance().is_finite());
            assert!(
                (layer.od_quad_start + layer.od_quad_end - layer.effective_distance()).abs() < 1e-8
            );
            assert!((layer.od_quad_start_fraction + layer.od_quad_end_fraction - 1.0).abs() < 1e-8);
            for point in [layer.entrance, layer.exit] {
                let weight_sum: f64 = point.interpolation.iter().map(|weight| weight.weight).sum();
                assert!((weight_sum - 1.0).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn refracted_ground_path_remains_surface_blocked() {
        let tracer = StructuredRayTracer2D::new(grid());
        let profile = refractive_profile(1.0);
        let ray = Ray::new(Vec3::new(0.0, 0.0, 40.0), Vec3::new(0.1, 0.0, -1.0));
        let traced = tracer.trace(
            ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert!(traced.ground_is_hit);
        assert!(traced.layers[0].exit.event.boundaries.contains_surface());
        assert!(
            traced
                .layers
                .iter()
                .all(|layer| layer.exit.position.z > 0.0)
        );
    }

    #[test]
    fn rotated_basis_preserves_refracted_geometry() {
        let profile = refractive_profile(1.0);
        let canonical_grid = grid();
        let rotated_basis = AngularBasis::new(Vec3::UNIT_X, -Vec3::UNIT_Z).unwrap();
        let rotated_grid = StructuredGrid2D::new_with_basis(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![-0.5, 0.0, 0.5],
            InterpolationMethod::Linear,
            rotated_basis,
        )
        .unwrap();
        let canonical_ray = Ray::new(Vec3::new(0.0, 0.0, 40.0), Vec3::new(0.375, 0.0, -1.0));
        let rotate = |vector: Vec3| {
            vector.x * rotated_basis.reference_x()
                + vector.y * rotated_basis.reference_y()
                + vector.z * rotated_basis.reference_z()
        };
        let rotated_ray = Ray::new(
            rotate(canonical_ray.origin),
            rotate(canonical_ray.direction),
        );

        let canonical = StructuredRayTracer2D::new(canonical_grid).trace(
            canonical_ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );
        let rotated = StructuredRayTracer2D::new(rotated_grid).trace(
            rotated_ray,
            TraceOptions2D {
                solar: None,
                refraction: Some(&profile),
            },
        );

        assert_eq!(canonical.ground_is_hit, rotated.ground_is_hit);
        assert_eq!(canonical.layers.len(), rotated.layers.len());
        assert_close(canonical.tangent_radius, rotated.tangent_radius);
        for (canonical, rotated) in canonical.layers.iter().zip(&rotated.layers) {
            assert_eq!(canonical.cell, rotated.cell);
            assert_close(canonical.effective_distance(), rotated.effective_distance());
        }
    }

    #[test]
    fn refracted_solar_metadata_matches_endpoint_geometry() {
        let tracer = StructuredRayTracer2D::new(grid());
        let profile = refractive_profile(1.0);
        let ray = Ray::new(Vec3::new(0.0, 0.0, 40.0), Vec3::new(0.375, 0.0, -1.0));
        let traced = tracer.trace(
            ray,
            TraceOptions2D {
                solar: Some(SolarContext::new(Vec3::UNIT_Z, GeometryKind::Spherical)),
                refraction: Some(&profile),
            },
        );

        assert!(!traced.layers.is_empty());
        for layer in &traced.layers {
            let (entrance_csz, entrance_saz) = calculate_csz_saz(
                Vec3::UNIT_Z,
                layer.entrance.position,
                layer.average_look_away,
                GeometryKind::Spherical,
            );
            let (exit_csz, exit_saz) = calculate_csz_saz(
                Vec3::UNIT_Z,
                layer.exit.position,
                layer.average_look_away,
                GeometryKind::Spherical,
            );
            assert_close(layer.cos_sza_entrance, entrance_csz);
            assert_close(layer.saz_entrance, entrance_saz);
            assert_close(layer.cos_sza_exit, exit_csz);
            assert_close(layer.saz_exit, exit_saz);
        }
    }

    #[test]
    fn sampled_refracted_rays_preserve_layer_invariants() {
        let tracer = StructuredRayTracer2D::new(grid());
        let profile = refractive_profile(1.0);
        let mut nonempty_traces = 0;

        for observer_angle in [-0.4, -0.1, 0.1, 0.4] {
            for look_x in [-0.8, -0.2, 0.2, 0.8] {
                let ray = Ray::new(
                    35.0 * radial_direction(observer_angle),
                    Vec3::new(look_x, 0.0, -1.0),
                );
                let traced = tracer.trace(
                    ray,
                    TraceOptions2D {
                        solar: None,
                        refraction: Some(&profile),
                    },
                );

                if !traced.layers.is_empty() {
                    nonempty_traces += 1;
                    assert!(!traced.is_straight);
                }

                for layer in &traced.layers {
                    assert!(layer.geometric_distance.is_finite());
                    assert!(layer.geometric_distance > 0.0);
                    assert!(layer.effective_distance().is_finite());
                    assert!(layer.effective_distance() > 0.0);
                    assert!(layer.cell.is_some());
                    assert!(
                        (layer.od_quad_start + layer.od_quad_end - layer.effective_distance())
                            .abs()
                            < 1e-8
                    );
                    let midpoint_direction = (layer.entrance.position.normalized()
                        + layer.exit.position.normalized())
                    .normalized();
                    let midpoint_radius =
                        (layer.entrance.position.norm() + layer.exit.position.norm()) / 2.0;
                    assert_eq!(
                        layer.cell,
                        tracer
                            .grid()
                            .locate_cell(midpoint_radius * midpoint_direction)
                    );
                    for point in [layer.entrance, layer.exit] {
                        let weight_sum: f64 =
                            point.interpolation.iter().map(|weight| weight.weight).sum();
                        assert!((weight_sum - 1.0).abs() < 1e-8);
                    }
                }
                for adjacent in traced.layers.windows(2) {
                    assert_eq!(adjacent[0].entrance.position, adjacent[1].exit.position);
                }
            }
        }
        assert!(nonempty_traces > 0);
    }

    #[test]
    fn sampled_rays_preserve_layer_invariants() {
        let tracer = StructuredRayTracer2D::new(grid());
        for observer_angle in [-0.75, -0.4, -0.1, 0.1, 0.4, 0.75] {
            for look_x in [-0.8, -0.2, 0.2, 0.8] {
                let origin = 35.0 * radial_direction(observer_angle);
                let ray = Ray::new(origin, Vec3::new(look_x, 0.0, -1.0));
                let traced = tracer.trace(ray, TraceOptions2D::default());

                for layer in &traced.layers {
                    assert!(layer.geometric_distance > 0.0);
                    assert!(layer.geometric_distance.is_finite());
                    let midpoint = (layer.entrance.position + layer.exit.position) / 2.0;
                    assert_eq!(layer.cell, tracer.grid().locate_cell(midpoint));
                    for point in [layer.entrance, layer.exit] {
                        let weight_sum: f64 =
                            point.interpolation.iter().map(|weight| weight.weight).sum();
                        assert!((weight_sum - 1.0).abs() < 1e-8);
                    }
                }
            }
        }
    }
}
