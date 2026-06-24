use super::{
    BoundaryTag, CellId, GeometryKind, InterpolationMethod, InterpolationWeight, LayerType, Ray,
    RefractiveProfile, SolarContext, TraceOptions, TracePoint, TracedRay, Vec3, VerticalGrid1D,
    VerticalRayTracer,
};
use crate::raytracer::trace::NADIR_VIEWING_CUTOFF;

const EARTH_RADIUS: f64 = 10.0;
const ALTITUDES: [f64; 4] = [0.0, 10.0, 20.0, 30.0];
const EPSILON: f64 = 1e-8;

fn tracer(geometry: GeometryKind, interpolation: InterpolationMethod) -> VerticalRayTracer {
    let grid = VerticalGrid1D::new(EARTH_RADIUS, ALTITUDES.to_vec(), interpolation, geometry)
        .expect("test grid should be valid");
    VerticalRayTracer::new(grid)
}

fn spherical_tracer() -> VerticalRayTracer {
    tracer(GeometryKind::Spherical, InterpolationMethod::Linear)
}

fn plane_tracer() -> VerticalRayTracer {
    tracer(GeometryKind::PlaneParallel, InterpolationMethod::Linear)
}

fn radial_ray(observer_altitude: f64, up: bool, geometry: GeometryKind) -> Ray {
    let z = EARTH_RADIUS + observer_altitude;
    let direction = if up { Vec3::UNIT_Z } else { -Vec3::UNIT_Z };
    match geometry {
        GeometryKind::Spherical => Ray::new(Vec3::new(0.0, 0.0, z), direction),
        GeometryKind::PlaneParallel => Ray::new(Vec3::new(4.0, -2.0, z), direction),
    }
}

fn spherical_limb_ray(observer_altitude: f64, tangent_altitude: f64) -> Ray {
    let observer_radius = EARTH_RADIUS + observer_altitude;
    let tangent_radius = EARTH_RADIUS + tangent_altitude;
    let ratio = tangent_radius / observer_radius;
    let dz = -(1.0 - ratio * ratio).sqrt();
    Ray::new(
        Vec3::new(0.0, 0.0, observer_radius),
        Vec3::new(ratio, 0.0, dz),
    )
}

fn spherical_downward_ray(observer_altitude: f64, cos_viewing_zenith: f64) -> Ray {
    let radius = EARTH_RADIUS + observer_altitude;
    let sin_viewing_zenith = (1.0 - cos_viewing_zenith * cos_viewing_zenith).sqrt();
    Ray::new(
        Vec3::new(0.0, 0.0, radius),
        Vec3::new(sin_viewing_zenith, 0.0, -cos_viewing_zenith),
    )
}

fn unity_profile() -> RefractiveProfile {
    RefractiveProfile::new(EARTH_RADIUS, ALTITUDES.to_vec(), vec![1.0; ALTITUDES.len()])
        .expect("unity profile should be valid")
}

fn atmosphere_like_profile() -> RefractiveProfile {
    RefractiveProfile::new(
        EARTH_RADIUS,
        ALTITUDES.to_vec(),
        vec![1.00030, 1.00018, 1.00008, 1.00001],
    )
    .expect("atmosphere-like profile should be valid")
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < EPSILON,
        "actual {actual} != expected {expected}"
    );
}

fn assert_layer_altitudes(traced: &TracedRay, expected: &[(f64, f64)]) {
    assert_eq!(
        traced.layers.len(),
        expected.len(),
        "unexpected layer count"
    );
    for (layer, (entrance_altitude, exit_altitude)) in traced.layers.iter().zip(expected) {
        assert_close(layer.entrance.altitude, *entrance_altitude);
        assert_close(layer.exit.altitude, *exit_altitude);
    }
}

fn assert_layer_types(traced: &TracedRay, expected: &[LayerType]) {
    assert_eq!(
        traced.layers.len(),
        expected.len(),
        "unexpected layer count"
    );
    for (layer, layer_type) in traced.layers.iter().zip(expected) {
        assert_eq!(layer.layer_type, *layer_type);
    }
}

fn stencil_weights(point: TracePoint) -> Vec<InterpolationWeight> {
    point.interpolation.iter().collect()
}

fn assert_exact_boundary(point: TracePoint, expected_altitude_index: usize) {
    assert!(point.on_exact_vertical_boundary);
    let weights = stencil_weights(point);
    assert_eq!(weights.len(), 1);
    assert_eq!(weights[0].index, expected_altitude_index);
    assert_close(weights[0].weight, 1.0);
}

fn assert_linear_stencil(point: TracePoint, expected: &[(usize, f64)]) {
    let weights = stencil_weights(point);
    assert_eq!(weights.len(), expected.len());
    for (weight, (index, expected_weight)) in weights.iter().zip(expected) {
        assert_eq!(weight.index, *index);
        assert_close(weight.weight, *expected_weight);
    }
}

#[test]
fn spherical_outside_ground_viewing_matches_cpp_layer_order() {
    let traced = spherical_tracer().trace(
        radial_ray(40.0, false, GeometryKind::Spherical),
        TraceOptions::default(),
    );

    assert!(traced.ground_is_hit);
    assert_close(traced.tangent_radius, 0.0);
    assert_layer_altitudes(&traced, &[(10.0, 0.0), (20.0, 10.0), (30.0, 20.0)]);
    assert_layer_types(
        &traced,
        &[
            LayerType::Complete,
            LayerType::Complete,
            LayerType::Complete,
        ],
    );
}

#[test]
fn spherical_outside_limb_viewing_matches_cpp_layer_order() {
    let traced = spherical_tracer().trace(spherical_limb_ray(40.0, 15.0), TraceOptions::default());

    assert!(!traced.ground_is_hit);
    assert_close(traced.tangent_radius, 25.0);
    assert_layer_altitudes(
        &traced,
        &[(20.0, 30.0), (15.0, 20.0), (20.0, 15.0), (30.0, 20.0)],
    );
    assert_layer_types(
        &traced,
        &[
            LayerType::Complete,
            LayerType::Tangent,
            LayerType::Tangent,
            LayerType::Complete,
        ],
    );
}

#[test]
fn spherical_outside_limb_tangent_on_boundary_is_not_duplicated() {
    let traced = spherical_tracer().trace(spherical_limb_ray(40.0, 10.0), TraceOptions::default());

    assert!(!traced.ground_is_hit);
    assert_layer_altitudes(
        &traced,
        &[(20.0, 30.0), (10.0, 20.0), (20.0, 10.0), (30.0, 20.0)],
    );
    assert_layer_types(
        &traced,
        &[
            LayerType::Complete,
            LayerType::Tangent,
            LayerType::Tangent,
            LayerType::Complete,
        ],
    );
    assert!(
        traced.layers[1]
            .entrance
            .event
            .boundaries
            .contains(BoundaryTag::Altitude { index: 1 })
    );
}

#[test]
fn spherical_inside_looking_up_matches_cpp_layer_order() {
    let traced = spherical_tracer().trace(
        radial_ray(15.0, true, GeometryKind::Spherical),
        TraceOptions::default(),
    );

    assert!(!traced.ground_is_hit);
    assert_layer_altitudes(&traced, &[(20.0, 30.0), (15.0, 20.0)]);
    assert_layer_types(&traced, &[LayerType::Complete, LayerType::Partial]);
}

#[test]
fn spherical_inside_ground_viewing_matches_cpp_layer_order() {
    let traced = spherical_tracer().trace(
        radial_ray(15.0, false, GeometryKind::Spherical),
        TraceOptions::default(),
    );

    assert!(traced.ground_is_hit);
    assert_layer_altitudes(&traced, &[(10.0, 0.0), (15.0, 10.0)]);
    assert_layer_types(&traced, &[LayerType::Complete, LayerType::Partial]);
}

#[test]
fn spherical_inside_limb_viewing_matches_cpp_layer_order() {
    let traced = spherical_tracer().trace(spherical_limb_ray(25.0, 5.0), TraceOptions::default());

    assert!(!traced.ground_is_hit);
    assert_close(traced.tangent_radius, 15.0);
    assert_layer_altitudes(
        &traced,
        &[
            (20.0, 30.0),
            (10.0, 20.0),
            (5.0, 10.0),
            (10.0, 5.0),
            (20.0, 10.0),
            (25.0, 20.0),
        ],
    );
    assert_layer_types(
        &traced,
        &[
            LayerType::Complete,
            LayerType::Complete,
            LayerType::Tangent,
            LayerType::Tangent,
            LayerType::Complete,
            LayerType::Partial,
        ],
    );
}

#[test]
fn spherical_observer_outside_looking_up_is_empty() {
    let traced = spherical_tracer().trace(
        radial_ray(40.0, true, GeometryKind::Spherical),
        TraceOptions::default(),
    );

    assert!(!traced.ground_is_hit);
    assert!(traced.layers.is_empty());
}

#[test]
fn plane_parallel_outside_ground_viewing_matches_cpp_layer_order() {
    let traced = plane_tracer().trace(
        radial_ray(40.0, false, GeometryKind::PlaneParallel),
        TraceOptions::default(),
    );

    assert!(traced.ground_is_hit);
    assert!(traced.tangent_radius.is_nan());
    assert_layer_altitudes(&traced, &[(10.0, 0.0), (20.0, 10.0), (30.0, 20.0)]);
    assert_layer_types(
        &traced,
        &[
            LayerType::Complete,
            LayerType::Complete,
            LayerType::Complete,
        ],
    );
}

#[test]
fn plane_parallel_inside_looking_up_matches_cpp_layer_order() {
    let traced = plane_tracer().trace(
        radial_ray(15.0, true, GeometryKind::PlaneParallel),
        TraceOptions::default(),
    );

    assert!(!traced.ground_is_hit);
    assert_layer_altitudes(&traced, &[(20.0, 30.0), (15.0, 20.0)]);
    assert_layer_types(&traced, &[LayerType::Complete, LayerType::Partial]);
}

#[test]
fn plane_parallel_inside_ground_viewing_matches_cpp_layer_order() {
    let traced = plane_tracer().trace(
        radial_ray(15.0, false, GeometryKind::PlaneParallel),
        TraceOptions::default(),
    );

    assert!(traced.ground_is_hit);
    assert_layer_altitudes(&traced, &[(10.0, 0.0), (15.0, 10.0)]);
    assert_layer_types(&traced, &[LayerType::Complete, LayerType::Partial]);
}

#[test]
fn plane_parallel_observer_outside_looking_up_is_empty() {
    let traced = plane_tracer().trace(
        radial_ray(40.0, true, GeometryKind::PlaneParallel),
        TraceOptions::default(),
    );

    assert!(!traced.ground_is_hit);
    assert!(traced.layers.is_empty());
}

#[test]
fn partial_layer_interpolation_metadata_matches_cpp_linear_locations() {
    let traced = spherical_tracer().trace(
        radial_ray(15.0, true, GeometryKind::Spherical),
        TraceOptions::default(),
    );

    let partial = traced
        .layers
        .last()
        .expect("inside-up ray has partial layer");
    assert_eq!(partial.cell, Some(CellId::AltitudeLayer(1)));
    assert!(!partial.entrance.on_exact_vertical_boundary);
    assert_linear_stencil(partial.entrance, &[(1, 0.5), (2, 0.5)]);
    assert_exact_boundary(partial.exit, 2);
}

#[test]
fn complete_layer_endpoint_metadata_matches_cpp_exact_locations() {
    let traced = spherical_tracer().trace(
        radial_ray(40.0, false, GeometryKind::Spherical),
        TraceOptions::default(),
    );

    let ground_layer = traced.layers.first().expect("ground layer exists");
    assert_eq!(ground_layer.cell, Some(CellId::AltitudeLayer(0)));
    assert_exact_boundary(ground_layer.entrance, 1);
    assert_exact_boundary(ground_layer.exit, 0);
    assert!(ground_layer.exit.event.boundaries.contains_surface());
}

#[test]
fn shell_interpolation_uses_equal_od_quadrature_halves() {
    let traced = tracer(GeometryKind::Spherical, InterpolationMethod::Shell)
        .trace(spherical_limb_ray(40.0, 15.0), TraceOptions::default());

    for layer in &traced.layers {
        assert_close(layer.od_quad_start, layer.effective_distance() / 2.0);
        assert_close(layer.od_quad_end, layer.effective_distance() / 2.0);
        assert_close(layer.od_quad_start_fraction, 0.5);
        assert_close(layer.od_quad_end_fraction, 0.5);
    }
}

#[test]
fn lower_interpolation_assigns_od_to_lower_altitude_endpoint() {
    let traced = tracer(GeometryKind::Spherical, InterpolationMethod::Lower).trace(
        radial_ray(40.0, false, GeometryKind::Spherical),
        TraceOptions::default(),
    );

    let ground_layer = traced.layers.first().expect("ground layer exists");
    assert_close(ground_layer.od_quad_start, 0.0);
    assert_close(ground_layer.od_quad_end, ground_layer.effective_distance());
    assert_close(ground_layer.od_quad_start_fraction, 0.5);
    assert_close(ground_layer.od_quad_end_fraction, 0.5);
}

#[test]
fn trace_into_reuses_result_and_scratch_without_leaking_previous_layers() {
    let tracer = spherical_tracer();
    let mut result = tracer.trace(
        radial_ray(40.0, false, GeometryKind::Spherical),
        TraceOptions::default(),
    );
    let mut scratch = super::TraceScratch::with_capacity(tracer.primitives().len());

    tracer.trace_into(
        radial_ray(40.0, true, GeometryKind::Spherical),
        &mut result,
        &mut scratch,
        TraceOptions::default(),
    );

    assert!(!result.ground_is_hit);
    assert!(result.layers.is_empty());
}

#[test]
fn solar_options_populate_layer_solar_parameters() {
    let traced = plane_tracer().trace(
        radial_ray(40.0, false, GeometryKind::PlaneParallel),
        TraceOptions {
            solar: Some(SolarContext::new(Vec3::UNIT_Z, GeometryKind::PlaneParallel)),
            ..TraceOptions::default()
        },
    );

    for layer in &traced.layers {
        assert_close(layer.cos_sza_entrance, 1.0);
        assert_close(layer.cos_sza_exit, 1.0);
        assert!(layer.saz_entrance.is_finite());
        assert!(layer.saz_exit.is_finite());
    }
}

#[test]
fn refracted_spherical_outside_limb_uses_lower_tangent_radius() {
    let tracer = spherical_tracer();
    let profile = atmosphere_like_profile();
    let ray = spherical_limb_ray(40.0, 15.0);

    let traced = tracer.trace(
        ray,
        TraceOptions {
            refraction: Some(&profile),
            ..TraceOptions::default()
        },
    );

    assert!(!traced.is_straight);
    assert!(traced.tangent_radius < ray.spherical_tangent_radius());
    assert!(!traced.ground_is_hit);
    assert!(
        traced
            .layers
            .iter()
            .any(|layer| layer.curvature_factor != 1.0)
    );
}

#[test]
fn refracted_spherical_inside_limb_has_finite_curvature() {
    let tracer = spherical_tracer();
    let profile = atmosphere_like_profile();
    let traced = tracer.trace(
        spherical_limb_ray(25.0, 5.0),
        TraceOptions {
            refraction: Some(&profile),
            ..TraceOptions::default()
        },
    );

    assert!(!traced.is_straight);
    assert!(!traced.ground_is_hit);
    assert!(!traced.layers.is_empty());
    for layer in &traced.layers {
        assert!(layer.geometric_distance.is_finite());
        assert!(layer.curvature_factor.is_finite());
        assert!(layer.curvature_factor > 0.0);
    }
}

#[test]
fn refracted_spherical_ground_viewing_remains_ground_hit() {
    let tracer = spherical_tracer();
    let profile = atmosphere_like_profile();
    let ray = spherical_downward_ray(40.0, 0.999);
    let traced = tracer.trace(
        ray,
        TraceOptions {
            refraction: Some(&profile),
            ..TraceOptions::default()
        },
    );

    assert!(!traced.is_straight);
    assert!(traced.ground_is_hit);
    assert_close(
        traced
            .layers
            .first()
            .expect("ground layer exists")
            .exit
            .altitude,
        0.0,
    );
}

#[test]
fn refraction_is_ignored_above_nadir_cutoff() {
    let tracer = spherical_tracer();
    let profile = atmosphere_like_profile();
    let cos_viewing_zenith = (NADIR_VIEWING_CUTOFF + 1.0) / 2.0;
    let ray = spherical_downward_ray(40.0, cos_viewing_zenith);
    let traced = tracer.trace(
        ray,
        TraceOptions {
            refraction: Some(&profile),
            ..TraceOptions::default()
        },
    );

    assert!(traced.is_straight);
    assert_close(traced.tangent_radius, ray.spherical_tangent_radius());
}

#[test]
fn refraction_is_enabled_below_nadir_cutoff() {
    let tracer = spherical_tracer();
    let profile = atmosphere_like_profile();
    let cos_viewing_zenith = NADIR_VIEWING_CUTOFF - 1e-4;
    let ray = spherical_downward_ray(40.0, cos_viewing_zenith);
    let traced = tracer.trace(
        ray,
        TraceOptions {
            refraction: Some(&profile),
            ..TraceOptions::default()
        },
    );

    assert!(!traced.is_straight);
    assert!(traced.tangent_radius < ray.spherical_tangent_radius());
}

#[test]
fn unity_refractive_index_matches_straight_spherical_geometry() {
    let tracer = spherical_tracer();
    let profile = unity_profile();
    let ray = spherical_limb_ray(40.0, 15.0);
    let straight = tracer.trace(ray, TraceOptions::default());
    let unity = tracer.trace(
        ray,
        TraceOptions {
            refraction: Some(&profile),
            ..TraceOptions::default()
        },
    );

    assert!(!unity.is_straight);
    assert_eq!(straight.layers.len(), unity.layers.len());
    for (straight_layer, unity_layer) in straight.layers.iter().zip(&unity.layers) {
        assert_close(
            straight_layer.entrance.altitude,
            unity_layer.entrance.altitude,
        );
        assert_close(straight_layer.exit.altitude, unity_layer.exit.altitude);
        assert!(
            (straight_layer.effective_distance() - unity_layer.effective_distance()).abs() < 1e-2
        );
        assert!(unity_layer.curvature_factor.is_finite());
        assert!(unity_layer.curvature_factor > 0.0);
    }
}
