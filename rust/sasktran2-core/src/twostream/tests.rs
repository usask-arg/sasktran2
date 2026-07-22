use super::*;

fn geometry(n: usize) -> Geometry {
    let mut chapman = vec![0.0; n * n];
    for boundary in 0..n {
        for layer in 0..=boundary {
            chapman[boundary * n + layer] = 2.0;
        }
    }
    Geometry::new(vec![1.0; n], chapman, 0.5)
}

fn atmosphere(n: usize, nw: usize, thermal: bool) -> AtmosphereBatch {
    let levels = (n + 1) * nw;
    AtmosphereBatch {
        num_wavelengths: nw,
        extinction: (0..levels).map(|i| 0.01 + i as f64 * 1.0e-5).collect(),
        single_scatter_albedo: vec![0.8; levels],
        first_legendre: vec![0.4; levels],
        emission: thermal.then(|| {
            (0..=n)
                .flat_map(|level| (0..nw).map(move |_| 2.0 - 0.02 * level as f64))
                .collect()
        }),
        surface_albedo: vec![0.2; nw],
        surface_emission: thermal.then(|| vec![1.5; nw]),
        solar_irradiance: (!thermal).then(|| vec![1.0; nw]),
    }
}

fn spherical_geometry(n: usize) -> SphericalGeometry {
    assert_eq!(n, 3);
    let sza_grid = vec![0.35, 0.75];
    let columns = sza_grid
        .iter()
        .map(|&cosine| {
            let mut column = geometry(n);
            column.solar_cosine = cosine;
            column
        })
        .collect();
    let segment_layers = vec![0, 1, 1, 0, 2, 1, 0, 1];
    let segment_lengths = [0.8, 1.1, 1.1, 0.8, 1.3, 1.0, 0.7, 2.4];
    let mut od_offsets = Vec::with_capacity(segment_layers.len() + 1);
    let mut od_indices = Vec::with_capacity(2 * segment_layers.len());
    let mut od_weights = Vec::with_capacity(2 * segment_layers.len());
    od_offsets.push(0);
    for (&layer, &distance) in segment_layers.iter().zip(&segment_lengths) {
        od_indices.extend([layer, layer + 1]);
        // These are the already path-integrated endpoint weights. The final
        // tangent segment includes a synthetic 1.2 curvature factor.
        od_weights.extend([0.5 * distance, 0.5 * distance]);
        od_offsets.push(od_indices.len());
    }
    let rays = SphericalRayGeometry::new(
        n + 1,
        vec![0, 4, 7, 8],
        vec![false, true, false],
        vec![0.45, 0.62, 0.55],
        segment_layers,
        vec![0.5, 0.4, 0.6, 0.5, 0.8, 0.5, 0.3, 0.5],
        vec![-0.75, -0.18, 0.18, 0.75, 0.82, 0.68, 0.55, 0.0],
        vec![0.2, 0.5, 0.7, 0.9, 0.4, 0.3, 0.2, 1.1],
        vec![0.4, 0.48, 0.57, 0.68, 0.62, 0.55, 0.46, 0.52],
        od_offsets,
        od_indices,
        od_weights,
    )
    .unwrap();
    SphericalGeometry::new(sza_grid, columns, rays).unwrap()
}

#[test]
fn spherical_paths_support_signed_cosines_and_curvature_scaled_od() {
    let spherical = spherical_geometry(3);
    assert!(
        spherical
            .rays
            .segment_cosines
            .iter()
            .any(|value| *value < 0.0)
    );
    assert!(
        spherical
            .rays
            .segment_cosines
            .iter()
            .any(|value| *value > 0.0)
    );
    assert!(spherical.rays.segment_cosines.contains(&0.0));
    let tangent = 7;
    let start = spherical.rays.od_offsets[tangent];
    let end = spherical.rays.od_offsets[tangent + 1];
    assert!((spherical.rays.od_weights[start..end].iter().sum::<f64>() - 2.4).abs() < 1e-14);
}

#[test]
fn spherical_geometry_compiles_straight_and_refracted_traced_rays() {
    use crate::raytracer::{
        GeometryKind, InterpolationMethod, Ray, RefractiveProfile, TraceOptions, Vec3,
        VerticalGrid1D, VerticalRayTracer,
    };

    let earth_radius = 10.0;
    let altitudes = vec![0.0, 10.0, 20.0, 30.0];
    let grid = VerticalGrid1D::new(
        earth_radius,
        altitudes.clone(),
        InterpolationMethod::Linear,
        GeometryKind::Spherical,
    )
    .unwrap();
    let tracer = VerticalRayTracer::new(grid);
    let observer_radius: f64 = 50.0;
    let tangent_radius: f64 = 25.0;
    let ray = Ray::new(
        Vec3::new(0.0, 0.0, observer_radius),
        Vec3::new(
            tangent_radius / observer_radius,
            0.0,
            -(1.0 - (tangent_radius / observer_radius).powi(2)).sqrt(),
        ),
    );
    let straight = tracer.trace(ray, TraceOptions::default());
    let unity = RefractiveProfile::new(earth_radius, altitudes.clone(), vec![1.0; altitudes.len()])
        .unwrap();
    let unity_traced = tracer.trace(
        ray,
        TraceOptions {
            refraction: Some(&unity),
            ..TraceOptions::default()
        },
    );
    let refractive = RefractiveProfile::new(
        earth_radius,
        altitudes.clone(),
        vec![1.00030, 1.00018, 1.00008, 1.00001],
    )
    .unwrap();
    let refracted = tracer.trace(
        ray,
        TraceOptions {
            refraction: Some(&refractive),
            ..TraceOptions::default()
        },
    );

    let straight_compiled =
        SphericalRayGeometry::from_traced_rays(&altitudes, std::slice::from_ref(&straight))
            .unwrap();
    let unity_compiled =
        SphericalRayGeometry::from_traced_rays(&altitudes, std::slice::from_ref(&unity_traced))
            .unwrap();
    assert_eq!(straight_compiled.od_indices, unity_compiled.od_indices);
    for segment in 0..straight.layers.len() {
        let straight_start = straight_compiled.od_offsets[segment];
        let straight_end = straight_compiled.od_offsets[segment + 1];
        let unity_start = unity_compiled.od_offsets[segment];
        let unity_end = unity_compiled.od_offsets[segment + 1];
        let straight_tau: f64 = straight_compiled.od_weights[straight_start..straight_end]
            .iter()
            .sum();
        let unity_tau: f64 = unity_compiled.od_weights[unity_start..unity_end]
            .iter()
            .sum();
        assert!((straight_tau - unity_tau).abs() < 1e-2);
    }

    let refracted_compiled =
        SphericalRayGeometry::from_traced_rays(&altitudes, std::slice::from_ref(&refracted))
            .unwrap();
    assert!(
        refracted
            .layers
            .iter()
            .any(|layer| layer.curvature_factor != 1.0)
    );
    for (segment, layer) in refracted.layers.iter().enumerate() {
        let start = refracted_compiled.od_offsets[segment];
        let end = refracted_compiled.od_offsets[segment + 1];
        let constant_extinction_tau: f64 = refracted_compiled.od_weights[start..end].iter().sum();
        let expected = layer.geometric_distance * layer.curvature_factor;
        assert!((constant_extinction_tau - expected).abs() < 1e-9 * expected.max(1.0));
    }
}

#[test]
fn spherical_explicit_adjoint_matches_central_differences() {
    let n = 3;
    let nw = 3;
    let spherical = spherical_geometry(n);
    let view = 1;
    let level = 2;
    let wave = 1;

    for mode in [SourceMode::Solar, SourceMode::Thermal] {
        let solver = TwoStreamSolver::new(spherical.columns[0].clone(), mode)
            .unwrap()
            .with_execution_policy(ExecutionPolicy::Serial);
        let mut atmosphere = atmosphere(n, nw, mode == SourceMode::Thermal);
        let (forward, jacobians) = solver
            .solve_spherical_atmosphere_with_jacobians(
                &atmosphere,
                &spherical,
                &mut Workspace::new(),
            )
            .unwrap();
        assert!(forward.values.iter().all(|value| value.is_finite()));
        let output_index = view * nw + wave;
        let level_index = (view * (n + 1) + level) * nw + wave;
        let surface_index = view * nw + wave;

        macro_rules! check_field {
            ($field:expr, $analytic:expr, $name:literal) => {{
                let original = $field;
                let step = 1.0e-5 * original.abs().max(1.0e-3);
                $field = original + step;
                let plus = solver
                    .solve_spherical_atmosphere(&atmosphere, &spherical, &mut Workspace::new())
                    .unwrap()
                    .values[output_index];
                $field = original - step;
                let minus = solver
                    .solve_spherical_atmosphere(&atmosphere, &spherical, &mut Workspace::new())
                    .unwrap()
                    .values[output_index];
                $field = original;
                let numeric = (plus - minus) / (2.0 * step);
                let analytic = $analytic;
                let tolerance = 2.0e-6 * numeric.abs().max(1.0);
                assert!(
                    (analytic - numeric).abs() < tolerance,
                    "{mode:?} {} mismatch: analytic={analytic}, numeric={numeric}",
                    $name,
                );
            }};
        }

        check_field!(
            atmosphere.extinction[level * nw + wave],
            jacobians.extinction[level_index],
            "extinction"
        );
        check_field!(
            atmosphere.single_scatter_albedo[level * nw + wave],
            jacobians.single_scatter_albedo[level_index],
            "SSA"
        );
        check_field!(
            atmosphere.first_legendre[level * nw + wave],
            jacobians.first_legendre[level_index],
            "first Legendre"
        );
        check_field!(
            atmosphere.surface_albedo[wave],
            jacobians.surface_albedo[surface_index],
            "surface albedo"
        );
        if mode == SourceMode::Thermal {
            check_field!(
                atmosphere.emission.as_mut().unwrap()[level * nw + wave],
                jacobians.emission.as_ref().unwrap()[level_index],
                "emission"
            );
            check_field!(
                atmosphere.surface_emission.as_mut().unwrap()[wave],
                jacobians.surface_emission.as_ref().unwrap()[surface_index],
                "surface emission"
            );
        }
    }
}

#[test]
fn spherical_single_sza_column_is_supported() {
    let n = 3;
    let nw = 5;
    let multi = spherical_geometry(n);
    let single = SphericalGeometry::new(
        vec![multi.sza_grid[0]],
        vec![multi.columns[0].clone()],
        multi.rays.clone(),
    )
    .unwrap();
    for mode in [SourceMode::Solar, SourceMode::Thermal] {
        let solver = TwoStreamSolver::new(single.columns[0].clone(), mode)
            .unwrap()
            .with_execution_policy(ExecutionPolicy::Serial);
        let result = solver
            .solve_spherical_atmosphere(
                &atmosphere(n, nw, mode == SourceMode::Thermal),
                &single,
                &mut Workspace::new(),
            )
            .unwrap();
        assert!(result.values.iter().all(|value| value.is_finite()));
    }
}

#[test]
fn batch_matches_independent_wavelengths() {
    let n = 4;
    let nw = 5;
    let solver = TwoStreamSolver::new(geometry(n), SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let inputs = solver.prepare(&atmosphere(n, nw, false)).unwrap();
    let views = [
        View {
            cosine: 0.7,
            relative_azimuth: 0.0,
        },
        View {
            cosine: 0.4,
            relative_azimuth: 1.2,
        },
    ];
    let output = solver
        .solve(&inputs, &views, &mut Workspace::new())
        .unwrap();
    assert_eq!(output.values.len(), nw * views.len());
    assert!(output.values.iter().all(|value| value.is_finite()));
    for wave in 1..nw {
        for view in 0..views.len() {
            assert!((output.values[view * nw + wave] - output.values[view * nw]).abs() < 1.0e-3);
        }
    }
}

#[test]
fn fused_vjp_returns_forward_radiance() {
    let n = 4;
    let nw = 5;
    let solver = TwoStreamSolver::new(geometry(n), SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let inputs = solver.prepare(&atmosphere(n, nw, false)).unwrap();
    let views = [
        View {
            cosine: 0.7,
            relative_azimuth: 0.0,
        },
        View {
            cosine: 0.4,
            relative_azimuth: 1.2,
        },
    ];
    let mut workspace = Workspace::new();
    let expected = solver.solve(&inputs, &views, &mut workspace).unwrap();
    let (actual, gradient) = solver
        .solve_with_vjp(
            &inputs,
            &views,
            &vec![1.0; views.len() * nw],
            &mut workspace,
        )
        .unwrap();
    for (actual, expected) in actual.values.iter().zip(expected.values) {
        assert!((actual - expected).abs() < 1.0e-14);
    }
    assert!(gradient.optical_depth.iter().all(|value| value.is_finite()));
}

fn assert_wide_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        let tolerance = 2.0e-11 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() < tolerance,
            "gradient mismatch at {index}: explicit={actual}, tape={expected}"
        );
    }
}

fn assert_atmosphere_close(actual: &AtmosphereAdjoints, expected: &AtmosphereAdjoints) {
    assert_wide_close(&actual.extinction, &expected.extinction);
    assert_wide_close(
        &actual.single_scatter_albedo,
        &expected.single_scatter_albedo,
    );
    assert_wide_close(&actual.first_legendre, &expected.first_legendre);
    assert_wide_close(&actual.surface_albedo, &expected.surface_albedo);
    match (&actual.emission, &expected.emission) {
        (Some(actual), Some(expected)) => assert_wide_close(actual, expected),
        (None, None) => {}
        _ => panic!("emission adjoint shape mismatch"),
    }
    match (&actual.surface_emission, &expected.surface_emission) {
        (Some(actual), Some(expected)) => assert_wide_close(actual, expected),
        (None, None) => {}
        _ => panic!("surface emission adjoint shape mismatch"),
    }
}

#[test]
fn shared_forward_jacobians_match_one_hot_vjps() {
    let n = 4;
    let nw = 17;
    let views = [
        View {
            cosine: 0.72,
            relative_azimuth: 0.1,
        },
        View {
            cosine: 0.43,
            relative_azimuth: 1.1,
        },
    ];

    for mode in [SourceMode::Solar, SourceMode::Thermal] {
        for execution in [ExecutionPolicy::Serial, ExecutionPolicy::Rayon] {
            let atmosphere = atmosphere(n, nw, mode == SourceMode::Thermal);
            let solver = TwoStreamSolver::new(geometry(n), mode)
                .unwrap()
                .with_execution_policy(execution);
            let forward = solver
                .solve_atmosphere(&atmosphere, &views, &mut Workspace::new())
                .unwrap();
            let (radiance, jacobians) = solver
                .solve_atmosphere_with_jacobians(&atmosphere, &views, &mut Workspace::new())
                .unwrap();
            assert_wide_close(&radiance.values, &forward.values);
            assert_eq!(jacobians.num_views, views.len());
            assert_eq!(jacobians.num_levels, n + 1);
            assert_eq!(jacobians.num_wavelengths, nw);

            let level_values = (n + 1) * nw;
            for view in 0..views.len() {
                let mut cotangent = vec![0.0; views.len() * nw];
                cotangent[view * nw..(view + 1) * nw].fill(1.0);
                let (_, expected) = solver
                    .solve_atmosphere_with_vjp(
                        &atmosphere,
                        &views,
                        &cotangent,
                        &mut Workspace::new(),
                    )
                    .unwrap();
                let levels = view * level_values..(view + 1) * level_values;
                let surface = view * nw..(view + 1) * nw;
                assert_wide_close(&jacobians.extinction[levels.clone()], &expected.extinction);
                assert_wide_close(
                    &jacobians.single_scatter_albedo[levels.clone()],
                    &expected.single_scatter_albedo,
                );
                assert_wide_close(
                    &jacobians.first_legendre[levels.clone()],
                    &expected.first_legendre,
                );
                assert_wide_close(
                    &jacobians.surface_albedo[surface.clone()],
                    &expected.surface_albedo,
                );
                if mode == SourceMode::Thermal {
                    assert_wide_close(
                        &jacobians.emission.as_ref().unwrap()[levels],
                        expected.emission.as_ref().unwrap(),
                    );
                    assert_wide_close(
                        &jacobians.surface_emission.as_ref().unwrap()[surface],
                        expected.surface_emission.as_ref().unwrap(),
                    );
                }
            }
        }
    }
}

#[test]
fn atmosphere_jacobians_match_end_to_end_finite_differences() {
    let n = 3;
    let nw = 5;
    let views = [
        View {
            cosine: 0.71,
            relative_azimuth: 0.3,
        },
        View {
            cosine: 0.44,
            relative_azimuth: 1.0,
        },
    ];
    let view = 1;
    let level = 2;
    let wave = 3;

    for mode in [SourceMode::Solar, SourceMode::Thermal] {
        let solver = TwoStreamSolver::new(geometry(n), mode)
            .unwrap()
            .with_execution_policy(ExecutionPolicy::Serial);
        let mut atmosphere = atmosphere(n, nw, mode == SourceMode::Thermal);
        let (_, jacobians) = solver
            .solve_atmosphere_with_jacobians(&atmosphere, &views, &mut Workspace::new())
            .unwrap();
        let jacobian_index = (view * (n + 1) + level) * nw + wave;
        let surface_index = view * nw + wave;
        let output_index = view * nw + wave;

        macro_rules! check_field {
            ($field:expr, $analytic:expr) => {{
                let original = $field;
                let step = 1.0e-5 * original.abs().max(1.0e-3);
                $field = original + step;
                let plus = solver
                    .solve_atmosphere(&atmosphere, &views, &mut Workspace::new())
                    .unwrap()
                    .values[output_index];
                $field = original - step;
                let minus = solver
                    .solve_atmosphere(&atmosphere, &views, &mut Workspace::new())
                    .unwrap()
                    .values[output_index];
                $field = original;
                let numeric = (plus - minus) / (2.0 * step);
                let analytic = $analytic;
                let tolerance = 5.0e-8 * numeric.abs().max(1.0);
                assert!(
                    (analytic - numeric).abs() < tolerance,
                    "{mode:?} Jacobian mismatch: analytic={analytic}, numeric={numeric}"
                );
            }};
        }

        check_field!(
            atmosphere.extinction[level * nw + wave],
            jacobians.extinction[jacobian_index]
        );
        check_field!(
            atmosphere.single_scatter_albedo[level * nw + wave],
            jacobians.single_scatter_albedo[jacobian_index]
        );
        check_field!(
            atmosphere.first_legendre[level * nw + wave],
            jacobians.first_legendre[jacobian_index]
        );
        check_field!(
            atmosphere.surface_albedo[wave],
            jacobians.surface_albedo[surface_index]
        );
        if mode == SourceMode::Thermal {
            check_field!(
                atmosphere.emission.as_mut().unwrap()[level * nw + wave],
                jacobians.emission.as_ref().unwrap()[jacobian_index]
            );
            check_field!(
                atmosphere.surface_emission.as_mut().unwrap()[wave],
                jacobians.surface_emission.as_ref().unwrap()[surface_index]
            );
        }
    }
}

#[test]
fn explicit_adjoint_matches_tape_for_all_fields() {
    let nw = 5;
    let views = [
        View {
            cosine: 0.72,
            relative_azimuth: 0.25,
        },
        View {
            cosine: 0.43,
            relative_azimuth: 1.1,
        },
    ];
    let cotangent: Vec<_> = (0..views.len() * nw)
        .map(|index| 0.25 + 0.03 * index as f64)
        .collect();
    for n in [1, 4] {
        for mode in [SourceMode::Solar, SourceMode::Thermal] {
            let solver = TwoStreamSolver::new(geometry(n), mode)
                .unwrap()
                .with_execution_policy(ExecutionPolicy::Serial);
            let atmosphere = atmosphere(n, nw, mode == SourceMode::Thermal);
            let inputs = solver.prepare(&atmosphere).unwrap();
            let mut workspace = Workspace::new();
            let (radiance, explicit) = solver
                .solve_with_vjp(&inputs, &views, &cotangent, &mut workspace)
                .unwrap();
            let (tape_radiance, tape) = super::reverse::solve_vjp(
                solver.geometry(),
                mode,
                &inputs,
                &views,
                &cotangent,
                &mut workspace.reverse,
            );
            assert_wide_close(&radiance.values, &tape_radiance.values);
            assert_wide_close(&explicit.optical_depth, &tape.optical_depth);
            assert_wide_close(&explicit.single_scatter_albedo, &tape.single_scatter_albedo);
            assert_wide_close(&explicit.first_legendre, &tape.first_legendre);
            assert_wide_close(&explicit.surface_albedo, &tape.surface_albedo);
            match mode {
                SourceMode::Solar => {
                    assert_wide_close(
                        explicit.transmission.as_ref().unwrap(),
                        tape.transmission.as_ref().unwrap(),
                    );
                    assert_wide_close(
                        explicit.average_secant.as_ref().unwrap(),
                        tape.average_secant.as_ref().unwrap(),
                    );
                }
                SourceMode::Thermal => {
                    assert_wide_close(
                        explicit.thermal_b0.as_ref().unwrap(),
                        tape.thermal_b0.as_ref().unwrap(),
                    );
                    assert_wide_close(
                        explicit.thermal_b1.as_ref().unwrap(),
                        tape.thermal_b1.as_ref().unwrap(),
                    );
                    assert_wide_close(
                        explicit.surface_emission.as_ref().unwrap(),
                        tape.surface_emission.as_ref().unwrap(),
                    );
                }
            }

            let mapped = solver
                .map_adjoint_to_atmosphere(&atmosphere, &inputs, &explicit)
                .unwrap();
            let (direct_radiance, direct) = solver
                .solve_with_atmosphere_vjp(&atmosphere, &inputs, &views, &cotangent, &mut workspace)
                .unwrap();
            assert_wide_close(&direct_radiance.values, &radiance.values);
            assert_atmosphere_close(&direct, &mapped);
        }
    }
}

#[test]
fn thermal_forward_and_vjp_are_finite() {
    let n = 3;
    let nw = 2;
    let solver = TwoStreamSolver::new(geometry(n), SourceMode::Thermal)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let atmosphere = atmosphere(n, nw, true);
    let inputs = solver.prepare(&atmosphere).unwrap();
    let views = [View {
        cosine: 0.6,
        relative_azimuth: 0.0,
    }];
    let mut workspace = Workspace::new();
    let output = solver.solve(&inputs, &views, &mut workspace).unwrap();
    assert!(output.values.iter().all(|value| value.is_finite()));
    let gradients = solver
        .vjp_layers(&inputs, &views, &[1.0, 1.0], &mut workspace)
        .unwrap();
    assert!(
        gradients
            .optical_depth
            .iter()
            .all(|value| value.is_finite())
    );
    let mapped = solver
        .map_adjoint_to_atmosphere(&atmosphere, &inputs, &gradients)
        .unwrap();
    assert!(mapped.extinction.iter().all(|value| value.is_finite()));
}

#[test]
fn parallel_matches_serial() {
    let n = 4;
    let nw = 17;
    let geo = geometry(n);
    let atmosphere = atmosphere(n, nw, false);
    let serial = TwoStreamSolver::new(geo.clone(), SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let parallel = TwoStreamSolver::new(geo, SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Rayon);
    let inputs = serial.prepare(&atmosphere).unwrap();
    let parallel_inputs = parallel.prepare(&atmosphere).unwrap();
    assert_wide_close(&parallel_inputs.optical_depth, &inputs.optical_depth);
    assert_wide_close(
        &parallel_inputs.single_scatter_albedo,
        &inputs.single_scatter_albedo,
    );
    assert_wide_close(&parallel_inputs.first_legendre, &inputs.first_legendre);
    assert_wide_close(
        parallel_inputs.transmission.as_ref().unwrap(),
        inputs.transmission.as_ref().unwrap(),
    );
    assert_wide_close(
        parallel_inputs.average_secant.as_ref().unwrap(),
        inputs.average_secant.as_ref().unwrap(),
    );
    let views = [View {
        cosine: 0.55,
        relative_azimuth: 0.4,
    }];
    let a = serial
        .solve(&inputs, &views, &mut Workspace::new())
        .unwrap();
    let b = parallel
        .solve(&inputs, &views, &mut Workspace::new())
        .unwrap();
    assert_eq!(a.values, b.values);

    let cotangent = vec![1.0; nw];
    let (serial_radiance, serial_gradient) = serial
        .solve_with_vjp(&inputs, &views, &cotangent, &mut Workspace::new())
        .unwrap();
    let (parallel_radiance, parallel_gradient) = parallel
        .solve_with_vjp(&inputs, &views, &cotangent, &mut Workspace::new())
        .unwrap();
    assert_wide_close(&parallel_radiance.values, &serial_radiance.values);
    assert_wide_close(
        &parallel_gradient.optical_depth,
        &serial_gradient.optical_depth,
    );
    assert_wide_close(
        &parallel_gradient.single_scatter_albedo,
        &serial_gradient.single_scatter_albedo,
    );
    assert_wide_close(
        &parallel_gradient.first_legendre,
        &serial_gradient.first_legendre,
    );
    assert_wide_close(
        parallel_gradient.transmission.as_ref().unwrap(),
        serial_gradient.transmission.as_ref().unwrap(),
    );
    assert_wide_close(
        parallel_gradient.average_secant.as_ref().unwrap(),
        serial_gradient.average_secant.as_ref().unwrap(),
    );
    assert_wide_close(
        &parallel_gradient.surface_albedo,
        &serial_gradient.surface_albedo,
    );

    let serial_atmosphere_gradient = serial
        .map_adjoint_to_atmosphere(&atmosphere, &inputs, &serial_gradient)
        .unwrap();
    let parallel_atmosphere_gradient = parallel
        .map_adjoint_to_atmosphere(&atmosphere, &inputs, &parallel_gradient)
        .unwrap();
    assert_wide_close(
        &parallel_atmosphere_gradient.extinction,
        &serial_atmosphere_gradient.extinction,
    );
    assert_wide_close(
        &parallel_atmosphere_gradient.single_scatter_albedo,
        &serial_atmosphere_gradient.single_scatter_albedo,
    );
    assert_wide_close(
        &parallel_atmosphere_gradient.first_legendre,
        &serial_atmosphere_gradient.first_legendre,
    );
    assert_wide_close(
        &parallel_atmosphere_gradient.surface_albedo,
        &serial_atmosphere_gradient.surface_albedo,
    );

    let (serial_fused_radiance, serial_fused_atmosphere) = serial
        .solve_with_atmosphere_vjp(
            &atmosphere,
            &inputs,
            &views,
            &cotangent,
            &mut Workspace::new(),
        )
        .unwrap();
    let (parallel_fused_radiance, parallel_fused_atmosphere) = parallel
        .solve_with_atmosphere_vjp(
            &atmosphere,
            &inputs,
            &views,
            &cotangent,
            &mut Workspace::new(),
        )
        .unwrap();
    assert_wide_close(&serial_fused_radiance.values, &serial_radiance.values);
    assert_wide_close(&parallel_fused_radiance.values, &serial_radiance.values);
    assert_atmosphere_close(&serial_fused_atmosphere, &serial_atmosphere_gradient);
    assert_atmosphere_close(&parallel_fused_atmosphere, &serial_atmosphere_gradient);

    let (serial_unprepared_radiance, serial_unprepared_atmosphere) = serial
        .solve_atmosphere_with_vjp(&atmosphere, &views, &cotangent, &mut Workspace::new())
        .unwrap();
    let (parallel_unprepared_radiance, parallel_unprepared_atmosphere) = parallel
        .solve_atmosphere_with_vjp(&atmosphere, &views, &cotangent, &mut Workspace::new())
        .unwrap();
    assert_wide_close(&serial_unprepared_radiance.values, &serial_radiance.values);
    assert_wide_close(
        &parallel_unprepared_radiance.values,
        &serial_radiance.values,
    );
    assert_atmosphere_close(&serial_unprepared_atmosphere, &serial_atmosphere_gradient);
    assert_atmosphere_close(&parallel_unprepared_atmosphere, &serial_atmosphere_gradient);
}

#[test]
fn parallel_thermal_atmosphere_vjp_matches_serial() {
    let n = 4;
    let nw = 17;
    let geo = geometry(n);
    let atmosphere = atmosphere(n, nw, true);
    let serial = TwoStreamSolver::new(geo.clone(), SourceMode::Thermal)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let parallel = TwoStreamSolver::new(geo, SourceMode::Thermal)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Rayon);
    let inputs = serial.prepare(&atmosphere).unwrap();
    let views = [View {
        cosine: 0.55,
        relative_azimuth: 0.4,
    }];
    let cotangent = vec![1.0; nw];
    let (_, serial_layer) = serial
        .solve_with_vjp(&inputs, &views, &cotangent, &mut Workspace::new())
        .unwrap();
    let (_, parallel_layer) = parallel
        .solve_with_vjp(&inputs, &views, &cotangent, &mut Workspace::new())
        .unwrap();
    let serial_atmosphere = serial
        .map_adjoint_to_atmosphere(&atmosphere, &inputs, &serial_layer)
        .unwrap();
    let parallel_atmosphere = parallel
        .map_adjoint_to_atmosphere(&atmosphere, &inputs, &parallel_layer)
        .unwrap();

    assert_wide_close(
        &parallel_atmosphere.extinction,
        &serial_atmosphere.extinction,
    );
    assert_wide_close(
        &parallel_atmosphere.single_scatter_albedo,
        &serial_atmosphere.single_scatter_albedo,
    );
    assert_wide_close(
        &parallel_atmosphere.first_legendre,
        &serial_atmosphere.first_legendre,
    );
    assert_wide_close(
        parallel_atmosphere.emission.as_ref().unwrap(),
        serial_atmosphere.emission.as_ref().unwrap(),
    );
    assert_wide_close(
        &parallel_atmosphere.surface_albedo,
        &serial_atmosphere.surface_albedo,
    );
    assert_wide_close(
        parallel_atmosphere.surface_emission.as_ref().unwrap(),
        serial_atmosphere.surface_emission.as_ref().unwrap(),
    );

    let (serial_fused_radiance, serial_fused_atmosphere) = serial
        .solve_with_atmosphere_vjp(
            &atmosphere,
            &inputs,
            &views,
            &cotangent,
            &mut Workspace::new(),
        )
        .unwrap();
    let (parallel_fused_radiance, parallel_fused_atmosphere) = parallel
        .solve_with_atmosphere_vjp(
            &atmosphere,
            &inputs,
            &views,
            &cotangent,
            &mut Workspace::new(),
        )
        .unwrap();
    assert_wide_close(
        &parallel_fused_radiance.values,
        &serial_fused_radiance.values,
    );
    assert_atmosphere_close(&serial_fused_atmosphere, &serial_atmosphere);
    assert_atmosphere_close(&parallel_fused_atmosphere, &serial_atmosphere);

    let (serial_unprepared_radiance, serial_unprepared_atmosphere) = serial
        .solve_atmosphere_with_vjp(&atmosphere, &views, &cotangent, &mut Workspace::new())
        .unwrap();
    let (parallel_unprepared_radiance, parallel_unprepared_atmosphere) = parallel
        .solve_atmosphere_with_vjp(&atmosphere, &views, &cotangent, &mut Workspace::new())
        .unwrap();
    assert_wide_close(
        &serial_unprepared_radiance.values,
        &serial_fused_radiance.values,
    );
    assert_wide_close(
        &parallel_unprepared_radiance.values,
        &serial_fused_radiance.values,
    );
    assert_atmosphere_close(&serial_unprepared_atmosphere, &serial_atmosphere);
    assert_atmosphere_close(&parallel_unprepared_atmosphere, &serial_atmosphere);
}

#[test]
fn matches_low_level_regression_oracles() {
    let n = 3;
    let geo = Geometry::new(vec![1.0; n], vec![0.0; n * n], 0.5);
    let view = [View {
        cosine: 0.6,
        relative_azimuth: 0.4,
    }];
    let base = LayerInputs {
        num_layers: n,
        num_wavelengths: 1,
        optical_depth: vec![0.01, 0.02, 0.03],
        single_scatter_albedo: vec![0.8, 0.75, 0.7],
        first_legendre: vec![0.4, 0.3, 0.2],
        transmission: Some(vec![
            1.0,
            (-0.02_f64).exp(),
            (-0.06_f64).exp(),
            (-0.12_f64).exp(),
        ]),
        average_secant: Some(vec![2.0; n]),
        thermal_b0: None,
        thermal_b1: None,
        surface_albedo: vec![0.2],
        surface_emission: None,
    };
    let solar = TwoStreamSolver::new(geo.clone(), SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial)
        .solve(&base, &view, &mut Workspace::new())
        .unwrap();
    assert!((solar.values[0] - 0.002_798_163_771_751_667_6).abs() < 1.0e-14);

    let thermal_inputs = LayerInputs {
        transmission: None,
        average_secant: None,
        thermal_b0: Some(vec![2.0, 1.8, 1.6]),
        thermal_b1: Some(vec![0.1, 0.15, 0.2]),
        surface_emission: Some(vec![1.5]),
        ..base
    };
    let thermal = TwoStreamSolver::new(geo, SourceMode::Thermal)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial)
        .solve(&thermal_inputs, &view, &mut Workspace::new())
        .unwrap();
    assert!((thermal.values[0] - 1.479_020_565_824_959_3).abs() < 1.0e-13);
}

#[test]
fn pure_absorption_thermal_matches_formal_solution() {
    let optical_depth = vec![0.15, 0.2, 0.1];
    let thermal_b0 = vec![2.0, 1.8, 1.65];
    let thermal_b1 = vec![0.2, -0.1, 0.05];
    let surface_emission = 1.4;
    let view_cosine = 0.63;
    let inputs = LayerInputs {
        num_layers: optical_depth.len(),
        num_wavelengths: 1,
        optical_depth: optical_depth.clone(),
        single_scatter_albedo: vec![0.0; optical_depth.len()],
        first_legendre: vec![0.0; optical_depth.len()],
        transmission: None,
        average_secant: None,
        thermal_b0: Some(thermal_b0.clone()),
        thermal_b1: Some(thermal_b1.clone()),
        surface_albedo: vec![0.0],
        surface_emission: Some(vec![surface_emission]),
    };
    let view = [View {
        cosine: view_cosine,
        relative_azimuth: 0.0,
    }];
    let solver = TwoStreamSolver::new(
        Geometry::new(vec![1.0; optical_depth.len()], vec![0.0; 9], 0.5),
        SourceMode::Thermal,
    )
    .unwrap()
    .with_execution_policy(ExecutionPolicy::Serial);
    let actual = solver
        .solve(&inputs, &view, &mut Workspace::new())
        .unwrap()
        .values[0];

    let mut attenuation = 1.0;
    let mut expected = 0.0;
    for layer in 0..optical_depth.len() {
        let exponent = thermal_b1[layer] + 1.0 / view_cosine;
        expected +=
            attenuation * thermal_b0[layer] * (1.0 - (-exponent * optical_depth[layer]).exp())
                / (1.0 + thermal_b1[layer] * view_cosine);
        attenuation *= (-optical_depth[layer] / view_cosine).exp();
    }
    expected += attenuation * surface_emission;
    assert!((actual - expected).abs() < 2.0e-13);
}

#[test]
fn transparent_thermal_atmosphere_is_finite_and_transmits_surface() {
    let n = 3;
    let nw = 5;
    let solver = TwoStreamSolver::new(geometry(n), SourceMode::Thermal)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let mut atmosphere = atmosphere(n, nw, true);
    atmosphere.extinction.fill(0.0);
    let views = [
        View {
            cosine: 0.72,
            relative_azimuth: 0.0,
        },
        View {
            cosine: 0.43,
            relative_azimuth: 1.0,
        },
    ];
    let cotangent = vec![1.0; views.len() * nw];
    let (radiance, adjoint) = solver
        .solve_atmosphere_with_vjp(&atmosphere, &views, &cotangent, &mut Workspace::new())
        .unwrap();

    for (view, values) in radiance.values.chunks_exact(nw).enumerate() {
        for (wave, &value) in values.iter().enumerate() {
            assert!(
                (value - atmosphere.surface_emission.as_ref().unwrap()[wave]).abs() < 1.0e-13,
                "transparent thermal radiance mismatch for view {view}, wavelength {wave}"
            );
        }
    }
    assert!(radiance.values.iter().all(|value| value.is_finite()));
    assert!(adjoint.extinction.iter().all(|value| value.is_finite()));
    assert!(
        adjoint
            .single_scatter_albedo
            .iter()
            .all(|value| value.is_finite())
    );
    assert!(
        adjoint
            .emission
            .as_ref()
            .unwrap()
            .iter()
            .all(|value| value.is_finite())
    );
}

#[test]
fn solar_atmosphere_mapping_matches_end_to_end_difference() {
    let n = 2;
    let nw = 1;
    let solver = TwoStreamSolver::new(geometry(n), SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let mut atmosphere = atmosphere(n, nw, false);
    let views = [View {
        cosine: 0.63,
        relative_azimuth: 0.2,
    }];
    let inputs = solver.prepare(&atmosphere).unwrap();
    let mut workspace = Workspace::new();
    let layer = solver
        .vjp_layers(&inputs, &views, &[1.0], &mut workspace)
        .unwrap();
    let mapped = solver
        .map_adjoint_to_atmosphere(&atmosphere, &inputs, &layer)
        .unwrap();

    let index = 1;
    let original = atmosphere.extinction[index];
    let h = 1.0e-7;
    atmosphere.extinction[index] = original + h;
    let plus = solver
        .solve(
            &solver.prepare(&atmosphere).unwrap(),
            &views,
            &mut workspace,
        )
        .unwrap()
        .values[0];
    atmosphere.extinction[index] = original - h;
    let minus = solver
        .solve(
            &solver.prepare(&atmosphere).unwrap(),
            &views,
            &mut workspace,
        )
        .unwrap()
        .values[0];
    let numerical = (plus - minus) / (2.0 * h);
    assert!((mapped.extinction[index] - numerical).abs() < 2.0e-7);
}

#[test]
fn zero_cotangent_produces_zero_gradient() {
    let solver = TwoStreamSolver::new(geometry(2), SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let inputs = solver.prepare(&atmosphere(2, 1, false)).unwrap();
    let views = [View {
        cosine: 0.7,
        relative_azimuth: 0.0,
    }];
    let gradient = solver
        .vjp_layers(&inputs, &views, &[0.0], &mut Workspace::new())
        .unwrap();
    assert!(
        gradient
            .optical_depth
            .iter()
            .chain(&gradient.single_scatter_albedo)
            .chain(&gradient.first_legendre)
            .all(|value| *value == 0.0)
    );
}

fn numerical_layer_derivative(
    solver: &TwoStreamSolver,
    inputs: &LayerInputs,
    views: &[View],
    original: f64,
    mut set: impl FnMut(&mut LayerInputs, f64),
) -> f64 {
    let h = 1.0e-6 * original.abs().max(1.0);
    let mut plus = inputs.clone();
    set(&mut plus, original + h);
    let plus = solver
        .solve(&plus, views, &mut Workspace::new())
        .unwrap()
        .values[0];
    let mut minus = inputs.clone();
    set(&mut minus, original - h);
    let minus = solver
        .solve(&minus, views, &mut Workspace::new())
        .unwrap()
        .values[0];
    (plus - minus) / (2.0 * h)
}

fn assert_derivative(analytic: f64, numerical: f64) {
    let tolerance = 2.0e-7 * analytic.abs().max(numerical.abs()).max(1.0);
    assert!(
        (analytic - numerical).abs() < tolerance,
        "analytic={analytic:e}, numerical={numerical:e}"
    );
}

#[test]
fn solar_reverse_mode_matches_layer_finite_differences() {
    let solver = TwoStreamSolver::new(geometry(3), SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let inputs = solver.prepare(&atmosphere(3, 1, false)).unwrap();
    let views = [View {
        cosine: 0.57,
        relative_azimuth: 0.7,
    }];
    let gradient = solver
        .vjp_layers(&inputs, &views, &[1.0], &mut Workspace::new())
        .unwrap();
    let index = 1;
    assert_derivative(
        gradient.optical_depth[index],
        numerical_layer_derivative(
            &solver,
            &inputs,
            &views,
            inputs.optical_depth[index],
            |i, v| i.optical_depth[index] = v,
        ),
    );
    assert_derivative(
        gradient.single_scatter_albedo[index],
        numerical_layer_derivative(
            &solver,
            &inputs,
            &views,
            inputs.single_scatter_albedo[index],
            |i, v| i.single_scatter_albedo[index] = v,
        ),
    );
    assert_derivative(
        gradient.first_legendre[index],
        numerical_layer_derivative(
            &solver,
            &inputs,
            &views,
            inputs.first_legendre[index],
            |i, v| i.first_legendre[index] = v,
        ),
    );
    assert_derivative(
        gradient.transmission.as_ref().unwrap()[2],
        numerical_layer_derivative(
            &solver,
            &inputs,
            &views,
            inputs.transmission.as_ref().unwrap()[2],
            |i, v| i.transmission.as_mut().unwrap()[2] = v,
        ),
    );
    assert_derivative(
        gradient.average_secant.as_ref().unwrap()[index],
        numerical_layer_derivative(
            &solver,
            &inputs,
            &views,
            inputs.average_secant.as_ref().unwrap()[index],
            |i, v| i.average_secant.as_mut().unwrap()[index] = v,
        ),
    );
    assert_derivative(
        gradient.surface_albedo[0],
        numerical_layer_derivative(
            &solver,
            &inputs,
            &views,
            inputs.surface_albedo[0],
            |i, v| i.surface_albedo[0] = v,
        ),
    );
}

#[test]
fn thermal_reverse_mode_matches_layer_finite_differences() {
    let solver = TwoStreamSolver::new(geometry(3), SourceMode::Thermal)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Serial);
    let inputs = solver.prepare(&atmosphere(3, 1, true)).unwrap();
    let views = [View {
        cosine: 0.61,
        relative_azimuth: 0.0,
    }];
    let gradient = solver
        .vjp_layers(&inputs, &views, &[1.0], &mut Workspace::new())
        .unwrap();
    let index = 1;
    for (analytic, original, field) in [
        (
            gradient.optical_depth[index],
            inputs.optical_depth[index],
            0,
        ),
        (
            gradient.single_scatter_albedo[index],
            inputs.single_scatter_albedo[index],
            1,
        ),
        (
            gradient.first_legendre[index],
            inputs.first_legendre[index],
            2,
        ),
        (
            gradient.thermal_b0.as_ref().unwrap()[index],
            inputs.thermal_b0.as_ref().unwrap()[index],
            3,
        ),
        (
            gradient.thermal_b1.as_ref().unwrap()[index],
            inputs.thermal_b1.as_ref().unwrap()[index],
            4,
        ),
    ] {
        let numerical =
            numerical_layer_derivative(&solver, &inputs, &views, original, |input, value| {
                match field {
                    0 => input.optical_depth[index] = value,
                    1 => input.single_scatter_albedo[index] = value,
                    2 => input.first_legendre[index] = value,
                    3 => input.thermal_b0.as_mut().unwrap()[index] = value,
                    _ => input.thermal_b1.as_mut().unwrap()[index] = value,
                }
            });
        assert_derivative(analytic, numerical);
    }
    assert_derivative(
        gradient.surface_emission.as_ref().unwrap()[0],
        numerical_layer_derivative(
            &solver,
            &inputs,
            &views,
            inputs.surface_emission.as_ref().unwrap()[0],
            |i, v| i.surface_emission.as_mut().unwrap()[0] = v,
        ),
    );
}
