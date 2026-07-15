use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rayon::ThreadPoolBuilder;
use sasktran2_core::twostream::{
    AtmosphereBatch, ExecutionPolicy, Geometry, SourceMode, TwoStreamSolver, View, Workspace,
};

fn fixture(layers: usize, wavelengths: usize) -> (Geometry, AtmosphereBatch) {
    let mut chapman = vec![0.0; layers * layers];
    for boundary in 0..layers {
        for layer in 0..=boundary {
            chapman[boundary * layers + layer] = 1.0 / 0.6;
        }
    }
    let levels = (layers + 1) * wavelengths;
    (
        Geometry::new(vec![1.0; layers], chapman, 0.6),
        AtmosphereBatch {
            num_wavelengths: wavelengths,
            extinction: vec![0.01; levels],
            single_scatter_albedo: vec![0.8; levels],
            first_legendre: vec![0.5; levels],
            emission: None,
            surface_albedo: vec![0.3; wavelengths],
            surface_emission: None,
            solar_irradiance: Some(vec![1.0; wavelengths]),
        },
    )
}

fn benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("twostream_solar_forward");
    for layers in [3, 20, 64] {
        for wavelengths in [1, 4, 5, 64, 1024] {
            let (geometry, atmosphere) = fixture(layers, wavelengths);
            let solver = TwoStreamSolver::new(geometry, SourceMode::Solar)
                .unwrap()
                .with_execution_policy(ExecutionPolicy::Serial);
            let inputs = solver.prepare(&atmosphere).unwrap();
            let views = [View {
                cosine: 0.6,
                relative_azimuth: 0.3,
            }];
            let mut workspace = Workspace::new();
            group.throughput(Throughput::Elements((layers * wavelengths) as u64));
            group.bench_with_input(
                BenchmarkId::new(format!("{layers}_layers"), wavelengths),
                &wavelengths,
                |b, _| b.iter(|| solver.solve(&inputs, &views, &mut workspace).unwrap()),
            );
        }
    }
    group.finish();

    let mut group = c.benchmark_group("twostream_solar_fused_atmosphere_vjp");
    for layers in [3, 20] {
        for wavelengths in [4, 64, 1024] {
            let (geometry, atmosphere) = fixture(layers, wavelengths);
            let solver = TwoStreamSolver::new(geometry, SourceMode::Solar)
                .unwrap()
                .with_execution_policy(ExecutionPolicy::Serial);
            let inputs = solver.prepare(&atmosphere).unwrap();
            let views = [View {
                cosine: 0.6,
                relative_azimuth: 0.3,
            }];
            let cotangent = vec![1.0; wavelengths];
            let mut workspace = Workspace::new();
            group.throughput(Throughput::Elements((layers * wavelengths) as u64));
            group.bench_with_input(
                BenchmarkId::new(format!("{layers}_layers"), wavelengths),
                &wavelengths,
                |b, _| {
                    b.iter(|| {
                        solver
                            .solve_with_atmosphere_vjp(
                                &atmosphere,
                                &inputs,
                                &views,
                                &cotangent,
                                &mut workspace,
                            )
                            .unwrap()
                    })
                },
            );
        }
    }
    group.finish();

    let layers = 20;
    let wavelengths = 8192;
    let (geometry, atmosphere) = fixture(layers, wavelengths);
    let solver = TwoStreamSolver::new(geometry, SourceMode::Solar)
        .unwrap()
        .with_execution_policy(ExecutionPolicy::Rayon);
    let inputs = solver.prepare(&atmosphere).unwrap();
    let views = [View {
        cosine: 0.6,
        relative_azimuth: 0.3,
    }];
    let cotangent = vec![1.0; wavelengths];
    let mut group = c.benchmark_group("twostream_solar_fused_vjp_threads");
    group.sample_size(30);
    group.throughput(Throughput::Elements((layers * wavelengths) as u64));
    for threads in [1, 2, 4, 8, 12] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut workspace = Workspace::new();
        group.bench_with_input(BenchmarkId::from_parameter(threads), &threads, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    solver
                        .solve_with_vjp(&inputs, &views, &cotangent, &mut workspace)
                        .unwrap()
                })
            })
        });
    }
    group.finish();

    let mut group = c.benchmark_group("twostream_solar_fused_atmosphere_vjp_threads");
    group.sample_size(30);
    group.throughput(Throughput::Elements((layers * wavelengths) as u64));
    for threads in [1, 2, 4, 8, 12] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut workspace = Workspace::new();
        group.bench_with_input(BenchmarkId::from_parameter(threads), &threads, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    solver
                        .solve_with_atmosphere_vjp(
                            &atmosphere,
                            &inputs,
                            &views,
                            &cotangent,
                            &mut workspace,
                        )
                        .unwrap()
                })
            })
        });
    }
    group.finish();

    let mut group = c.benchmark_group("twostream_solar_atmosphere_vjp_fusion_control");
    group.sample_size(20);
    group.throughput(Throughput::Elements((layers * wavelengths) as u64));
    for threads in [1, 8, 12] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut two_stage_workspace = Workspace::new();
        group.bench_with_input(BenchmarkId::new("two_stage", threads), &threads, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    let (radiance, gradient) = solver
                        .solve_with_vjp(&inputs, &views, &cotangent, &mut two_stage_workspace)
                        .unwrap();
                    let atmosphere_gradient = solver
                        .map_adjoint_to_atmosphere(&atmosphere, &inputs, &gradient)
                        .unwrap();
                    (radiance, atmosphere_gradient)
                })
            })
        });
        let mut fused_workspace = Workspace::new();
        group.bench_with_input(BenchmarkId::new("fused", threads), &threads, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    solver
                        .solve_with_atmosphere_vjp(
                            &atmosphere,
                            &inputs,
                            &views,
                            &cotangent,
                            &mut fused_workspace,
                        )
                        .unwrap()
                })
            })
        });
    }
    group.finish();

    let mut group = c.benchmark_group("twostream_solar_prepare_threads");
    group.sample_size(30);
    group.throughput(Throughput::Elements((layers * wavelengths) as u64));
    for threads in [1, 8, 12] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(threads), &threads, |b, _| {
            b.iter(|| pool.install(|| solver.prepare(&atmosphere).unwrap()))
        });
    }
    group.finish();

    let mut group = c.benchmark_group("twostream_solar_end_to_end_vjp_threads");
    group.sample_size(30);
    group.throughput(Throughput::Elements((layers * wavelengths) as u64));
    for threads in [1, 8, 12] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut workspace = Workspace::new();
        group.bench_with_input(BenchmarkId::from_parameter(threads), &threads, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    let inputs = solver.prepare(&atmosphere).unwrap();
                    solver
                        .solve_with_atmosphere_vjp(
                            &atmosphere,
                            &inputs,
                            &views,
                            &cotangent,
                            &mut workspace,
                        )
                        .unwrap()
                })
            })
        });
    }
    group.finish();

    let mut group = c.benchmark_group("twostream_solar_direct_atmosphere_vjp_threads");
    group.sample_size(30);
    group.throughput(Throughput::Elements((layers * wavelengths) as u64));
    for threads in [1, 8, 12] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut workspace = Workspace::new();
        group.bench_with_input(BenchmarkId::from_parameter(threads), &threads, |b, _| {
            b.iter(|| {
                pool.install(|| {
                    solver
                        .solve_atmosphere_with_vjp(&atmosphere, &views, &cotangent, &mut workspace)
                        .unwrap()
                })
            })
        });
    }
    group.finish();
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
