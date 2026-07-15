use ndarray::s;
use sasktran2_rs::bindings::config::{EmissionSource, ThreadingLib, TwoStreamBackend};
use sasktran2_rs::bindings::prelude::*;

fn make_atmosphere(nlevel: usize, nwavel: usize, thermal: bool) -> Atmosphere {
    let mut atmosphere = Atmosphere::new(nwavel, nlevel, 4, true, thermal, Stokes::Stokes1, None);
    for level in 0..nlevel {
        for wave in 0..nwavel {
            atmosphere.storage.total_extinction[[level, wave]] =
                8.0e-6 + 2.0e-7 * level as f64 + 1.0e-8 * wave as f64;
            atmosphere.storage.ssa[[level, wave]] = 0.72 + 0.001 * level as f64;
            atmosphere.storage.leg_coeff[[0, level, wave]] = 1.0;
            atmosphere.storage.leg_coeff[[1, level, wave]] = 0.32 + 0.002 * level as f64;
            atmosphere.storage.emission_source[[level, wave]] =
                1.5 + 0.03 * level as f64 + 0.001 * wave as f64;
        }
    }
    atmosphere.storage.solar_irradiance.fill(1.1);
    atmosphere.surface.brdf_args.fill(0.23);
    atmosphere.surface.emission.fill(1.4);

    {
        let mapping = atmosphere
            .storage
            .get_derivative_mapping("wf_extinction")
            .unwrap();
        mapping.d_extinction().fill(1.0);
        mapping.d_ssa().fill(0.0);
    }
    {
        let mapping = atmosphere.storage.get_derivative_mapping("wf_ssa").unwrap();
        mapping.d_ssa().fill(1.0);
        mapping.d_extinction().fill(0.0);
    }
    {
        let mapping = atmosphere.storage.get_derivative_mapping("wf_b1").unwrap();
        mapping.d_leg_coeff().slice_mut(s![1, .., ..]).fill(1.0);
        mapping.d_extinction().fill(0.0);
        mapping.d_ssa().fill(0.0);
        mapping.scat_factor().fill(1.0);
    }
    if thermal {
        let mapping = atmosphere
            .storage
            .get_derivative_mapping("wf_emission")
            .unwrap();
        mapping.d_emission().fill(1.0);
        mapping.d_extinction().fill(0.0);
        mapping.d_ssa().fill(0.0);
    }
    atmosphere.storage.finalize_scattering_derivatives();

    {
        let mapping = atmosphere
            .surface
            .get_derivative_mapping("wf_albedo")
            .unwrap();
        mapping.d_brdf().fill(1.0);
    }
    if thermal {
        let mapping = atmosphere
            .surface
            .get_derivative_mapping("wf_surface_emission")
            .unwrap();
        mapping.d_emission().fill(1.0);
    }
    atmosphere
}

fn config(backend: TwoStreamBackend, thermal: bool) -> Config {
    let mut config = Config::new();
    config.with_num_threads(1).unwrap();
    config.with_num_streams(2).unwrap();
    config.with_two_stream_backend(backend).unwrap();
    config
        .with_single_scatter_source(SingleScatterSource::None)
        .unwrap();
    if thermal {
        config
            .with_emission_source(EmissionSource::TwoStream)
            .unwrap();
    } else {
        config
            .with_multiple_scatter_source(MultipleScatterSource::TwoStream)
            .unwrap();
    }
    config
}

fn assert_close(actual: &[f64], expected: &[f64], relative_tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
        let tolerance = relative_tolerance * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tolerance,
            "mismatch at {index}: Rust={actual}, C++={expected}, tolerance={tolerance}"
        );
    }
}

fn geometry_and_views(nlevel: usize) -> (Geometry1D, ViewingGeometry) {
    let geometry = Geometry1D::new(
        0.6,
        0.0,
        6_371_000.0,
        (0..nlevel).map(|level| level as f64 * 5_000.0).collect(),
        InterpolationMethod::Linear,
        GeometryType::PlaneParallel,
    );
    let mut viewing = ViewingGeometry::new();
    viewing.add_ground_viewing_solar(0.6, 0.2, 200_000.0, 0.72);
    viewing.add_ground_viewing_solar(0.6, 1.1, 200_000.0, 0.43);
    (geometry, viewing)
}

#[test]
fn rust_twostream_engine_matches_cpp_radiances_and_jacobians() -> Result<()> {
    let nlevel = 6;
    let nwavel = 17;
    let (geometry, viewing) = geometry_and_views(nlevel);

    for thermal in [false, true] {
        let mut atmosphere = make_atmosphere(nlevel, nwavel, thermal);
        let cpp_config = config(TwoStreamBackend::Cpp, thermal);
        let rust_config = config(TwoStreamBackend::Rust, thermal);
        let cpp_engine = Engine::new(&cpp_config, &geometry, &viewing)?;
        let rust_engine = Engine::new(&rust_config, &geometry, &viewing)?;

        let cpp = cpp_engine.calculate_radiance(&atmosphere)?;
        let rust = rust_engine.calculate_radiance(&atmosphere)?;
        assert_close(
            rust.radiance.as_slice().unwrap(),
            cpp.radiance.as_slice().unwrap(),
            1.0e-10,
        );
        for (name, cpp_derivative) in &cpp.d_radiance {
            // The existing C++ thermal backprop has known inconsistencies
            // with end-to-end finite differences for these fields. They are
            // checked against finite differences below instead.
            if thermal && (name == "wf_emission" || name == "wf_b1" || name == "wf_extinction") {
                continue;
            }
            assert_close(
                rust.d_radiance[name].as_slice().unwrap(),
                cpp_derivative.as_slice().unwrap(),
                3.0e-8,
            );
        }
        for (name, cpp_derivative) in &cpp.d_radiance_surf {
            // The C++ surface-albedo adjoint does not match a perturbation of
            // the Lambertian albedo, while the Rust explicit adjoint does.
            if name == "wf_albedo" || (thermal && name == "wf_surface_emission") {
                continue;
            }
            assert_close(
                rust.d_radiance_surf[name].as_slice().unwrap(),
                cpp_derivative.as_slice().unwrap(),
                3.0e-8,
            );
        }

        // Reusing the Rust engine must refresh the batch cache.
        let original_extinction = atmosphere.storage.total_extinction[[2, 3]];
        atmosphere.storage.total_extinction[[2, 3]] *= 1.02;
        let cpp_updated = cpp_engine.calculate_radiance(&atmosphere)?;
        let rust_updated = rust_engine.calculate_radiance(&atmosphere)?;
        assert_close(
            rust_updated.radiance.as_slice().unwrap(),
            cpp_updated.radiance.as_slice().unwrap(),
            1.0e-10,
        );
        atmosphere.storage.total_extinction[[2, 3]] = original_extinction;

        // Verify the engine-level indexing and derivative mapping against
        // end-to-end finite differences, including the fields where the C++
        // adjoint is not a reliable oracle.
        let level = 2;
        let wave = 3;
        let los = 1;
        macro_rules! check_numeric {
            ($field:expr, $analytic:expr, $name:literal) => {{
                let original = $field;
                let step = 1.0e-5 * original.abs().max(1.0e-3);
                $field = original + step;
                let plus = rust_engine.calculate_radiance(&atmosphere)?.radiance[[wave, los, 0]];
                $field = original - step;
                let minus = rust_engine.calculate_radiance(&atmosphere)?.radiance[[wave, los, 0]];
                $field = original;
                let numeric = (plus - minus) / (2.0 * step);
                let analytic = $analytic;
                let tolerance = 2.0e-6 * numeric.abs().max(1.0);
                assert!(
                    (analytic - numeric).abs() < tolerance,
                    "{} (thermal={}): analytic={}, numeric={}",
                    $name,
                    thermal,
                    analytic,
                    numeric
                );
            }};
        }
        check_numeric!(
            atmosphere.storage.total_extinction[[level, wave]],
            rust.d_radiance["wf_extinction"][[level, wave, los, 0]],
            "extinction"
        );
        check_numeric!(
            atmosphere.storage.ssa[[level, wave]],
            rust.d_radiance["wf_ssa"][[level, wave, los, 0]],
            "single-scatter albedo"
        );
        check_numeric!(
            atmosphere.storage.leg_coeff[[1, level, wave]],
            rust.d_radiance["wf_b1"][[level, wave, los, 0]],
            "first Legendre coefficient"
        );
        check_numeric!(
            atmosphere.surface.brdf_args[[0, wave]],
            rust.d_radiance_surf["wf_albedo"][[wave, los, 0]],
            "surface albedo"
        );
        if thermal {
            check_numeric!(
                atmosphere.storage.emission_source[[level, wave]],
                rust.d_radiance["wf_emission"][[level, wave, los, 0]],
                "emission"
            );
            check_numeric!(
                atmosphere.surface.emission[wave],
                rust.d_radiance_surf["wf_surface_emission"][[wave, los, 0]],
                "surface emission"
            );
        }
    }
    Ok(())
}

#[test]
fn rust_twostream_multithread_matches_serial() -> Result<()> {
    let nlevel = 8;
    let atmosphere = make_atmosphere(nlevel, 65, true);
    let (geometry, viewing) = geometry_and_views(nlevel);
    let serial_config = config(TwoStreamBackend::Rust, true);
    let mut parallel_config = config(TwoStreamBackend::Rust, true);
    parallel_config.with_num_threads(4)?;
    parallel_config.with_threading_lib(ThreadingLib::Rayon)?;
    let serial =
        Engine::new(&serial_config, &geometry, &viewing)?.calculate_radiance(&atmosphere)?;
    let parallel =
        Engine::new(&parallel_config, &geometry, &viewing)?.calculate_radiance(&atmosphere)?;

    assert_close(
        parallel.radiance.as_slice().unwrap(),
        serial.radiance.as_slice().unwrap(),
        1.0e-12,
    );
    for (name, serial_derivative) in &serial.d_radiance {
        assert_close(
            parallel.d_radiance[name].as_slice().unwrap(),
            serial_derivative.as_slice().unwrap(),
            1.0e-11,
        );
    }
    for (name, serial_derivative) in &serial.d_radiance_surf {
        assert_close(
            parallel.d_radiance_surf[name].as_slice().unwrap(),
            serial_derivative.as_slice().unwrap(),
            1.0e-11,
        );
    }
    Ok(())
}

#[test]
fn rust_solar_and_thermal_twostream_sources_coexist() -> Result<()> {
    let nlevel = 6;
    let atmosphere = make_atmosphere(nlevel, 17, true);
    let (geometry, viewing) = geometry_and_views(nlevel);
    let mut cpp_config = config(TwoStreamBackend::Cpp, true);
    cpp_config.with_multiple_scatter_source(MultipleScatterSource::TwoStream)?;
    let mut rust_config = config(TwoStreamBackend::Rust, true);
    rust_config.with_multiple_scatter_source(MultipleScatterSource::TwoStream)?;
    let cpp = Engine::new(&cpp_config, &geometry, &viewing)?.calculate_radiance(&atmosphere)?;
    let rust = Engine::new(&rust_config, &geometry, &viewing)?.calculate_radiance(&atmosphere)?;
    assert_close(
        rust.radiance.as_slice().unwrap(),
        cpp.radiance.as_slice().unwrap(),
        1.0e-10,
    );
    Ok(())
}

#[test]
#[ignore = "release-mode performance benchmark"]
fn benchmark_twostream_engine_backends() -> Result<()> {
    use std::hint::black_box;
    use std::time::Instant;

    let nlevel = 21;
    let atmosphere = make_atmosphere(nlevel, 8192, false);
    let (geometry, viewing) = geometry_and_views(nlevel);
    for threads in [1, 4] {
        for backend in [TwoStreamBackend::Cpp, TwoStreamBackend::Rust] {
            let mut benchmark_config = config(backend, false);
            benchmark_config.with_num_threads(threads)?;
            benchmark_config.with_threading_lib(ThreadingLib::Rayon)?;
            let engine = Engine::new(&benchmark_config, &geometry, &viewing)?;
            black_box(engine.calculate_radiance(&atmosphere)?);
            let mut samples = Vec::with_capacity(5);
            for _ in 0..5 {
                let start = Instant::now();
                let output = engine.calculate_radiance(&atmosphere)?;
                black_box(output.radiance[[0, 0, 0]]);
                samples.push(start.elapsed().as_secs_f64() * 1.0e3);
            }
            samples.sort_by(f64::total_cmp);
            eprintln!(
                "backend={backend:?}, threads={threads}, median_ms={:.3}",
                samples[samples.len() / 2]
            );
        }
    }
    Ok(())
}
