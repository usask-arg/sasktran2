use std::pin::Pin;

use anyhow::{Result, anyhow};

use super::{
    AtmosphereBatch, AtmosphereJacobians, ExecutionPolicy, Geometry, RadianceBatch, SourceMode,
    TwoStreamSolver, View, Workspace,
};

#[cxx::bridge(namespace = "sasktran2::rust::twostream")]
pub mod ffi {
    extern "Rust" {
        type RustTwoStreamSource;

        fn new_rust_twostream_source(
            layer_thickness: &[f64],
            chapman_factors: &[f64],
            solar_cosine: f64,
            source_mode: i32,
            num_threads: usize,
        ) -> Result<Box<RustTwoStreamSource>>;

        fn set_views(
            source: Pin<&mut RustTwoStreamSource>,
            viewing_cosines: &[f64],
            relative_azimuths: &[f64],
        ) -> Result<()>;

        #[allow(clippy::too_many_arguments)]
        fn solve(
            source: Pin<&mut RustTwoStreamSource>,
            num_wavelengths: usize,
            num_levels: usize,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            num_legendre: usize,
            emission: &[f64],
            surface_albedo: &[f64],
            surface_emission: &[f64],
            solar_irradiance: &[f64],
            calculate_jacobians: bool,
        ) -> Result<()>;

        fn radiance(source: &RustTwoStreamSource) -> &[f64];
        fn extinction_jacobian(source: &RustTwoStreamSource) -> &[f64];
        fn ssa_jacobian(source: &RustTwoStreamSource) -> &[f64];
        fn b1_jacobian(source: &RustTwoStreamSource) -> &[f64];
        fn emission_jacobian(source: &RustTwoStreamSource) -> &[f64];
        fn surface_albedo_jacobian(source: &RustTwoStreamSource) -> &[f64];
        fn surface_emission_jacobian(source: &RustTwoStreamSource) -> &[f64];
    }
}

pub struct RustTwoStreamSource {
    solver: TwoStreamSolver,
    views: Vec<View>,
    workspace: Workspace,
    pool: rayon::ThreadPool,
    radiance: Option<RadianceBatch>,
    jacobians: Option<AtmosphereJacobians>,
}

fn new_rust_twostream_source(
    layer_thickness: &[f64],
    chapman_factors: &[f64],
    solar_cosine: f64,
    source_mode: i32,
    num_threads: usize,
) -> Result<Box<RustTwoStreamSource>> {
    let mode = match source_mode {
        0 => SourceMode::Solar,
        1 => SourceMode::Thermal,
        _ => return Err(anyhow!("invalid Rust two-stream source mode")),
    };
    let geometry = Geometry::new(
        layer_thickness.to_vec(),
        chapman_factors.to_vec(),
        solar_cosine,
    );
    let solver = TwoStreamSolver::new(geometry, mode)?.with_execution_policy(if num_threads == 1 {
        ExecutionPolicy::Serial
    } else {
        ExecutionPolicy::Rayon
    });
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .thread_name(|index| format!("sasktran2-twostream-{index}"))
        .build()?;
    Ok(Box::new(RustTwoStreamSource {
        solver,
        views: Vec::new(),
        workspace: Workspace::new(),
        pool,
        radiance: None,
        jacobians: None,
    }))
}

fn set_views(
    mut source: Pin<&mut RustTwoStreamSource>,
    viewing_cosines: &[f64],
    relative_azimuths: &[f64],
) -> Result<()> {
    if viewing_cosines.len() != relative_azimuths.len() {
        return Err(anyhow!(
            "view cosine and azimuth arrays have different lengths"
        ));
    }
    source.views = viewing_cosines
        .iter()
        .zip(relative_azimuths)
        .map(|(&cosine, &relative_azimuth)| View {
            cosine,
            relative_azimuth,
        })
        .collect();
    source.radiance = None;
    source.jacobians = None;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn solve(
    mut source: Pin<&mut RustTwoStreamSource>,
    num_wavelengths: usize,
    num_levels: usize,
    extinction: &[f64],
    single_scatter_albedo: &[f64],
    legendre_coefficients: &[f64],
    num_legendre: usize,
    emission: &[f64],
    surface_albedo: &[f64],
    surface_emission: &[f64],
    solar_irradiance: &[f64],
    calculate_jacobians: bool,
) -> Result<()> {
    if num_levels != source.solver.geometry().num_layers() + 1 {
        return Err(anyhow!("atmosphere and geometry level counts differ"));
    }
    let level_values = num_levels
        .checked_mul(num_wavelengths)
        .ok_or_else(|| anyhow!("atmosphere is too large"))?;
    if extinction.len() != level_values || single_scatter_albedo.len() != level_values {
        return Err(anyhow!("atmospheric optical arrays have the wrong shape"));
    }
    if num_legendre < 2
        || legendre_coefficients.len()
            != num_legendre
                .checked_mul(level_values)
                .ok_or_else(|| anyhow!("Legendre array is too large"))?
    {
        return Err(anyhow!("Legendre coefficient array has the wrong shape"));
    }
    if surface_albedo.len() != num_wavelengths {
        return Err(anyhow!("surface albedo array has the wrong shape"));
    }

    let mode = source.solver.source_mode();
    if mode == SourceMode::Solar && solar_irradiance.len() != num_wavelengths {
        return Err(anyhow!("solar irradiance array has the wrong shape"));
    }
    if mode == SourceMode::Thermal
        && (emission.len() != level_values || surface_emission.len() != num_wavelengths)
    {
        return Err(anyhow!("thermal emission arrays have the wrong shape"));
    }

    // Eigen stores the C++ level-by-wavelength matrices wavelength-major in
    // raw memory.  Reverse the bottom-up C++ level order while transposing to
    // the Rust `[top-down level, wavelength]` layout once per atmosphere.
    let mut rust_extinction = vec![0.0; level_values];
    let mut rust_ssa = vec![0.0; level_values];
    let mut rust_b1 = vec![0.0; level_values];
    let mut rust_emission = (mode == SourceMode::Thermal).then(|| vec![0.0; level_values]);
    for rust_level in 0..num_levels {
        let cpp_level = num_levels - 1 - rust_level;
        for wave in 0..num_wavelengths {
            let rust_index = rust_level * num_wavelengths + wave;
            let cpp_index = wave * num_levels + cpp_level;
            rust_extinction[rust_index] = extinction[cpp_index];
            rust_ssa[rust_index] = single_scatter_albedo[cpp_index];
            rust_b1[rust_index] = legendre_coefficients
                [wave * num_levels * num_legendre + cpp_level * num_legendre + 1];
            if let Some(rust_emission) = &mut rust_emission {
                rust_emission[rust_index] = emission[cpp_index];
            }
        }
    }
    let atmosphere = AtmosphereBatch {
        num_wavelengths,
        extinction: rust_extinction,
        single_scatter_albedo: rust_ssa,
        first_legendre: rust_b1,
        emission: rust_emission,
        surface_albedo: surface_albedo.to_vec(),
        surface_emission: (mode == SourceMode::Thermal).then(|| surface_emission.to_vec()),
        solar_irradiance: (mode == SourceMode::Solar).then(|| solar_irradiance.to_vec()),
    };

    let this = source.as_mut().get_mut();
    let solver = &this.solver;
    let views = &this.views;
    let workspace = &mut this.workspace;
    let result = this.pool.install(|| {
        if calculate_jacobians {
            solver
                .solve_atmosphere_with_jacobians(&atmosphere, views, workspace)
                .map(|(radiance, jacobians)| (radiance, Some(jacobians)))
        } else {
            solver
                .solve_atmosphere(&atmosphere, views, workspace)
                .map(|radiance| (radiance, None))
        }
    })?;
    this.radiance = Some(result.0);
    this.jacobians = result.1;
    Ok(())
}

fn radiance(source: &RustTwoStreamSource) -> &[f64] {
    source
        .radiance
        .as_ref()
        .map_or(&[], |radiance| radiance.values.as_slice())
}

fn extinction_jacobian(source: &RustTwoStreamSource) -> &[f64] {
    source
        .jacobians
        .as_ref()
        .map_or(&[], |jacobians| jacobians.extinction.as_slice())
}

fn ssa_jacobian(source: &RustTwoStreamSource) -> &[f64] {
    source
        .jacobians
        .as_ref()
        .map_or(&[], |jacobians| jacobians.single_scatter_albedo.as_slice())
}

fn b1_jacobian(source: &RustTwoStreamSource) -> &[f64] {
    source
        .jacobians
        .as_ref()
        .map_or(&[], |jacobians| jacobians.first_legendre.as_slice())
}

fn emission_jacobian(source: &RustTwoStreamSource) -> &[f64] {
    source
        .jacobians
        .as_ref()
        .and_then(|jacobians| jacobians.emission.as_deref())
        .unwrap_or(&[])
}

fn surface_albedo_jacobian(source: &RustTwoStreamSource) -> &[f64] {
    source
        .jacobians
        .as_ref()
        .map_or(&[], |jacobians| jacobians.surface_albedo.as_slice())
}

fn surface_emission_jacobian(source: &RustTwoStreamSource) -> &[f64] {
    source
        .jacobians
        .as_ref()
        .and_then(|jacobians| jacobians.surface_emission.as_deref())
        .unwrap_or(&[])
}
