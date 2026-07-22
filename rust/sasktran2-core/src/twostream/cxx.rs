use std::{pin::Pin, sync::Arc};

use anyhow::{Result, anyhow};

use super::{
    AtmosphereBatch, AtmosphereJacobians, ExecutionPolicy, Geometry, RadianceBatch, SourceMode,
    SphericalGeometry, SphericalRayGeometry, TwoStreamSolver, View, Workspace,
};

#[cxx::bridge(namespace = "sasktran2::rust::twostream")]
pub mod ffi {
    extern "Rust" {
        type RustTwoStreamSource;
        type RustSphericalGeometry;

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
        fn new_rust_spherical_geometry(
            source: &RustTwoStreamSource,
            sza_grid: &[f64],
            chapman_factors: &[f64],
            ray_offsets: &[usize],
            ground_hit: &[u8],
            ground_cos_sza: &[f64],
            segment_layers: &[usize],
            segment_fractions: &[f64],
            segment_cosines: &[f64],
            segment_relative_azimuths: &[f64],
            segment_cos_sza: &[f64],
            od_offsets: &[usize],
            od_indices: &[usize],
            od_weights: &[f64],
        ) -> Result<Box<RustSphericalGeometry>>;

        fn set_spherical_geometry(
            source: Pin<&mut RustTwoStreamSource>,
            geometry: &RustSphericalGeometry,
        );

        #[allow(clippy::too_many_arguments)]
        fn solve(
            source: Pin<&mut RustTwoStreamSource>,
            num_wavelengths: usize,
            num_levels: usize,
            extinction: &[f64],
            single_scatter_albedo: &[f64],
            legendre_coefficients: &[f64],
            num_legendre: usize,
            delta_m_fraction: &[f64],
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
    spherical: Option<Arc<SphericalGeometry>>,
    atmosphere: AtmosphereBatch,
    workspace: Workspace,
    pool: Option<rayon::ThreadPool>,
    radiance: Option<RadianceBatch>,
    jacobians: Option<AtmosphereJacobians>,
}

pub struct RustSphericalGeometry {
    geometry: Arc<SphericalGeometry>,
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
    let use_rayon = num_threads > 1;
    let solver = TwoStreamSolver::new(geometry, mode)?.with_execution_policy(if use_rayon {
        ExecutionPolicy::Rayon
    } else {
        ExecutionPolicy::Serial
    });
    let pool = use_rayon
        .then(|| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .thread_name(|index| format!("sasktran2-twostream-{index}"))
                .build()
        })
        .transpose()?;
    Ok(Box::new(RustTwoStreamSource {
        solver,
        views: Vec::new(),
        spherical: None,
        atmosphere: AtmosphereBatch {
            num_wavelengths: 0,
            extinction: Vec::new(),
            single_scatter_albedo: Vec::new(),
            first_legendre: Vec::new(),
            emission: (mode == SourceMode::Thermal).then(Vec::new),
            surface_albedo: Vec::new(),
            surface_emission: (mode == SourceMode::Thermal).then(Vec::new),
            solar_irradiance: (mode == SourceMode::Solar).then(Vec::new),
        },
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
    source.spherical = None;
    source.radiance = None;
    source.jacobians = None;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn new_rust_spherical_geometry(
    source: &RustTwoStreamSource,
    sza_grid: &[f64],
    chapman_factors: &[f64],
    ray_offsets: &[usize],
    ground_hit: &[u8],
    ground_cos_sza: &[f64],
    segment_layers: &[usize],
    segment_fractions: &[f64],
    segment_cosines: &[f64],
    segment_relative_azimuths: &[f64],
    segment_cos_sza: &[f64],
    od_offsets: &[usize],
    od_indices: &[usize],
    od_weights: &[f64],
) -> Result<Box<RustSphericalGeometry>> {
    let num_levels = source.solver.geometry().num_layers() + 1;
    let num_layers = num_levels - 1;
    let expected_chapman = sza_grid
        .len()
        .checked_mul(num_layers)
        .and_then(|value| value.checked_mul(num_layers))
        .ok_or_else(|| anyhow!("spherical Chapman geometry is too large"))?;
    if chapman_factors.len() != expected_chapman {
        return Err(anyhow!("spherical Chapman factors have the wrong shape"));
    }
    let layer_thickness = source.solver.geometry().layer_thickness.clone();
    let mut columns = Vec::with_capacity(sza_grid.len());
    for (column, &solar_cosine) in sza_grid.iter().enumerate() {
        let start = column * num_layers * num_layers;
        columns.push(Geometry::new(
            layer_thickness.clone(),
            chapman_factors[start..start + num_layers * num_layers].to_vec(),
            solar_cosine,
        ));
    }
    let rays = SphericalRayGeometry::new(
        num_levels,
        ray_offsets.to_vec(),
        ground_hit.iter().map(|&value| value != 0).collect(),
        ground_cos_sza.to_vec(),
        segment_layers.to_vec(),
        segment_fractions.to_vec(),
        segment_cosines.to_vec(),
        segment_relative_azimuths.to_vec(),
        segment_cos_sza.to_vec(),
        od_offsets.to_vec(),
        od_indices.to_vec(),
        od_weights.to_vec(),
    )?;
    let spherical = SphericalGeometry::new(sza_grid.to_vec(), columns, rays)?;
    Ok(Box::new(RustSphericalGeometry {
        geometry: Arc::new(spherical),
    }))
}

fn set_spherical_geometry(
    mut source: Pin<&mut RustTwoStreamSource>,
    geometry: &RustSphericalGeometry,
) {
    let this = source.as_mut().get_mut();
    this.views.clear();
    this.spherical = Some(Arc::clone(&geometry.geometry));
    this.radiance = None;
    this.jacobians = None;
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
    delta_m_fraction: &[f64],
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
    if delta_m_fraction.len() != level_values {
        return Err(anyhow!("delta-M fraction array has the wrong shape"));
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
    // raw memory. Reverse the bottom-up C++ level order while transposing to
    // the Rust `[top-down level, wavelength]` layout. The worker owns these
    // buffers so repeated wavelength blocks reuse their allocations.
    let this = source.as_mut().get_mut();
    let atmosphere = &mut this.atmosphere;
    atmosphere.num_wavelengths = num_wavelengths;
    atmosphere.extinction.resize(level_values, 0.0);
    atmosphere.single_scatter_albedo.resize(level_values, 0.0);
    atmosphere.first_legendre.resize(level_values, 0.0);
    if let Some(emission) = &mut atmosphere.emission {
        emission.resize(level_values, 0.0);
    }
    for rust_level in 0..num_levels {
        let cpp_level = num_levels - 1 - rust_level;
        for wave in 0..num_wavelengths {
            let rust_index = rust_level * num_wavelengths + wave;
            let cpp_index = wave * num_levels + cpp_level;
            atmosphere.extinction[rust_index] = extinction[cpp_index];
            atmosphere.single_scatter_albedo[rust_index] = single_scatter_albedo[cpp_index];
            let delta_m = delta_m_fraction[cpp_index];
            let one_minus_delta_m = 1.0 - delta_m;
            if !delta_m.is_finite() || one_minus_delta_m == 0.0 {
                return Err(anyhow!("delta-M fraction contains an invalid value"));
            }
            let scaled_b1 = legendre_coefficients
                [wave * num_levels * num_legendre + cpp_level * num_legendre + 1];
            // Atmosphere::apply_delta_m_scaling performs the common half
            // transform. Complete the multiple-scatter phase transform used
            // by the discrete-ordinate layer preparation for moment l = 1.
            atmosphere.first_legendre[rust_index] = scaled_b1 - 3.0 * delta_m / one_minus_delta_m;
            if let Some(rust_emission) = &mut atmosphere.emission {
                rust_emission[rust_index] = emission[cpp_index];
            }
        }
    }
    atmosphere.surface_albedo.resize(num_wavelengths, 0.0);
    atmosphere.surface_albedo.copy_from_slice(surface_albedo);
    if let Some(target) = &mut atmosphere.surface_emission {
        target.resize(num_wavelengths, 0.0);
        target.copy_from_slice(surface_emission);
    }
    if let Some(target) = &mut atmosphere.solar_irradiance {
        target.resize(num_wavelengths, 0.0);
        target.copy_from_slice(solar_irradiance);
    }

    let solver = &this.solver;
    let views = &this.views;
    let spherical = this.spherical.as_ref();
    let workspace = &mut this.workspace;
    let mut solve = || match (spherical, calculate_jacobians) {
        (Some(geometry), true) => solver
            .solve_spherical_atmosphere_with_jacobians(atmosphere, geometry, workspace)
            .map(|(radiance, jacobians)| (radiance, Some(jacobians))),
        (Some(geometry), false) => solver
            .solve_spherical_atmosphere(atmosphere, geometry, workspace)
            .map(|radiance| (radiance, None)),
        (None, true) => solver
            .solve_atmosphere_with_jacobians(atmosphere, views, workspace)
            .map(|(radiance, jacobians)| (radiance, Some(jacobians))),
        (None, false) => solver
            .solve_atmosphere(atmosphere, views, workspace)
            .map(|radiance| (radiance, None)),
    };
    let result = if let Some(pool) = &this.pool {
        pool.install(solve)
    } else {
        solve()
    }?;
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
