#![allow(clippy::too_many_arguments)]

use super::{
    AngularBasis, CellId, GeometryKind, InterpolationMethod, Layer, LayerType, Ray,
    RefractiveProfile, SolarContext, StructuredGrid2D, StructuredRayTracer2D, TraceOptions,
    TraceOptions2D, TraceScratch, TraceScratch2D, TracedRay, Vec3, VerticalGrid1D,
    VerticalRayTracer,
};
use anyhow::{Result, anyhow};
use cxx::CxxVector;
use std::{cell::RefCell, pin::Pin};

#[cxx::bridge(namespace = "sasktran2::rust::raytracer")]
pub mod ffi {
    struct RustTraceSummary {
        ground_is_hit: bool,
        is_straight: bool,
        tangent_radius: f64,
        num_layers: usize,
    }

    struct RustTraceLayer {
        layer_type: i32,
        entrance_x: f64,
        entrance_y: f64,
        entrance_z: f64,
        exit_x: f64,
        exit_y: f64,
        exit_z: f64,
        entrance_altitude: f64,
        exit_altitude: f64,
        entrance_on_exact_altitude: bool,
        exit_on_exact_altitude: bool,
        entrance_lower_alt_index: i32,
        exit_lower_alt_index: i32,
        layer_distance: f64,
        curvature_factor: f64,
        od_quad_start: f64,
        od_quad_end: f64,
        od_quad_start_fraction: f64,
        od_quad_end_fraction: f64,
        cos_sza_entrance: f64,
        cos_sza_exit: f64,
        saz_entrance: f64,
        saz_exit: f64,
        cell_altitude_index: i32,
        cell_horizontal_index: i32,
    }

    struct RustTimingSummary {
        checksum: f64,
        num_traces: usize,
        total_layers: usize,
    }

    unsafe extern "C++" {
        include!("sasktran2/raytracing_rust_bridge.h");

        type CppTraceResult;
        type CppTraceResult2D;

        fn prepare_trace_result(result: Pin<&mut CppTraceResult>, summary: &RustTraceSummary);
        fn set_trace_layer(result: Pin<&mut CppTraceResult>, index: usize, layer: &RustTraceLayer);
        fn prepare_trace_result_2d(result: Pin<&mut CppTraceResult2D>, summary: &RustTraceSummary);
        fn set_trace_layer_2d(
            result: Pin<&mut CppTraceResult2D>,
            index: usize,
            layer: &RustTraceLayer,
        );
    }

    extern "Rust" {
        type RustVerticalTracer;
        type RustStructuredTracer2D;

        fn new_vertical_tracer(
            earth_radius: f64,
            altitudes: &CxxVector<f64>,
            refractive_index: &CxxVector<f64>,
            geometry_type: i32,
            interpolation_method: i32,
            sun_x: f64,
            sun_y: f64,
            sun_z: f64,
        ) -> Result<Box<RustVerticalTracer>>;

        fn trace_vertical_ray_with_tracer(
            tracer: &RustVerticalTracer,
            origin_x: f64,
            origin_y: f64,
            origin_z: f64,
            look_x: f64,
            look_y: f64,
            look_z: f64,
            include_refraction: bool,
            layers: Pin<&mut CxxVector<RustTraceLayer>>,
        ) -> Result<RustTraceSummary>;

        fn trace_vertical_ray_into_cpp_result(
            tracer: &RustVerticalTracer,
            origin_x: f64,
            origin_y: f64,
            origin_z: f64,
            look_x: f64,
            look_y: f64,
            look_z: f64,
            include_refraction: bool,
            result: Pin<&mut CppTraceResult>,
        ) -> Result<RustTraceSummary>;

        fn trace_vertical_ray(
            earth_radius: f64,
            altitudes: &CxxVector<f64>,
            refractive_index: &CxxVector<f64>,
            geometry_type: i32,
            interpolation_method: i32,
            sun_x: f64,
            sun_y: f64,
            sun_z: f64,
            origin_x: f64,
            origin_y: f64,
            origin_z: f64,
            look_x: f64,
            look_y: f64,
            look_z: f64,
            include_refraction: bool,
            layers: Pin<&mut CxxVector<RustTraceLayer>>,
        ) -> Result<RustTraceSummary>;

        fn trace_vertical_ray_batch_checksum(
            earth_radius: f64,
            altitudes: &CxxVector<f64>,
            refractive_index: &CxxVector<f64>,
            geometry_type: i32,
            interpolation_method: i32,
            sun_x: f64,
            sun_y: f64,
            sun_z: f64,
            origins: &CxxVector<f64>,
            looks: &CxxVector<f64>,
            include_refraction: bool,
            iterations: usize,
        ) -> Result<RustTimingSummary>;

        fn new_structured_tracer_2d(
            earth_radius: f64,
            altitudes: &CxxVector<f64>,
            horizontal_angles: &CxxVector<f64>,
            interpolation_method: i32,
            reference_x_x: f64,
            reference_x_y: f64,
            reference_x_z: f64,
            reference_z_x: f64,
            reference_z_y: f64,
            reference_z_z: f64,
            sun_x: f64,
            sun_y: f64,
            sun_z: f64,
        ) -> Result<Box<RustStructuredTracer2D>>;

        fn trace_structured_ray_2d_into_cpp_result(
            tracer: &RustStructuredTracer2D,
            origin_x: f64,
            origin_y: f64,
            origin_z: f64,
            look_x: f64,
            look_y: f64,
            look_z: f64,
            refractive_index: &[f64],
            result: Pin<&mut CppTraceResult2D>,
        ) -> Result<RustTraceSummary>;
    }
}

pub struct RustVerticalTracer {
    tracer: VerticalRayTracer,
    profile: Option<RefractiveProfile>,
    geometry: GeometryKind,
    solar: SolarContext,
}

pub struct RustStructuredTracer2D {
    tracer: StructuredRayTracer2D,
    solar: SolarContext,
}

struct ReusableTraceStorage {
    result: TracedRay,
    scratch: TraceScratch,
}

struct ReusableTraceStorage2D {
    result: TracedRay,
    scratch: TraceScratch2D,
}

impl ReusableTraceStorage {
    fn new(
        geometry: GeometryKind,
        interpolation: InterpolationMethod,
        primitive_capacity: usize,
    ) -> Self {
        Self {
            result: TracedRay::new(
                Ray::new(Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)),
                geometry,
                interpolation,
            ),
            scratch: TraceScratch::with_capacity(primitive_capacity),
        }
    }
}

impl ReusableTraceStorage2D {
    fn new(interpolation: InterpolationMethod, primitive_capacity: usize) -> Self {
        Self {
            result: TracedRay::new(
                Ray::new(Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)),
                GeometryKind::Spherical,
                interpolation,
            ),
            scratch: TraceScratch2D::with_capacity(primitive_capacity),
        }
    }
}

thread_local! {
    static TRACE_STORAGE: RefCell<Option<ReusableTraceStorage>> = const { RefCell::new(None) };
    static TRACE_STORAGE_2D: RefCell<Option<ReusableTraceStorage2D>> = const { RefCell::new(None) };
}

fn new_vertical_tracer(
    earth_radius: f64,
    altitudes: &CxxVector<f64>,
    refractive_index: &CxxVector<f64>,
    geometry_type: i32,
    interpolation_method: i32,
    sun_x: f64,
    sun_y: f64,
    sun_z: f64,
) -> Result<Box<RustVerticalTracer>> {
    let altitudes: Vec<f64> = altitudes.iter().copied().collect();
    let refractive_index: Vec<f64> = refractive_index.iter().copied().collect();
    let geometry = geometry_from_i32(geometry_type)?;
    let interpolation = interpolation_from_i32(interpolation_method)?;
    let grid = VerticalGrid1D::new(earth_radius, altitudes.clone(), interpolation, geometry)?;
    let tracer = VerticalRayTracer::new(grid);
    let profile = if geometry == GeometryKind::Spherical {
        Some(RefractiveProfile::new_with_interpolation(
            earth_radius,
            altitudes,
            refractive_index,
            interpolation,
        )?)
    } else {
        None
    };

    Ok(Box::new(RustVerticalTracer {
        tracer,
        profile,
        geometry,
        solar: SolarContext::new(Vec3::new(sun_x, sun_y, sun_z), geometry),
    }))
}

fn new_structured_tracer_2d(
    earth_radius: f64,
    altitudes: &CxxVector<f64>,
    horizontal_angles: &CxxVector<f64>,
    interpolation_method: i32,
    reference_x_x: f64,
    reference_x_y: f64,
    reference_x_z: f64,
    reference_z_x: f64,
    reference_z_y: f64,
    reference_z_z: f64,
    sun_x: f64,
    sun_y: f64,
    sun_z: f64,
) -> Result<Box<RustStructuredTracer2D>> {
    let interpolation = interpolation_from_i32(interpolation_method)?;
    let basis = AngularBasis::new(
        Vec3::new(reference_z_x, reference_z_y, reference_z_z),
        Vec3::new(reference_x_x, reference_x_y, reference_x_z),
    )?;
    let grid = StructuredGrid2D::new_with_basis(
        earth_radius,
        altitudes.iter().copied().collect(),
        horizontal_angles.iter().copied().collect(),
        interpolation,
        basis,
    )?;

    Ok(Box::new(RustStructuredTracer2D {
        tracer: StructuredRayTracer2D::new(grid),
        solar: SolarContext::new(Vec3::new(sun_x, sun_y, sun_z), GeometryKind::Spherical),
    }))
}

fn trace_structured_ray_2d_into_cpp_result(
    tracer: &RustStructuredTracer2D,
    origin_x: f64,
    origin_y: f64,
    origin_z: f64,
    look_x: f64,
    look_y: f64,
    look_z: f64,
    refractive_index: &[f64],
    mut result: Pin<&mut ffi::CppTraceResult2D>,
) -> Result<ffi::RustTraceSummary> {
    let profile = if refractive_index.is_empty() {
        None
    } else {
        Some(RefractiveProfile::new_with_interpolation(
            tracer.tracer.grid().earth_radius(),
            tracer.tracer.grid().altitudes().to_vec(),
            refractive_index.to_vec(),
            tracer.tracer.grid().altitude_interpolation(),
        )?)
    };
    let ray = Ray::new(
        Vec3::new(origin_x, origin_y, origin_z),
        Vec3::new(look_x, look_y, look_z),
    );

    Ok(TRACE_STORAGE_2D.with(|storage| {
        let mut storage = storage.borrow_mut();
        let storage = storage.get_or_insert_with(|| {
            ReusableTraceStorage2D::new(
                tracer.tracer.grid().altitude_interpolation(),
                tracer.tracer.primitives().len(),
            )
        });
        tracer.tracer.trace_into(
            ray,
            &mut storage.result,
            &mut storage.scratch,
            TraceOptions2D {
                solar: Some(tracer.solar),
                refraction: profile.as_ref(),
            },
        );

        let summary = summary_to_ffi(&storage.result);
        ffi::prepare_trace_result_2d(result.as_mut(), &summary);
        for (index, layer) in storage.result.layers.iter().enumerate() {
            ffi::set_trace_layer_2d(result.as_mut(), index, &layer_to_ffi(layer));
        }
        summary
    }))
}

fn trace_vertical_ray_with_tracer(
    tracer: &RustVerticalTracer,
    origin_x: f64,
    origin_y: f64,
    origin_z: f64,
    look_x: f64,
    look_y: f64,
    look_z: f64,
    include_refraction: bool,
    mut layers: Pin<&mut CxxVector<ffi::RustTraceLayer>>,
) -> Result<ffi::RustTraceSummary> {
    let ray = Ray::new(
        Vec3::new(origin_x, origin_y, origin_z),
        Vec3::new(look_x, look_y, look_z),
    );
    let traced = trace_with_cached_tracer(tracer, ray, include_refraction);

    for layer in &traced.layers {
        layers.as_mut().push(layer_to_ffi(layer));
    }

    Ok(summary_to_ffi(&traced))
}

fn trace_vertical_ray_into_cpp_result(
    tracer: &RustVerticalTracer,
    origin_x: f64,
    origin_y: f64,
    origin_z: f64,
    look_x: f64,
    look_y: f64,
    look_z: f64,
    include_refraction: bool,
    mut result: Pin<&mut ffi::CppTraceResult>,
) -> Result<ffi::RustTraceSummary> {
    let ray = Ray::new(
        Vec3::new(origin_x, origin_y, origin_z),
        Vec3::new(look_x, look_y, look_z),
    );
    Ok(trace_with_thread_local_storage(
        tracer,
        ray,
        include_refraction,
        |storage| {
            let summary = summary_to_ffi(&storage.result);
            ffi::prepare_trace_result(result.as_mut(), &summary);
            for (index, layer) in storage.result.layers.iter().enumerate() {
                ffi::set_trace_layer(result.as_mut(), index, &layer_to_ffi(layer));
            }
            summary
        },
    ))
}

fn trace_with_cached_tracer(
    tracer: &RustVerticalTracer,
    ray: Ray,
    include_refraction: bool,
) -> TracedRay {
    tracer.tracer.trace(
        ray,
        TraceOptions {
            solar: Some(tracer.solar),
            refraction: if include_refraction && tracer.geometry == GeometryKind::Spherical {
                tracer.profile.as_ref()
            } else {
                None
            },
        },
    )
}

fn trace_with_thread_local_storage<T>(
    tracer: &RustVerticalTracer,
    ray: Ray,
    include_refraction: bool,
    with_traced: impl FnOnce(&mut ReusableTraceStorage) -> T,
) -> T {
    let refraction = if include_refraction && tracer.geometry == GeometryKind::Spherical {
        tracer.profile.as_ref()
    } else {
        None
    };
    TRACE_STORAGE.with(|storage| {
        let mut storage = storage.borrow_mut();
        let storage = storage.get_or_insert_with(|| {
            ReusableTraceStorage::new(
                tracer.geometry,
                tracer.tracer.grid().interpolation_method(),
                tracer.tracer.primitives().len(),
            )
        });
        tracer.tracer.trace_into(
            ray,
            &mut storage.result,
            &mut storage.scratch,
            TraceOptions {
                solar: Some(tracer.solar),
                refraction,
            },
        );
        with_traced(storage)
    })
}

fn trace_vertical_ray(
    earth_radius: f64,
    altitudes: &CxxVector<f64>,
    refractive_index: &CxxVector<f64>,
    geometry_type: i32,
    interpolation_method: i32,
    sun_x: f64,
    sun_y: f64,
    sun_z: f64,
    origin_x: f64,
    origin_y: f64,
    origin_z: f64,
    look_x: f64,
    look_y: f64,
    look_z: f64,
    include_refraction: bool,
    mut layers: Pin<&mut CxxVector<ffi::RustTraceLayer>>,
) -> Result<ffi::RustTraceSummary> {
    let altitudes: Vec<f64> = altitudes.iter().copied().collect();
    let refractive_index: Vec<f64> = refractive_index.iter().copied().collect();
    let geometry = geometry_from_i32(geometry_type)?;
    let interpolation = interpolation_from_i32(interpolation_method)?;

    let grid = VerticalGrid1D::new(earth_radius, altitudes.clone(), interpolation, geometry)?;
    let tracer = VerticalRayTracer::new(grid);
    let profile = if include_refraction && geometry == GeometryKind::Spherical {
        Some(RefractiveProfile::new_with_interpolation(
            earth_radius,
            altitudes,
            refractive_index,
            interpolation,
        )?)
    } else {
        None
    };

    let ray = Ray::new(
        Vec3::new(origin_x, origin_y, origin_z),
        Vec3::new(look_x, look_y, look_z),
    );
    let traced = tracer.trace(
        ray,
        TraceOptions {
            solar: Some(SolarContext::new(Vec3::new(sun_x, sun_y, sun_z), geometry)),
            refraction: profile.as_ref(),
        },
    );

    for layer in &traced.layers {
        layers.as_mut().push(layer_to_ffi(layer));
    }

    Ok(summary_to_ffi(&traced))
}

fn trace_vertical_ray_batch_checksum(
    earth_radius: f64,
    altitudes: &CxxVector<f64>,
    refractive_index: &CxxVector<f64>,
    geometry_type: i32,
    interpolation_method: i32,
    sun_x: f64,
    sun_y: f64,
    sun_z: f64,
    origins: &CxxVector<f64>,
    looks: &CxxVector<f64>,
    include_refraction: bool,
    iterations: usize,
) -> Result<ffi::RustTimingSummary> {
    let altitudes: Vec<f64> = altitudes.iter().copied().collect();
    let refractive_index: Vec<f64> = refractive_index.iter().copied().collect();
    let origins: Vec<f64> = origins.iter().copied().collect();
    let looks: Vec<f64> = looks.iter().copied().collect();

    if origins.len() != looks.len() || !origins.len().is_multiple_of(3) {
        return Err(anyhow!(
            "origins and looks must have the same length and contain xyz triples"
        ));
    }

    let geometry = geometry_from_i32(geometry_type)?;
    let interpolation = interpolation_from_i32(interpolation_method)?;
    let grid = VerticalGrid1D::new(earth_radius, altitudes.clone(), interpolation, geometry)?;
    let tracer = VerticalRayTracer::new(grid);
    let profile = if include_refraction && geometry == GeometryKind::Spherical {
        Some(RefractiveProfile::new_with_interpolation(
            earth_radius,
            altitudes,
            refractive_index,
            interpolation,
        )?)
    } else {
        None
    };

    let rays: Vec<Ray> = origins
        .chunks_exact(3)
        .zip(looks.chunks_exact(3))
        .map(|(origin, look)| {
            Ray::new(
                Vec3::new(origin[0], origin[1], origin[2]),
                Vec3::new(look[0], look[1], look[2]),
            )
        })
        .collect();

    if rays.is_empty() || iterations == 0 {
        return Ok(ffi::RustTimingSummary {
            checksum: 0.0,
            num_traces: 0,
            total_layers: 0,
        });
    }

    let mut result = TracedRay::new(rays[0], geometry, interpolation);
    let mut scratch = TraceScratch::with_capacity(tracer.primitives().len());
    let options = TraceOptions {
        solar: Some(SolarContext::new(Vec3::new(sun_x, sun_y, sun_z), geometry)),
        refraction: profile.as_ref(),
    };

    let mut checksum = 0.0;
    let mut total_layers = 0usize;
    for _ in 0..iterations {
        for &ray in &rays {
            tracer.trace_into(ray, &mut result, &mut scratch, options);
            checksum += traced_ray_checksum(&result);
            total_layers += result.layers.len();
        }
    }

    Ok(ffi::RustTimingSummary {
        checksum,
        num_traces: iterations * rays.len(),
        total_layers,
    })
}

fn traced_ray_checksum(traced: &TracedRay) -> f64 {
    let mut checksum = traced.tangent_radius * 1e-9;
    checksum += if traced.ground_is_hit { 7.0 } else { 11.0 };
    checksum += if traced.is_straight { 13.0 } else { 17.0 };
    checksum += traced.layers.len() as f64;
    for layer in &traced.layers {
        checksum += layer.geometric_distance * 1e-6;
        checksum += layer.curvature_factor;
        checksum += layer.od_quad_start * 1e-7;
        checksum += layer.od_quad_end * 1e-7;
        checksum += layer.cos_sza_entrance * 1e-3;
        checksum += layer.cos_sza_exit * 1e-3;
    }
    checksum
}

fn summary_to_ffi(traced: &TracedRay) -> ffi::RustTraceSummary {
    ffi::RustTraceSummary {
        ground_is_hit: traced.ground_is_hit,
        is_straight: traced.is_straight,
        tangent_radius: traced.tangent_radius,
        num_layers: traced.layers.len(),
    }
}

fn geometry_from_i32(value: i32) -> Result<GeometryKind> {
    match value {
        0 | 1 => Ok(GeometryKind::PlaneParallel),
        2 => Ok(GeometryKind::Spherical),
        _ => Err(anyhow!("unsupported geometry type {value}")),
    }
}

fn interpolation_from_i32(value: i32) -> Result<InterpolationMethod> {
    match value {
        0 => Ok(InterpolationMethod::Shell),
        1 => Ok(InterpolationMethod::Linear),
        2 => Ok(InterpolationMethod::Lower),
        _ => Err(anyhow!("unsupported interpolation method {value}")),
    }
}

fn layer_to_ffi(layer: &Layer) -> ffi::RustTraceLayer {
    let (cell_altitude_index, cell_horizontal_index) = match layer.cell {
        Some(CellId::Structured2D {
            altitude_index,
            horizontal_index,
        }) => (altitude_index as i32, horizontal_index as i32),
        Some(CellId::AltitudeLayer(altitude_index)) => (altitude_index as i32, -1),
        Some(CellId::Unstructured(_)) => (-1, -1),
        None => (-1, -1),
    };
    ffi::RustTraceLayer {
        layer_type: match layer.layer_type {
            LayerType::Complete => 0,
            LayerType::Partial => 1,
            LayerType::Tangent => 2,
        },
        entrance_x: layer.entrance.position.x,
        entrance_y: layer.entrance.position.y,
        entrance_z: layer.entrance.position.z,
        exit_x: layer.exit.position.x,
        exit_y: layer.exit.position.y,
        exit_z: layer.exit.position.z,
        entrance_altitude: layer.entrance.altitude,
        exit_altitude: layer.exit.altitude,
        entrance_on_exact_altitude: layer.entrance.on_exact_vertical_boundary,
        exit_on_exact_altitude: layer.exit.on_exact_vertical_boundary,
        entrance_lower_alt_index: lower_altitude_index(layer.entrance.interpolation),
        exit_lower_alt_index: lower_altitude_index(layer.exit.interpolation),
        layer_distance: layer.geometric_distance,
        curvature_factor: layer.curvature_factor,
        od_quad_start: layer.od_quad_start,
        od_quad_end: layer.od_quad_end,
        od_quad_start_fraction: layer.od_quad_start_fraction,
        od_quad_end_fraction: layer.od_quad_end_fraction,
        cos_sza_entrance: layer.cos_sza_entrance,
        cos_sza_exit: layer.cos_sza_exit,
        saz_entrance: layer.saz_entrance,
        saz_exit: layer.saz_exit,
        cell_altitude_index,
        cell_horizontal_index,
    }
}

fn lower_altitude_index(interpolation: super::InterpolationStencil) -> i32 {
    interpolation
        .iter()
        .next()
        .map(|weight| weight.index as i32)
        .unwrap_or(-1)
}
