use crate::threading;

use super::atmosphere::Atmosphere;
use super::common::openmp_support_enabled;
use super::config::{Config, ThreadingLib};
use super::geometry::{Geometry1D, Geometry2D};
use super::output::{JvpOutput, Output, VjpOutput};
use super::prelude::*;
use super::viewing_geometry::ViewingGeometry;
use ndarray::{Array1, Array2, Array3, ArrayView2};
use rayon::current_thread_index;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use sasktran2_sys::ffi;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum LinearizationMode {
    Jacobian = 0,
    Jvp = 1,
    Vjp = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum LinearizationBackend {
    Unavailable = 0,
    StreamingJacobian = 1,
    Native = 2,
}

pub enum EngineGeometry<'a> {
    OneDimensional(&'a Geometry1D),
    TwoDimensional(&'a Geometry2D),
}

/// Wrapper around the C++ SASKTRAN2 Engine object. The referenced inputs must
/// outlive the engine because the C++ implementation retains their pointers.
pub struct Engine<'a> {
    pub engine: *mut ffi::Engine,
    pub config: &'a Config,
    pub geometry: EngineGeometry<'a>,
    pub viewing_geometry: &'a ViewingGeometry,
}

// Newtype wrapper so you own the impl.
pub struct SafeFFIEngine(pub *mut ffi::Engine);

// “Hey Rust, trust me, this is thread‐safe.”
unsafe impl Send for SafeFFIEngine {}
unsafe impl Sync for SafeFFIEngine {}

// Newtype wrapper so you own the impl.
pub struct SafeFFIOutput(pub *mut ffi::OutputC);

// “Hey Rust, trust me, this is thread‐safe.”
unsafe impl Send for SafeFFIOutput {}
unsafe impl Sync for SafeFFIOutput {}

pub struct SafeFFIJvpOutput(pub *mut ffi::OutputJVP);

unsafe impl Send for SafeFFIJvpOutput {}
unsafe impl Sync for SafeFFIJvpOutput {}

pub struct SafeFFIVjpOutput(pub *mut ffi::OutputVJP);

unsafe impl Send for SafeFFIVjpOutput {}
unsafe impl Sync for SafeFFIVjpOutput {}

fn rayon_for_each_worker<F>(num_items: usize, num_threads: usize, function: F) -> Result<()>
where
    F: Fn(usize, i32) -> Result<()> + Send + Sync,
{
    let thread_pool = threading::thread_pool()?;
    thread_pool.install(|| {
        (0..num_items).into_par_iter().try_for_each(|item| {
            let thread_idx = current_thread_index()
                .ok_or_else(|| anyhow::anyhow!("Rayon worker index is unavailable"))?
                as i32;
            if thread_idx >= num_threads as i32 {
                return Err(anyhow::anyhow!("Thread index out of bounds"));
            }
            function(item, thread_idx)
        })
    })
}

fn safe_calc_block_thread(
    engine: &SafeFFIEngine,
    output: &SafeFFIOutput,
    wavelength_start: i32,
    wavelength_count: i32,
    thread_idx: i32,
) -> Result<()> {
    let result = unsafe {
        ffi::sk_engine_calculate_radiance_block_thread(
            engine.0,
            output.0,
            wavelength_start,
            wavelength_count,
            thread_idx,
        )
    };
    if result == 0 {
        Ok(())
    } else {
        Err(anyhow::anyhow!(
            "Failed to calculate wavelength block: {}",
            result
        ))
    }
}

fn safe_calc_jvp_wavelength_thread(
    engine: &SafeFFIEngine,
    output: &SafeFFIJvpOutput,
    wavelength: i32,
    thread_idx: i32,
) -> Result<()> {
    let result = unsafe {
        ffi::sk_engine_calculate_jvp_wavelength_thread(engine.0, output.0, wavelength, thread_idx)
    };
    if result == 0 {
        Ok(())
    } else {
        Err(anyhow::anyhow!(
            "Failed to calculate JVP wavelength: {}",
            result
        ))
    }
}

fn safe_calc_vjp_block_thread(
    engine: &SafeFFIEngine,
    output: &SafeFFIVjpOutput,
    wavelength_start: i32,
    wavelength_count: i32,
    thread_idx: i32,
) -> Result<()> {
    let result = unsafe {
        ffi::sk_engine_calculate_vjp_block_thread(
            engine.0,
            output.0,
            wavelength_start,
            wavelength_count,
            thread_idx,
        )
    };
    if result == 0 {
        Ok(())
    } else {
        Err(anyhow::anyhow!(
            "Failed to calculate VJP wavelength block: {}",
            result
        ))
    }
}

struct RadianceBlockTask {
    engine: SafeFFIEngine,
    output: SafeFFIOutput,
    wavelength_start: usize,
    wavelength_count: usize,
}

struct JvpWavelengthTask {
    engine: SafeFFIEngine,
    output: SafeFFIJvpOutput,
    wavelength: usize,
}

struct VjpBlockTask {
    engine: SafeFFIEngine,
    output: SafeFFIVjpOutput,
    wavelength_start: usize,
    wavelength_count: usize,
}

struct RefractiveProfileTask {
    engine: SafeFFIEngine,
    profile_index: usize,
}

/// Execute already-initialized radiance outputs over one shared pool of
/// `(engine, wavelength block)` tasks. This is intended for composite engines
/// whose children share a wavelength-threaded configuration.
pub fn calculate_prepared_radiances(
    engines: &[&Engine<'_>],
    outputs: &[Output],
    num_wavel: usize,
    num_threads: usize,
) -> Result<()> {
    if engines.len() != outputs.len() {
        return Err(anyhow::anyhow!(
            "Prepared radiance engine/output counts do not match"
        ));
    }
    crate::threading::set_num_threads(num_threads)?;
    let mut tasks = Vec::new();
    for (engine, output) in engines.iter().zip(outputs) {
        let batch_size = engine.effective_wavelength_batch_size(num_wavel)?;
        for wavelength_start in (0..num_wavel).step_by(batch_size) {
            tasks.push(RadianceBlockTask {
                engine: SafeFFIEngine(engine.engine),
                output: SafeFFIOutput(output.output),
                wavelength_start,
                wavelength_count: (num_wavel - wavelength_start).min(batch_size),
            });
        }
    }
    rayon_for_each_worker(tasks.len(), num_threads, |task_index, thread_idx| {
        let task = &tasks[task_index];
        safe_calc_block_thread(
            &task.engine,
            &task.output,
            task.wavelength_start as i32,
            task.wavelength_count as i32,
            thread_idx,
        )
    })
}

/// Execute already-initialized native JVP outputs over one shared pool of
/// `(engine, wavelength)` tasks.
pub fn calculate_prepared_jvps(
    engines: &[&Engine<'_>],
    outputs: &[JvpOutput],
    num_wavel: usize,
    num_threads: usize,
) -> Result<()> {
    if engines.len() != outputs.len() {
        return Err(anyhow::anyhow!(
            "Prepared JVP engine/output counts do not match"
        ));
    }
    crate::threading::set_num_threads(num_threads)?;
    let mut tasks = Vec::with_capacity(engines.len() * num_wavel);
    for (engine, output) in engines.iter().zip(outputs) {
        for wavelength in 0..num_wavel {
            tasks.push(JvpWavelengthTask {
                engine: SafeFFIEngine(engine.engine),
                output: SafeFFIJvpOutput(output.output),
                wavelength,
            });
        }
    }
    rayon_for_each_worker(tasks.len(), num_threads, |task_index, thread_idx| {
        let task = &tasks[task_index];
        safe_calc_jvp_wavelength_thread(
            &task.engine,
            &task.output,
            task.wavelength as i32,
            thread_idx,
        )
    })
}

/// Execute already-initialized native VJP outputs over one shared pool of
/// `(engine, wavelength block)` tasks, then reduce every output.
pub fn calculate_prepared_vjps(
    engines: &[&Engine<'_>],
    outputs: &[VjpOutput],
    num_wavel: usize,
    num_threads: usize,
) -> Result<()> {
    if engines.len() != outputs.len() {
        return Err(anyhow::anyhow!(
            "Prepared VJP engine/output counts do not match"
        ));
    }
    crate::threading::set_num_threads(num_threads)?;
    let mut tasks = Vec::new();
    for (engine, output) in engines.iter().zip(outputs) {
        let batch_size = engine.effective_wavelength_batch_size(num_wavel)?;
        for wavelength_start in (0..num_wavel).step_by(batch_size) {
            tasks.push(VjpBlockTask {
                engine: SafeFFIEngine(engine.engine),
                output: SafeFFIVjpOutput(output.output),
                wavelength_start,
                wavelength_count: (num_wavel - wavelength_start).min(batch_size),
            });
        }
    }
    let worker_result =
        rayon_for_each_worker(tasks.len(), num_threads, |task_index, thread_idx| {
            let task = &tasks[task_index];
            safe_calc_vjp_block_thread(
                &task.engine,
                &task.output,
                task.wavelength_start as i32,
                task.wavelength_count as i32,
                thread_idx,
            )
        });
    let mut finalize_result = Ok(());
    for output in outputs {
        let result = unsafe { ffi::sk_output_vjp_finalize(output.output) };
        if result != 0 && finalize_result.is_ok() {
            finalize_result = Err(anyhow::anyhow!(
                "Failed to finalize prepared VJP: {}",
                result
            ));
        }
    }
    worker_result?;
    finalize_result
}

/// Refresh independent structured-2D engines on the shared Rayon pool.
///
/// Every engine must be unique because the underlying C++ call mutates its
/// geometry-dependent state. Per-engine results are retained so a composite
/// caller can update the cache metadata for every successful refresh even if
/// another group rejects its profile.
pub fn set_2d_refractive_profiles_parallel(
    engines: &[&Engine<'_>],
    profiles: &[Array2<f64>],
    num_threads: usize,
) -> Result<Vec<Result<()>>> {
    if engines.len() != profiles.len() {
        return Err(anyhow::anyhow!(
            "Refractive-profile engine/profile counts do not match"
        ));
    }
    if num_threads == 0 {
        return Err(anyhow::anyhow!(
            "Refractive-profile thread count must be positive"
        ));
    }

    let mut engine_addresses = HashSet::with_capacity(engines.len());
    let mut tasks = Vec::with_capacity(engines.len());
    for (profile_index, (engine, profile)) in engines.iter().zip(profiles).enumerate() {
        if !matches!(&engine.geometry, EngineGeometry::TwoDimensional(_)) {
            return Err(anyhow::anyhow!(
                "Per-ray refractive profiles require Geometry2D engines"
            ));
        }
        if !engine_addresses.insert(engine.engine as usize) {
            return Err(anyhow::anyhow!(
                "Parallel refractive-profile refresh requires unique engines"
            ));
        }
        if !profile.is_standard_layout() {
            return Err(anyhow::anyhow!(
                "Parallel refractive profiles must use standard row-major layout"
            ));
        }
        if profile.nrows() > i32::MAX as usize || profile.ncols() > i32::MAX as usize {
            return Err(anyhow::anyhow!(
                "Parallel refractive-profile dimensions exceed the C API range"
            ));
        }
        tasks.push(RefractiveProfileTask {
            engine: SafeFFIEngine(engine.engine),
            profile_index,
        });
    }

    crate::threading::set_num_threads(num_threads)?;
    let thread_pool = threading::thread_pool()?;
    Ok(thread_pool.install(|| {
        tasks
            .into_par_iter()
            .map(|task| {
                let profile = &profiles[task.profile_index];
                let result = unsafe {
                    ffi::sk_engine_set_2d_refractive_profiles(
                        task.engine.0,
                        profile.as_ptr(),
                        profile.nrows() as i32,
                        profile.ncols() as i32,
                    )
                };
                if result == 0 {
                    Ok(())
                } else {
                    Err(anyhow::anyhow!(
                        "Failed to set Geometry2D refractive profiles for engine {}: {}",
                        task.profile_index,
                        result
                    ))
                }
            })
            .collect()
    }))
}

impl<'a> Engine<'a> {
    fn use_rayon_wavelength_threading(&self) -> Result<bool> {
        Ok(
            (self.config.threading_lib() == ThreadingLib::Rayon || !openmp_support_enabled())
                && self.config.num_threads()? > 1
                && self.config.threading_model()?
                    == crate::bindings::config::ThreadingModel::Wavelength,
        )
    }

    /// Creates a new engine object
    pub fn new(
        config: &'a Config,
        geometry: &'a Geometry1D,
        viewing_geometry: &'a ViewingGeometry,
    ) -> Result<Self> {
        crate::threading::set_num_threads(config.num_threads().unwrap_or(1))?;

        // Create the ffi engine first
        let engine = unsafe {
            ffi::sk_engine_create(
                config.config,
                geometry.geometry,
                viewing_geometry.viewing_geometry,
            )
        };

        // Check for null pointer
        if engine.is_null() {
            return Err(anyhow::anyhow!("Failed to create Engine"));
        }

        Ok(Engine {
            engine,
            config,
            geometry: EngineGeometry::OneDimensional(geometry),
            viewing_geometry,
        })
    }

    /// Creates an engine on a structured 2D geometry.
    pub fn new_2d(
        config: &'a Config,
        geometry: &'a Geometry2D,
        viewing_geometry: &'a ViewingGeometry,
    ) -> Result<Self> {
        crate::threading::set_num_threads(config.num_threads().unwrap_or(1))?;

        let engine = unsafe {
            ffi::sk_engine_create_2d(
                config.config,
                geometry.geometry,
                viewing_geometry.viewing_geometry,
            )
        };

        if engine.is_null() {
            return Err(anyhow::anyhow!("Failed to create Geometry2D Engine"));
        }

        Ok(Engine {
            engine,
            config,
            geometry: EngineGeometry::TwoDimensional(geometry),
            viewing_geometry,
        })
    }

    pub fn supports_linearization(&self, mode: LinearizationMode) -> Result<bool> {
        let mut supported = 0i32;
        let result = unsafe {
            ffi::sk_engine_supports_linearization(self.engine, mode as i32, &mut supported)
        };
        if result == 0 {
            Ok(supported != 0)
        } else {
            Err(anyhow::anyhow!(
                "Failed to query linearization support: {}",
                result
            ))
        }
    }

    /// Refresh a structured-2D engine with one altitude-only refractive
    /// profile per line of sight. The C++ engine retains its allocations and
    /// rebuilds only geometry-dependent state.
    pub fn set_2d_refractive_profiles(&mut self, profiles: ArrayView2<'_, f64>) -> Result<()> {
        if !matches!(self.geometry, EngineGeometry::TwoDimensional(_)) {
            return Err(anyhow::anyhow!(
                "Per-ray refractive profiles require a Geometry2D engine"
            ));
        }
        let profiles = profiles.as_standard_layout();
        let shape = profiles.shape();
        let result = unsafe {
            ffi::sk_engine_set_2d_refractive_profiles(
                self.engine,
                profiles.as_ptr(),
                shape[0] as i32,
                shape[1] as i32,
            )
        };
        if result == 0 {
            Ok(())
        } else {
            Err(anyhow::anyhow!(
                "Failed to set Geometry2D refractive profiles: {}",
                result
            ))
        }
    }

    pub fn surface_interpolation_weights(&self) -> Result<Array2<f64>> {
        let num_horizontal = match self.geometry {
            EngineGeometry::TwoDimensional(geometry) => geometry.location_shape()?.0,
            EngineGeometry::OneDimensional(_) => {
                return Err(anyhow::anyhow!(
                    "Surface interpolation weights require a Geometry2D engine"
                ));
            }
        };
        let num_rays = self.viewing_geometry.num_rays()?;
        let mut weights = Array2::zeros((num_rays, num_horizontal));
        let result = unsafe {
            ffi::sk_engine_get_2d_surface_interpolation_weights(
                self.engine,
                weights.as_mut_ptr(),
                num_rays as i32,
                num_horizontal as i32,
            )
        };
        if result == 0 {
            Ok(weights)
        } else {
            Err(anyhow::anyhow!(
                "Failed to obtain Geometry2D surface interpolation weights: {}",
                result
            ))
        }
    }

    pub fn horizontal_edge_usage(&self) -> Result<Array2<i32>> {
        match self.geometry {
            EngineGeometry::TwoDimensional(_) => {}
            EngineGeometry::OneDimensional(_) => {
                return Err(anyhow::anyhow!(
                    "Horizontal edge usage requires a Geometry2D engine"
                ));
            }
        }
        let num_rays = self.viewing_geometry.num_rays()?;
        let mut usage = Array2::zeros((num_rays, 2));
        let result = unsafe {
            ffi::sk_engine_get_2d_horizontal_edge_usage(
                self.engine,
                usage.as_mut_ptr(),
                num_rays as i32,
            )
        };
        if result == 0 {
            Ok(usage)
        } else {
            Err(anyhow::anyhow!(
                "Failed to obtain Geometry2D horizontal edge usage: {}",
                result
            ))
        }
    }

    pub fn linearization_backend(&self, mode: LinearizationMode) -> Result<LinearizationBackend> {
        let mut backend = 0i32;
        let result =
            unsafe { ffi::sk_engine_linearization_backend(self.engine, mode as i32, &mut backend) };
        if result != 0 {
            return Err(anyhow::anyhow!(
                "Failed to select linearization backend: {}",
                result
            ));
        }
        match backend {
            0 => Ok(LinearizationBackend::Unavailable),
            1 => Ok(LinearizationBackend::StreamingJacobian),
            2 => Ok(LinearizationBackend::Native),
            _ => Err(anyhow::anyhow!(
                "Unknown linearization backend: {}",
                backend
            )),
        }
    }

    pub fn calculate_radiance(&self, atmosphere: &Atmosphere) -> Result<Output> {
        let derivative_names = atmosphere
            .storage
            .derivative_mapping_names()
            .map_err(|e| anyhow::anyhow!(e))?;
        let surface_names = atmosphere
            .surface
            .derivative_mapping_names()
            .map_err(|e| anyhow::anyhow!(e))?;
        self.calculate_radiance_with_mappings(atmosphere, &derivative_names, &surface_names)
    }

    pub fn effective_wavelength_batch_size(&self, num_wavel: usize) -> Result<usize> {
        let batch_size = unsafe {
            ffi::sk_engine_effective_wavelength_batch_size(self.engine, num_wavel as i32)
        };
        if batch_size < 1 {
            Err(anyhow::anyhow!(
                "Failed to determine wavelength batch size: {}",
                batch_size
            ))
        } else {
            Ok(batch_size as usize)
        }
    }

    /// Initialize a primal-only output without executing wavelength blocks.
    pub fn initialize_radiance_only(&self, atmosphere: &Atmosphere) -> Result<Output> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;
        let output = Output::new(
            atmosphere.num_wavel(),
            self.viewing_geometry.num_rays()?,
            self.viewing_geometry.num_flux_observers()?,
            self.config.num_flux_types()?,
            self.config.num_stokes()?,
        );
        let result = unsafe {
            ffi::sk_engine_calculate_radiance(self.engine, atmosphere.atmosphere, output.output, 1)
        };
        if result != 0 {
            return Err(anyhow::anyhow!("Failed to initialize radiance: {}", result));
        }
        Ok(output)
    }

    /// Initialize a native JVP output without executing wavelength tasks.
    pub fn initialize_jvp(
        &self,
        atmosphere: &Atmosphere,
        derivative_tangents: &HashMap<String, Array1<f64>>,
        surface_tangents: &HashMap<String, Array1<f64>>,
    ) -> Result<JvpOutput> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;
        let mut output = JvpOutput::try_new(
            atmosphere.num_wavel(),
            self.viewing_geometry.num_rays()?,
            self.config.num_stokes()?,
        )
        .map_err(anyhow::Error::msg)?;
        for (name, tangent) in derivative_tangents {
            output
                .with_derivative_tangent(name, tangent)
                .map_err(anyhow::Error::msg)?;
        }
        for (name, tangent) in surface_tangents {
            output
                .with_surface_tangent(name, tangent)
                .map_err(anyhow::Error::msg)?;
        }
        let result = unsafe {
            ffi::sk_engine_initialize_jvp(self.engine, atmosphere.atmosphere, output.output)
        };
        if result != 0 {
            return Err(anyhow::anyhow!("Failed to initialize JVP: {}", result));
        }
        Ok(output)
    }

    /// Initialize a native VJP output without executing wavelength blocks.
    pub fn initialize_vjp(
        &self,
        atmosphere: &Atmosphere,
        cotangent: &Array3<f64>,
        derivative_sizes: &HashMap<String, usize>,
        surface_sizes: &HashMap<String, usize>,
    ) -> Result<VjpOutput> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;
        let expected = (
            atmosphere.num_wavel(),
            self.viewing_geometry.num_rays()?,
            self.config.num_stokes()?,
        );
        if cotangent.dim() != expected {
            return Err(anyhow::anyhow!(
                "Radiance cotangent shape {:?} does not match {:?}",
                cotangent.dim(),
                expected
            ));
        }
        let mut output = VjpOutput::try_new(cotangent).map_err(anyhow::Error::msg)?;
        for (name, size) in derivative_sizes {
            output
                .with_derivative_gradient(name, *size)
                .map_err(anyhow::Error::msg)?;
        }
        for (name, size) in surface_sizes {
            output
                .with_surface_gradient(name, *size)
                .map_err(anyhow::Error::msg)?;
        }
        let result = unsafe {
            ffi::sk_engine_initialize_vjp(self.engine, atmosphere.atmosphere, output.output)
        };
        if result != 0 {
            return Err(anyhow::anyhow!("Failed to initialize VJP: {}", result));
        }
        Ok(output)
    }

    /// Calculate radiance while allocating only the requested derivative
    /// outputs. An empty pair of mapping collections is a true radiance-only
    /// calculation even when the atmosphere contains derivative mappings.
    pub fn calculate_radiance_with_mappings(
        &self,
        atmosphere: &Atmosphere,
        derivative_names: &[String],
        surface_names: &[String],
    ) -> Result<Output> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;

        let num_stokes = self.config.num_stokes()?;
        let num_los = self.viewing_geometry.num_rays()?;
        let num_flux = self.viewing_geometry.num_flux_observers()?;
        let num_flux_types = self.config.num_flux_types()?;
        let num_wavel = atmosphere.num_wavel();

        let mut output = Output::new(num_wavel, num_los, num_flux, num_flux_types, num_stokes);

        // Assign the memory for the derivatives
        for deriv_name in derivative_names {
            let mapping = atmosphere
                .storage
                .get_derivative_mapping(deriv_name)
                .map_err(|e| anyhow::anyhow!(e))?;
            let num_deriv_output = mapping.num_output();

            output.with_derivative(deriv_name, num_deriv_output);
        }

        for deriv_name in surface_names {
            // Resolve the name here so an invalid selection fails before the
            // native calculation instead of creating an unattached output.
            atmosphere
                .surface
                .get_derivative_mapping(deriv_name)
                .map_err(|e| anyhow::anyhow!(e))?;
            output.with_surface_derivative(deriv_name);
        }

        // Rayon partitions wavelengths, so only use it with wavelength
        // threading. Source threading remains inside the C++ engine.
        let use_rayon_threading = self.use_rayon_wavelength_threading()?;

        if !use_rayon_threading {
            let result = unsafe {
                ffi::sk_engine_calculate_radiance(
                    self.engine,
                    atmosphere.atmosphere,
                    output.output,
                    0,
                )
            };
            if result != 0 {
                return Err(anyhow::anyhow!("Failed to calculate radiance: {}", result));
            }
        } else {
            // Use the thread pool to calculate the radiance
            let result = unsafe {
                ffi::sk_engine_calculate_radiance(
                    self.engine,
                    atmosphere.atmosphere,
                    output.output,
                    1,
                )
            };
            if result != 0 {
                return Err(anyhow::anyhow!("Failed to calculate radiance: {}", result));
            }

            // To determine the min length, start with how many wavelengths we have per thread
            let num_threads = self.config.num_threads()?;
            let wavel_per_thread = (num_wavel / num_threads).min(1);

            // And then take the sqrt
            let min_length = (num_wavel as f64 / wavel_per_thread as f64).sqrt() as usize;

            let min_length = if num_threads >= num_wavel {
                1
            } else {
                min_length
            };

            let safe_engine = SafeFFIEngine(self.engine);
            let safe_output = SafeFFIOutput(output.output);

            let thread_pool = threading::thread_pool()?;
            let batch_size = unsafe {
                ffi::sk_engine_effective_wavelength_batch_size(self.engine, num_wavel as i32)
            };
            if batch_size < 1 {
                return Err(anyhow::anyhow!(
                    "Failed to determine wavelength batch size: {}",
                    batch_size
                ));
            }

            if batch_size > 1 {
                let batch_size = batch_size as usize;
                let num_batches = num_wavel.div_ceil(batch_size);
                thread_pool.install(|| {
                    (0..num_batches)
                        .into_par_iter()
                        .try_for_each(|batch_index| {
                            let thread_idx = current_thread_index().unwrap() as i32;
                            if thread_idx >= num_threads as i32 {
                                return Err(anyhow::anyhow!("Thread index out of bounds"));
                            }
                            let wavelength_start = batch_index * batch_size;
                            let wavelength_count = (num_wavel - wavelength_start).min(batch_size);
                            safe_calc_block_thread(
                                &safe_engine,
                                &safe_output,
                                wavelength_start as i32,
                                wavelength_count as i32,
                                thread_idx,
                            )
                        })
                })?;
            } else {
                thread_pool.install(|| {
                    (0..num_wavel)
                        .into_par_iter()
                        .with_min_len(min_length)
                        .try_for_each(|w| {
                            let thread_idx = current_thread_index().unwrap() as i32;
                            if thread_idx >= num_threads as i32 {
                                return Err(anyhow::anyhow!("Thread index out of bounds"));
                            }
                            safe_calc_block_thread(
                                &safe_engine,
                                &safe_output,
                                w as i32,
                                1,
                                thread_idx,
                            )
                        })
                })?;
            }
        }

        Ok(output)
    }

    pub fn calculate_jvp(
        &self,
        atmosphere: &Atmosphere,
        derivative_tangents: &HashMap<String, Array1<f64>>,
        surface_tangents: &HashMap<String, Array1<f64>>,
    ) -> Result<JvpOutput> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;
        let mut output = JvpOutput::try_new(
            atmosphere.num_wavel(),
            self.viewing_geometry.num_rays()?,
            self.config.num_stokes()?,
        )
        .map_err(anyhow::Error::msg)?;
        for (name, tangent) in derivative_tangents {
            output
                .with_derivative_tangent(name, tangent)
                .map_err(anyhow::Error::msg)?;
        }
        for (name, tangent) in surface_tangents {
            output
                .with_surface_tangent(name, tangent)
                .map_err(anyhow::Error::msg)?;
        }
        let use_rayon_threading = self.use_rayon_wavelength_threading()?
            && self.linearization_backend(LinearizationMode::Jvp)? == LinearizationBackend::Native;
        if !use_rayon_threading {
            let result = unsafe {
                ffi::sk_engine_calculate_jvp(self.engine, atmosphere.atmosphere, output.output)
            };
            if result != 0 {
                return Err(anyhow::anyhow!("Failed to calculate JVP: {}", result));
            }
        } else {
            let result = unsafe {
                ffi::sk_engine_initialize_jvp(self.engine, atmosphere.atmosphere, output.output)
            };
            if result != 0 {
                return Err(anyhow::anyhow!("Failed to initialize JVP: {}", result));
            }
            let engine = SafeFFIEngine(self.engine);
            let output = SafeFFIJvpOutput(output.output);
            rayon_for_each_worker(
                atmosphere.num_wavel(),
                self.config.num_threads()?,
                |wavelength, thread_idx| {
                    safe_calc_jvp_wavelength_thread(&engine, &output, wavelength as i32, thread_idx)
                },
            )?;
        }
        Ok(output)
    }

    pub fn calculate_vjp(
        &self,
        atmosphere: &Atmosphere,
        cotangent: &Array3<f64>,
        derivative_sizes: &HashMap<String, usize>,
        surface_sizes: &HashMap<String, usize>,
    ) -> Result<VjpOutput> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;
        let expected = (
            atmosphere.num_wavel(),
            self.viewing_geometry.num_rays()?,
            self.config.num_stokes()?,
        );
        if cotangent.dim() != expected {
            return Err(anyhow::anyhow!(
                "Radiance cotangent shape {:?} does not match {:?}",
                cotangent.dim(),
                expected
            ));
        }
        let mut output = VjpOutput::try_new(cotangent).map_err(anyhow::Error::msg)?;
        for (name, size) in derivative_sizes {
            output
                .with_derivative_gradient(name, *size)
                .map_err(anyhow::Error::msg)?;
        }
        for (name, size) in surface_sizes {
            output
                .with_surface_gradient(name, *size)
                .map_err(anyhow::Error::msg)?;
        }
        let use_rayon_threading = self.use_rayon_wavelength_threading()?
            && self.linearization_backend(LinearizationMode::Vjp)? == LinearizationBackend::Native;
        if !use_rayon_threading {
            let result = unsafe {
                ffi::sk_engine_calculate_vjp(self.engine, atmosphere.atmosphere, output.output)
            };
            if result != 0 {
                return Err(anyhow::anyhow!("Failed to calculate VJP: {}", result));
            }
        } else {
            let result = unsafe {
                ffi::sk_engine_initialize_vjp(self.engine, atmosphere.atmosphere, output.output)
            };
            if result != 0 {
                return Err(anyhow::anyhow!("Failed to initialize VJP: {}", result));
            }
            let num_wavel = atmosphere.num_wavel();
            let batch_size = unsafe {
                ffi::sk_engine_effective_wavelength_batch_size(self.engine, num_wavel as i32)
            };
            if batch_size < 1 {
                return Err(anyhow::anyhow!(
                    "Failed to determine VJP wavelength batch size: {}",
                    batch_size
                ));
            }
            let batch_size = batch_size as usize;
            let num_batches = num_wavel.div_ceil(batch_size);
            let engine = SafeFFIEngine(self.engine);
            let safe_output = SafeFFIVjpOutput(output.output);
            let worker_result = rayon_for_each_worker(
                num_batches,
                self.config.num_threads()?,
                |batch_index, thread_idx| {
                    let start = batch_index * batch_size;
                    let count = (num_wavel - start).min(batch_size);
                    safe_calc_vjp_block_thread(
                        &engine,
                        &safe_output,
                        start as i32,
                        count as i32,
                        thread_idx,
                    )
                },
            );
            let finalize_result = unsafe { ffi::sk_output_vjp_finalize(output.output) };
            worker_result?;
            if finalize_result != 0 {
                return Err(anyhow::anyhow!(
                    "Failed to finalize VJP: {}",
                    finalize_result
                ));
            }
        }
        Ok(output)
    }
}
impl<'a> Drop for Engine<'a> {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_engine_destroy(self.engine);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::config::{
        EmissionSource, MultipleScatterSource, OccultationSource, SingleScatterSource,
    };
    use super::super::geometry::Geometry2D;
    use super::super::geometry::GeometryType;
    use super::super::geometry::InterpolationMethod;

    use super::*;
    use ndarray::ShapeBuilder;

    #[test]
    fn test_engine() {
        let config = Config::new();
        let cos_sza = 0.5;
        let saa = 0.0;
        let earth_radius = 6371000.0;
        let grid_values = vec![1.0, 2.0, 3.0];
        let interp_method = InterpolationMethod::Linear;
        let geotype = GeometryType::Spherical;

        let geometry = Geometry1D::new(
            cos_sza,
            saa,
            earth_radius,
            grid_values,
            interp_method,
            geotype,
        );
        let viewing_geometry = ViewingGeometry::new();

        let engine = Engine::new(&config, &geometry, &viewing_geometry).unwrap();

        assert!(!engine.engine.is_null());
        assert!(
            engine
                .supports_linearization(LinearizationMode::Jacobian)
                .unwrap()
        );
        assert!(
            engine
                .supports_linearization(LinearizationMode::Jvp)
                .unwrap()
        );
        assert!(
            engine
                .supports_linearization(LinearizationMode::Vjp)
                .unwrap()
        );
        assert_eq!(
            engine
                .linearization_backend(LinearizationMode::Jacobian)
                .unwrap(),
            LinearizationBackend::Native
        );
        assert_eq!(
            engine
                .linearization_backend(LinearizationMode::Jvp)
                .unwrap(),
            LinearizationBackend::Native
        );
        assert_eq!(
            engine
                .linearization_backend(LinearizationMode::Vjp)
                .unwrap(),
            LinearizationBackend::Native
        );
    }

    #[test]
    fn test_native_products_reject_derivative_free_atmosphere() {
        let config = Config::new();
        let geometry = Geometry1D::new(
            0.5,
            0.0,
            6_371_000.0,
            vec![0.0, 10_000.0, 20_000.0],
            InterpolationMethod::Linear,
            GeometryType::Spherical,
        );
        let mut viewing_geometry = ViewingGeometry::new();
        viewing_geometry.add_ground_viewing_solar(0.5, 0.0, 100_000.0, 0.5);
        let engine = Engine::new(&config, &geometry, &viewing_geometry).unwrap();
        let atmosphere = Atmosphere::new(1, 3, 3, false, false, Stokes::Stokes1, None);

        let derivative_tangents: HashMap<String, Array1<f64>> = HashMap::new();
        let surface_tangents: HashMap<String, Array1<f64>> = HashMap::new();
        assert!(
            engine
                .calculate_jvp(&atmosphere, &derivative_tangents, &surface_tangents)
                .is_err()
        );

        let cotangent = Array3::zeros((1, 1, 1));
        let derivative_sizes: HashMap<String, usize> = HashMap::new();
        let surface_sizes: HashMap<String, usize> = HashMap::new();
        assert!(
            engine
                .calculate_vjp(&atmosphere, &cotangent, &derivative_sizes, &surface_sizes)
                .is_err()
        );
    }

    #[test]
    fn test_native_products_reject_incorrect_parameter_size() {
        let config = Config::new();
        let geometry = Geometry1D::new(
            0.5,
            0.0,
            6_371_000.0,
            vec![0.0, 10_000.0, 20_000.0],
            InterpolationMethod::Linear,
            GeometryType::Spherical,
        );
        let mut viewing_geometry = ViewingGeometry::new();
        viewing_geometry.add_ground_viewing_solar(0.5, 0.0, 100_000.0, 0.5);
        let engine = Engine::new(&config, &geometry, &viewing_geometry).unwrap();
        let mut atmosphere = Atmosphere::new(
            1,
            3,
            config.num_singlescatter_moments().unwrap(),
            true,
            false,
            Stokes::Stokes1,
            None,
        );
        atmosphere.storage.total_extinction.fill(1e-5);
        atmosphere.storage.ssa.fill(0.9);
        for location in 0..atmosphere.num_location() {
            atmosphere.storage.leg_coeff[[0, location, 0]] = 1.0;
        }
        atmosphere
            .storage
            .get_derivative_mapping("volume")
            .unwrap()
            .d_extinction()
            .fill(1.0);

        let derivative_tangents = HashMap::from([(
            "volume".to_string(),
            Array1::zeros(atmosphere.num_location() - 1),
        )]);
        let surface_tangents: HashMap<String, Array1<f64>> = HashMap::new();
        assert!(
            engine
                .calculate_jvp(&atmosphere, &derivative_tangents, &surface_tangents)
                .is_err()
        );

        let cotangent = Array3::ones((1, 1, 1));
        let derivative_sizes =
            HashMap::from([("volume".to_string(), atmosphere.num_location() - 1)]);
        let surface_sizes: HashMap<String, usize> = HashMap::new();
        assert!(
            engine
                .calculate_vjp(&atmosphere, &cotangent, &derivative_sizes, &surface_sizes)
                .is_err()
        );
    }

    #[test]
    fn test_native_products_accept_nonstandard_ndarray_layouts() {
        let config = Config::new();
        let geometry = Geometry1D::new(
            0.5,
            0.0,
            6_371_000.0,
            vec![0.0, 10_000.0, 20_000.0],
            InterpolationMethod::Linear,
            GeometryType::Spherical,
        );
        let mut viewing_geometry = ViewingGeometry::new();
        viewing_geometry.add_ground_viewing_solar(0.5, 0.1, 100_000.0, 0.5);
        viewing_geometry.add_tangent_altitude_solar(5_000.0, -0.2, 100_000.0, 0.5);
        let engine = Engine::new(&config, &geometry, &viewing_geometry).unwrap();
        let mut atmosphere = Atmosphere::new(
            2,
            3,
            config.num_singlescatter_moments().unwrap(),
            true,
            false,
            Stokes::Stokes1,
            None,
        );
        atmosphere.storage.total_extinction.fill(1e-5);
        atmosphere.storage.ssa.fill(0.9);
        for wavelength in 0..atmosphere.num_wavel() {
            for location in 0..atmosphere.num_location() {
                atmosphere.storage.leg_coeff[[0, location, wavelength]] = 1.0;
            }
        }
        let mapping = atmosphere.storage.get_derivative_mapping("volume").unwrap();
        for ((location, wavelength), value) in mapping.d_extinction().indexed_iter_mut() {
            *value = 0.5 + location as f64 + 0.25 * wavelength as f64;
        }

        let standard_tangent = ndarray::array![1.0, 2.0, 4.0];
        let strided_tangent = ndarray::array![4.0, 2.0, 1.0].slice_move(ndarray::s![..;-1]);
        assert_eq!(standard_tangent, strided_tangent);
        assert!(!strided_tangent.is_standard_layout());

        let no_surface_tangents: HashMap<String, Array1<f64>> = HashMap::new();
        let standard_jvp = engine
            .calculate_jvp(
                &atmosphere,
                &HashMap::from([("volume".to_string(), standard_tangent)]),
                &no_surface_tangents,
            )
            .unwrap();
        let strided_jvp = engine
            .calculate_jvp(
                &atmosphere,
                &HashMap::from([("volume".to_string(), strided_tangent)]),
                &no_surface_tangents,
            )
            .unwrap();
        assert!(standard_jvp.jvp.iter().any(|value| value.abs() > 1e-15));
        for (standard, strided) in standard_jvp.jvp.iter().zip(strided_jvp.jvp.iter()) {
            assert!((standard - strided).abs() <= 1e-12 * (1.0 + standard.abs()));
        }

        let standard_cotangent = Array3::from_shape_fn((2, 2, 1), |(wavelength, los, _)| {
            1.0 + 2.0 * wavelength as f64 + 0.5 * los as f64
        });
        let mut fortran_cotangent = Array3::zeros((2, 2, 1).f());
        for ((wavelength, los, stokes), value) in fortran_cotangent.indexed_iter_mut() {
            *value = standard_cotangent[[wavelength, los, stokes]];
        }
        assert_eq!(standard_cotangent, fortran_cotangent);
        assert!(!fortran_cotangent.is_standard_layout());

        let derivative_sizes = HashMap::from([("volume".to_string(), 3)]);
        let no_surface_sizes: HashMap<String, usize> = HashMap::new();
        let standard_vjp = engine
            .calculate_vjp(
                &atmosphere,
                &standard_cotangent,
                &derivative_sizes,
                &no_surface_sizes,
            )
            .unwrap();
        let fortran_vjp = engine
            .calculate_vjp(
                &atmosphere,
                &fortran_cotangent,
                &derivative_sizes,
                &no_surface_sizes,
            )
            .unwrap();
        let standard_gradient = &standard_vjp.derivative_gradients["volume"];
        let fortran_gradient = &fortran_vjp.derivative_gradients["volume"];
        assert!(standard_gradient.iter().any(|value| value.abs() > 1e-15));
        for (standard, fortran) in standard_gradient.iter().zip(fortran_gradient.iter()) {
            assert!((standard - fortran).abs() <= 1e-12 * (1.0 + standard.abs()));
        }
    }

    #[test]
    fn test_native_product_engine_rejects_stokes_mismatches() {
        let config = Config::new();
        let geometry = Geometry1D::new(
            0.5,
            0.0,
            6_371_000.0,
            vec![0.0, 10_000.0, 20_000.0],
            InterpolationMethod::Linear,
            GeometryType::Spherical,
        );
        let mut viewing_geometry = ViewingGeometry::new();
        viewing_geometry.add_ground_viewing_solar(0.5, 0.0, 100_000.0, 0.5);
        let engine = Engine::new(&config, &geometry, &viewing_geometry).unwrap();
        let atmosphere = Atmosphere::new(1, 3, 3, true, false, Stokes::Stokes1, None);

        let jvp = JvpOutput::try_new(1, 1, 3).unwrap();
        assert_eq!(
            unsafe {
                ffi::sk_engine_calculate_jvp(engine.engine, atmosphere.atmosphere, jvp.output)
            },
            -2
        );

        let cotangent = Array3::ones((1, 1, 3));
        let vjp = VjpOutput::try_new(&cotangent).unwrap();
        assert_eq!(
            unsafe {
                ffi::sk_engine_calculate_vjp(engine.engine, atmosphere.atmosphere, vjp.output)
            },
            -2
        );

        let atmosphere_stokes3 = Atmosphere::new(1, 3, 3, true, false, Stokes::Stokes3, None);
        let jvp_stokes1 = JvpOutput::try_new(1, 1, 1).unwrap();
        assert_eq!(
            unsafe {
                ffi::sk_engine_calculate_jvp(
                    engine.engine,
                    atmosphere_stokes3.atmosphere,
                    jvp_stokes1.output,
                )
            },
            -2
        );
    }

    #[test]
    fn test_geometry2d_engine_configuration_boundary() {
        let mut config = Config::new();
        let geometry = Geometry2D::new(
            0.5,
            0.0,
            6_371_000.0,
            vec![0.0, 10_000.0, 30_000.0],
            vec![-0.5, 0.0, 0.5],
            InterpolationMethod::Linear,
        )
        .unwrap();
        let mut viewing_geometry = ViewingGeometry::new();
        viewing_geometry.add_tangent_altitude_solar(15_000.0, 0.0, 100_000.0, 0.5);

        let engine = Engine::new_2d(&config, &geometry, &viewing_geometry).unwrap();
        assert!(!engine.engine.is_null());
        drop(engine);

        config
            .with_single_scatter_source(SingleScatterSource::SolarTable)
            .unwrap();
        let engine = Engine::new_2d(&config, &geometry, &viewing_geometry).unwrap();
        assert!(!engine.engine.is_null());
        drop(engine);

        config
            .with_multiple_scatter_source(MultipleScatterSource::SuccessiveOrders)
            .unwrap();
        let engine = Engine::new_2d(&config, &geometry, &viewing_geometry).unwrap();
        assert!(!engine.engine.is_null());
        drop(engine);

        config
            .with_multiple_scatter_source(MultipleScatterSource::SuccessiveOrdersLegacy)
            .unwrap();
        assert!(Engine::new_2d(&config, &geometry, &viewing_geometry).is_err());

        config
            .with_multiple_scatter_source(MultipleScatterSource::None)
            .unwrap()
            .with_single_scatter_source(SingleScatterSource::DiscreteOrdinates)
            .unwrap();
        assert!(Engine::new_2d(&config, &geometry, &viewing_geometry).is_err());

        config
            .with_single_scatter_source(SingleScatterSource::None)
            .unwrap()
            .with_multiple_scatter_source(MultipleScatterSource::None)
            .unwrap()
            .with_emission_source(EmissionSource::None)
            .unwrap()
            .with_occultation_source(OccultationSource::Standard)
            .unwrap();

        let engine = Engine::new_2d(&config, &geometry, &viewing_geometry).unwrap();
        assert!(!engine.engine.is_null());
    }
}
