use crate::threading;

use super::atmosphere::Atmosphere;
use super::common::openmp_support_enabled;
use super::config::{Config, ThreadingLib};
use super::geometry::{Geometry1D, Geometry2D};
use super::output::{JacobianVjpOutput, JvpOutput, Output, VjpOutput};
use super::prelude::*;
use super::viewing_geometry::ViewingGeometry;
use ndarray::{Array1, Array3};
use rayon::current_thread_index;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sasktran2_sys::ffi;
use std::collections::HashMap;

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

pub struct SafeFFIJacobianVjpOutput(pub *mut ffi::OutputJacobianVJP);

unsafe impl Send for SafeFFIJacobianVjpOutput {}
unsafe impl Sync for SafeFFIJacobianVjpOutput {}

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

fn safe_calc_vjp_wavelength_thread(
    engine: &SafeFFIEngine,
    output: &SafeFFIVjpOutput,
    wavelength: i32,
    thread_idx: i32,
) -> Result<()> {
    let result = unsafe {
        ffi::sk_engine_calculate_vjp_wavelength_thread(engine.0, output.0, wavelength, thread_idx)
    };
    if result == 0 {
        Ok(())
    } else {
        Err(anyhow::anyhow!(
            "Failed to calculate VJP wavelength: {}",
            result
        ))
    }
}

fn safe_calc_jacobian_vjp_wavelength_thread(
    engine: &SafeFFIEngine,
    output: &SafeFFIJacobianVjpOutput,
    wavelength: i32,
    thread_idx: i32,
) -> Result<()> {
    let result = unsafe {
        ffi::sk_engine_calculate_jacobian_vjp_wavelength_thread(
            engine.0, output.0, wavelength, thread_idx,
        )
    };
    if result == 0 {
        Ok(())
    } else {
        Err(anyhow::anyhow!(
            "Failed to calculate VJP Jacobian wavelength: {}",
            result
        ))
    }
}

impl<'a> Engine<'a> {
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

    /// Creates a transmission-only engine on a structured 2D geometry.
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
            return Err(anyhow::anyhow!(
                "Failed to create transmission-only Geometry2D Engine"
            ));
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

    pub fn retained_forward_state_bytes(&self) -> Result<usize> {
        let mut bytes = 0usize;
        let result =
            unsafe { ffi::sk_engine_retained_forward_state_bytes(self.engine, &mut bytes) };
        if result == 0 {
            Ok(bytes)
        } else {
            Err(anyhow::anyhow!(
                "Failed to query retained forward state: {}",
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

    fn use_rayon_wavelength_threading(&self) -> Result<bool> {
        Ok(
            (self.config.threading_lib() == ThreadingLib::Rayon || !openmp_support_enabled())
                && self.config.num_threads()? > 1
                && self.config.threading_model()?
                    == crate::bindings::config::ThreadingModel::Wavelength,
        )
    }

    pub fn calculate_radiance(&self, atmosphere: &Atmosphere) -> Result<Output> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;

        let num_stokes = self.config.num_stokes()?;
        let num_los = self.viewing_geometry.num_rays()?;
        let num_flux = self.viewing_geometry.num_flux_observers()?;
        let num_flux_types = self.config.num_flux_types()?;
        let num_wavel = atmosphere.num_wavel();

        let mut output = Output::new(num_wavel, num_los, num_flux, num_flux_types, num_stokes);

        let deriv_names = atmosphere
            .storage
            .derivative_mapping_names()
            .map_err(|e| anyhow::anyhow!(e))?;

        // Assign the memory for the derivatives
        for deriv_name in deriv_names.iter() {
            let mapping = atmosphere
                .storage
                .get_derivative_mapping(deriv_name)
                .map_err(|e| anyhow::anyhow!(e))?;
            let num_deriv_output = mapping.num_output();

            output.with_derivative(deriv_name, num_deriv_output);
        }

        let deriv_names = atmosphere
            .surface
            .derivative_mapping_names()
            .map_err(|e| anyhow::anyhow!(e))?;
        for deriv_name in deriv_names.iter() {
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

            let num_threads = self.config.num_threads()?;

            let safe_engine = SafeFFIEngine(self.engine);
            let safe_output = SafeFFIOutput(output.output);

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
                rayon_for_each_worker(num_batches, num_threads, |batch_index, thread_idx| {
                    let wavelength_start = batch_index * batch_size;
                    let wavelength_count = (num_wavel - wavelength_start).min(batch_size);
                    safe_calc_block_thread(
                        &safe_engine,
                        &safe_output,
                        wavelength_start as i32,
                        wavelength_count as i32,
                        thread_idx,
                    )
                })?;
            } else {
                rayon_for_each_worker(num_wavel, num_threads, |wavelength, thread_idx| {
                    safe_calc_block_thread(
                        &safe_engine,
                        &safe_output,
                        wavelength as i32,
                        1,
                        thread_idx,
                    )
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
        let mut output = JvpOutput::new(
            atmosphere.num_wavel(),
            self.viewing_geometry.num_rays()?,
            self.config.num_stokes()?,
        );
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
            let num_threads = self.config.num_threads()?;
            let initialize_result = unsafe {
                ffi::sk_engine_initialize_jvp(self.engine, atmosphere.atmosphere, output.output)
            };
            if initialize_result != 0 {
                return Err(anyhow::anyhow!(
                    "Failed to initialize JVP: {}",
                    initialize_result
                ));
            }

            let safe_engine = SafeFFIEngine(self.engine);
            let safe_output = SafeFFIJvpOutput(output.output);
            let worker_result = rayon_for_each_worker(
                atmosphere.num_wavel(),
                num_threads,
                |wavelength, thread_idx| {
                    safe_calc_jvp_wavelength_thread(
                        &safe_engine,
                        &safe_output,
                        wavelength as i32,
                        thread_idx,
                    )
                },
            );
            let finalize_result = unsafe { ffi::sk_engine_finalize_jvp(self.engine) };
            worker_result?;
            if finalize_result != 0 {
                return Err(anyhow::anyhow!(
                    "Failed to finalize JVP: {}",
                    finalize_result
                ));
            }
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
        let mut output = VjpOutput::new(cotangent);
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
            let num_threads = self.config.num_threads()?;
            let initialize_result = unsafe {
                ffi::sk_engine_initialize_vjp(self.engine, atmosphere.atmosphere, output.output)
            };
            if initialize_result != 0 {
                return Err(anyhow::anyhow!(
                    "Failed to initialize VJP: {}",
                    initialize_result
                ));
            }

            let safe_engine = SafeFFIEngine(self.engine);
            let safe_output = SafeFFIVjpOutput(output.output);
            rayon_for_each_worker(
                atmosphere.num_wavel(),
                num_threads,
                |wavelength, thread_idx| {
                    safe_calc_vjp_wavelength_thread(
                        &safe_engine,
                        &safe_output,
                        wavelength as i32,
                        thread_idx,
                    )
                },
            )?;
            let finalize_result = unsafe { ffi::sk_output_vjp_finalize(output.output) };
            if finalize_result != 0 {
                return Err(anyhow::anyhow!(
                    "Failed to finalize VJP: {}",
                    finalize_result
                ));
            }
        }
        Ok(output)
    }

    pub fn calculate_jacobian_vjp(
        &self,
        atmosphere: &Atmosphere,
        derivative_sizes: &HashMap<String, usize>,
        surface_sizes: &HashMap<String, usize>,
    ) -> Result<JacobianVjpOutput> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;
        let mut output = JacobianVjpOutput::new(
            atmosphere.num_wavel(),
            self.viewing_geometry.num_rays()?,
            self.config.num_stokes()?,
        );
        for (name, size) in derivative_sizes {
            output
                .with_derivative_jacobian(name, *size)
                .map_err(anyhow::Error::msg)?;
        }
        for (name, size) in surface_sizes {
            output
                .with_surface_jacobian(name, *size)
                .map_err(anyhow::Error::msg)?;
        }

        if !self.use_rayon_wavelength_threading()? {
            let result = unsafe {
                ffi::sk_engine_calculate_jacobian_vjp(
                    self.engine,
                    atmosphere.atmosphere,
                    output.output,
                )
            };
            if result != 0 {
                return Err(anyhow::anyhow!(
                    "Failed to materialize native VJP Jacobian: {}",
                    result
                ));
            }
        } else {
            let num_threads = self.config.num_threads()?;
            let initialize_result = unsafe {
                ffi::sk_engine_initialize_jacobian_vjp(
                    self.engine,
                    atmosphere.atmosphere,
                    output.output,
                )
            };
            if initialize_result != 0 {
                return Err(anyhow::anyhow!(
                    "Failed to initialize native VJP Jacobian: {}",
                    initialize_result
                ));
            }

            let safe_engine = SafeFFIEngine(self.engine);
            let safe_output = SafeFFIJacobianVjpOutput(output.output);
            rayon_for_each_worker(
                atmosphere.num_wavel(),
                num_threads,
                |wavelength, thread_idx| {
                    safe_calc_jacobian_vjp_wavelength_thread(
                        &safe_engine,
                        &safe_output,
                        wavelength as i32,
                        thread_idx,
                    )
                },
            )?;
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
