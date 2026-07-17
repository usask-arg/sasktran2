use crate::threading;

use super::atmosphere::Atmosphere;
use super::common::openmp_support_enabled;
use super::config::{Config, ThreadingLib};
use super::geometry::{Geometry1D, Geometry2D};
use super::output::Output;
use super::prelude::*;
use super::viewing_geometry::ViewingGeometry;
use rayon::current_thread_index;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use sasktran2_sys::ffi;

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
        let use_rayon_threading = (self.config.threading_lib() == ThreadingLib::Rayon
            || !openmp_support_enabled())
            && self.config.num_threads()? > 1
            && self.config.threading_model()?
                == crate::bindings::config::ThreadingModel::Wavelength;

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
