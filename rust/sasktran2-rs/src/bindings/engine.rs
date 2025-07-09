use crate::threading;

use super::atmosphere::Atmosphere;
use super::common::openmp_support_enabled;
use super::config::{Config, ThreadingLib};
use super::geometry::Geometry1D;
use super::output::Output;
use super::prelude::*;
use super::viewing_geometry::ViewingGeometry;
use rayon::current_thread_index;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use sasktran2_sys::ffi;

/// Wrapper around the c++ sasktran2 Engine object.  Note that we also keep references to the
/// rust Config, Geometry1D and ViewingGeometry objects which is why the lifetime is required
pub struct Engine<'a> {
    pub engine: *mut ffi::Engine,
    pub config: &'a Config,
    pub geometry: &'a Geometry1D,
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

// Newtype wrapper so you own the impl.
pub struct SafeFFIAtmosphere(pub *mut ffi::Atmosphere);

// “Hey Rust, trust me, this is thread‐safe.”
unsafe impl Send for SafeFFIAtmosphere {}
unsafe impl Sync for SafeFFIAtmosphere {}

/// Allows us to call the c++ thread-safe function from multiple threads
fn safe_calc_thread(
    engine: &SafeFFIEngine,
    atmosphere: &SafeFFIAtmosphere,
    output: &SafeFFIOutput,
    wavel: i32,
    thread_idx: i32,
) {
    unsafe {
        ffi::sk_engine_calculate_radiance_thread(
            engine.0,
            atmosphere.0,
            output.0,
            wavel,
            thread_idx,
        );
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
        Ok(Engine {
            engine: unsafe {
                ffi::sk_engine_create(
                    config.config,
                    geometry.geometry,
                    viewing_geometry.viewing_geometry,
                )
            },
            config,
            geometry,
            viewing_geometry,
        })
    }

    pub fn calculate_radiance(&self, atmosphere: &Atmosphere) -> Result<Output> {
        crate::threading::set_num_threads(self.config.num_threads()?)?;

        let num_stokes = self.config.num_stokes()?;
        let num_los = self.viewing_geometry.num_rays()?;
        let num_wavel = atmosphere.num_wavel();

        let mut output = Output::new(num_wavel, num_los, num_stokes);

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

        // We use rayon threading either when the user explicitly enables it, or when
        // openMP support is not enabled.  Also only use it if the number of threads is greater than 1.
        let use_rayon_threading = (self.config.threading_lib() == ThreadingLib::Rayon
            && self.config.num_threads()? > 1)
            || (!openmp_support_enabled() && self.config.num_threads()? > 1)
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
            let safe_atmosphere = SafeFFIAtmosphere(atmosphere.atmosphere);

            let thread_pool = threading::thread_pool()?;

            thread_pool.install(|| {
                (0..num_wavel)
                    .into_par_iter()
                    .with_min_len(min_length)
                    .for_each(|w| {
                        let thread_idx = current_thread_index().unwrap() as i32;
                        if thread_idx > num_threads as i32 {
                            panic!("Thread index out of bounds");
                        }
                        let wavel = w as i32;
                        safe_calc_thread(
                            &safe_engine,
                            &safe_atmosphere,
                            &safe_output,
                            wavel,
                            thread_idx,
                        );
                    })
            })
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
}
