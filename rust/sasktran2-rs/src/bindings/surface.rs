use super::brdf::{IsCBRDF, Lambertian};
use super::deriv_mapping::SurfaceDerivativeMapping;
use super::prelude::*;
use ndarray::*;
use sasktran2_sys::ffi;

/// C++ wrapper around the Surface object
pub struct Surface {
    pub surface: *mut ffi::Surface,
    pub emission: Array1<f64>,
    pub brdf_args: Array2<f64>,
}

impl Surface {
    pub fn new(num_wavel: usize, num_stokes: usize) -> Self {
        let mut emission = Array1::zeros(num_wavel);
        let mut_ptr = emission.as_mut_ptr();
        let mut surf = Surface {
            surface: unsafe {
                ffi::sk_surface_create(num_wavel as i32, num_stokes as i32, mut_ptr)
            },
            emission,
            brdf_args: Array2::zeros((0, num_wavel).f()),
        };

        // Set the default brdf to be lambertian
        let brdf = Lambertian::new(num_stokes);
        surf.set_brdf(&brdf).unwrap();

        surf
    }

    pub fn set_brdf(&mut self, brdf: &dyn IsCBRDF) -> Result<()> {
        let num_brdf_args = brdf.num_args()?;
        self.brdf_args = Array2::zeros((num_brdf_args, self.emission.len()).f());
        let mut_ptr = self.brdf_args.as_mut_ptr();

        let result = unsafe { ffi::sk_surface_set_brdf(self.surface, brdf.as_cbrdf(), mut_ptr) };

        if result != 0 {
            return Err(anyhow!("Failed to set BRDF for surface: {}", result));
        }

        Ok(())
    }

    pub fn get_derivative_mapping(&self, name: &str) -> Result<SurfaceDerivativeMapping, String> {
        let mut mapping: *mut ffi::SurfaceDerivativeMapping = std::ptr::null_mut();
        let c_name = std::ffi::CString::new(name).unwrap();

        let result = unsafe {
            ffi::sk_surface_get_derivative_mapping(self.surface, c_name.as_ptr(), &mut mapping)
        };

        if result != 0 {
            return Err("Failed to get derivative mapping".to_string());
        }

        Ok(SurfaceDerivativeMapping::new(mapping))
    }

    pub fn derivative_mapping_names(&self) -> Result<Vec<String>, String> {
        let mut names: Vec<String> = Vec::new();

        let mut num_mappings: i32 = 0;
        let result =
            unsafe { ffi::sk_surface_get_num_derivative_mappings(self.surface, &mut num_mappings) };
        if result != 0 {
            return Err("Failed to get number of derivative mappings".to_string());
        }

        for i in 0..num_mappings {
            let mut name_ptr: *const std::ffi::c_char = std::ptr::null();
            let result = unsafe {
                ffi::sk_surface_get_derivative_mapping_name(self.surface, i, &mut name_ptr)
            };
            if result != 0 {
                return Err("Failed to get derivative mapping name".to_string());
            }

            let c_str = unsafe { std::ffi::CStr::from_ptr(name_ptr) };
            let name = c_str.to_string_lossy().into_owned();
            names.push(name);
        }

        Ok(names)
    }

    pub fn set_zero(&mut self) -> Result<()> {
        let result = unsafe { ffi::sk_surface_set_zero(self.surface) };
        if result != 0 {
            return Err(anyhow!("Failed to set surface to zero: {}", result));
        }
        Ok(())
    }
}

impl Drop for Surface {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_surface_destroy(self.surface);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_surface() {
        let num_wavel = 10;
        let num_stokes = 3;
        let _surface = Surface::new(num_wavel, num_stokes);
    }

    #[test]
    fn test_surface_set_brdf() {
        let num_wavel = 5;
        let num_stokes = 1;
        let mut surface = Surface::new(num_wavel, num_stokes);

        let brdf = Lambertian::new(num_stokes);
        let result = surface.set_brdf(&brdf);
        assert!(result.is_ok());
    }

    #[test]
    fn test_surface_set_zero() {
        let num_wavel = 3;
        let num_stokes = 1;
        let mut surface = Surface::new(num_wavel, num_stokes);

        let result = surface.set_zero();
        assert!(result.is_ok());
    }

    #[test]
    fn test_surface_derivative_mapping_names() {
        let num_wavel = 2;
        let num_stokes = 1;
        let surface = Surface::new(num_wavel, num_stokes);

        let names = surface.derivative_mapping_names().unwrap();
        assert!(names.is_empty());
    }
}
