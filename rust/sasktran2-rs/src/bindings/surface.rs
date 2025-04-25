use sasktran2_sys::ffi;
use super::prelude::*;
use super::brdf::{IsCBRDF, Lambertian};
use ndarray::*;

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
            surface: unsafe { ffi::sk_surface_create(num_wavel as i32, num_stokes as i32, mut_ptr) },
            emission: emission,
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

        let result = unsafe {
            ffi::sk_surface_set_brdf(self.surface, brdf.as_cbrdf(), mut_ptr)
        };

        if result != 0 {
            return Err(anyhow!("Failed to set BRDF for surface: {}", result));
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
        let mut surface = Surface::new(num_wavel, num_stokes);
    }
}
