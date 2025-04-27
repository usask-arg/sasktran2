use crate::prelude::*;
use sasktran2_sys::ffi;

#[repr(i32)]
pub enum InterpolationMethod {
    Shell = 0,
    Linear = 1,
    Lower = 2,
}

pub enum GeometryType {
    PlaneParallel = 0,
    PseudoSpherical = 1,
    Spherical = 2,
    Ellipsoidal = 3,
}

pub struct Geometry1D {
    pub geometry: *mut ffi::Geometry1D,
}

impl Geometry1D {
    pub fn new(
        cos_sza: f64,
        saa: f64,
        earth_radius: f64,
        grid_values: Vec<f64>,
        interp_method: InterpolationMethod,
        geotype: GeometryType,
    ) -> Self {
        let ngrid_values = grid_values.len() as i32;
        let geometry = unsafe {
            ffi::sk_geometry1d_create(
                cos_sza,
                saa,
                earth_radius,
                grid_values.clone().as_mut_ptr(),
                ngrid_values,
                interp_method as i32,
                geotype as i32,
            )
        };
        Geometry1D { geometry }
    }

    fn get_num_altitudes(&self) -> i32 {
        unsafe { ffi::sk_geometry1d_get_num_altitudes(self.geometry) }
    }

    pub fn altitudes_m(&self) -> Result<Array1<f64>> {
        let mut altitudes = vec![0.0; self.get_num_altitudes() as usize];
        let altitudes_ptr = altitudes.as_mut_ptr();
        let result = unsafe { ffi::sk_geometry1d_get_altitudes(self.geometry, altitudes_ptr) };
        if result != 0 {
            return Err(anyhow::anyhow!("Failed to get altitudes"));
        }
        Ok(Array1::from(altitudes))
    }

    pub fn refractive_index_mut(&self) -> Result<ArrayViewMut1<f64>> {
        let mut refractive_index = vec![0.0; self.get_num_altitudes() as usize];
        let mut refractive_index_ptr = refractive_index.as_mut_ptr();
        let result = unsafe {
            ffi::sk_geometry1d_get_refractive_index_ptr(self.geometry, &mut refractive_index_ptr)
        };
        if result != 0 {
            return Err(anyhow::anyhow!("Failed to get refractive index"));
        }
        let num_alt = self.get_num_altitudes() as usize;
        unsafe { Ok(ArrayViewMut1::from_shape_ptr(num_alt, refractive_index_ptr)) }
    }
}

impl Drop for Geometry1D {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_geometry1d_destroy(self.geometry);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geometry1d() {
        let cos_sza = 0.5;
        let saa = 0.0;
        let earth_radius = 6371000.0;
        let grid_values = vec![7.0, 8.0, 9.0];

        let geometry = Geometry1D::new(
            cos_sza,
            saa,
            earth_radius,
            grid_values,
            InterpolationMethod::Linear,
            GeometryType::Spherical,
        );

        assert!(!geometry.geometry.is_null(), "Geometry should not be null");

        let num_altitudes = geometry.get_num_altitudes();
        assert_eq!(num_altitudes, 3, "Number of altitudes should be 3");
    }
}
