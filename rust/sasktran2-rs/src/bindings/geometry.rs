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

/// Wrapper around the c++ Geometry1D object
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

    pub fn refractive_index_mut(&self) -> Result<ArrayViewMut1<'_, f64>> {
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

    #[test]
    fn test_geometry1d_altitudes() {
        let cos_sza = 0.5;
        let saa = 0.0;
        let earth_radius = 6371000.0;
        let grid_values = vec![7000.0, 8000.0, 9000.0];

        let geometry = Geometry1D::new(
            cos_sza,
            saa,
            earth_radius,
            grid_values.clone(),
            InterpolationMethod::Linear,
            GeometryType::Spherical,
        );

        let altitudes = geometry.altitudes_m().unwrap();
        assert_eq!(
            altitudes.len(),
            grid_values.len(),
            "Altitude array length should match input grid values"
        );

        for (i, &val) in grid_values.iter().enumerate() {
            assert_eq!(
                altitudes[i], val,
                "Altitude values should match input grid values"
            );
        }
    }

    #[test]
    fn test_geometry1d_refractive_index() {
        let cos_sza = 0.6;
        let saa = 45.0;
        let earth_radius = 6371000.0;
        let grid_values = vec![10000.0, 20000.0, 30000.0, 40000.0];

        let geometry = Geometry1D::new(
            cos_sza,
            saa,
            earth_radius,
            grid_values,
            InterpolationMethod::Shell,
            GeometryType::PseudoSpherical,
        );

        let mut ref_index = geometry.refractive_index_mut().unwrap();
        assert_eq!(
            ref_index.len(),
            4,
            "Refractive index array length should match grid size"
        );

        // Default values should be zeros
        for val in ref_index.iter() {
            assert_eq!(*val, 1.0, "Default refractive index should be 1.0");
        }

        // Test we can modify the refractive index
        ref_index[0] = 1.0;
        ref_index[1] = 1.1;
        assert_eq!(ref_index[0], 1.0);
        assert_eq!(ref_index[1], 1.1);
    }

    #[test]
    fn test_different_interpolation_methods() {
        let earth_radius = 6371000.0;
        let grid_values = vec![5000.0, 10000.0];

        let geom_linear = Geometry1D::new(
            0.5,
            0.0,
            earth_radius,
            grid_values.clone(),
            InterpolationMethod::Linear,
            GeometryType::Spherical,
        );

        let geom_shell = Geometry1D::new(
            0.5,
            0.0,
            earth_radius,
            grid_values.clone(),
            InterpolationMethod::Shell,
            GeometryType::Spherical,
        );

        let geom_lower = Geometry1D::new(
            0.5,
            0.0,
            earth_radius,
            grid_values,
            InterpolationMethod::Lower,
            GeometryType::Spherical,
        );

        assert!(!geom_linear.geometry.is_null());
        assert!(!geom_shell.geometry.is_null());
        assert!(!geom_lower.geometry.is_null());
    }
}
