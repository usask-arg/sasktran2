use crate::ffi;

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
