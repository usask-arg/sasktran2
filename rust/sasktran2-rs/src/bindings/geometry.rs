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

/// Owning wrapper around the C++ structured Geometry2D object.
pub struct Geometry2D {
    pub geometry: *mut ffi::Geometry2D,
}

impl Geometry2D {
    pub fn new(
        cos_sza: f64,
        saa: f64,
        earth_radius: f64,
        altitude_grid_values: Vec<f64>,
        horizontal_angle_grid_values: Vec<f64>,
        altitude_interp_method: InterpolationMethod,
    ) -> Result<Self> {
        let num_altitudes = i32::try_from(altitude_grid_values.len())
            .map_err(|_| anyhow!("Altitude grid is too large"))?;
        let num_horizontal_locations = i32::try_from(horizontal_angle_grid_values.len())
            .map_err(|_| anyhow!("Horizontal angle grid is too large"))?;
        let geometry = unsafe {
            ffi::sk_geometry2d_create(
                cos_sza,
                saa,
                earth_radius,
                altitude_grid_values.as_ptr(),
                num_altitudes,
                horizontal_angle_grid_values.as_ptr(),
                num_horizontal_locations,
                altitude_interp_method as i32,
            )
        };
        if geometry.is_null() {
            return Err(anyhow!("Failed to create Geometry2D"));
        }
        Ok(Self { geometry })
    }

    pub fn location_shape(&self) -> Result<(usize, usize)> {
        let mut num_horizontal_locations = 0;
        let mut num_altitudes = 0;
        let result = unsafe {
            ffi::sk_geometry2d_get_location_shape(
                self.geometry,
                &mut num_horizontal_locations,
                &mut num_altitudes,
            )
        };
        if result != 0 {
            return Err(anyhow!("Failed to get Geometry2D location shape"));
        }
        Ok((num_horizontal_locations as usize, num_altitudes as usize))
    }

    pub fn altitudes_m(&self) -> Result<Array1<f64>> {
        let (_, num_altitudes) = self.location_shape()?;
        let mut altitudes = vec![0.0; num_altitudes];
        let result =
            unsafe { ffi::sk_geometry2d_get_altitudes(self.geometry, altitudes.as_mut_ptr()) };
        if result != 0 {
            return Err(anyhow!("Failed to get Geometry2D altitudes"));
        }
        Ok(Array1::from(altitudes))
    }

    pub fn horizontal_angles(&self) -> Result<Array1<f64>> {
        let (num_horizontal_locations, _) = self.location_shape()?;
        let mut horizontal_angles = vec![0.0; num_horizontal_locations];
        let result = unsafe {
            ffi::sk_geometry2d_get_horizontal_angles(self.geometry, horizontal_angles.as_mut_ptr())
        };
        if result != 0 {
            return Err(anyhow!("Failed to get Geometry2D horizontal angles"));
        }
        Ok(Array1::from(horizontal_angles))
    }

    pub fn location_index(&self, altitude_index: usize, horizontal_index: usize) -> Result<usize> {
        let altitude_index =
            i32::try_from(altitude_index).map_err(|_| anyhow!("Altitude index is too large"))?;
        let horizontal_index = i32::try_from(horizontal_index)
            .map_err(|_| anyhow!("Horizontal index is too large"))?;
        let mut location_index = 0;
        let result = unsafe {
            ffi::sk_geometry2d_get_location_index(
                self.geometry,
                altitude_index,
                horizontal_index,
                &mut location_index,
            )
        };
        if result != 0 {
            return Err(anyhow!("Geometry2D location index is out of range"));
        }
        Ok(location_index as usize)
    }
}

impl Drop for Geometry2D {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_geometry2d_destroy(self.geometry);
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

    #[test]
    fn test_geometry2d_grids_shape_and_indexing() {
        let altitudes = vec![0.0, 10_000.0, 20_000.0];
        let horizontal_angles = vec![-0.2, 0.0, 0.3, 0.6];
        let geometry = Geometry2D::new(
            0.5,
            0.1,
            6_371_000.0,
            altitudes.clone(),
            horizontal_angles.clone(),
            InterpolationMethod::Linear,
        )
        .unwrap();

        assert_eq!(geometry.location_shape().unwrap(), (4, 3));
        assert_eq!(geometry.altitudes_m().unwrap().to_vec(), altitudes);
        assert_eq!(
            geometry.horizontal_angles().unwrap().to_vec(),
            horizontal_angles
        );
        assert_eq!(geometry.location_index(2, 3).unwrap(), 11);
        assert!(geometry.location_index(3, 0).is_err());
        assert!(geometry.location_index(0, 4).is_err());
    }

    #[test]
    fn test_geometry2d_rejects_invalid_grids() {
        assert!(
            Geometry2D::new(
                0.5,
                0.1,
                6_371_000.0,
                vec![0.0],
                vec![-0.2, 0.2],
                InterpolationMethod::Linear,
            )
            .is_err()
        );
        assert!(
            Geometry2D::new(
                0.5,
                0.1,
                6_371_000.0,
                vec![0.0, 10_000.0],
                vec![0.2, -0.2],
                InterpolationMethod::Linear,
            )
            .is_err()
        );
    }
}
