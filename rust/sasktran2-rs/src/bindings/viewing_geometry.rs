use super::prelude::*;
use sasktran2_sys::ffi;

/// Wrapper around the c++ viewing geometry object
/// This works a little differently than the c++ api in that we provide functions to add rays
/// instead of having separate types passed in
pub struct ViewingGeometry {
    pub viewing_geometry: *mut ffi::ViewingGeometry,
}

impl Default for ViewingGeometry {
    fn default() -> Self {
        ViewingGeometry::new()
    }
}

impl ViewingGeometry {
    pub fn new() -> Self {
        ViewingGeometry {
            viewing_geometry: unsafe { ffi::sk_viewing_geometry_create() },
        }
    }

    pub fn add_ground_viewing_solar(
        &mut self,
        cos_sza: f64,
        relative_azimuth_angle: f64,
        observer_altitude: f64,
        cos_viewing_zenith: f64,
    ) {
        unsafe {
            ffi::sk_viewing_geometry_add_ground_viewing_solar(
                self.viewing_geometry,
                cos_sza,
                relative_azimuth_angle,
                observer_altitude,
                cos_viewing_zenith,
            );
        }
    }

    pub fn add_tangent_altitude_solar(
        &mut self,
        tangent_altitude_m: f64,
        relative_azimuth_angle: f64,
        observer_altitude: f64,
        cos_sza: f64,
    ) {
        unsafe {
            ffi::sk_viewing_geometry_add_tangent_altitude_solar(
                self.viewing_geometry,
                tangent_altitude_m,
                relative_azimuth_angle,
                observer_altitude,
                cos_sza,
            );
        }
    }

    pub fn add_solar_angles_observer_location(
        &mut self,
        cos_sza: f64,
        relative_azimuth_angle: f64,
        cos_viewing_zenith: f64,
        observer_altitude: f64,
    ) {
        unsafe {
            ffi::sk_viewing_geometry_add_solar_angles_observer_location(
                self.viewing_geometry,
                cos_sza,
                relative_azimuth_angle,
                cos_viewing_zenith,
                observer_altitude,
            );
        }
    }

    pub fn num_rays(&self) -> Result<usize> {
        let mut num_rays = 0i32;
        let error_code =
            unsafe { ffi::sk_viewing_geometry_num_rays(self.viewing_geometry, &mut num_rays) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of rays: error code {}",
                error_code
            ))
        } else {
            Ok(num_rays as usize)
        }
    }
}

impl Drop for ViewingGeometry {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_viewing_geometry_destroy(self.viewing_geometry);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_viewing_geometry() {
        let mut viewing_geometry = ViewingGeometry::new();
        let cos_sza = 0.5;
        let relative_azimuth_angle = 0.0;
        let observer_altitude = 6371000.0;
        let cos_viewing_zenith = 1.0;

        viewing_geometry.add_ground_viewing_solar(
            cos_sza,
            relative_azimuth_angle,
            observer_altitude,
            cos_viewing_zenith,
        );
    }

    #[test]
    fn test_viewing_geometry_tangent_altitude() {
        let mut viewing_geometry = ViewingGeometry::new();
        let tangent_altitude_m = 30000.0;
        let relative_azimuth_angle = 0.0;
        let observer_altitude = 700000.0;
        let cos_sza = 0.6;

        viewing_geometry.add_tangent_altitude_solar(
            tangent_altitude_m,
            relative_azimuth_angle,
            observer_altitude,
            cos_sza,
        );

        assert_eq!(viewing_geometry.num_rays().unwrap(), 1);
    }

    #[test]
    fn test_viewing_geometry_solar_angles_observer() {
        let mut viewing_geometry = ViewingGeometry::new();
        let cos_sza = 0.7;
        let relative_azimuth_angle = 45.0;
        let cos_viewing_zenith = 0.8;
        let observer_altitude = 8000000.0;

        viewing_geometry.add_solar_angles_observer_location(
            cos_sza,
            relative_azimuth_angle,
            cos_viewing_zenith,
            observer_altitude,
        );

        assert_eq!(viewing_geometry.num_rays().unwrap(), 1);
    }

    #[test]
    fn test_multiple_rays() {
        let mut viewing_geometry = ViewingGeometry::new();

        // Add three different ray types
        viewing_geometry.add_ground_viewing_solar(0.5, 0.0, 6371000.0, 1.0);
        viewing_geometry.add_tangent_altitude_solar(30000.0, 90.0, 700000.0, 0.6);
        viewing_geometry.add_solar_angles_observer_location(0.7, 180.0, 0.8, 8000000.0);

        assert_eq!(viewing_geometry.num_rays().unwrap(), 3);
    }

    #[test]
    fn test_default_creation() {
        let viewing_geometry = ViewingGeometry::default();
        assert_eq!(viewing_geometry.num_rays().unwrap(), 0);
    }
}
