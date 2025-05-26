use super::prelude::*;
use sasktran2_sys::ffi;

/// Wrapper around the c++ geodetic object
pub struct Geodetic {
    pub geodetic: *mut ffi::Geodetic,
}

impl Drop for Geodetic {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_geodetic_destroy(self.geodetic);
        }
    }
}

impl Geodetic {
    pub fn new(equatorial_radius: f64, flattening_factor: f64) -> Result<Self> {
        let geodetic = unsafe { ffi::sk_geodetic_create(equatorial_radius, flattening_factor) };
        if geodetic.is_null() {
            return Err(anyhow!("Failed to create Geodetic object"));
        }
        Ok(Self { geodetic })
    }

    pub fn altitude(&self) -> Result<f64> {
        let mut altitude = 0.0;
        let result = unsafe { ffi::sk_geodetic_get_altitude(self.geodetic, &mut altitude) };
        if result != 0 {
            return Err(anyhow!("Failed to get altitude"));
        }
        Ok(altitude)
    }

    pub fn latitude(&self) -> Result<f64> {
        let mut latitude = 0.0;
        let result = unsafe { ffi::sk_geodetic_get_latitude(self.geodetic, &mut latitude) };
        if result != 0 {
            return Err(anyhow!("Failed to get latitude"));
        }
        Ok(latitude)
    }

    pub fn longitude(&self) -> Result<f64> {
        let mut longitude = 0.0;
        let result = unsafe { ffi::sk_geodetic_get_longitude(self.geodetic, &mut longitude) };
        if result != 0 {
            return Err(anyhow!("Failed to get longitude"));
        }
        Ok(longitude)
    }

    pub fn location(&self) -> Result<[f64; 3]> {
        let mut loc: [f64; 3] = [0.0; 3];
        let result = unsafe {
            ffi::sk_geodetic_get_location(self.geodetic, &mut loc[0], &mut loc[1], &mut loc[2])
        };
        if result != 0 {
            return Err(anyhow!("Failed to get location"));
        }
        Ok(loc)
    }

    pub fn local_south(&self) -> Result<[f64; 3]> {
        let mut south: [f64; 3] = [0.0; 3];
        let result = unsafe {
            ffi::sk_geodetic_get_local_south(
                self.geodetic,
                &mut south[0],
                &mut south[1],
                &mut south[2],
            )
        };
        if result != 0 {
            return Err(anyhow!("Failed to get local south"));
        }
        Ok(south)
    }

    pub fn local_up(&self) -> Result<[f64; 3]> {
        let mut up: [f64; 3] = [0.0; 3];
        let result = unsafe {
            ffi::sk_geodetic_get_local_up(self.geodetic, &mut up[0], &mut up[1], &mut up[2])
        };
        if result != 0 {
            return Err(anyhow!("Failed to get local up"));
        }
        Ok(up)
    }

    pub fn local_west(&self) -> Result<[f64; 3]> {
        let mut west: [f64; 3] = [0.0; 3];
        let result = unsafe {
            ffi::sk_geodetic_get_local_west(self.geodetic, &mut west[0], &mut west[1], &mut west[2])
        };
        if result != 0 {
            return Err(anyhow!("Failed to get local west"));
        }
        Ok(west)
    }

    pub fn altitude_intercepts(
        &self,
        altitude: f64,
        observer: [f64; 3],
        look_vector: [f64; 3],
    ) -> Result<([f64; 3], [f64; 3])> {
        let mut intercept_1: [f64; 3] = [0.0; 3];
        let mut intercept_2: [f64; 3] = [0.0; 3];
        let result = unsafe {
            ffi::sk_geodetic_get_altitude_intercepts(
                self.geodetic,
                altitude,
                observer[0],
                observer[1],
                observer[2],
                look_vector[0],
                look_vector[1],
                look_vector[2],
                &mut intercept_1[0],
                &mut intercept_1[1],
                &mut intercept_1[2],
                &mut intercept_2[0],
                &mut intercept_2[1],
                &mut intercept_2[2],
            )
        };
        if result != 0 {
            return Err(anyhow!("Failed to get altitude intercepts"));
        }
        Ok((intercept_1, intercept_2))
    }

    pub fn from_lat_lon_alt(&mut self, latitude: f64, longitude: f64, altitude: f64) -> Result<()> {
        let result = unsafe {
            ffi::sk_geodetic_from_lat_lon_altitude(self.geodetic, latitude, longitude, altitude)
        };
        if result != 0 {
            return Err(anyhow!("Failed to set latitude, longitude and altitude"));
        }
        Ok(())
    }

    pub fn from_tangent_altitude(
        &mut self,
        altitude: f64,
        observer: [f64; 3],
        boresight: [f64; 3],
    ) -> Result<[f64; 3]> {
        let mut look_vector: [f64; 3] = [0.0; 3];
        let result = unsafe {
            ffi::sk_geodetic_from_tangent_altitude(
                self.geodetic,
                altitude,
                observer[0],
                observer[1],
                observer[2],
                boresight[0],
                boresight[1],
                boresight[2],
                &mut look_vector[0],
                &mut look_vector[1],
                &mut look_vector[2],
            )
        };
        if result != 0 {
            return Err(anyhow!("Failed to set tangent altitude"));
        }
        Ok(look_vector)
    }

    pub fn from_tangent_point(&mut self, observer: [f64; 3], look_vector: [f64; 3]) -> Result<()> {
        let result = unsafe {
            ffi::sk_geodetic_from_tangent_point(
                self.geodetic,
                observer[0],
                observer[1],
                observer[2],
                look_vector[0],
                look_vector[1],
                look_vector[2],
            )
        };
        if result != 0 {
            return Err(anyhow!("Failed to set tangent point"));
        }
        Ok(())
    }

    pub fn from_xyz(&mut self, x: f64, y: f64, z: f64) -> Result<()> {
        let result = unsafe { ffi::sk_geodetic_from_xyz(self.geodetic, x, y, z) };
        if result != 0 {
            return Err(anyhow!("Failed to set XYZ coordinates"));
        }
        Ok(())
    }

    pub fn is_valid(&self) -> Result<bool> {
        let mut is_valid = 0;
        let result = unsafe { ffi::sk_geodetic_is_valid(self.geodetic, &mut is_valid) };
        if result != 0 {
            return Err(anyhow!("Failed to check validity"));
        }
        Ok(is_valid != 0)
    }

    pub fn osculating_spheroid(&mut self) -> Result<(f64, [f64; 3])> {
        let mut radius = 0.0;
        let mut location: [f64; 3] = [0.0; 3];
        let result = unsafe {
            ffi::sk_geodetic_get_osculating_spheroid(
                self.geodetic,
                &mut radius,
                &mut location[0],
                &mut location[1],
                &mut location[2],
            )
        };
        if result != 0 {
            return Err(anyhow!("Failed to get osculating spheroid"));
        }
        Ok((radius, location))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geodetic_construction() {
        let _ = Geodetic::new(6378137.0, 1.0 / 298.257223563).unwrap();
    }

    #[test]
    fn test_location_methods() {
        let mut geodetic = Geodetic::new(6378137.0, 1.0 / 298.257223563).unwrap();

        // Set location using lat/lon/alt
        geodetic
            .from_lat_lon_alt(45.0_f64.to_radians(), 90.0_f64.to_radians(), 100.0)
            .unwrap();

        // Test retrieving values
        let latitude = geodetic.latitude().unwrap();
        let longitude = geodetic.longitude().unwrap();
        let altitude = geodetic.altitude().unwrap();

        assert!((latitude - 45.0_f64.to_radians()).abs() < 1e-10);
        assert!((longitude - 90.0_f64.to_radians()).abs() < 1e-10);
        assert!((altitude - 100.0).abs() < 1e-10);

        // Test location retrieval
        let location = geodetic.location().unwrap();
        assert!(location[0] != 0.0 && location[1] != 0.0 && location[2] != 0.0);

        // Test local vectors
        let up = geodetic.local_up().unwrap();
        let south = geodetic.local_south().unwrap();
        let west = geodetic.local_west().unwrap();

        // Verify vectors are unit vectors
        let up_norm = (up[0] * up[0] + up[1] * up[1] + up[2] * up[2]).sqrt();
        let south_norm = (south[0] * south[0] + south[1] * south[1] + south[2] * south[2]).sqrt();
        let west_norm = (west[0] * west[0] + west[1] * west[1] + west[2] * west[2]).sqrt();

        assert!((up_norm - 1.0).abs() < 1e-10);
        assert!((south_norm - 1.0).abs() < 1e-10);
        assert!((west_norm - 1.0).abs() < 1e-10);

        // Test validity
        assert!(geodetic.is_valid().unwrap());

        // Test osculating spheroid
        let (radius, center) = geodetic.osculating_spheroid().unwrap();
        assert!(radius > 0.0);
        assert!(center[0] != 0.0 || center[1] != 0.0 || center[2] != 0.0);
    }
}
