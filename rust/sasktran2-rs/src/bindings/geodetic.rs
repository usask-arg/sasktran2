use super::prelude::*;
use sasktran2_sys::ffi;

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
}
