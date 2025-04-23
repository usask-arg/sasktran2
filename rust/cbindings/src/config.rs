use crate::ffi;
use crate::prelude::*;

#[repr(i32)]
pub enum MultipleScatterSource {
    DiscreteOrdinates = 0,
    SuccessiveOrders = 1,
    TwoStream = 2,
    None = 3,
}

#[repr(i32)]
pub enum SingleScatterSource {
    Exact = 0,
    SolarTable = 1,
    DiscreteOrdinates = 2,
    None = 3,
}

pub struct Config {
    pub config: *mut ffi::Config,
}

impl Config {
    pub fn new() -> Self {
        Config {
            config: unsafe { ffi::sk_config_create() },
        }
    }

    pub fn num_stokes(&self) -> Result<usize> {
        let mut num_stokes = 0i32;
        let error_code = unsafe { ffi::sk_config_get_num_stokes(self.config, &mut num_stokes) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of stokes: error code {}",
                error_code
            ))
        } else {
            Ok(num_stokes as usize)
        }
    }

    pub fn with_num_stokes(self, num_stokes: usize) -> Result<Self> {
        let error_code = unsafe { ffi::sk_config_set_num_stokes(self.config, num_stokes as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of stokes: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn with_multiple_scatter_source(self, source: MultipleScatterSource) -> Result<Self> {
        let error_code =
            unsafe { ffi::sk_config_set_multiple_scatter_source(self.config, source as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting multiple scatter source: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn with_single_scatter_source(self, source: SingleScatterSource) -> Result<Self> {
        let error_code =
            unsafe { ffi::sk_config_set_single_scatter_source(self.config, source as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting single scatter source: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn with_num_streams(self, num_streams: usize) -> Result<Self> {
        let error_code = unsafe { ffi::sk_config_set_num_streams(self.config, num_streams as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of streams: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }
}

impl Drop for Config {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_config_destroy(self.config);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config() {
        let mut config = Config::new();
    }
}
