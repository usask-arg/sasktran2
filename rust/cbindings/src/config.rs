use crate::ffi;
use crate::prelude::*;

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
            Err(anyhow!("Error getting number of stokes: error code {}", error_code))
        } else {
            Ok(num_stokes as usize)
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
