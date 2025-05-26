use super::prelude::*;
use sasktran2_sys::ffi;

/// Wrapper around the c++ BRDF object, should not be used directly, instead prefer
/// the Specific implementations of the BRDF such as Lambertian, MODIS, etc.
#[allow(clippy::upper_case_acronyms)]
struct BRDF {
    pub brdf: *mut ffi::BRDF,
}

impl Drop for BRDF {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_brdf_destroy(self.brdf);
        }
    }
}

/// Trait indicating that the object is a c++ BRDF object and can be passed to the
/// c++ code
pub trait IsCBRDF {
    /// The raw ffi object
    fn as_cbrdf(&self) -> *mut ffi::BRDF;

    /// The number of arguments for the BRDF, needed to allocate memory
    /// e.g. a Lambertian BRDF has 1 argument, the albedo, other BRDFs may have more
    fn num_args(&self) -> Result<usize> {
        let mut num_args = 0i32;
        let error_code = unsafe { ffi::sk_brdf_get_num_args(self.as_cbrdf(), &mut num_args) };

        if error_code != 0 {
            return Err(anyhow!(
                "Failed to get number of arguments for BRDF: {}",
                error_code
            ));
        }

        Ok(num_args as usize)
    }
}

/// A lambertian BRDF, implemented on the c++ side
pub struct Lambertian {
    internal: BRDF,
}

impl Lambertian {
    pub fn new(nstokes: usize) -> Self {
        let brdf = unsafe { ffi::sk_brdf_create_lambertian(nstokes as i32) };
        if brdf.is_null() {
            panic!("Failed to create Lambertian BRDF");
        }
        Lambertian {
            internal: BRDF { brdf },
        }
    }
}

impl IsCBRDF for Lambertian {
    fn as_cbrdf(&self) -> *mut ffi::BRDF {
        self.internal.brdf
    }
}

/// A Kokhanovsky BRDF, implemented on the c++ side
pub struct SnowKokhanovsky {
    internal: BRDF,
}

impl SnowKokhanovsky {
    pub fn new(nstokes: usize) -> Self {
        let brdf = unsafe { ffi::sk_brdf_create_kokhanovsky(nstokes as i32) };
        if brdf.is_null() {
            panic!("Failed to create Kokhanovsky BRDF");
        }
        SnowKokhanovsky {
            internal: BRDF { brdf },
        }
    }
}
impl IsCBRDF for SnowKokhanovsky {
    fn as_cbrdf(&self) -> *mut ffi::BRDF {
        self.internal.brdf
    }
}

/// A MODIS BRDF, implemented on the c++ side
pub struct MODIS {
    internal: BRDF,
}

impl MODIS {
    pub fn new(nstokes: usize) -> Self {
        let brdf = unsafe { ffi::sk_brdf_create_modis(nstokes as i32) };
        if brdf.is_null() {
            panic!("Failed to create MODIS BRDF");
        }
        MODIS {
            internal: BRDF { brdf },
        }
    }
}

impl IsCBRDF for MODIS {
    fn as_cbrdf(&self) -> *mut ffi::BRDF {
        self.internal.brdf
    }
}
