use ndarray::*;
use sasktran2_sys::ffi;

/// Wrapper around the c++ DerivativeMapping object, just for internal use
pub struct DerivativeMapping {
    pub mapping: *mut ffi::DerivativeMapping,
}

impl DerivativeMapping {
    pub fn new(mapping: *mut ffi::DerivativeMapping) -> Self {
        DerivativeMapping { mapping }
    }

    pub fn d_ssa(&self) -> ArrayViewMut2<'_, f64> {
        let mut d_ssa: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_deriv_mapping_get_d_ssa(self.mapping, &mut d_ssa);
            ArrayViewMut2::from_shape_ptr((self.num_location(), self.num_wavel()).f(), d_ssa)
        }
    }

    pub fn d_extinction(&self) -> ArrayViewMut2<'_, f64> {
        let mut d_extinction: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_deriv_mapping_get_d_extinction(self.mapping, &mut d_extinction);
            ArrayViewMut2::from_shape_ptr((self.num_location(), self.num_wavel()).f(), d_extinction)
        }
    }

    pub fn d_emission(&self) -> ArrayViewMut2<'_, f64> {
        let mut d_emission: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_deriv_mapping_get_d_emission(self.mapping, &mut d_emission);
            ArrayViewMut2::from_shape_ptr((self.num_location(), self.num_wavel()).f(), d_emission)
        }
    }

    pub fn scat_factor(&self) -> ArrayViewMut2<'_, f64> {
        let mut scat_factor: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_deriv_mapping_get_scat_factor(self.mapping, &mut scat_factor);
            ArrayViewMut2::from_shape_ptr((self.num_location(), self.num_wavel()).f(), scat_factor)
        }
    }

    pub fn d_leg_coeff(&self) -> ArrayViewMut3<'_, f64> {
        let mut d_leg_coeff: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_deriv_mapping_get_d_legendre(self.mapping, &mut d_leg_coeff);
            ArrayViewMut3::from_shape_ptr(
                (self.num_legendre(), self.num_location(), self.num_wavel()).f(),
                d_leg_coeff,
            )
        }
    }

    fn num_location(&self) -> usize {
        let mut num_location: i32 = 0;
        unsafe {
            ffi::sk_deriv_mapping_get_num_location(self.mapping, &mut num_location);
            num_location as usize
        }
    }

    fn num_wavel(&self) -> usize {
        let mut num_wavel: i32 = 0;
        unsafe {
            ffi::sk_deriv_mapping_get_num_wavel(self.mapping, &mut num_wavel);
            num_wavel as usize
        }
    }

    fn num_legendre(&self) -> usize {
        let mut num_legendre: i32 = 0;
        unsafe {
            ffi::sk_deriv_mapping_get_num_legendre(self.mapping, &mut num_legendre);
            num_legendre as usize
        }
    }

    pub fn num_output(&self) -> usize {
        let mut num_output: i32 = 0;
        unsafe {
            ffi::sk_deriv_mapping_get_num_output(self.mapping, &mut num_output);
            num_output as usize
        }
    }

    pub fn get_assign_name(&self) -> String {
        let mut name: *const libc::c_char = std::ptr::null();
        unsafe {
            ffi::sk_deriv_mapping_get_assign_name(self.mapping, &mut name);
            let c_str = std::ffi::CStr::from_ptr(name);
            c_str.to_string_lossy().into_owned()
        }
    }

    pub fn set_assign_name(&mut self, name: &str) {
        let c_name = std::ffi::CString::new(name).unwrap();
        unsafe {
            ffi::sk_deriv_mapping_set_assign_name(self.mapping, c_name.as_ptr());
        }
    }

    pub fn set_interp_dim(&mut self, name: &str) {
        let c_name = std::ffi::CString::new(name).unwrap();
        unsafe {
            ffi::sk_deriv_mapping_set_interp_dim(self.mapping, c_name.as_ptr());
        }
    }

    pub fn get_interp_dim(&self) -> String {
        let mut interp_dim: *const libc::c_char = std::ptr::null();
        unsafe {
            ffi::sk_deriv_mapping_get_interp_dim(self.mapping, &mut interp_dim);
            let c_str = std::ffi::CStr::from_ptr(interp_dim);
            c_str.to_string_lossy().into_owned()
        }
    }

    pub fn set_interpolator(&mut self, interpolator: &mut Array2<f64>) {
        let (dim1, dim2) = interpolator.dim();

        // Ensure that the interpolator is fortran ordered
        let mut f_interpolator = Array2::zeros((dim1, dim2).f());
        f_interpolator.assign(interpolator);

        // The deriv mapping copies the data, so it's okay to use the local array
        unsafe {
            ffi::sk_deriv_mapping_set_interpolator(
                self.mapping,
                f_interpolator.as_mut_ptr(),
                dim1 as i32,
                dim2 as i32,
            );
        }
    }

    pub fn get_interpolator(&self) -> ArrayView2<'_, f64> {
        let mut interpolator: *mut f64 = std::ptr::null_mut();
        let mut dim1: i32 = 0;
        let mut dim2: i32 = 0;
        unsafe {
            ffi::sk_deriv_mapping_get_interpolator(
                self.mapping,
                &mut interpolator,
                &mut dim1,
                &mut dim2,
            );
            ArrayView2::from_shape_ptr((dim1 as usize, dim2 as usize).f(), interpolator)
        }
    }

    pub fn set_log_radiance_space(&mut self, log_radiance_space: bool) {
        let log_radiance_space = if log_radiance_space { 1 } else { 0 };
        unsafe {
            ffi::sk_deriv_mapping_set_log_radiance_space(self.mapping, log_radiance_space);
        }
    }
}

impl Drop for DerivativeMapping {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_deriv_mapping_destroy(self.mapping);
        }
    }
}

pub struct SurfaceDerivativeMapping {
    pub mapping: *mut ffi::SurfaceDerivativeMapping,
}

impl Drop for SurfaceDerivativeMapping {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_surface_deriv_mapping_destroy(self.mapping);
        }
    }
}

impl SurfaceDerivativeMapping {
    pub fn new(mapping: *mut ffi::SurfaceDerivativeMapping) -> Self {
        SurfaceDerivativeMapping { mapping }
    }

    pub fn set_zero(&mut self) {
        unsafe {
            ffi::sk_surface_deriv_mapping_set_zero(self.mapping);
        }
    }

    pub fn get_interp_dim(&self) -> String {
        let mut interp_dim: *const libc::c_char = std::ptr::null();
        unsafe {
            ffi::sk_surface_deriv_mapping_get_interp_dim(self.mapping, &mut interp_dim);
            let c_str = std::ffi::CStr::from_ptr(interp_dim);
            c_str.to_string_lossy().into_owned()
        }
    }

    pub fn set_interp_dim(&mut self, name: &str) {
        let c_name = std::ffi::CString::new(name).unwrap();
        unsafe {
            ffi::sk_surface_deriv_mapping_set_interp_dim(self.mapping, c_name.as_ptr());
        }
    }

    pub fn set_interpolator(&mut self, interpolator: &mut Array2<f64>) {
        let (dim1, dim2) = interpolator.dim();

        // Ensure that the interpolator is fortran ordered
        let mut f_interpolator = Array2::zeros((dim1, dim2).f());
        f_interpolator.assign(interpolator);
        // The deriv mapping copies the data, so it's okay to use the local array
        unsafe {
            ffi::sk_surface_deriv_mapping_set_interpolator(
                self.mapping,
                f_interpolator.as_mut_ptr(),
                dim1 as i32,
                dim2 as i32,
            );
        }
    }

    pub fn get_interpolator(&self) -> ArrayView2<'_, f64> {
        let mut interpolator: *mut f64 = std::ptr::null_mut();
        let mut dim1: i32 = 0;
        let mut dim2: i32 = 0;
        unsafe {
            ffi::sk_surface_deriv_mapping_get_interpolator(
                self.mapping,
                &mut interpolator,
                &mut dim1,
                &mut dim2,
            );
            ArrayView2::from_shape_ptr((dim1 as usize, dim2 as usize).f(), interpolator)
        }
    }

    fn num_brdf_args(&self) -> usize {
        let mut num_brdf_args: i32 = 0;
        unsafe {
            ffi::sk_surface_deriv_mapping_get_num_brdf_args(self.mapping, &mut num_brdf_args);
            num_brdf_args as usize
        }
    }

    fn num_wavel(&self) -> usize {
        let mut num_wavel: i32 = 0;
        unsafe {
            ffi::sk_surface_deriv_mapping_get_num_wavel(self.mapping, &mut num_wavel);
            num_wavel as usize
        }
    }

    pub fn d_emission(&self) -> ArrayViewMut1<'_, f64> {
        let mut d_emission: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_surface_deriv_mapping_get_d_emission(self.mapping, &mut d_emission);
            ArrayViewMut1::from_shape_ptr((self.num_wavel()).f(), d_emission)
        }
    }

    pub fn d_brdf(&self) -> ArrayViewMut2<'_, f64> {
        let mut d_brdf: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_surface_deriv_mapping_get_d_brdf(self.mapping, &mut d_brdf);
            ArrayViewMut2::from_shape_ptr((self.num_wavel(), self.num_brdf_args()).f(), d_brdf)
        }
    }
}
