use crate::ffi;
use crate::prelude::*;
use ndarray::*;

pub struct DerivativeMapping {
    pub mapping: *mut ffi::DerivativeMapping,
}

impl DerivativeMapping {
    pub fn new(mapping: *mut ffi::DerivativeMapping) -> Self {
        DerivativeMapping { mapping: mapping }
    }

    pub fn d_ssa(&self) -> ArrayViewMut2<f64> {
        let mut d_ssa: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_deriv_mapping_get_d_ssa(self.mapping, &mut d_ssa);
            ArrayViewMut2::from_shape_ptr((self.num_location(), self.num_wavel()).f(), d_ssa)
        }
    }

    pub fn d_extinction(&self) -> ArrayViewMut2<f64> {
        let mut d_extinction: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_deriv_mapping_get_d_extinction(self.mapping, &mut d_extinction);
            ArrayViewMut2::from_shape_ptr(
                (self.num_location(), self.num_legendre()).f(),
                d_extinction,
            )
        }
    }

    pub fn scat_factor(&self) -> ArrayViewMut2<f64> {
        let mut scat_factor: *mut f64 = std::ptr::null_mut();
        unsafe {
            ffi::sk_deriv_mapping_get_scat_factor(self.mapping, &mut scat_factor);
            ArrayViewMut2::from_shape_ptr(
                (self.num_location(), self.num_legendre()).f(),
                scat_factor,
            )
        }
    }

    pub fn d_leg_coeff(&self) -> ArrayViewMut3<f64> {
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
}

impl Drop for DerivativeMapping {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_deriv_mapping_destroy(self.mapping);
        }
    }
}
