use std::collections::HashMap;
use std::ffi::CString;

use ndarray::*;
use sasktran2_sys::ffi;

pub struct Output {
    pub output: *mut ffi::OutputC,
    pub radiance: Array3<f64>,
    num_wavel: usize,
    num_los: usize,
    num_stokes: usize,
    pub d_radiance: HashMap<String, Array4<f64>>,
    pub d_radiance_surf: HashMap<String, Array3<f64>>,
}

impl Output {
    pub fn new(num_wavel: usize, num_los: usize, num_stokes: usize) -> Self {
        let mut radiance = Array3::<f64>::zeros((num_wavel, num_los, num_stokes));

        let num_radiance = num_wavel * num_los * num_stokes;
        let radiance_ptr = radiance.as_mut_ptr();
        let output =
            unsafe { ffi::sk_output_create(radiance_ptr, num_radiance as i32, num_stokes as i32) };

        Output {
            output,
            radiance,
            num_wavel,
            num_los,
            num_stokes,
            d_radiance: HashMap::new(),
            d_radiance_surf: HashMap::new(),
        }
    }

    pub fn with_derivative(&mut self, deriv_name: &str, num_deriv_output: usize) -> &mut Self {
        let mut d_radiance_internal = Array4::<f64>::zeros((
            num_deriv_output,
            self.num_wavel,
            self.num_los,
            self.num_stokes,
        ));
        let d_radiance_ptr = d_radiance_internal.as_mut_ptr();

        let nrad = (self.num_wavel * self.num_los) as i32;

        let c_deriv_name = CString::new(deriv_name).unwrap();

        let result = unsafe {
            ffi::sk_output_assign_derivative_memory(
                self.output,
                c_deriv_name.as_ptr(),
                d_radiance_ptr,
                nrad,
                self.num_stokes as i32,
                num_deriv_output as i32,
            )
        };

        self.d_radiance
            .insert(deriv_name.to_string(), d_radiance_internal);

        if result != 0 {
            panic!("Error assigning derivative memory");
        }

        self
    }

    pub fn with_surface_derivative(&mut self, deriv_name: &str) -> &mut Self {
        let mut d_radiance_internal =
            Array3::<f64>::zeros((self.num_wavel, self.num_los, self.num_stokes));
        let d_radiance_ptr = d_radiance_internal.as_mut_ptr();

        let nrad = (self.num_wavel * self.num_los) as i32;

        let c_deriv_name = CString::new(deriv_name).unwrap();

        let result = unsafe {
            ffi::sk_output_assign_surface_derivative_memory(
                self.output,
                c_deriv_name.as_ptr(),
                d_radiance_ptr,
                nrad,
                self.num_stokes as i32,
            )
        };

        self.d_radiance_surf
            .insert(deriv_name.to_string(), d_radiance_internal);

        if result != 0 {
            panic!("Error assigning surface derivative memory");
        }

        self
    }

    pub fn los_optical_depth(&self) -> Array2<f64> {
        let mut internal: *mut f64 = std::ptr::null_mut();
        let internal_view = unsafe {
            ffi::sk_output_get_los_optical_depth(self.output, &mut internal);
            ArrayView2::from_shape_ptr((self.num_wavel, self.num_los).f(), internal)
        };

        let mut output = Array2::<f64>::zeros((self.num_wavel, self.num_los));
        output.assign(&internal_view);

        output
    }
}

impl Drop for Output {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_output_destroy(self.output);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_output() {
        let _output = Output::new(10, 10, 3);
    }
}
