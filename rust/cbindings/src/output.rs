use crate::ffi;
use ndarray::*;

pub struct Output {
    pub output: *mut ffi::OutputC,
    pub radiance: Array3<f64>,
}

impl Output {
    pub fn new(num_wavel: usize, num_los: usize, num_stokes: usize) -> Self {
        let mut radiance = Array3::<f64>::zeros((num_wavel, num_los, num_stokes));

        let num_radiance = num_wavel * num_los * num_stokes;
        let radiance_ptr = radiance.as_mut_ptr();
        let output = unsafe { ffi::sk_output_create(radiance_ptr, num_radiance as i32) };

        Output {
            output: output,
            radiance: radiance,
        }
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
        let mut output = Output::new(10, 10, 3);
    }
}
