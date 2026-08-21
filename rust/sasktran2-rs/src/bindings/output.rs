use std::collections::HashMap;
use std::ffi::CString;

use ndarray::*;
use sasktran2_sys::ffi;

pub struct JvpOutput {
    pub output: *mut ffi::OutputJVP,
    pub radiance: Array3<f64>,
    pub jvp: Array3<f64>,
}

impl JvpOutput {
    pub fn new(num_wavel: usize, num_los: usize, num_stokes: usize) -> Self {
        Self::try_new(num_wavel, num_los, num_stokes).expect("Failed to create native JVP output")
    }

    pub fn try_new(num_wavel: usize, num_los: usize, num_stokes: usize) -> Result<Self, String> {
        let num_radiance = num_wavel
            .checked_mul(num_los)
            .ok_or_else(|| "JVP output dimensions overflow usize".to_string())?;
        let num_values = num_radiance
            .checked_mul(num_stokes)
            .ok_or_else(|| "JVP output dimensions overflow usize".to_string())?;
        i32::try_from(num_values)
            .map_err(|_| "JVP output size exceeds the C API limit".to_string())?;
        let num_radiance_c = i32::try_from(num_radiance)
            .map_err(|_| "JVP output dimensions exceed the C API limit".to_string())?;
        let num_stokes_c = i32::try_from(num_stokes)
            .map_err(|_| "JVP Stokes dimension exceeds the C API limit".to_string())?;

        let mut radiance = Array3::<f64>::zeros((num_wavel, num_los, num_stokes));
        let mut jvp = Array3::<f64>::zeros((num_wavel, num_los, num_stokes));
        let output = unsafe {
            ffi::sk_output_jvp_create(
                radiance.as_mut_ptr(),
                jvp.as_mut_ptr(),
                num_radiance_c,
                num_stokes_c,
            )
        };
        if output.is_null() {
            return Err("Failed to create native JVP output".to_string());
        }
        Ok(Self {
            output,
            radiance,
            jvp,
        })
    }

    pub fn with_derivative_tangent(
        &mut self,
        name: &str,
        tangent: &Array1<f64>,
    ) -> Result<(), String> {
        let name = CString::new(name).map_err(|err| err.to_string())?;
        let nparam = i32::try_from(tangent.len())
            .map_err(|_| "JVP tangent size exceeds the C API limit".to_string())?;
        // The C API accepts a pointer to a dense vector and has no stride
        // argument. Keep the temporary alive until C++ has copied the tangent.
        let tangent = tangent.as_standard_layout();
        let result = unsafe {
            ffi::sk_output_jvp_assign_derivative_tangent(
                self.output,
                name.as_ptr(),
                tangent.as_ptr(),
                nparam,
            )
        };
        if result == 0 {
            Ok(())
        } else {
            Err(format!("Failed to register JVP tangent: {result}"))
        }
    }

    pub fn with_surface_tangent(
        &mut self,
        name: &str,
        tangent: &Array1<f64>,
    ) -> Result<(), String> {
        let name = CString::new(name).map_err(|err| err.to_string())?;
        let nparam = i32::try_from(tangent.len())
            .map_err(|_| "Surface JVP tangent size exceeds the C API limit".to_string())?;
        let tangent = tangent.as_standard_layout();
        let result = unsafe {
            ffi::sk_output_jvp_assign_surface_tangent(
                self.output,
                name.as_ptr(),
                tangent.as_ptr(),
                nparam,
            )
        };
        if result == 0 {
            Ok(())
        } else {
            Err(format!("Failed to register surface JVP tangent: {result}"))
        }
    }
}

impl Drop for JvpOutput {
    fn drop(&mut self) {
        unsafe { ffi::sk_output_jvp_destroy(self.output) }
    }
}

pub struct VjpOutput {
    pub output: *mut ffi::OutputVJP,
    pub radiance: Array3<f64>,
    pub derivative_gradients: HashMap<String, Array1<f64>>,
    pub surface_gradients: HashMap<String, Array1<f64>>,
    _cotangent: Array3<f64>,
}

impl VjpOutput {
    pub fn new(cotangent: &Array3<f64>) -> Self {
        Self::try_new(cotangent).expect("Failed to create native VJP output")
    }

    pub fn try_new(cotangent: &Array3<f64>) -> Result<Self, String> {
        let (num_wavel, num_los, num_stokes) = cotangent.dim();
        let num_radiance = num_wavel
            .checked_mul(num_los)
            .ok_or_else(|| "VJP output dimensions overflow usize".to_string())?;
        let num_values = num_radiance
            .checked_mul(num_stokes)
            .ok_or_else(|| "VJP output dimensions overflow usize".to_string())?;
        i32::try_from(num_values)
            .map_err(|_| "VJP output size exceeds the C API limit".to_string())?;
        let num_radiance_c = i32::try_from(num_radiance)
            .map_err(|_| "VJP output dimensions exceed the C API limit".to_string())?;
        let num_stokes_c = i32::try_from(num_stokes)
            .map_err(|_| "VJP Stokes dimension exceeds the C API limit".to_string())?;

        // The C++ output holds a non-owning map for the duration of its
        // lifetime and assumes standard ndarray order, so normalize the
        // layout and keep the owned cotangent alongside it.
        let cotangent = cotangent.as_standard_layout().into_owned();
        let mut radiance = Array3::<f64>::zeros(cotangent.raw_dim());
        let output = unsafe {
            ffi::sk_output_vjp_create(
                radiance.as_mut_ptr(),
                cotangent.as_ptr(),
                num_radiance_c,
                num_stokes_c,
            )
        };
        if output.is_null() {
            return Err("Failed to create native VJP output".to_string());
        }
        Ok(Self {
            output,
            radiance,
            derivative_gradients: HashMap::new(),
            surface_gradients: HashMap::new(),
            _cotangent: cotangent,
        })
    }

    pub fn with_derivative_gradient(&mut self, name: &str, size: usize) -> Result<(), String> {
        let size_c = i32::try_from(size)
            .map_err(|_| "VJP gradient size exceeds the C API limit".to_string())?;
        self.derivative_gradients
            .insert(name.to_string(), Array1::zeros(size));
        let gradient = self.derivative_gradients.get_mut(name).unwrap();
        let name_c = CString::new(name).map_err(|err| err.to_string())?;
        let result = unsafe {
            ffi::sk_output_vjp_assign_derivative_gradient(
                self.output,
                name_c.as_ptr(),
                gradient.as_mut_ptr(),
                size_c,
            )
        };
        if result == 0 {
            Ok(())
        } else {
            Err(format!("Failed to register VJP gradient: {result}"))
        }
    }

    pub fn with_surface_gradient(&mut self, name: &str, size: usize) -> Result<(), String> {
        let size_c = i32::try_from(size)
            .map_err(|_| "Surface VJP gradient size exceeds the C API limit".to_string())?;
        self.surface_gradients
            .insert(name.to_string(), Array1::zeros(size));
        let gradient = self.surface_gradients.get_mut(name).unwrap();
        let name_c = CString::new(name).map_err(|err| err.to_string())?;
        let result = unsafe {
            ffi::sk_output_vjp_assign_surface_gradient(
                self.output,
                name_c.as_ptr(),
                gradient.as_mut_ptr(),
                size_c,
            )
        };
        if result == 0 {
            Ok(())
        } else {
            Err(format!("Failed to register surface VJP gradient: {result}"))
        }
    }
}

impl Drop for VjpOutput {
    fn drop(&mut self) {
        unsafe { ffi::sk_output_vjp_destroy(self.output) }
    }
}

///  Wrapper around the C++ Output object
/// This is typically only constructed internally by the Engine, and then used by the user
pub struct Output {
    pub output: *mut ffi::OutputC,
    pub radiance: Array3<f64>,
    num_wavel: usize,
    num_los: usize,
    num_stokes: usize,
    num_flux_obs: usize,
    num_flux_types: usize,
    pub d_radiance: HashMap<String, Array4<f64>>,
    pub d_radiance_surf: HashMap<String, Array3<f64>>,

    pub d_flux: HashMap<String, Array4<f64>>,
    pub d_flux_surf: HashMap<String, Array3<f64>>,

    pub flux: Array3<f64>,
    los_optical_depth_override: Option<Array2<f64>>,
}

impl Output {
    pub fn new(
        num_wavel: usize,
        num_los: usize,
        num_flux_obs: usize,
        num_flux_types: usize,
        num_stokes: usize,
    ) -> Self {
        let mut radiance = Array3::<f64>::zeros((num_wavel, num_los, num_stokes));

        // TODO: Get this from the config, based on if upwelling/downwelling fluxes are requested
        let mut flux = Array3::<f64>::zeros((num_flux_types, num_wavel, num_flux_obs));

        let num_radiance = num_wavel * num_los * num_stokes;
        let num_flux = num_wavel * num_flux_obs * num_flux_types;
        let radiance_ptr = radiance.as_mut_ptr();
        let flux_ptr = flux.as_mut_ptr();
        let output = unsafe {
            ffi::sk_output_create(
                radiance_ptr,
                num_radiance as i32,
                num_stokes as i32,
                flux_ptr,
                num_flux as i32,
            )
        };

        Output {
            output,
            radiance,
            num_wavel,
            num_los,
            num_stokes,
            num_flux_obs,
            num_flux_types,
            d_radiance: HashMap::new(),
            d_radiance_surf: HashMap::new(),
            d_flux: HashMap::new(),
            d_flux_surf: HashMap::new(),
            flux,
            los_optical_depth_override: None,
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

        // And the flux derivative
        let num_flux = (self.num_flux_obs * self.num_wavel * self.num_flux_types) as i32;
        if num_flux > 0 {
            let mut d_flux_internal = Array4::<f64>::zeros((
                num_deriv_output,
                self.num_flux_types,
                self.num_wavel,
                self.num_flux_obs,
            ));

            let result = unsafe {
                ffi::sk_output_assign_flux_derivative_memory(
                    self.output,
                    c_deriv_name.as_ptr(),
                    d_flux_internal.as_mut_ptr(),
                    num_flux,
                    num_deriv_output as i32,
                )
            };

            self.d_flux.insert(deriv_name.to_string(), d_flux_internal);

            if result != 0 {
                panic!("Error assigning flux derivative memory");
            }
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

        let nflux = (self.num_wavel * self.num_flux_obs * self.num_flux_types) as i32;
        if nflux > 0 {
            let mut d_flux_internal =
                Array3::<f64>::zeros((self.num_flux_types, self.num_wavel, self.num_flux_obs));
            let d_flux_ptr = d_flux_internal.as_mut_ptr();

            let result = unsafe {
                ffi::sk_output_assign_surface_flux_derivative_memory(
                    self.output,
                    c_deriv_name.as_ptr(),
                    d_flux_ptr,
                    nflux,
                )
            };

            self.d_flux_surf
                .insert(deriv_name.to_string(), d_flux_internal);

            if result != 0 {
                panic!("Error assigning surface flux derivative memory");
            }
        }

        self
    }

    pub fn set_los_optical_depth(&mut self, values: Array2<f64>) -> Result<(), String> {
        if values.dim() != (self.num_wavel, self.num_los) {
            return Err(format!(
                "LOS optical-depth shape {:?} does not match ({}, {})",
                values.dim(),
                self.num_wavel,
                self.num_los
            ));
        }
        self.los_optical_depth_override = Some(values);
        Ok(())
    }

    pub fn los_optical_depth(&self) -> Array2<f64> {
        if let Some(values) = &self.los_optical_depth_override {
            return values.clone();
        }
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
        let _output = Output::new(10, 10, 0, 2, 3);
    }

    #[test]
    fn test_output_with_derivative() {
        let mut output = Output::new(10, 10, 0, 2, 3);
        output.with_derivative("test_deriv", 5);
        assert!(output.d_radiance.contains_key("test_deriv"));
        assert_eq!(output.d_radiance["test_deriv"].shape(), &[5, 10, 10, 3]);
    }

    #[test]
    fn test_output_with_surface_derivative() {
        let mut output = Output::new(10, 10, 0, 2, 3);
        output.with_surface_derivative("test_surf_deriv");
        assert!(output.d_radiance_surf.contains_key("test_surf_deriv"));
        assert_eq!(
            output.d_radiance_surf["test_surf_deriv"].shape(),
            &[10, 10, 3]
        );
    }

    #[test]
    fn test_output_dimensions() {
        let output = Output::new(5, 8, 0, 2, 2);
        assert_eq!(output.radiance.shape(), &[5, 8, 2]);
    }

    #[test]
    fn los_optical_depth_override_validates_and_owns_its_shape() {
        let mut output = Output::new(2, 3, 0, 0, 1);
        assert!(output.set_los_optical_depth(Array2::zeros((3, 2))).is_err());

        let expected = Array2::from_shape_vec((2, 3), (0..6).map(f64::from).collect()).unwrap();
        output.set_los_optical_depth(expected.clone()).unwrap();

        assert_eq!(output.los_optical_depth(), expected);
    }

    #[test]
    fn native_product_output_constructors_reject_invalid_descriptors() {
        let mut radiance = [0.0];
        let mut product = [0.0];

        unsafe {
            assert!(
                ffi::sk_output_jvp_create(std::ptr::null_mut(), product.as_mut_ptr(), 1, 1,)
                    .is_null()
            );
            assert!(
                ffi::sk_output_jvp_create(radiance.as_mut_ptr(), product.as_mut_ptr(), -1, 1,)
                    .is_null()
            );
            assert!(
                ffi::sk_output_jvp_create(radiance.as_mut_ptr(), product.as_mut_ptr(), 1, 2,)
                    .is_null()
            );
            assert!(
                ffi::sk_output_vjp_create(radiance.as_mut_ptr(), product.as_ptr(), i32::MAX, 3,)
                    .is_null()
            );
        }

        assert!(JvpOutput::try_new(1, 1, 2).is_err());
        assert!(JvpOutput::try_new(i32::MAX as usize, 1, 3).is_err());
        assert!(VjpOutput::try_new(&Array3::zeros((1, 1, 2))).is_err());
    }

    #[test]
    fn native_product_registration_rejects_null_inputs() {
        let mut radiance = [0.0];
        let mut product = [0.0];
        let name = CString::new("parameter").unwrap();

        unsafe {
            let jvp = ffi::sk_output_jvp_create(radiance.as_mut_ptr(), product.as_mut_ptr(), 1, 1);
            assert!(!jvp.is_null());
            assert_eq!(
                ffi::sk_output_jvp_assign_derivative_tangent(
                    jvp,
                    std::ptr::null(),
                    product.as_ptr(),
                    1,
                ),
                -1
            );
            assert_eq!(
                ffi::sk_output_jvp_assign_derivative_tangent(
                    jvp,
                    name.as_ptr(),
                    std::ptr::null(),
                    1,
                ),
                -1
            );
            ffi::sk_output_jvp_destroy(jvp);

            let vjp = ffi::sk_output_vjp_create(radiance.as_mut_ptr(), product.as_ptr(), 1, 1);
            assert!(!vjp.is_null());
            assert_eq!(
                ffi::sk_output_vjp_assign_derivative_gradient(
                    vjp,
                    std::ptr::null(),
                    product.as_mut_ptr(),
                    1,
                ),
                -1
            );
            assert_eq!(
                ffi::sk_output_vjp_assign_derivative_gradient(
                    vjp,
                    name.as_ptr(),
                    std::ptr::null_mut(),
                    1,
                ),
                -1
            );
            ffi::sk_output_vjp_destroy(vjp);
        }
    }
}
