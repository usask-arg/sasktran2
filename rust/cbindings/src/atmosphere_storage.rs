use crate::deriv_mapping::DerivativeMapping;
use crate::ffi;
use crate::prelude::*;
use ndarray::*;

pub struct AtmosphereStorage {
    pub storage: *mut ffi::AtmosphereStorage,
    pub ssa: Array2<f64>,
    pub total_extinction: Array2<f64>,
    pub emission_source: Array2<f64>,
    pub f: Array2<f64>,
    pub leg_coeff: Array3<f64>,
    pub d_leg_coeff: Array4<f64>,
    pub d_f: Array3<f64>,
    pub solar_irradiance: Array1<f64>,
}

impl AtmosphereStorage {
    pub fn new(
        num_wavel: usize,
        num_location: usize,
        num_legendre: usize,
        num_deriv: usize,
        stokes: Stokes,
    ) -> Self {
        let mut ssa = Array2::<f64>::zeros((num_location, num_wavel).f());
        let mut total_extinction = Array2::<f64>::zeros((num_location, num_wavel).f());
        let mut emission_source = Array2::<f64>::zeros((num_location, num_wavel).f());
        let mut f = Array2::<f64>::zeros((num_location, num_wavel).f());
        let mut leg_coeff = Array3::<f64>::zeros(
            (
                num_legendre * stokes.num_legendre(),
                num_location,
                num_wavel,
            )
                .f(),
        );

        let mut d_leg_coeff = Array4::<f64>::zeros(
            (
                num_legendre * stokes.num_legendre(),
                num_location,
                num_wavel,
                num_deriv,
            )
                .f(),
        );
        let mut d_f = Array3::<f64>::zeros((num_location, num_wavel, num_deriv).f());
        let mut solar_irradiance = Array1::<f64>::zeros(num_wavel);

        solar_irradiance.fill(1.0);

        AtmosphereStorage {
            storage: unsafe {
                ffi::sk_atmosphere_storage_create(
                    num_location as i32,
                    num_wavel as i32,
                    (num_legendre * stokes.num_legendre()) as i32,
                    stokes.num_stokes() as i32, // nstokes
                    num_deriv as i32,
                    ssa.as_mut_ptr(),
                    total_extinction.as_mut_ptr(),
                    emission_source.as_mut_ptr(),
                    f.as_mut_ptr(),
                    leg_coeff.as_mut_ptr(),
                    d_leg_coeff.as_mut_ptr(),
                    d_f.as_mut_ptr(),
                    solar_irradiance.as_mut_ptr(),
                )
            },
            ssa,
            total_extinction,
            emission_source,
            f,
            leg_coeff,
            d_leg_coeff,
            d_f,
            solar_irradiance,
        }
    }

    pub fn get_derivative_mapping(&self, name: &str) -> Result<DerivativeMapping, String> {
        let mut mapping: *mut ffi::DerivativeMapping = std::ptr::null_mut();
        let c_name = std::ffi::CString::new(name).unwrap();

        let result = unsafe {
            ffi::sk_atmosphere_storage_get_derivative_mapping(
                self.storage,
                c_name.as_ptr(),
                &mut mapping,
            )
        };

        if result != 0 {
            return Err("Failed to get derivative mapping".to_string());
        }

        Ok(DerivativeMapping::new(mapping))
    }
}

impl Drop for AtmosphereStorage {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_atmosphere_storage_destroy(self.storage);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atmosphere_storage() {
        let num_wavel = 10;
        let num_location = 5;
        let num_legendre = 3;
        let num_deriv = 2;

        let storage = AtmosphereStorage::new(
            num_wavel,
            num_location,
            num_legendre,
            num_deriv,
            Stokes::Stokes1,
        );

        let mapping = storage.get_derivative_mapping("wf_test").unwrap();

        let d_ssa = mapping.d_ssa();

        assert_eq!(storage.ssa.shape(), &[num_location, num_wavel]);
        assert_eq!(storage.total_extinction.shape(), &[num_location, num_wavel]);
        assert_eq!(storage.emission_source.shape(), &[num_location, num_wavel]);
        assert_eq!(storage.f.shape(), &[num_location, num_wavel]);
        assert_eq!(
            storage.leg_coeff.shape(),
            &[num_legendre, num_location, num_wavel]
        );
        assert_eq!(
            storage.d_leg_coeff.shape(),
            &[num_legendre, num_location, num_wavel, num_deriv]
        );
        assert_eq!(storage.d_f.shape(), &[num_location, num_wavel, num_deriv]);
        assert_eq!(storage.solar_irradiance.shape(), &[num_wavel]);
    }
}
