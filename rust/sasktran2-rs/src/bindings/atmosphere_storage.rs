use super::deriv_mapping::DerivativeMapping;
use super::prelude::*;
use ndarray::*;
use sasktran2_sys::ffi;

/// Binding of the sasktran2 AtmosphereStorage object
/// Memory for the storage is allocated on the rust side and views passed into the C++ side
/// Note that the storage arrays are in Fortran order to match the C++ side
pub struct AtmosphereStorage {
    pub storage: *mut ffi::AtmosphereStorage,
    pub ssa: Array2<f64>,
    pub total_extinction: Array2<f64>,
    pub emission_source: Array2<f64>,
    pub leg_coeff: Array3<f64>,
    pub solar_irradiance: Array1<f64>,
}

impl AtmosphereStorage {
    /// Creates a new AtmosphereStorage object with the appropriate memory allocated
    pub fn new(num_wavel: usize, num_location: usize, num_legendre: usize, stokes: Stokes) -> Self {
        let mut ssa = Array2::<f64>::zeros((num_location, num_wavel).f());
        let mut total_extinction = Array2::<f64>::zeros((num_location, num_wavel).f());
        let mut emission_source = Array2::<f64>::zeros((num_location, num_wavel).f());
        let mut leg_coeff = Array3::<f64>::zeros(
            (
                num_legendre * stokes.num_legendre(),
                num_location,
                num_wavel,
            )
                .f(),
        );

        let mut solar_irradiance = Array1::<f64>::zeros(num_wavel);

        solar_irradiance.fill(1.0);

        AtmosphereStorage {
            storage: unsafe {
                ffi::sk_atmosphere_storage_create(
                    num_location as i32,
                    num_wavel as i32,
                    (num_legendre * stokes.num_legendre()) as i32,
                    stokes.num_stokes() as i32, // nstokes
                    ssa.as_mut_ptr(),
                    total_extinction.as_mut_ptr(),
                    emission_source.as_mut_ptr(),
                    leg_coeff.as_mut_ptr(),
                    solar_irradiance.as_mut_ptr(),
                )
            },
            ssa,
            total_extinction,
            emission_source,
            leg_coeff,
            solar_irradiance,
        }
    }

    /// Returns back the derivative mapping for name if it exists, if not creates a new one
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

    /// Returns back a vector of all the derivative mapping names
    pub fn derivative_mapping_names(&self) -> Result<Vec<String>, String> {
        let mut names: Vec<String> = Vec::new();

        let mut num_mappings: i32 = 0;
        let result = unsafe {
            ffi::sk_atmosphere_storage_get_num_derivative_mappings(self.storage, &mut num_mappings)
        };

        if result != 0 {
            return Err("Failed to get number of derivative mappings".to_string());
        }

        for i in 0..num_mappings {
            let mut name_ptr: *const std::ffi::c_char = std::ptr::null();
            let result = unsafe {
                ffi::sk_atmosphere_storage_get_derivative_mapping_name(
                    self.storage,
                    i,
                    &mut name_ptr,
                )
            };
            if result != 0 {
                return Err("Failed to get derivative mapping name".to_string());
            }

            let c_str = unsafe { std::ffi::CStr::from_ptr(name_ptr) };
            let name = c_str.to_string_lossy().into_owned();
            names.push(name);
        }

        Ok(names)
    }

    /// Called after ssa is initialized to scattering extinction, and legendre need to be normalized.
    /// This divides the legendre coefficients by the ssa and then normalizes the ssa by the total extinction
    pub fn normalize_by_extinctions(&mut self) {
        // Start by dividing the leg_coeff by the ssa
        // Since these are fortran ordered we iterate over the last axis first
        let thread_pool = crate::threading::thread_pool().unwrap();

        thread_pool.install(|| {
            Zip::from(self.leg_coeff.axis_iter_mut(Axis(2)))
                .and(self.ssa.axis_iter(Axis(1)))
                .par_for_each(|mut leg_coeff, ssa| {
                    Zip::from(ssa)
                        .and(leg_coeff.axis_iter_mut(Axis(1)))
                        .for_each(|ssa, mut leg_coeff| {
                            let v = match ssa {
                                0.0 => 1.0,
                                _ => *ssa,
                            };
                            leg_coeff.mapv_inplace(|l| l / v);
                        });
                });

            // Now we need to divide the ssa by the total extinction
            Zip::from(self.ssa.axis_iter_mut(Axis(1)))
                .and(self.total_extinction.axis_iter(Axis(1)))
                .par_for_each(|ssa, total_extinction| {
                    Zip::from(ssa)
                        .and(total_extinction)
                        .for_each(|ssa, total_extinction| {
                            let v = *total_extinction;
                            *ssa /= v;

                            if *ssa > 1.0 {
                                *ssa = 1.0;
                            }
                        });
                });
        });
    }

    /// Finalizes the scattering derivatives for the storage, done on the c++ side
    pub fn finalize_scattering_derivatives(&mut self) {
        unsafe {
            ffi::sk_atmosphere_storage_finalize_scattering_derivatives(self.storage);
        }
    }

    /// Sets the storage to zero, done on the c++ side
    pub fn set_zero(&mut self) {
        unsafe {
            ffi::sk_atmosphere_storage_set_zero(self.storage);
        }
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
    fn test_atmosphere_storage_new() {
        let num_wavel = 10;
        let num_location = 5;
        let num_legendre = 3;

        let storage =
            AtmosphereStorage::new(num_wavel, num_location, num_legendre, Stokes::Stokes1);

        assert_eq!(storage.ssa.shape(), &[num_location, num_wavel]);
        assert_eq!(storage.total_extinction.shape(), &[num_location, num_wavel]);
        assert_eq!(storage.emission_source.shape(), &[num_location, num_wavel]);
        assert_eq!(
            storage.leg_coeff.shape(),
            &[
                num_legendre * Stokes::Stokes1.num_legendre(),
                num_location,
                num_wavel
            ]
        );
        assert_eq!(storage.solar_irradiance.shape(), &[num_wavel]);
    }

    #[test]
    fn test_atmosphere_storage_set_zero() {
        let num_wavel = 10;
        let num_location = 5;
        let num_legendre = 3;

        let mut storage =
            AtmosphereStorage::new(num_wavel, num_location, num_legendre, Stokes::Stokes1);

        storage.ssa.fill(1.0);
        storage.total_extinction.fill(1.0);
        storage.emission_source.fill(1.0);
        storage.leg_coeff.fill(1.0);

        storage.set_zero();

        assert_eq!(storage.ssa.sum(), 0.0);
        assert_eq!(storage.total_extinction.sum(), 0.0);
        assert_eq!(storage.emission_source.sum(), 0.0);
        assert_eq!(storage.leg_coeff.sum(), 0.0);
    }

    #[test]
    fn test_atmosphere_storage() {
        let num_wavel = 10;
        let num_location = 5;
        let num_legendre = 3;

        let storage =
            AtmosphereStorage::new(num_wavel, num_location, num_legendre, Stokes::Stokes1);

        let _mapping = storage.get_derivative_mapping("wf_test").unwrap();

        println!("mapping_names {:?}", storage.derivative_mapping_names());

        assert_eq!(storage.ssa.shape(), &[num_location, num_wavel]);
        assert_eq!(storage.total_extinction.shape(), &[num_location, num_wavel]);
        assert_eq!(storage.emission_source.shape(), &[num_location, num_wavel]);
        assert_eq!(
            storage.leg_coeff.shape(),
            &[num_legendre, num_location, num_wavel]
        );
        assert_eq!(storage.solar_irradiance.shape(), &[num_wavel]);
    }
}
