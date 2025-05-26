use super::atmosphere_storage::AtmosphereStorage;
use super::prelude::*;
use super::surface::Surface;
use anyhow::{Result, anyhow};
use sasktran2_sys::ffi;

/// Wrapping of the c++ sasktran2 Atmosphere object, which consists of the AtmosphereStorage and
/// Surface objects
pub struct Atmosphere {
    pub atmosphere: *mut ffi::Atmosphere,
    pub storage: AtmosphereStorage,
    pub surface: Surface,
    nstokes: usize,
}
impl Atmosphere {
    /// Creates a new Atmosphere object with the appropriate memory allocated
    pub fn new(
        num_wavel: usize,
        num_location: usize,
        num_legendre: usize,
        calc_derivatives: bool,
        calc_emission_derivatives: bool,
        stokes: Stokes,
    ) -> Self {
        let storage = AtmosphereStorage::new(num_wavel, num_location, num_legendre, stokes);
        let surface = Surface::new(num_wavel, stokes.num_stokes());

        // convert calc_derivatives to 0 for false
        // and 1 for true
        let calc_derivatives = if calc_derivatives { 1 } else { 0 };
        let calc_emission_derivatives = if calc_emission_derivatives { 1 } else { 0 };

        let atmosphere = unsafe {
            ffi::sk_atmosphere_create(
                storage.storage,
                surface.surface,
                calc_derivatives,
                calc_emission_derivatives,
            )
        };

        Atmosphere {
            atmosphere,
            storage,
            surface,
            nstokes: stokes.num_stokes(),
        }
    }

    /// Number of wavelengths in the atmosphere
    pub fn num_wavel(&self) -> usize {
        self.storage.ssa.dim().1
    }

    /// Number of geometry grid points in the atmosphere
    pub fn num_location(&self) -> usize {
        self.storage.ssa.dim().0
    }

    /// Applies delta_m scaling to the atmosphere, applied on the c++ side
    pub fn apply_delta_m_scaling(&mut self, order: usize) -> Result<()> {
        let result =
            unsafe { ffi::sk_atmosphere_apply_delta_m_scaling(self.atmosphere, order as i32) };
        if result != 0 {
            return Err(anyhow!("Error applying delta m scaling: {}", result));
        }
        Ok(())
    }

    /// Number of stokes parameters
    pub fn num_stokes(&self) -> usize {
        self.nstokes
    }
}
impl Drop for Atmosphere {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_atmosphere_destroy(self.atmosphere);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atmosphere() {
        let num_wavel = 10;
        let num_location = 5;
        let num_legendre = 3;
        let calc_derivatives = true;

        let atmosphere = Atmosphere::new(
            num_wavel,
            num_location,
            num_legendre,
            calc_derivatives,
            false,
            Stokes::Stokes1,
        );

        assert_eq!(atmosphere.storage.ssa.shape(), &[num_location, num_wavel]);

        // check that atmosphere storage is not null
        assert!(!atmosphere.atmosphere.is_null());
    }

    #[test]
    fn test_apply_delta_m_scaling() {
        let num_wavel = 10;
        let num_location = 5;
        let num_legendre = 3;
        let calc_derivatives = true;

        let mut atmosphere = Atmosphere::new(
            num_wavel,
            num_location,
            num_legendre,
            calc_derivatives,
            false,
            Stokes::Stokes1,
        );

        // apply delta m scaling
        let result = atmosphere.apply_delta_m_scaling(2);
        assert!(result.is_ok());
    }
}
