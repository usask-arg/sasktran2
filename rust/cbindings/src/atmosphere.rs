use crate::atmosphere_storage::AtmosphereStorage;
use crate::ffi;
use crate::surface::Surface;

pub struct Atmosphere {
    pub atmosphere: *mut ffi::Atmosphere,
    pub storage: AtmosphereStorage,
    pub surface: Surface,
}
impl Atmosphere {
    pub fn new(
        num_wavel: usize,
        num_location: usize,
        num_legendre: usize,
        num_deriv: usize,
        calc_derivatives: bool,
    ) -> Self {
        let storage = AtmosphereStorage::new(num_wavel, num_location, num_legendre, num_deriv);
        let surface = Surface::new(num_wavel, 1); // Assuming 1 for NSTOKES

        // convert calc_derivatives to 0 for false
        // and 1 for true
        let calc_derivatives = if calc_derivatives { 1 } else { 0 };

        let atmosphere = unsafe {
            ffi::sk_atmosphere_create(storage.storage, surface.surface, calc_derivatives)
        };

        Atmosphere {
            atmosphere,
            storage,
            surface,
        }
    }

    pub fn num_wavel(&self) -> usize {
        self.storage.ssa.dim().1
    }

    pub fn num_location(&self) -> usize {
        self.storage.ssa.dim().0
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
        let num_deriv = 2;
        let calc_derivatives = true;

        let mut atmosphere = Atmosphere::new(
            num_wavel,
            num_location,
            num_legendre,
            num_deriv,
            calc_derivatives,
        );

        assert_eq!(atmosphere.storage.ssa.shape(), &[num_location, num_wavel]);

        // check that atmosphere storage is not null
        assert!(!atmosphere.atmosphere.is_null());
    }
}
