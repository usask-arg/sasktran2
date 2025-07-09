use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::prelude::*;

const PLANCK: f64 = 6.62607015e-34; // J s = kg m^2 / s
const SPEED_OF_LIGHT: f64 = 299792458.0; // m/s
const K_BOLTZMANN: f64 = 1.380649e-23; // J/K = m^2 kg / (s^2 K)

fn planck_blackbody_radiance(temperature_k: f64, wavelength_nm: f64) -> f64 {
    let wavelength_m: f64 = wavelength_nm * 1.0e-9;
    {
        (2.0 * PLANCK * SPEED_OF_LIGHT.powi(2) / wavelength_m.powi(5))
            / ((PLANCK * SPEED_OF_LIGHT / (wavelength_m * K_BOLTZMANN * temperature_k)).exp() - 1.0)
            * 1.0e-9
    }
}

fn d_planck_blackbody_radiance_d_temperature(temperature_k: f64, wavelength_nm: f64) -> f64 {
    let wavelength_m: f64 = wavelength_nm * 1.0e-9;
    let exponent: f64 =
        (PLANCK * SPEED_OF_LIGHT / (wavelength_m * K_BOLTZMANN * temperature_k)).exp();
    {
        (2.0 * PLANCK.powi(2) * SPEED_OF_LIGHT.powi(3) * exponent)
            / (wavelength_m.powi(6)
                * K_BOLTZMANN
                * temperature_k.powi(2)
                * (exponent - 1.0).powi(2))
            * 1.0e-9
    }
}

pub struct ThermalEmission;

impl Default for ThermalEmission {
    fn default() -> Self {
        ThermalEmission::new()
    }
}

impl ThermalEmission {
    pub fn new() -> Self {
        ThermalEmission
    }
}

impl Constituent for ThermalEmission {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (inputs, outputs) = storage.split_inputs_outputs();
        let mut outputs = outputs.mut_view();

        let wavelengths_nm = inputs.wavelengths_nm().unwrap();
        let temperature_k = inputs.temperature_k().unwrap();

        let emission_source_array = &mut outputs.emission_source;

        let thread_pool = crate::threading::thread_pool()?;

        thread_pool.install(|| {
            Zip::from(emission_source_array.columns_mut())
                .and(wavelengths_nm)
                .par_for_each(|emission_col, wav| {
                    Zip::from(emission_col)
                        .and(temperature_k)
                        .for_each(|emission, tem| {
                            *emission += planck_blackbody_radiance(*tem, *wav);
                        })
                });
        });

        Ok(())
    }

    fn register_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()> {
        let (inputs, _, derivative_generator) = storage.split_inputs_outputs_deriv();

        let wavelengths_nm = inputs.wavelengths_nm().unwrap();
        let temperature_k = inputs.temperature_k().unwrap();

        let thread_pool = crate::threading::thread_pool()?;

        let full_name = "wf_".to_owned() + constituent_name + "_temperature_k";
        let mut deriv = derivative_generator.get_derivative_mapping(&full_name);
        let mut deriv_view = deriv.mut_view();

        // TODO: do extinction and ssa derivatives need to be explicitly set to zero?
        thread_pool.install(|| {
            Zip::from(deriv_view.d_emission.columns_mut())
                .and(wavelengths_nm)
                .par_for_each(|d_emis_col, wav| {
                    Zip::from(d_emis_col)
                        .and(temperature_k)
                        .for_each(|d_emis, tem| {
                            *d_emis += d_planck_blackbody_radiance_d_temperature(*tem, *wav);
                        })
                });
        });

        deriv.set_interp_dim("altitude");
        // TODO: set interpolator?
        deriv.set_assign_name("wf_temperature_k");

        Ok(())
    }
}
