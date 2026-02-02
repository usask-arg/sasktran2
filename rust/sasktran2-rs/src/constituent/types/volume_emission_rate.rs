use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::interpolation::grid1d::Grid1D;
use crate::interpolation::linear::{Interp1Weights, linear_interpolating_matrix};
use crate::prelude::*;

pub struct MonochromaticVolumeEmissionRate {
    pub ver: Array1<f64>,
    pub altitudes: Array1<f64>,
    wavelength_nm: f64,
    interp_mode: crate::interpolation::OutOfBoundsMode,
}

impl MonochromaticVolumeEmissionRate {
    pub fn new(altitudes: Array1<f64>, ver: Array1<f64>, wavelength_nm: f64) -> Self {
        MonochromaticVolumeEmissionRate {
            altitudes,
            ver,
            wavelength_nm,
            interp_mode: crate::interpolation::OutOfBoundsMode::Zero,
        }
    }

    pub fn with_interp_mode(mut self, interp_mode: crate::interpolation::OutOfBoundsMode) -> Self {
        self.interp_mode = interp_mode;
        self
    }
}

impl Constituent for MonochromaticVolumeEmissionRate {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;
        let (inputs, outputs) = storage.split_inputs_outputs();

        let altitudes_m = inputs.altitude_m();

        let interp_matrix = linear_interpolating_matrix(
            &self.altitudes,
            &altitudes_m,
            crate::interpolation::OutOfBoundsMode::Zero,
        );
        let interp_ver = interp_matrix.dot(&self.ver);

        let mut emission = outputs.mut_view().emission_source;
        let emission_wavenumber = 1e7 / self.wavelength_nm;

        if inputs.finite_resolution_mode() {
            let wavenumber_left = inputs.wavenumbers_cminv_left().ok_or(anyhow!(
                "Must be in Finite resolution mode to add volume emission rate constituent"
            ))?;
            let wavenumber_right = inputs.wavenumbers_cminv_right().ok_or(anyhow!(
                "Must be in Finite resolution mode to add volume emission rate constituent"
            ))?;

            Zip::from(emission.axis_iter_mut(Axis(1)))
                .and(wavenumber_left)
                .and(wavenumber_right)
                .for_each(|mut emission_slice, &wn_left, &wn_right| {
                    let min_wvnum = wn_left.min(wn_right);
                    let max_wvnum = wn_left.max(wn_right);

                    if emission_wavenumber >= min_wvnum && emission_wavenumber < max_wvnum {
                        // multiply by width in nm
                        let width = 1e7 / min_wvnum - 1e7 / max_wvnum;
                        Zip::from(&mut emission_slice).and(&interp_ver).for_each(
                            |emission_val, &ver_val| {
                                *emission_val += ver_val / width / FOUR_PI;
                            },
                        );
                    }
                });
        } else {
            // set up a sorted wavelength grid
            let wavelength_nm = inputs.wavelengths_nm().ok_or(anyhow!(
                "Wavelengths must be set to add volume emission rate constituent"
            ))?;

            let mut sorted_wvlen = wavelength_nm.as_slice().unwrap().to_owned();
            sorted_wvlen.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let wavelen_grid = Grid1D::new(sorted_wvlen.into());

            // Find the interpolation index
            let interp_weights = wavelen_grid.interp1_weights(
                self.wavelength_nm,
                crate::interpolation::OutOfBoundsMode::Zero,
            );

            for interp_weight in interp_weights {
                let left_width = match interp_weight.0 {
                    0 => 0.0,
                    _ => (wavelen_grid.x[interp_weight.0] - wavelen_grid.x[interp_weight.0 - 1])
                        .abs(),
                };
                let right_width = match interp_weight.0 + 1 {
                    x if x >= wavelen_grid.x.len() => 0.0,
                    _ => (wavelen_grid.x[interp_weight.0 + 1] - wavelen_grid.x[interp_weight.0])
                        .abs(),
                };

                let width = 0.5 * (left_width + right_width);

                Zip::from(emission.column_mut(interp_weight.0))
                    .and(&interp_ver)
                    .for_each(|emission_val, &ver_val| {
                        *emission_val += interp_weight.1 * ver_val / width / FOUR_PI;
                    });
            }
        }

        Ok(())
    }

    fn register_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()> {
        const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;

        let (inputs, _, derivative_generator) = storage.split_inputs_outputs_deriv();

        let altitudes_m = inputs.altitude_m();

        let interp_matrix = linear_interpolating_matrix(
            &self.altitudes,
            &altitudes_m,
            crate::interpolation::OutOfBoundsMode::Zero,
        );

        let emission_wavenumber = 1e7 / self.wavelength_nm;

        let full_name = "wf_".to_owned() + constituent_name + "_ver";
        let mut deriv = derivative_generator.get_derivative_mapping(&full_name);

        deriv.set_interp_dim(format!("{}_altitude", constituent_name).as_str());
        // TODO: set interpolator?
        deriv.set_assign_name(full_name.as_str());
        deriv.set_interpolator(&interp_matrix);
        let deriv_view = deriv.mut_view();
        let mut d_emission = deriv_view.d_emission;

        if inputs.finite_resolution_mode() {
            let wavenumber_left = inputs.wavenumbers_cminv_left().ok_or(anyhow!(
                "Must be in Finite resolution mode to add volume emission rate constituent"
            ))?;
            let wavenumber_right = inputs.wavenumbers_cminv_right().ok_or(anyhow!(
                "Must be in Finite resolution mode to add volume emission rate constituent"
            ))?;

            Zip::from(d_emission.axis_iter_mut(Axis(1)))
                .and(wavenumber_left)
                .and(wavenumber_right)
                .for_each(|mut emission_slice, &wn_left, &wn_right| {
                    let min_wvnum = wn_left.min(wn_right);
                    let max_wvnum = wn_left.max(wn_right);

                    if emission_wavenumber >= min_wvnum && emission_wavenumber < max_wvnum {
                        // multiply by width in nm
                        let width = 1e7 / min_wvnum - 1e7 / max_wvnum;
                        Zip::from(&mut emission_slice).for_each(|emission_val| {
                            *emission_val += 1.0 / width / FOUR_PI;
                        });
                    }
                });
        } else {
            // set up a sorted wavelength grid
            let wavelength_nm = inputs.wavelengths_nm().ok_or(anyhow!(
                "Wavelengths must be set to add volume emission rate constituent"
            ))?;

            let mut sorted_wvlen = wavelength_nm.as_slice().unwrap().to_owned();
            sorted_wvlen.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let wavelen_grid = Grid1D::new(sorted_wvlen.into());

            // Find the interpolation index
            let interp_weights = wavelen_grid.interp1_weights(
                self.wavelength_nm,
                crate::interpolation::OutOfBoundsMode::Zero,
            );

            for interp_weight in interp_weights {
                let left_width = match interp_weight.0 {
                    0 => 0.0,
                    _ => (wavelen_grid.x[interp_weight.0] - wavelen_grid.x[interp_weight.0 - 1])
                        .abs(),
                };
                let right_width = match interp_weight.0 + 1 {
                    x if x >= wavelen_grid.x.len() => 0.0,
                    _ => (wavelen_grid.x[interp_weight.0 + 1] - wavelen_grid.x[interp_weight.0])
                        .abs(),
                };

                let width = 0.5 * (left_width + right_width);

                Zip::from(d_emission.column_mut(interp_weight.0)).for_each(|emission_val| {
                    *emission_val += interp_weight.1 / width / FOUR_PI;
                });
            }
        }

        Ok(())
    }
}
