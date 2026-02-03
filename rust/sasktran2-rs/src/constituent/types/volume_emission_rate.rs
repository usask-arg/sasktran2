use crate::atmosphere::types::SpectralGrid;
use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::interpolation::grid1d::Grid1D;
use crate::interpolation::linear::{Interp1Weights, linear_interpolating_matrix};
use crate::prelude::*;
use rebasis::basis::Delta;
use rebasis::grid::{Grid, mapping_matrix, MappingMatrix};

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

    fn mapping_values(&self, spectral_grid: &SpectralGrid) -> MappingMatrix {
        let internal_grid = match spectral_grid.coordinate() {
            crate::atmosphere::types::LineshapeCoordinate::WavelengthNm => {
                let delta = Delta::new( self.wavelength_nm );
                Grid::new( vec![ rebasis::basis::BasisType::Delta( delta ) ] )
            }
            crate::atmosphere::types::LineshapeCoordinate::WavenumberCminv => {
                let delta = Delta::new( 1e7 / self.wavelength_nm );
                Grid::new( vec![ rebasis::basis::BasisType::Delta( delta ) ] )
            }
        };

        // shape (spectral,)
        mapping_matrix(&internal_grid, spectral_grid.basis_grid())
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


        let mapping_vals = self.mapping_values(inputs.spectral_grid().ok_or(anyhow!(
            "Spectral grid must be set to add volume emission rate constituent"
        ))?);

        let mapping_vals = mapping_vals.matrix.column(0).to_owned();

        for (i, &map_val) in mapping_vals.iter().enumerate() {
            Zip::from(emission.column_mut(i))
                .and(&interp_ver)
                .for_each(|emission_val, &ver_val| {
                    *emission_val += map_val * ver_val / FOUR_PI;
                });
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

        let wavelength_nm = inputs.spectral_grid().ok_or(anyhow!(
            "Spectral grid must be set to add volume emission rate constituent"
        ))?.central_wavelengths_nm();

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

        Ok(())
    }
}
