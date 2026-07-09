use crate::atmosphere::types::{LineshapeCoordinate, SpectralGrid};
use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::interpolation::linear::linear_interpolating_matrix;
use crate::optical::line::OpticalLine;
use crate::optical::types::line_absorber::assign_normalized_doppler_line_shape;
use crate::prelude::*;
use rebasis::basis::Delta;
use rebasis::grid::{Grid, MappingMatrix, mapping_matrix};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MonochromaticLineShape {
    Delta,
    Doppler {
        emitter_molecular_weight_g_per_mol: f64,
    },
}

pub struct MonochromaticVolumeEmissionRate {
    pub ver: Array1<f64>,
    pub altitudes: Array1<f64>,
    wavelength_nm: f64,
    interp_mode: crate::interpolation::OutOfBoundsMode,
    line_shape: MonochromaticLineShape,
}

impl MonochromaticVolumeEmissionRate {
    pub fn new(altitudes: Array1<f64>, ver: Array1<f64>, wavelength_nm: f64) -> Self {
        MonochromaticVolumeEmissionRate {
            altitudes,
            ver,
            wavelength_nm,
            interp_mode: crate::interpolation::OutOfBoundsMode::Zero,
            line_shape: MonochromaticLineShape::Delta,
        }
    }

    pub fn with_interp_mode(mut self, interp_mode: crate::interpolation::OutOfBoundsMode) -> Self {
        self.interp_mode = interp_mode;
        self
    }

    pub fn with_line_shape(mut self, line_shape: MonochromaticLineShape) -> Self {
        self.line_shape = line_shape;
        self
    }

    fn mapping_values(&self, spectral_grid: &SpectralGrid) -> MappingMatrix {
        let internal_grid = match spectral_grid.coordinate() {
            crate::atmosphere::types::LineshapeCoordinate::WavelengthNm => {
                let delta = Delta::new(self.wavelength_nm);
                Grid::new(vec![rebasis::basis::BasisType::Delta(delta)])
            }
            crate::atmosphere::types::LineshapeCoordinate::WavenumberCminv => {
                let delta = Delta::new(1e7 / self.wavelength_nm);
                Grid::new(vec![rebasis::basis::BasisType::Delta(delta)])
            }
        };

        // shape (spectral,)
        mapping_matrix(&internal_grid, spectral_grid.basis_grid())
    }

    fn doppler_spectral_weights(
        &self,
        spectral_grid: &SpectralGrid,
        temperature_k: ArrayView1<'_, f64>,
        emitter_molecular_weight_g_per_mol: f64,
    ) -> Array2<f64> {
        let line_center_cminv = 1.0e7 / self.wavelength_nm;
        let wavenumber_cminv = spectral_grid.central_wavenumber_cminv();
        let mut spectral_weights = Array2::zeros((temperature_k.len(), wavenumber_cminv.len()));

        Zip::from(spectral_weights.rows_mut())
            .and(temperature_k)
            .for_each(|mut spectrum, &temperature| {
                let doppler_width_cminv = OpticalLine::doppler_width_cminv(
                    line_center_cminv,
                    temperature,
                    emitter_molecular_weight_g_per_mol,
                );

                let mut spectrum_owned = Array1::zeros(spectrum.len());
                assign_normalized_doppler_line_shape(
                    wavenumber_cminv.as_slice().unwrap(),
                    line_center_cminv,
                    doppler_width_cminv,
                    1.0,
                    spectrum_owned.as_slice_mut().unwrap(),
                );

                if spectral_grid.coordinate() == LineshapeCoordinate::WavelengthNm {
                    Zip::from(&mut spectrum_owned)
                        .and(spectral_grid.central_wavelengths_nm())
                        .for_each(|value, &wavelength_nm| {
                            *value *= 1.0e7 / wavelength_nm.powi(2);
                        });
                }

                spectrum.assign(&spectrum_owned);
            });

        spectral_weights
    }

    fn spectral_weights(
        &self,
        spectral_grid: &SpectralGrid,
        temperature_k: Option<ArrayView1<'_, f64>>,
        num_altitudes: usize,
    ) -> Result<Array2<f64>> {
        match self.line_shape {
            MonochromaticLineShape::Delta => {
                let mapping_vals = self.mapping_values(spectral_grid);
                let mapping_vals = mapping_vals.matrix.column(0).to_owned();
                Ok(Array2::from_shape_fn(
                    (num_altitudes, mapping_vals.len()),
                    |(_, wavelength_idx)| mapping_vals[wavelength_idx],
                ))
            }
            MonochromaticLineShape::Doppler {
                emitter_molecular_weight_g_per_mol,
            } => {
                let temperature_k = temperature_k.ok_or(anyhow!(
                    "Temperature must be set to add Doppler-broadened volume emission"
                ))?;
                Ok(self.doppler_spectral_weights(
                    spectral_grid,
                    temperature_k,
                    emitter_molecular_weight_g_per_mol,
                ))
            }
        }
    }
}

impl Constituent for MonochromaticVolumeEmissionRate {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;
        let (inputs, outputs) = storage.split_inputs_outputs();

        let altitudes_m = inputs.altitude_m();

        let interp_matrix =
            linear_interpolating_matrix(&self.altitudes, &altitudes_m, self.interp_mode);
        let interp_ver = interp_matrix.dot(&self.ver);

        let mut emission = outputs.mut_view().emission_source;

        let spectral_grid = inputs.spectral_grid().ok_or(anyhow!(
            "Spectral grid must be set to add volume emission rate constituent"
        ))?;
        let spectral_weights =
            self.spectral_weights(spectral_grid, inputs.temperature_k(), altitudes_m.len())?;

        Zip::from(&mut emission)
            .and(&spectral_weights)
            .and_broadcast(&interp_ver.insert_axis(Axis(1)))
            .for_each(|emission_val, &spectral_weight, &ver_val| {
                *emission_val += spectral_weight * ver_val / FOUR_PI;
            });

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

        let interp_matrix =
            linear_interpolating_matrix(&self.altitudes, &altitudes_m, self.interp_mode);

        let full_name = "wf_".to_owned() + constituent_name + "_ver";
        let mut deriv = derivative_generator.get_derivative_mapping(&full_name);

        deriv.set_interp_dim(format!("{}_altitude", constituent_name).as_str());
        // TODO: set interpolator?
        deriv.set_assign_name(full_name.as_str());
        deriv.set_interpolator(&interp_matrix);
        let deriv_view = deriv.mut_view();
        let mut d_emission = deriv_view.d_emission;

        let spectral_grid = inputs.spectral_grid().ok_or(anyhow!(
            "Spectral grid must be set to add volume emission rate constituent"
        ))?;
        let spectral_weights =
            self.spectral_weights(spectral_grid, inputs.temperature_k(), altitudes_m.len())?;

        Zip::from(&mut d_emission).and(&spectral_weights).for_each(
            |emission_val, &spectral_weight| {
                *emission_val += spectral_weight / FOUR_PI;
            },
        );

        Ok(())
    }
}
