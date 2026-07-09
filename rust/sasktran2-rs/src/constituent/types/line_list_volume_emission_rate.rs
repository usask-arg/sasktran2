use crate::atmosphere::types::{LineshapeCoordinate, SpectralGrid};
use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::interpolation::linear::linear_interpolating_matrix;
use crate::optical::line::OpticalLine;
use crate::optical::types::line_absorber::assign_normalized_doppler_line_shape;
use crate::prelude::*;

const O2_MOLECULAR_MASS_G_PER_MOL: f64 = 31.9988;

pub struct LineListVolumeEmissionRate {
    pub photon_ver: Array1<f64>,
    pub altitudes: Array1<f64>,
    pub wavelengths_nm: Array1<f64>,
    pub weights: Array2<f64>,
    interp_mode: crate::interpolation::OutOfBoundsMode,
}

impl LineListVolumeEmissionRate {
    pub fn new(
        altitudes: Array1<f64>,
        photon_ver: Array1<f64>,
        wavelengths_nm: Array1<f64>,
        weights: Array2<f64>,
    ) -> Result<Self> {
        if altitudes.len() != photon_ver.len() {
            return Err(anyhow!(
                "Altitude length ({}) must match photon VER length ({})",
                altitudes.len(),
                photon_ver.len()
            ));
        }

        if weights.nrows() != altitudes.len() {
            return Err(anyhow!(
                "Line-list weight row count ({}) must match altitude length ({})",
                weights.nrows(),
                altitudes.len()
            ));
        }

        if wavelengths_nm.len() != weights.ncols() {
            return Err(anyhow!(
                "Line wavelength length ({}) must match weight column count ({})",
                wavelengths_nm.len(),
                weights.ncols()
            ));
        }

        if wavelengths_nm.is_empty() {
            return Err(anyhow!("Line-list emission requires at least one line"));
        }

        if wavelengths_nm
            .iter()
            .any(|wavelength| !wavelength.is_finite() || *wavelength <= 0.0)
        {
            return Err(anyhow!(
                "Line wavelengths must all be positive finite values"
            ));
        }

        let weights = normalize_weight_rows(weights)?;

        Ok(Self {
            altitudes,
            photon_ver,
            wavelengths_nm,
            weights,
            interp_mode: crate::interpolation::OutOfBoundsMode::Zero,
        })
    }

    pub fn with_interp_mode(mut self, interp_mode: crate::interpolation::OutOfBoundsMode) -> Self {
        self.interp_mode = interp_mode;
        self
    }

    fn has_altitude_independent_weights(&self) -> bool {
        if self.weights.nrows() <= 1 {
            return true;
        }

        let first = self.weights.row(0);
        self.weights.rows().into_iter().all(|row| {
            row.iter()
                .zip(first.iter())
                .all(|(a, b)| (a - b).abs() < 1.0e-14)
        })
    }

    fn add_doppler_broadened_spectrum(
        &self,
        spectral_grid: &SpectralGrid,
        temperature_k: f64,
        line_areas: ArrayView1<'_, f64>,
        spectrum: &mut Array1<f64>,
    ) {
        let wavenumber_cminv = spectral_grid.central_wavenumber_cminv();

        for (line_idx, &line_area) in line_areas.iter().enumerate() {
            if line_area == 0.0 {
                continue;
            }

            let line_center_cminv = 1.0e7 / self.wavelengths_nm[line_idx];
            let doppler_width_cminv = OpticalLine::doppler_width_cminv(
                line_center_cminv,
                temperature_k,
                O2_MOLECULAR_MASS_G_PER_MOL,
            );

            assign_normalized_doppler_line_shape(
                wavenumber_cminv.as_slice().unwrap(),
                line_center_cminv,
                doppler_width_cminv,
                line_area,
                spectrum.as_slice_mut().unwrap(),
            );
        }

        if spectral_grid.coordinate() == LineshapeCoordinate::WavelengthNm {
            Zip::from(spectrum)
                .and(spectral_grid.central_wavelengths_nm())
                .for_each(|value, &wavelength_nm| {
                    *value *= 1.0e7 / wavelength_nm.powi(2);
                });
        }
    }

    fn broadened_spectral_emission(
        &self,
        spectral_grid: &SpectralGrid,
        temperature_k: ArrayView1<'_, f64>,
        line_areas_by_alt: ArrayView2<'_, f64>,
    ) -> Array2<f64> {
        let mut spectral_emission = Array2::zeros((
            line_areas_by_alt.nrows(),
            spectral_grid.central_wavenumber_cminv().len(),
        ));

        Zip::from(spectral_emission.rows_mut())
            .and(temperature_k)
            .and(line_areas_by_alt.rows())
            .for_each(|mut spectrum, &temperature, line_areas| {
                let mut spectrum_owned = Array1::zeros(spectrum.len());
                self.add_doppler_broadened_spectrum(
                    spectral_grid,
                    temperature,
                    line_areas,
                    &mut spectrum_owned,
                );
                spectrum.assign(&spectrum_owned);
            });

        spectral_emission
    }

    fn broadened_spectral_weights(
        &self,
        spectral_grid: &SpectralGrid,
        temperature_k: ArrayView1<'_, f64>,
        model_altitudes_m: ArrayView1<'_, f64>,
    ) -> Result<Array2<f64>> {
        let weights_by_alt = if self.has_altitude_independent_weights() {
            Array2::from_shape_fn(
                (model_altitudes_m.len(), self.weights.ncols()),
                |(_, line_idx)| self.weights[[0, line_idx]],
            )
        } else {
            if self.altitudes.len() != model_altitudes_m.len()
                || self
                    .altitudes
                    .iter()
                    .zip(model_altitudes_m.iter())
                    .any(|(lhs, rhs)| (lhs - rhs).abs() > 1.0e-9)
            {
                return Err(anyhow!(
                    "LineListVolumeEmissionRate derivatives with altitude-dependent line weights require the emission altitude grid to match the model altitude grid"
                ));
            }

            self.weights.clone()
        };

        Ok(self.broadened_spectral_emission(spectral_grid, temperature_k, weights_by_alt.view()))
    }
}

impl Constituent for LineListVolumeEmissionRate {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;
        let (inputs, outputs) = storage.split_inputs_outputs();

        let altitudes_m = inputs.altitude_m();
        let temperature_k = inputs.temperature_k().ok_or(anyhow!(
            "Temperature must be set to add Doppler-broadened line-list volume emission"
        ))?;
        let spectral_grid = inputs.spectral_grid().ok_or(anyhow!(
            "Spectral grid must be set to add line-list volume emission rate constituent"
        ))?;

        let interp_matrix =
            linear_interpolating_matrix(&self.altitudes, &altitudes_m, self.interp_mode);
        let mut line_ver = self.weights.clone();
        Zip::from(line_ver.rows_mut())
            .and(&self.photon_ver)
            .for_each(|mut row, &ver| {
                row *= ver;
            });
        let interp_line_ver = interp_matrix.dot(&line_ver);

        let spectral_emission =
            self.broadened_spectral_emission(spectral_grid, temperature_k, interp_line_ver.view());

        let mut emission = outputs.mut_view().emission_source;

        Zip::from(&mut emission).and(&spectral_emission).for_each(
            |emission_val, &spectral_emission| {
                *emission_val += spectral_emission / FOUR_PI;
            },
        );

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
        let temperature_k = inputs.temperature_k().ok_or(anyhow!(
            "Temperature must be set to register Doppler-broadened line-list volume emission derivatives"
        ))?;
        let spectral_grid = inputs.spectral_grid().ok_or(anyhow!(
            "Spectral grid must be set to add line-list volume emission rate constituent"
        ))?;

        let spectral_weights_by_alt =
            self.broadened_spectral_weights(spectral_grid, temperature_k, altitudes_m)?;

        let interp_matrix = if self.has_altitude_independent_weights() {
            linear_interpolating_matrix(&self.altitudes, &altitudes_m, self.interp_mode)
        } else {
            if self.altitudes.len() != altitudes_m.len()
                || self
                    .altitudes
                    .iter()
                    .zip(altitudes_m.iter())
                    .any(|(lhs, rhs)| (lhs - rhs).abs() > 1.0e-9)
            {
                return Err(anyhow!(
                    "LineListVolumeEmissionRate derivatives with altitude-dependent line weights require the emission altitude grid to match the model altitude grid"
                ));
            }

            Array2::eye(self.altitudes.len())
        };

        let full_name = "wf_".to_owned() + constituent_name + "_photon_ver";
        let mut deriv = derivative_generator.get_derivative_mapping(&full_name);

        deriv.set_interp_dim(format!("{}_altitude", constituent_name).as_str());
        deriv.set_assign_name(full_name.as_str());
        deriv.set_interpolator(&interp_matrix);

        let deriv_view = deriv.mut_view();
        let mut d_emission = deriv_view.d_emission;

        Zip::from(&mut d_emission)
            .and(&spectral_weights_by_alt)
            .for_each(|emission_val, &spectral_weight| {
                *emission_val += spectral_weight / FOUR_PI;
            });

        Ok(())
    }
}

fn normalize_weight_rows(mut weights: Array2<f64>) -> Result<Array2<f64>> {
    for mut row in weights.rows_mut() {
        let weight_sum: f64 = row.iter().sum();
        if !weight_sum.is_finite() || weight_sum <= 0.0 {
            return Err(anyhow!(
                "Each line-list weight row must sum to a positive finite value"
            ));
        }

        for weight in row.iter() {
            if !weight.is_finite() || *weight < 0.0 {
                return Err(anyhow!(
                    "Line-list weights must be non-negative finite values"
                ));
            }
        }

        row /= weight_sum;
    }

    Ok(weights)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trapz(y: ArrayView1<'_, f64>, x: ArrayView1<'_, f64>) -> f64 {
        y.iter()
            .zip(y.iter().skip(1))
            .zip(x.iter().zip(x.iter().skip(1)))
            .map(|((&y0, &y1), (&x0, &x1))| 0.5 * (y0 + y1) * (x1 - x0))
            .sum()
    }

    #[test]
    fn line_list_weights_are_normalized_on_construction() {
        let ver = LineListVolumeEmissionRate::new(
            array![0.0, 1.0],
            array![10.0, 20.0],
            array![760.0, 761.0],
            array![[2.0, 6.0], [1.0, 3.0]],
        )
        .unwrap();

        assert!((ver.weights.row(0).sum() - 1.0).abs() < 1.0e-12);
        assert!((ver.weights.row(1).sum() - 1.0).abs() < 1.0e-12);
        assert!((ver.weights[[0, 0]] - 0.25).abs() < 1.0e-12);
        assert!((ver.weights[[0, 1]] - 0.75).abs() < 1.0e-12);
    }

    #[test]
    fn doppler_broadened_profile_conserves_line_area_on_wavelength_grid() {
        let ver = LineListVolumeEmissionRate::new(
            array![0.0],
            array![10.0],
            array![760.0],
            array![[1.0]],
        )
        .unwrap();

        let wavelengths_nm = Array1::from_iter((0..20001).map(|i| 759.0 + i as f64 * 0.0001));
        let spectral_grid = SpectralGrid::from_wavelengths_nm(&wavelengths_nm).unwrap();
        let mut spectrum = Array1::zeros(wavelengths_nm.len());

        ver.add_doppler_broadened_spectrum(
            &spectral_grid,
            200.0,
            array![1.0].view(),
            &mut spectrum,
        );

        let area = trapz(spectrum.view(), wavelengths_nm.view());
        assert!((area - 1.0).abs() < 1.0e-4);
        assert!(spectrum.iter().filter(|&&value| value > 0.0).count() > 1);
    }

    #[test]
    fn altitude_dependent_weights_change_spectral_shape() {
        let ver = LineListVolumeEmissionRate::new(
            array![0.0, 1.0],
            array![4.0 * std::f64::consts::PI, 8.0 * std::f64::consts::PI],
            array![760.0, 761.0],
            array![[1.0, 0.0], [0.0, 1.0]],
        )
        .unwrap();

        let wavelengths_nm = Array1::from_iter((0..20001).map(|i| 759.5 + i as f64 * 0.0001));
        let spectral_grid = SpectralGrid::from_wavelengths_nm(&wavelengths_nm).unwrap();

        let mut line_ver = ver.weights.clone();
        Zip::from(line_ver.rows_mut())
            .and(&ver.photon_ver)
            .for_each(|mut row, &photon_ver| row *= photon_ver);
        let spectral_emission = ver.broadened_spectral_emission(
            &spectral_grid,
            array![200.0, 200.0].view(),
            line_ver.view(),
        );

        let row0_peak_idx = spectral_emission
            .row(0)
            .iter()
            .enumerate()
            .max_by(|(_, lhs), (_, rhs)| lhs.partial_cmp(rhs).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();
        let row1_peak_idx = spectral_emission
            .row(1)
            .iter()
            .enumerate()
            .max_by(|(_, lhs), (_, rhs)| lhs.partial_cmp(rhs).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();

        assert!((wavelengths_nm[row0_peak_idx] - 760.0).abs() < 0.001);
        assert!((wavelengths_nm[row1_peak_idx] - 761.0).abs() < 0.001);
        assert!(
            (trapz(spectral_emission.row(0), wavelengths_nm.view()) - 4.0 * std::f64::consts::PI)
                .abs()
                < 1.0e-3
        );
        assert!(
            (trapz(spectral_emission.row(1), wavelengths_nm.view()) - 8.0 * std::f64::consts::PI)
                .abs()
                < 1.0e-3
        );
    }
}
