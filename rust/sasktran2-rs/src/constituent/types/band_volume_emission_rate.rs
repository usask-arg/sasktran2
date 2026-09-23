//! An independently specified band photon VER, distributed among lines at the
//! current atmospheric temperature. Population conversion belongs upstream.
use crate::atmosphere::types::LineshapeCoordinate;
use crate::atmosphere::*;
use crate::bindings::config::SpectralGridMode;
use crate::constituent::traits::Constituent;
use crate::interpolation::{OutOfBoundsMode, linear::linear_interpolating_matrix};
use crate::math::errorfunctions::optimized::SQRT_PI;
use crate::optical::line::shape::{LineShape, LineShapeDirection, assign_with_derivative};
use crate::optical::line::{AdjustedLineParameters, OpticalLine, OpticalLineDB};
use crate::optical::types::line_absorber::assign_normalized_doppler_line_shape;
use crate::photchem::emission::*;
use crate::prelude::*;

const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;
const O2_MASS: f64 = 31.9988;

/// One band per constituent gives a separate VER retrieval parameter, even
/// when several vibrational bands overlap spectrally.
#[derive(Clone)]
pub struct BandVolumeEmissionRate {
    pub altitudes: Array1<f64>,
    pub photon_ver: Array1<f64>,
    pub band: EmissionBand,
    pub line_weight_model: AEmissionLineWeightModel,
    interp_mode: OutOfBoundsMode,
}

impl BandVolumeEmissionRate {
    pub fn new(
        altitudes: Array1<f64>,
        photon_ver: Array1<f64>,
        band: EmissionBand,
        line_weight_model: AEmissionLineWeightModel,
    ) -> Result<Self> {
        anyhow::ensure!(
            !altitudes.is_empty(),
            "Emission altitude grid cannot be empty"
        );
        anyhow::ensure!(
            altitudes.iter().all(|a| a.is_finite())
                && altitudes
                    .iter()
                    .zip(altitudes.iter().skip(1))
                    .all(|(a, b)| b > a),
            "Emission altitudes must be finite and strictly increasing"
        );
        validate_photon_ver(photon_ver.view(), altitudes.len())?;
        anyhow::ensure!(!band.lines.is_empty(), "Emission band must contain lines");
        let first = &band.lines[0];
        anyhow::ensure!(
            band.lines.iter().all(|line| line.upper_vibrational_state
                == first.upper_vibrational_state
                && line.lower_vibrational_state == first.lower_vibrational_state),
            "A band VER must describe a single vibrational transition"
        );
        Ok(Self {
            altitudes,
            photon_ver,
            band,
            line_weight_model,
            interp_mode: OutOfBoundsMode::Zero,
        })
    }

    pub fn with_interp_mode(mut self, mode: OutOfBoundsMode) -> Self {
        self.interp_mode = mode;
        self
    }

    pub fn line_weights(&self, temperature: ArrayView1<'_, f64>) -> Result<Array2<f64>> {
        oxygen_a_band_lte_line_weights(&self.band, temperature, self.line_weight_model)
    }

    fn interpolation(&self, inputs: &impl StorageInputs) -> Array2<f64> {
        linear_interpolating_matrix(&self.altitudes, &inputs.altitude_m(), self.interp_mode)
    }

    fn spectrum(
        &self,
        inputs: &impl StorageInputs,
        derivative: bool,
        ver: Option<ArrayView1<'_, f64>>,
    ) -> Result<(Array2<f64>, Option<Array2<f64>>)> {
        let temperature = inputs
            .temperature_k()
            .ok_or_else(|| anyhow!("Temperature must be set for band emission"))?;
        let (mut weights, mut d_weights) = if derivative {
            let (weights, derivatives) =
                oxygen_a_band_lte_line_weights_with_temperature_derivative(
                    &self.band,
                    temperature,
                    self.line_weight_model,
                )?;
            (weights, Some(derivatives))
        } else {
            (self.line_weights(temperature)?, None)
        };
        if let Some(ver) = ver {
            Zip::from(weights.rows_mut())
                .and(ver)
                .for_each(|mut row, &v| row *= v);
            if let Some(d_weights) = d_weights.as_mut() {
                Zip::from(d_weights.rows_mut())
                    .and(ver)
                    .for_each(|mut row, &v| row *= v);
            }
        }
        doppler_spectrum(
            inputs,
            self.band.wavelengths_nm().view(),
            weights.view(),
            d_weights.as_ref().map(|d| d.view()),
        )
    }

    /// Population wrappers request T at fixed populations/Einstein coefficients;
    /// direct VER constituents additionally expose their own VER mapping.
    pub fn register_emission_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        name: &str,
        register_ver: bool,
    ) -> Result<()> {
        let (inputs, _, generator) = storage.split_inputs_outputs_deriv();
        let temperature_derivative = inputs.calculate_temperature_derivative();
        if !register_ver && !temperature_derivative {
            return Ok(());
        }
        validate_photon_ver(self.photon_ver.view(), self.altitudes.len())?;
        let interpolation = self.interpolation(inputs);
        let ver = interpolation.dot(&self.photon_ver);
        let (spectrum, d_spectrum) = self.spectrum(
            inputs,
            temperature_derivative,
            if register_ver { None } else { Some(ver.view()) },
        )?;
        if register_ver {
            let full_name = format!("wf_{name}_photon_ver");
            let mut mapping = generator.get_derivative_mapping(&full_name);
            mapping.set_assign_name(&full_name);
            mapping.set_interp_dim(&format!("{name}_altitude"));
            mapping.set_interpolator(&interpolation);
            mapping.mut_view().d_emission.assign(&(&spectrum / FOUR_PI));
        }
        if let Some(mut d_spectrum) = d_spectrum {
            Zip::from(d_spectrum.rows_mut())
                .and(&ver)
                .for_each(|mut row, &v| {
                    row *= if register_ver {
                        v / FOUR_PI
                    } else {
                        1.0 / FOUR_PI
                    }
                });
            let full_name = format!("wf_{name}_temperature_k");
            let mut mapping = generator.get_derivative_mapping(&full_name);
            mapping.set_assign_name("wf_temperature_k");
            mapping.set_interp_dim("altitude");
            mapping.mut_view().d_emission.assign(&d_spectrum);
        }
        Ok(())
    }
}

impl Constituent for BandVolumeEmissionRate {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (inputs, outputs) = storage.split_inputs_outputs();
        validate_photon_ver(self.photon_ver.view(), self.altitudes.len())?;
        let ver = self.interpolation(inputs).dot(&self.photon_ver);
        let (mut spectrum, _) = self.spectrum(inputs, false, Some(ver.view()))?;
        spectrum /= FOUR_PI;
        outputs.mut_view().emission_source += &spectrum;
        Ok(())
    }

    fn register_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        name: &str,
    ) -> Result<()> {
        self.register_emission_derivatives(storage, name, true)
    }
}

pub fn validate_photon_ver(ver: ArrayView1<'_, f64>, size: usize) -> Result<()> {
    anyhow::ensure!(
        ver.len() == size,
        "Photon VER length must match emission altitudes"
    );
    anyhow::ensure!(
        ver.iter().all(|v| v.is_finite() && *v >= 0.0),
        "Photon VER must be finite and non-negative"
    );
    Ok(())
}

/// Select the same O2 line lists as the population interface, but separate the
/// A-band 0-0 and 1-1 transitions so their total VERs can be retrieved independently.
pub fn oxygen_emission_band(db: &OpticalLineDB, transition: &str) -> Result<EmissionBand> {
    let (mut band, upper, lower, a) = match transition {
        "0-0" => (
            EmissionBand::oxygen_a_band_from_hitran(db)?,
            "O2(b)",
            "O2(X)",
            O2_B0_X0_EINSTEIN_A_S,
        ),
        "1-1" => (
            EmissionBand::oxygen_a_band_from_hitran(db)?,
            "O2(b, v=1)",
            "O2(X, v=1)",
            O2_B1_X1_EINSTEIN_A_S,
        ),
        "1-0" => (
            EmissionBand::oxygen_b_band_from_hitran(db)?
                .ok_or_else(|| anyhow!("No O2 1-0 emission lines in database"))?,
            "O2(b, v=1)",
            "O2(X)",
            O2_B1_X0_EINSTEIN_A_S,
        ),
        _ => {
            return Err(anyhow!(
                "Unsupported O2 band '{transition}'; expected '0-0', '1-1', or '1-0'"
            ));
        }
    };
    band.lines.retain(|line| {
        line.upper_vibrational_state == upper && line.lower_vibrational_state == lower
    });
    anyhow::ensure!(
        !band.lines.is_empty(),
        "No O2 {transition} emission lines in database"
    );
    EmissionBand::new(transition, upper, lower, a, band.lines)
}

/// Values use the original Gaussian kernel; derivative arrays and directional
/// kernels are only used when d_areas is supplied. A zero d_areas includes just
/// Doppler broadening, whereas band weights also contribute an amplitude term.
pub(crate) fn doppler_spectrum(
    inputs: &impl StorageInputs,
    wavelengths: ArrayView1<'_, f64>,
    areas: ArrayView2<'_, f64>,
    d_areas: Option<ArrayView2<'_, f64>>,
) -> Result<(Array2<f64>, Option<Array2<f64>>)> {
    let temperature = inputs
        .temperature_k()
        .ok_or_else(|| anyhow!("Temperature must be set for line emission"))?;
    anyhow::ensure!(
        temperature.iter().all(|t| t.is_finite() && *t > 0.0),
        "Emission temperatures must be positive and finite"
    );
    let integrated =
        inputs.spectral_integration_mode() == SpectralGridMode::AtmosphereIntegratedLineShape;
    let grid = if integrated {
        inputs.fine_spectral_grid()
    } else {
        inputs.spectral_grid()
    }
    .ok_or_else(|| anyhow!("Spectral grid must be set for line emission"))?;
    let wavenumbers = grid.central_wavenumber_cminv();
    let mut values = Array2::zeros((temperature.len(), wavenumbers.len()));
    let mut derivatives = d_areas.map(|_| Array2::zeros(values.raw_dim()));
    for (alt, &t) in temperature.iter().enumerate() {
        let mut row = values.row_mut(alt);
        for (i, &wavelength) in wavelengths.iter().enumerate() {
            let area = areas[[alt, i]];
            let center = 1.0e7 / wavelength;
            let width = OpticalLine::doppler_width_cminv(center, t, O2_MASS);
            if let (Some(d_areas), Some(derivatives)) = (d_areas, derivatives.as_mut()) {
                let d_area = d_areas[[alt, i]];
                if area == 0.0 && d_area == 0.0 {
                    continue;
                }
                let line = AdjustedLineParameters {
                    line_center: center,
                    doppler_width: width,
                    line_intensity_re: area / (SQRT_PI * width),
                    ..Default::default()
                };
                let direction = LineShapeDirection {
                    doppler_width: width / (2.0 * t),
                    line_intensity_re: (d_area - area / (2.0 * t)) / (SQRT_PI * width),
                    ..Default::default()
                };
                assign_with_derivative(
                    wavenumbers.as_slice().unwrap(),
                    &line,
                    &direction,
                    LineShape::Gaussian,
                    row.as_slice_mut().unwrap(),
                    derivatives.row_mut(alt).as_slice_mut().unwrap(),
                );
            } else if area != 0.0 {
                assign_normalized_doppler_line_shape(
                    wavenumbers.as_slice().unwrap(),
                    center,
                    width,
                    area,
                    row.as_slice_mut().unwrap(),
                );
            }
        }
    }
    if grid.coordinate() == LineshapeCoordinate::WavelengthNm {
        for (i, &wavelength) in grid.central_wavelengths_nm().iter().enumerate() {
            let jacobian = 1.0e7 / wavelength.powi(2);
            values.column_mut(i).mapv_inplace(|v| v * jacobian);
            if let Some(d) = derivatives.as_mut() {
                d.column_mut(i).mapv_inplace(|v| v * jacobian);
            }
        }
    }
    if integrated {
        let mapping = inputs
            .spectral_mapping_matrix()
            .ok_or_else(|| anyhow!("Missing emission spectral mapping"))?;
        values = mapping.dot(values.view());
        derivatives = derivatives.map(|d| mapping.dot(d.view()));
    }
    Ok((values, derivatives))
}
