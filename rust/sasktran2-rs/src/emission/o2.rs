//! O2 spectroscopy and rotational populations for the generic band VER source.

use super::{BandEmissionModel, BandLineData};
use crate::atmosphere::StorageInputs;
use crate::optical::line::OpticalLineDB;
use crate::photchem::emission::{
    AEmissionLineWeightModel, EmissionBand, O2_B0_X0_EINSTEIN_A_S, O2_B1_X0_EINSTEIN_A_S,
    O2_B1_X1_EINSTEIN_A_S, oxygen_a_band_lte_line_weights,
    oxygen_a_band_lte_line_weights_with_temperature_derivative,
};
use crate::prelude::*;

const O2_MASS_AMU: f64 = 31.9988;

/// Rotational LTE at the atmosphere's kinetic temperature, with an independent
/// vibrational band VER. Both historical O2 line-weight conventions are retained.
#[derive(Clone)]
pub struct O2BandEmissionModel {
    lines: BandLineData,
    band: EmissionBand,
    line_weight_model: AEmissionLineWeightModel,
}

impl O2BandEmissionModel {
    pub fn new(band: EmissionBand, line_weight_model: AEmissionLineWeightModel) -> Result<Self> {
        let first = band
            .lines
            .first()
            .ok_or_else(|| anyhow!("Emission band must contain lines"))?;
        anyhow::ensure!(
            band.lines.iter().all(|line| line.upper_vibrational_state
                == first.upper_vibrational_state
                && line.lower_vibrational_state == first.lower_vibrational_state),
            "An O2 band VER must describe a single vibrational transition"
        );
        Ok(Self {
            lines: BandLineData::new(&band.name, band.wavelengths_nm(), O2_MASS_AMU)?,
            band,
            line_weight_model,
        })
    }

    pub fn line_weights_at_temperature(
        &self,
        temperature_k: ArrayView1<'_, f64>,
    ) -> Result<Array2<f64>> {
        oxygen_a_band_lte_line_weights(&self.band, temperature_k, self.line_weight_model)
    }
}

impl BandEmissionModel for O2BandEmissionModel {
    fn line_data(&self) -> &BandLineData {
        &self.lines
    }

    fn line_weights(&self, inputs: &impl StorageInputs) -> Result<Array2<f64>> {
        self.line_weights_at_temperature(
            inputs
                .temperature_k()
                .ok_or_else(|| anyhow!("Temperature must be set for O2 band emission"))?,
        )
    }

    fn line_weights_with_temperature_derivative(
        &self,
        inputs: &impl StorageInputs,
    ) -> Result<(Array2<f64>, Array2<f64>)> {
        oxygen_a_band_lte_line_weights_with_temperature_derivative(
            &self.band,
            inputs
                .temperature_k()
                .ok_or_else(|| anyhow!("Temperature must be set for O2 band emission"))?,
            self.line_weight_model,
        )
    }
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
