use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::constituent::types::line_list_volume_emission_rate::LineListVolumeEmissionRate;
use crate::constituent::types::volume_emission_rate::MonochromaticVolumeEmissionRate;
use crate::optical::line::OpticalLineDB;
use crate::photchem::emission::{
    AEmissionLineWeightModel, EmissionBand, oxygen_a_band_line_list_weights_from_populations,
    oxygen_b_band_line_list_weights_from_populations,
};
use crate::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PopulationEmissionSpecies {
    O2,
}

impl std::str::FromStr for PopulationEmissionSpecies {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "O2" | "o2" => Ok(Self::O2),
            other => Err(anyhow!(
                "Unknown population emission species '{}'. Only 'O2' is currently supported",
                other
            )),
        }
    }
}

pub struct PopulationEmissionProfiles {
    pub altitudes_m: Array1<f64>,
    pub temperature_k: Array1<f64>,
    pub populations: HashMap<String, Array1<f64>>,
}

impl PopulationEmissionProfiles {
    pub fn new(
        altitudes_m: Array1<f64>,
        temperature_k: Array1<f64>,
        populations: HashMap<String, Array1<f64>>,
    ) -> Result<Self> {
        if altitudes_m.is_empty() {
            return Err(anyhow!("Population emission altitude grid cannot be empty"));
        }

        if altitudes_m.len() != temperature_k.len() {
            return Err(anyhow!(
                "Altitude length ({}) must match temperature length ({})",
                altitudes_m.len(),
                temperature_k.len()
            ));
        }

        for (name, population) in &populations {
            if population.len() != altitudes_m.len() {
                return Err(anyhow!(
                    "Population '{}' length ({}) must match altitude length ({})",
                    name,
                    population.len(),
                    altitudes_m.len()
                ));
            }
        }

        Ok(Self {
            altitudes_m,
            temperature_k,
            populations,
        })
    }

    fn required_population(&self, name: &str) -> Result<ArrayView1<'_, f64>> {
        self.populations
            .get(name)
            .map(|population| population.view())
            .ok_or_else(|| anyhow!("Required population '{}' is missing", name))
    }

    fn optional_population(&self, name: &str) -> Option<ArrayView1<'_, f64>> {
        self.populations
            .get(name)
            .map(|population| population.view())
    }
}

pub struct PopulationEmissionRate {
    pub line_list_emissions: Vec<LineListVolumeEmissionRate>,
    pub monochromatic_emissions: Vec<MonochromaticVolumeEmissionRate>,
}

impl PopulationEmissionRate {
    pub fn new(
        profiles: PopulationEmissionProfiles,
        species: &[PopulationEmissionSpecies],
        line_weight_model: AEmissionLineWeightModel,
        o2_hitran_db: Option<&OpticalLineDB>,
    ) -> Result<Self> {
        if species.is_empty() {
            return Err(anyhow!(
                "PopulationEmissionRate requires at least one emission species"
            ));
        }

        let mut line_list_emissions = Vec::new();
        let monochromatic_emissions = Vec::new();

        for species in species {
            match species {
                PopulationEmissionSpecies::O2 => {
                    let db = o2_hitran_db.ok_or_else(|| {
                        anyhow!("O2 population emission requires an O2 HITRAN line database")
                    })?;
                    line_list_emissions.push(o2_a_band_line_list_emission(
                        &profiles,
                        db,
                        line_weight_model,
                    )?);
                    if let Some(emission) =
                        o2_b_band_line_list_emission(&profiles, db, line_weight_model)?
                    {
                        line_list_emissions.push(emission);
                    }
                }
            }
        }

        Ok(Self {
            line_list_emissions,
            monochromatic_emissions,
        })
    }

    pub fn primary_line_list(&self) -> Option<&LineListVolumeEmissionRate> {
        self.line_list_emissions.first()
    }

    pub fn with_interp_mode(mut self, interp_mode: crate::interpolation::OutOfBoundsMode) -> Self {
        self.line_list_emissions = self
            .line_list_emissions
            .into_iter()
            .map(|emission| emission.with_interp_mode(interp_mode))
            .collect();
        self.monochromatic_emissions = self
            .monochromatic_emissions
            .into_iter()
            .map(|emission| emission.with_interp_mode(interp_mode))
            .collect();
        self
    }
}

impl Constituent for PopulationEmissionRate {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        for emission in &self.line_list_emissions {
            emission.add_to_atmosphere(storage)?;
        }

        for emission in &self.monochromatic_emissions {
            emission.add_to_atmosphere(storage)?;
        }

        Ok(())
    }

    fn register_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()> {
        for emission in &self.line_list_emissions {
            emission.register_derivatives(storage, constituent_name)?;
        }

        for emission in &self.monochromatic_emissions {
            emission.register_derivatives(storage, constituent_name)?;
        }

        Ok(())
    }
}

fn o2_a_band_line_list_emission(
    profiles: &PopulationEmissionProfiles,
    db: &OpticalLineDB,
    line_weight_model: AEmissionLineWeightModel,
) -> Result<LineListVolumeEmissionRate> {
    let band = EmissionBand::oxygen_a_band_from_hitran(db)?;
    let o2_b0 = profiles.required_population("O2(b)")?;
    let o2_b1 = profiles.optional_population("O2(b, v=1)");
    let o2_b2 = profiles.optional_population("O2(b, v=2)");

    let (photon_ver, weights) = oxygen_a_band_line_list_weights_from_populations(
        &band,
        profiles.temperature_k.view(),
        o2_b0,
        o2_b1,
        o2_b2,
        line_weight_model,
    )?;

    LineListVolumeEmissionRate::new(
        profiles.altitudes_m.clone(),
        photon_ver,
        band.wavelengths_nm(),
        weights,
    )
}

fn o2_b_band_line_list_emission(
    profiles: &PopulationEmissionProfiles,
    db: &OpticalLineDB,
    line_weight_model: AEmissionLineWeightModel,
) -> Result<Option<LineListVolumeEmissionRate>> {
    let Some(band) = EmissionBand::oxygen_b_band_from_hitran(db)? else {
        return Ok(None);
    };
    let o2_b1 = profiles.optional_population("O2(b, v=1)");

    let (photon_ver, weights) = oxygen_b_band_line_list_weights_from_populations(
        &band,
        profiles.temperature_k.view(),
        o2_b1,
        line_weight_model,
    )?;

    Ok(Some(LineListVolumeEmissionRate::new(
        profiles.altitudes_m.clone(),
        photon_ver,
        band.wavelengths_nm(),
        weights,
    )?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_optional_b1_and_b2_populations_contribute_zero() {
        let line = crate::photchem::emission::EmissionBandLine {
            wavelength_nm: 760.0,
            wavenumber_cminv: 1.0e7 / 760.0,
            line_intensity_296: 1.0,
            einstein_a_s: 1.0,
            isotope_id: 1,
            isotope_abundance: 1.0,
            lower_energy_cminv: 0.0,
            upper_energy_cminv: 1.0e7 / 760.0,
            upper_vibrational_state: "O2(b)".to_string(),
            lower_vibrational_state: "O2(X)".to_string(),
            upper_state_id: "iso=1 b 0 N'=1 J'=1".to_string(),
            lower_state_id: "X 0".to_string(),
            upper_statistical_weight: Some(3.0),
            lower_statistical_weight: Some(1.0),
            upper_branching_ratio: 0.0,
            relative_weight: 1.0,
        };
        let band = EmissionBand::new("test", "O2(b)", "O2(X)", 1.0, vec![line]).unwrap();

        let (photon_ver, weights) = oxygen_a_band_line_list_weights_from_populations(
            &band,
            array![200.0].view(),
            array![10.0].view(),
            None,
            None,
            AEmissionLineWeightModel::EinsteinABranching,
        )
        .unwrap();

        assert!(
            (photon_ver[0] - 10.0 * crate::photchem::emission::O2_B0_X0_EINSTEIN_A_S).abs()
                < 1.0e-12
        );
        assert!((weights[[0, 0]] - 1.0).abs() < 1.0e-12);
    }
}
