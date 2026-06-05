use crate::optical::line::{OpticalLine, OpticalLineDB};
use crate::prelude::*;

pub const OXYGEN_GREEN_LINE_WAVELENGTH_NM: f64 = 557.7;
pub const OXYGEN_GREEN_LINE_EINSTEIN_A_S: f64 = 1.26;
pub const MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_558_S: f64 = 1.18;
pub const MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_1S_S: f64 = 1.35;
pub const MCDADE_OXYGEN_GREEN_LINE_C0: f64 = 0.0;
pub const MCDADE_OXYGEN_GREEN_LINE_C1: f64 = 211.0;
pub const MCDADE_OXYGEN_GREEN_LINE_C2: f64 = 15.0;
pub const O2_A_BAND_CENTER_WAVELENGTH_NM: f64 = 762.0;
pub const O2_A_BAND_TOTAL_EINSTEIN_A_S: f64 = 7.58e-2;
pub const O2_B0_X0_EINSTEIN_A_S: f64 = 7.58e-2;
pub const O2_B1_X0_EINSTEIN_A_S: f64 = 7.0e-2;
pub const O2_B1_X1_EINSTEIN_A_S: f64 = 7.0e-2;
pub const O2_B2_X2_EINSTEIN_A_S: f64 = 5.4e-2;
const C2_K_CM: f64 = 1.4387769;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AEmissionLineWeightModel {
    EinsteinABranching,
    HitranLineStrength,
}

impl std::str::FromStr for AEmissionLineWeightModel {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "einstein_a_branching" => Ok(Self::EinsteinABranching),
            "hitran_line_strength" => Ok(Self::HitranLineStrength),
            other => Err(anyhow!(
                "Unknown A-band line weight model '{}'. Expected 'einstein_a_branching' or 'hitran_line_strength'",
                other
            )),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EmissionTransition {
    pub name: String,
    pub upper_state: String,
    pub lower_state: String,
    pub wavelength_nm: f64,
    pub einstein_a_s: f64,
}

impl EmissionTransition {
    pub fn new(
        name: impl Into<String>,
        upper_state: impl Into<String>,
        lower_state: impl Into<String>,
        wavelength_nm: f64,
        einstein_a_s: f64,
    ) -> Result<Self> {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(anyhow!(
                "Emission wavelength must be positive and finite, got {}",
                wavelength_nm
            ));
        }

        if !einstein_a_s.is_finite() || einstein_a_s < 0.0 {
            return Err(anyhow!(
                "Einstein-A coefficient must be non-negative and finite, got {}",
                einstein_a_s
            ));
        }

        Ok(Self {
            name: name.into(),
            upper_state: upper_state.into(),
            lower_state: lower_state.into(),
            wavelength_nm,
            einstein_a_s,
        })
    }

    pub fn photon_ver(&self, upper_population: ArrayView1<'_, f64>) -> Array1<f64> {
        upper_population.mapv(|population| population * self.einstein_a_s)
    }
}

pub fn oxygen_green_line_transition() -> EmissionTransition {
    EmissionTransition::new(
        "oxygen_green_5577",
        "O(1S)",
        "O(1D)",
        OXYGEN_GREEN_LINE_WAVELENGTH_NM,
        OXYGEN_GREEN_LINE_EINSTEIN_A_S,
    )
    .expect("hard-coded oxygen green line constants should be valid")
}

pub fn mcdade_oxygen_green_line_photon_ver(
    temperature_k: ArrayView1<'_, f64>,
    atomic_oxygen_density_m3: ArrayView1<'_, f64>,
    o2_density_m3: ArrayView1<'_, f64>,
    n2_density_m3: ArrayView1<'_, f64>,
) -> Result<Array1<f64>> {
    let n = temperature_k.len();
    if atomic_oxygen_density_m3.len() != n || o2_density_m3.len() != n || n2_density_m3.len() != n {
        return Err(anyhow!(
            "Temperature, O, O2, and N2 profiles must have the same length"
        ));
    }

    let mut photon_ver_m3 = Array1::zeros(n);
    for i in 0..n {
        let temperature = temperature_k[i];
        let atomic_o_m3 = atomic_oxygen_density_m3[i];
        let o2_m3 = o2_density_m3[i];
        let n2_m3 = n2_density_m3[i];

        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(anyhow!(
                "Temperature must be positive and finite, got {} at index {}",
                temperature,
                i
            ));
        }

        for (name, value) in [("O", atomic_o_m3), ("O2", o2_m3), ("N2", n2_m3)] {
            if !value.is_finite() || value < 0.0 {
                return Err(anyhow!(
                    "{} density must be non-negative and finite, got {} at index {}",
                    name,
                    value,
                    i
                ));
            }
        }

        let atomic_o_cm3 = atomic_o_m3 / 1.0e6;
        let o2_cm3 = o2_m3 / 1.0e6;
        let n2_cm3 = n2_m3 / 1.0e6;

        let denominator = MCDADE_OXYGEN_GREEN_LINE_C0
            + MCDADE_OXYGEN_GREEN_LINE_C1 * atomic_o_cm3
            + MCDADE_OXYGEN_GREEN_LINE_C2 * o2_cm3;
        if atomic_o_cm3 == 0.0 || denominator <= 0.0 {
            photon_ver_m3[i] = 0.0;
            continue;
        }

        let k1_cm6_s = 4.7e-33 * (300.0 / temperature).powi(2);
        let three_k5_cm3_s = 4.0e-12 * (-865.0 / temperature).exp();
        let radiative_branch = MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_558_S
            / (MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_1S_S + three_k5_cm3_s * o2_cm3);

        let photon_ver_cm3 = k1_cm6_s * atomic_o_cm3.powi(2) * (n2_cm3 + o2_cm3) * atomic_o_cm3
            / denominator
            * radiative_branch;

        photon_ver_m3[i] = photon_ver_cm3 * 1.0e6;
    }

    Ok(photon_ver_m3)
}

pub fn mcdade_oxygen_green_line_o1s_population(
    temperature_k: ArrayView1<'_, f64>,
    atomic_oxygen_density_m3: ArrayView1<'_, f64>,
    o2_density_m3: ArrayView1<'_, f64>,
    n2_density_m3: ArrayView1<'_, f64>,
) -> Result<Array1<f64>> {
    Ok(mcdade_oxygen_green_line_photon_ver(
        temperature_k,
        atomic_oxygen_density_m3,
        o2_density_m3,
        n2_density_m3,
    )?
    .mapv(|ver| ver / MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_558_S))
}

#[derive(Clone, Debug, PartialEq)]
pub struct EmissionBandLine {
    pub wavelength_nm: f64,
    pub wavenumber_cminv: f64,
    pub line_intensity_296: f64,
    pub einstein_a_s: f64,
    pub isotope_id: i32,
    pub isotope_abundance: f64,
    pub lower_energy_cminv: f64,
    pub upper_energy_cminv: f64,
    pub upper_vibrational_state: String,
    pub lower_vibrational_state: String,
    pub upper_state_id: String,
    pub lower_state_id: String,
    pub upper_statistical_weight: Option<f64>,
    pub lower_statistical_weight: Option<f64>,
    pub upper_branching_ratio: f64,
    pub relative_weight: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EmissionBand {
    pub name: String,
    pub upper_state: String,
    pub lower_state: String,
    pub total_einstein_a_s: f64,
    pub lines: Vec<EmissionBandLine>,
}

impl EmissionBand {
    pub fn new(
        name: impl Into<String>,
        upper_state: impl Into<String>,
        lower_state: impl Into<String>,
        total_einstein_a_s: f64,
        mut lines: Vec<EmissionBandLine>,
    ) -> Result<Self> {
        if lines.is_empty() {
            return Err(anyhow!("Emission band must contain at least one line"));
        }

        if !total_einstein_a_s.is_finite() || total_einstein_a_s < 0.0 {
            return Err(anyhow!(
                "Band Einstein-A coefficient must be non-negative and finite, got {}",
                total_einstein_a_s
            ));
        }

        normalize_band_line_weights(&mut lines)?;

        Ok(Self {
            name: name.into(),
            upper_state: upper_state.into(),
            lower_state: lower_state.into(),
            total_einstein_a_s,
            lines,
        })
    }

    pub fn from_hitran_lines(
        name: impl Into<String>,
        upper_state: impl Into<String>,
        lower_state: impl Into<String>,
        total_einstein_a_s: f64,
        db: &OpticalLineDB,
        min_wavelength_nm: f64,
        max_wavelength_nm: f64,
    ) -> Result<Self> {
        if min_wavelength_nm >= max_wavelength_nm {
            return Err(anyhow!(
                "Invalid band wavelength range: {} to {} nm",
                min_wavelength_nm,
                max_wavelength_nm
            ));
        }

        let mut lines: Vec<EmissionBandLine> = db
            .lines
            .iter()
            .filter_map(|line| emission_band_line_from_optical_line(line).transpose())
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .filter(|line| {
                line.wavelength_nm >= min_wavelength_nm && line.wavelength_nm <= max_wavelength_nm
            })
            .collect();

        lines.sort_by(|lhs, rhs| lhs.wavelength_nm.partial_cmp(&rhs.wavelength_nm).unwrap());

        Self::new(name, upper_state, lower_state, total_einstein_a_s, lines)
    }

    pub fn oxygen_a_band_from_hitran(db: &OpticalLineDB) -> Result<Self> {
        let mut lines: Vec<EmissionBandLine> = db
            .lines
            .iter()
            .filter(|line| {
                line_matches_o2_a_band_vibrational_sequence(line)
                    && line.wavelength_nm() >= 759.0
                    && line.wavelength_nm() <= 772.0
            })
            .filter_map(|line| emission_band_line_from_optical_line(line).transpose())
            .collect::<Result<Vec<_>>>()?;

        lines.sort_by(|lhs, rhs| lhs.wavelength_nm.partial_cmp(&rhs.wavelength_nm).unwrap());

        Self::new(
            "oxygen_a_band",
            "O2(b, v=0..1)",
            "O2(X, v=0..1)",
            O2_A_BAND_TOTAL_EINSTEIN_A_S,
            lines,
        )
    }

    pub fn oxygen_b_band_from_hitran(db: &OpticalLineDB) -> Result<Option<Self>> {
        let mut lines: Vec<EmissionBandLine> = db
            .lines
            .iter()
            .filter(|line| {
                line_matches_o2_b_band_vibrational_sequence(line)
                    && line.wavelength_nm() >= 680.0
                    && line.wavelength_nm() <= 700.0
            })
            .filter_map(|line| emission_band_line_from_optical_line(line).transpose())
            .collect::<Result<Vec<_>>>()?;

        if lines.is_empty() {
            return Ok(None);
        }

        lines.sort_by(|lhs, rhs| lhs.wavelength_nm.partial_cmp(&rhs.wavelength_nm).unwrap());

        Self::new(
            "oxygen_b_band",
            "O2(b, v=1)",
            "O2(X)",
            O2_B1_X0_EINSTEIN_A_S,
            lines,
        )
        .map(Some)
    }

    pub fn photon_ver(&self, upper_population: ArrayView1<'_, f64>) -> Array1<f64> {
        upper_population.mapv(|population| population * self.total_einstein_a_s)
    }

    pub fn wavelengths_nm(&self) -> Array1<f64> {
        Array1::from_iter(self.lines.iter().map(|line| line.wavelength_nm))
    }

    pub fn weights(&self) -> Array1<f64> {
        Array1::from_iter(self.lines.iter().map(|line| line.relative_weight))
    }
}

pub trait PopulationProvider {
    fn population(&self, state_name: &str) -> Option<ArrayView1<'_, f64>>;
}

pub struct PopulationMap<'a> {
    profiles: HashMap<String, ArrayView1<'a, f64>>,
}

impl<'a> PopulationMap<'a> {
    pub fn new(profiles: HashMap<String, ArrayView1<'a, f64>>) -> Self {
        Self { profiles }
    }
}

impl PopulationProvider for PopulationMap<'_> {
    fn population(&self, state_name: &str) -> Option<ArrayView1<'_, f64>> {
        self.profiles.get(state_name).map(|profile| profile.view())
    }
}

pub fn photon_ver_from_state_profile(
    state_profiles: &HashMap<String, Vec<f64>>,
    state_name: &str,
    einstein_a_s: f64,
) -> Result<Vec<f64>> {
    let profile = state_profiles
        .get(state_name)
        .ok_or_else(|| anyhow!("Population state '{}' is not available", state_name))?;

    Ok(profile
        .iter()
        .map(|population| population * einstein_a_s)
        .collect())
}

fn emission_band_line_from_optical_line(line: &OpticalLine) -> Result<Option<EmissionBandLine>> {
    let Some(einstein_a_s) = line.einstein_a else {
        return Ok(None);
    };

    if !einstein_a_s.is_finite() || einstein_a_s <= 0.0 {
        return Ok(None);
    }

    Ok(Some(EmissionBandLine {
        wavelength_nm: line.wavelength_nm(),
        wavenumber_cminv: line.line_center,
        line_intensity_296: line.line_intensity,
        einstein_a_s,
        isotope_id: line.iso_id,
        isotope_abundance: o2_hitran_isotope_abundance(line.iso_id),
        lower_energy_cminv: line.lower_energy,
        upper_energy_cminv: line.lower_energy + line.line_center,
        upper_vibrational_state: o2_vibrational_state_name(&line.upper_quanta),
        lower_vibrational_state: o2_vibrational_state_name(&line.lower_quanta),
        upper_state_id: resolved_upper_state_id(line),
        lower_state_id: format!("{} {}", line.lower_quanta, line.lower_local_quanta)
            .trim()
            .to_string(),
        upper_statistical_weight: line.upper_statistical_weight,
        lower_statistical_weight: line.lower_statistical_weight,
        upper_branching_ratio: 0.0,
        relative_weight: einstein_a_s * o2_hitran_isotope_abundance(line.iso_id),
    }))
}

fn resolved_upper_state_id(line: &OpticalLine) -> String {
    if let Some(upper_rotational_state) = o2_group6_upper_rotational_state_id(line) {
        return format!(
            "iso={} {} {}",
            line.iso_id, line.upper_quanta, upper_rotational_state
        )
        .trim()
        .to_string();
    }

    let local_quanta = if line.upper_local_quanta.trim().is_empty() {
        line.lower_local_quanta.trim()
    } else {
        line.upper_local_quanta.trim()
    };

    format!("iso={} {} {}", line.iso_id, line.upper_quanta, local_quanta)
        .trim()
        .to_string()
}

fn o2_hitran_isotope_abundance(iso_id: i32) -> f64 {
    match iso_id {
        1 => 0.995_261_6,
        2 => 0.003_991_41,
        3 => 0.000_742_235_2,
        _ => 0.0,
    }
}

fn o2_group6_upper_rotational_state_id(line: &OpticalLine) -> Option<String> {
    if line.mol_id != 7 {
        return None;
    }

    let local = line.lower_local_quanta.as_str();
    let offset = if local.len() >= 15 { 1 } else { 0 };
    if local.len() < offset + 14 {
        return None;
    }

    let n_branch = local.get(offset..offset + 1)?.chars().next()?;
    let lower_n: i32 = local.get(offset + 1..offset + 4)?.trim().parse().ok()?;
    let j_branch = local.get(offset + 4..offset + 5)?.chars().next()?;
    let lower_j: i32 = local.get(offset + 5..offset + 8)?.trim().parse().ok()?;

    let upper_n = lower_n + branch_delta(n_branch)?;
    let upper_j = lower_j + branch_delta(j_branch)?;
    if upper_n < 0 || upper_j < 0 {
        return None;
    }

    Some(format!("N'={upper_n} J'={upper_j}"))
}

fn branch_delta(branch: char) -> Option<i32> {
    match branch {
        'M' => Some(-4),
        'N' => Some(-3),
        'O' => Some(-2),
        'P' => Some(-1),
        'Q' => Some(0),
        'R' => Some(1),
        'S' => Some(2),
        'T' => Some(3),
        'U' => Some(4),
        _ => None,
    }
}

fn line_matches_o2_a_band_vibrational_sequence(line: &OpticalLine) -> bool {
    let upper_tokens: Vec<&str> = line.upper_quanta.split_whitespace().collect();
    let lower_tokens: Vec<&str> = line.lower_quanta.split_whitespace().collect();

    matches!(
        (upper_tokens.as_slice(), lower_tokens.as_slice()),
        (["b", upper_v], ["X", lower_v])
            if upper_v == lower_v && (*upper_v == "0" || *upper_v == "1")
    )
}

fn line_matches_o2_b_band_vibrational_sequence(line: &OpticalLine) -> bool {
    let upper_tokens: Vec<&str> = line.upper_quanta.split_whitespace().collect();
    let lower_tokens: Vec<&str> = line.lower_quanta.split_whitespace().collect();

    matches!(
        (upper_tokens.as_slice(), lower_tokens.as_slice()),
        (["b", "1"], ["X", "0"])
    )
}

fn o2_vibrational_state_name(quanta: &str) -> String {
    let tokens: Vec<&str> = quanta.split_whitespace().collect();
    match tokens.as_slice() {
        ["b", "0"] => "O2(b)".to_string(),
        ["X", "0"] => "O2(X)".to_string(),
        [electronic, vibrational] => {
            format!("O2({electronic}, v={vibrational})")
        }
        _ => quanta.trim().to_string(),
    }
}

pub fn oxygen_a_band_lte_line_weights(
    band: &EmissionBand,
    temperature_k: ArrayView1<'_, f64>,
    line_weight_model: AEmissionLineWeightModel,
) -> Result<Array2<f64>> {
    match line_weight_model {
        AEmissionLineWeightModel::EinsteinABranching => {
            oxygen_a_band_einstein_branching_line_weights(band, temperature_k)
        }
        AEmissionLineWeightModel::HitranLineStrength => {
            oxygen_a_band_hitran_line_strength_weights(band, temperature_k)
        }
    }
}

pub fn oxygen_a_band_line_list_weights_from_populations<'a>(
    band: &EmissionBand,
    temperature_k: ArrayView1<'_, f64>,
    o2_b0_population: ArrayView1<'a, f64>,
    o2_b1_population: Option<ArrayView1<'a, f64>>,
    o2_b2_population: Option<ArrayView1<'a, f64>>,
    line_weight_model: AEmissionLineWeightModel,
) -> Result<(Array1<f64>, Array2<f64>)> {
    validate_same_len(
        "temperature",
        temperature_k.len(),
        "O2(b)",
        o2_b0_population.len(),
    )?;
    if let Some(o2_b1_population) = o2_b1_population {
        validate_same_len(
            "temperature",
            temperature_k.len(),
            "O2(b, v=1)",
            o2_b1_population.len(),
        )?;
    }
    if let Some(o2_b2_population) = o2_b2_population {
        validate_same_len(
            "temperature",
            temperature_k.len(),
            "O2(b, v=2)",
            o2_b2_population.len(),
        )?;
    }

    let branches = [
        EmissionPopulationBranch {
            upper_vibrational_state: "O2(b)",
            population: Some(o2_b0_population),
            einstein_a_s: O2_B0_X0_EINSTEIN_A_S,
        },
        EmissionPopulationBranch {
            upper_vibrational_state: "O2(b, v=1)",
            population: o2_b1_population,
            einstein_a_s: O2_B1_X1_EINSTEIN_A_S,
        },
        EmissionPopulationBranch {
            upper_vibrational_state: "O2(b, v=2)",
            population: o2_b2_population,
            einstein_a_s: O2_B2_X2_EINSTEIN_A_S,
        },
    ];

    line_list_weights_from_population_branches(band, temperature_k, &branches, line_weight_model)
}

pub fn oxygen_b_band_line_list_weights_from_populations<'a>(
    band: &EmissionBand,
    temperature_k: ArrayView1<'_, f64>,
    o2_b1_population: Option<ArrayView1<'a, f64>>,
    line_weight_model: AEmissionLineWeightModel,
) -> Result<(Array1<f64>, Array2<f64>)> {
    if let Some(o2_b1_population) = o2_b1_population {
        validate_same_len(
            "temperature",
            temperature_k.len(),
            "O2(b, v=1)",
            o2_b1_population.len(),
        )?;
    }

    let branches = [EmissionPopulationBranch {
        upper_vibrational_state: "O2(b, v=1)",
        population: o2_b1_population,
        einstein_a_s: O2_B1_X0_EINSTEIN_A_S,
    }];

    line_list_weights_from_population_branches(band, temperature_k, &branches, line_weight_model)
}

struct EmissionPopulationBranch<'a> {
    upper_vibrational_state: &'a str,
    population: Option<ArrayView1<'a, f64>>,
    einstein_a_s: f64,
}

fn line_list_weights_from_population_branches(
    band: &EmissionBand,
    temperature_k: ArrayView1<'_, f64>,
    branches: &[EmissionPopulationBranch<'_>],
    line_weight_model: AEmissionLineWeightModel,
) -> Result<(Array1<f64>, Array2<f64>)> {
    let weights = oxygen_a_band_lte_line_weights(band, temperature_k, line_weight_model)?;
    let mut line_photon_ver = Array2::<f64>::zeros((temperature_k.len(), band.lines.len()));

    for (line_idx, line) in band.lines.iter().enumerate() {
        let state_ver = branches
            .iter()
            .find(|branch| branch.upper_vibrational_state == line.upper_vibrational_state)
            .and_then(|branch| {
                branch.population.map(|population| {
                    population
                        .iter()
                        .map(|population| population * branch.einstein_a_s)
                        .collect::<Array1<_>>()
                })
            })
            .unwrap_or_else(|| Array1::zeros(temperature_k.len()));

        Zip::from(line_photon_ver.column_mut(line_idx))
            .and(&state_ver)
            .and(weights.column(line_idx))
            .for_each(|line_ver, &state_ver, &weight| {
                *line_ver = state_ver * weight;
            });
    }

    let total_photon_ver = line_photon_ver.sum_axis(Axis(1));
    let fallback_weights = normalized_fallback_weights(band)?;
    let mut combined_weights = Array2::<f64>::zeros(line_photon_ver.raw_dim());

    Zip::from(combined_weights.rows_mut())
        .and(line_photon_ver.rows())
        .and(&total_photon_ver)
        .for_each(|mut row, line_ver, &total_ver| {
            if total_ver > 0.0 {
                row.assign(&(&line_ver / total_ver));
            } else {
                row.assign(&fallback_weights);
            }
        });

    Ok((total_photon_ver, combined_weights))
}

fn oxygen_a_band_einstein_branching_line_weights(
    band: &EmissionBand,
    temperature_k: ArrayView1<'_, f64>,
) -> Result<Array2<f64>> {
    let mut weights = Array2::<f64>::zeros((temperature_k.len(), band.lines.len()));

    for vibrational_state in unique_upper_vibrational_states(band) {
        let line_indices = line_indices_for_upper_vibrational_state(band, &vibrational_state);

        for (alt_idx, &temperature) in temperature_k.iter().enumerate() {
            let mut row_sum = 0.0;
            for &line_idx in &line_indices {
                let line = &band.lines[line_idx];
                let upper_g = line.upper_statistical_weight.unwrap_or(0.0);
                let upper_population_weight = line.isotope_abundance
                    * upper_g
                    * (-C2_K_CM * line.upper_energy_cminv / temperature).exp();
                let weight = upper_population_weight * line.upper_branching_ratio;
                weights[[alt_idx, line_idx]] = weight;
                row_sum += weight;
            }

            normalize_weight_row(&mut weights.row_mut(alt_idx), &line_indices, row_sum)?;
        }
    }

    Ok(weights)
}

fn oxygen_a_band_hitran_line_strength_weights(
    band: &EmissionBand,
    temperature_k: ArrayView1<'_, f64>,
) -> Result<Array2<f64>> {
    let mut weights = Array2::<f64>::zeros((temperature_k.len(), band.lines.len()));

    for vibrational_state in unique_upper_vibrational_states(band) {
        let line_indices = line_indices_for_upper_vibrational_state(band, &vibrational_state);

        for (alt_idx, &temperature) in temperature_k.iter().enumerate() {
            let mut max_log_weight = f64::NEG_INFINITY;
            for &line_idx in &line_indices {
                let line = &band.lines[line_idx];
                let log_weight = hitran_line_strength_emission_log_weight(line, temperature);
                weights[[alt_idx, line_idx]] = log_weight;
                max_log_weight = max_log_weight.max(log_weight);
            }

            let mut row_sum = 0.0;
            for &line_idx in &line_indices {
                let weight = (weights[[alt_idx, line_idx]] - max_log_weight).exp();
                weights[[alt_idx, line_idx]] = weight;
                row_sum += weight;
            }

            normalize_weight_row(&mut weights.row_mut(alt_idx), &line_indices, row_sum)?;
        }
    }

    Ok(weights)
}

fn hitran_line_strength_emission_log_weight(line: &EmissionBandLine, temperature: f64) -> f64 {
    line.line_intensity_296.ln()
        + (296.0 / temperature).ln()
        + C2_K_CM * line.lower_energy_cminv * (temperature - 296.0) / (temperature * 296.0)
        + 2.0 * line.wavenumber_cminv.ln()
        - C2_K_CM * line.wavenumber_cminv / temperature
}

fn unique_upper_vibrational_states(band: &EmissionBand) -> Vec<String> {
    let mut states = Vec::new();
    for line in &band.lines {
        if !states.contains(&line.upper_vibrational_state) {
            states.push(line.upper_vibrational_state.clone());
        }
    }
    states
}

fn line_indices_for_upper_vibrational_state(
    band: &EmissionBand,
    vibrational_state: &str,
) -> Vec<usize> {
    band.lines
        .iter()
        .enumerate()
        .filter_map(|(idx, line)| {
            (line.upper_vibrational_state == vibrational_state).then_some(idx)
        })
        .collect()
}

fn normalize_weight_row(
    row: &mut ndarray::ArrayViewMut1<'_, f64>,
    line_indices: &[usize],
    row_sum: f64,
) -> Result<()> {
    if !row_sum.is_finite() || row_sum <= 0.0 {
        return Err(anyhow!(
            "A-band line weights must contain at least one positive finite value"
        ));
    }

    for &line_idx in line_indices {
        row[line_idx] /= row_sum;
    }

    Ok(())
}

fn normalized_fallback_weights(band: &EmissionBand) -> Result<Array1<f64>> {
    let mut weights = band.weights();
    let sum = weights.sum();
    if !sum.is_finite() || sum <= 0.0 {
        return Err(anyhow!(
            "A-band fallback line weights must have a positive sum"
        ));
    }
    weights /= sum;
    Ok(weights)
}

fn validate_same_len(lhs_name: &str, lhs_len: usize, rhs_name: &str, rhs_len: usize) -> Result<()> {
    if lhs_len != rhs_len {
        return Err(anyhow!(
            "{} length ({}) must match {} length ({})",
            lhs_name,
            lhs_len,
            rhs_name,
            rhs_len
        ));
    }

    Ok(())
}

fn normalize_band_line_weights(lines: &mut [EmissionBandLine]) -> Result<()> {
    let mut upper_level_a_sums: HashMap<String, f64> = HashMap::new();
    let mut vibrational_band_weight_sums: HashMap<String, f64> = HashMap::new();
    for line in lines.iter() {
        *upper_level_a_sums
            .entry(line.upper_state_id.clone())
            .or_insert(0.0) += line.einstein_a_s;
        if line.relative_weight.is_finite() && line.relative_weight > 0.0 {
            *vibrational_band_weight_sums
                .entry(line.upper_vibrational_state.clone())
                .or_insert(0.0) += line.relative_weight;
        }
    }

    for line in lines {
        if !line.relative_weight.is_finite() || line.relative_weight < 0.0 {
            return Err(anyhow!(
                "Emission band line weight must be non-negative and finite, got {}",
                line.relative_weight
            ));
        }

        let upper_a_sum = upper_level_a_sums
            .get(&line.upper_state_id)
            .copied()
            .unwrap_or(0.0);
        if upper_a_sum > 0.0 {
            line.upper_branching_ratio = line.einstein_a_s / upper_a_sum;
        }

        let band_weight_sum = vibrational_band_weight_sums
            .get(&line.upper_vibrational_state)
            .copied()
            .unwrap_or(0.0);
        if band_weight_sum <= 0.0 {
            return Err(anyhow!(
                "Emission band line weights must contain at least one positive finite value for {}",
                line.upper_vibrational_state
            ));
        }

        line.relative_weight /= band_weight_sum;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transition_photon_ver_is_population_times_einstein_a() {
        let transition = EmissionTransition::new("test", "upper", "lower", 500.0, 2.0).unwrap();
        let population = array![1.0, 3.0, 5.0];

        let ver = transition.photon_ver(population.view());

        assert_eq!(ver, array![2.0, 6.0, 10.0]);
    }

    #[test]
    fn band_line_weights_are_normalized() {
        let band = EmissionBand::new(
            "test_band",
            "upper",
            "lower",
            1.0,
            vec![
                EmissionBandLine {
                    wavelength_nm: 760.0,
                    wavenumber_cminv: 1.0e7 / 760.0,
                    line_intensity_296: 1.0,
                    einstein_a_s: 2.0,
                    isotope_id: 1,
                    isotope_abundance: 1.0,
                    lower_energy_cminv: 0.0,
                    upper_energy_cminv: 1.0e7 / 760.0,
                    upper_vibrational_state: "O2(b)".to_string(),
                    lower_vibrational_state: "O2(X)".to_string(),
                    upper_state_id: "u1".to_string(),
                    lower_state_id: "l1".to_string(),
                    upper_statistical_weight: Some(3.0),
                    lower_statistical_weight: Some(1.0),
                    upper_branching_ratio: 0.0,
                    relative_weight: 2.0,
                },
                EmissionBandLine {
                    wavelength_nm: 761.0,
                    wavenumber_cminv: 1.0e7 / 761.0,
                    line_intensity_296: 1.0,
                    einstein_a_s: 6.0,
                    isotope_id: 1,
                    isotope_abundance: 1.0,
                    lower_energy_cminv: 0.0,
                    upper_energy_cminv: 1.0e7 / 761.0,
                    upper_vibrational_state: "O2(b)".to_string(),
                    lower_vibrational_state: "O2(X)".to_string(),
                    upper_state_id: "u2".to_string(),
                    lower_state_id: "l2".to_string(),
                    upper_statistical_weight: Some(5.0),
                    lower_statistical_weight: Some(3.0),
                    upper_branching_ratio: 0.0,
                    relative_weight: 6.0,
                },
            ],
        )
        .unwrap();

        let weight_sum: f64 = band.lines.iter().map(|line| line.relative_weight).sum();
        assert!((weight_sum - 1.0).abs() < 1.0e-12);
        assert!((band.lines[0].relative_weight - 0.25).abs() < 1.0e-12);
        assert!((band.lines[1].relative_weight - 0.75).abs() < 1.0e-12);
    }

    #[test]
    fn upper_branching_ratios_are_normalized_by_upper_level() {
        let band = EmissionBand::new(
            "test_band",
            "upper",
            "lower",
            1.0,
            vec![
                test_band_line(760.0, "upper_a", "lower_1", 2.0),
                test_band_line(761.0, "upper_a", "lower_2", 6.0),
                test_band_line(762.0, "upper_b", "lower_3", 5.0),
            ],
        )
        .unwrap();

        assert!((band.lines[0].upper_branching_ratio - 0.25).abs() < 1.0e-12);
        assert!((band.lines[1].upper_branching_ratio - 0.75).abs() < 1.0e-12);
        assert!((band.lines[2].upper_branching_ratio - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn o2_group6_lower_local_quanta_infers_upper_rotational_state() {
        let mut line = test_line(760.0, "b 0", "X 0", 1.0);
        line.upper_local_quanta = "".to_string();

        line.lower_local_quanta = "P 37P 37     d".to_string();
        assert_eq!(
            resolved_upper_state_id(&line),
            "iso=1 b 0 N'=36 J'=36".to_string()
        );

        line.lower_local_quanta = "N 19O 18     q".to_string();
        assert_eq!(
            resolved_upper_state_id(&line),
            "iso=1 b 0 N'=16 J'=16".to_string()
        );

        line.lower_local_quanta = "T  3S  4     q".to_string();
        assert_eq!(
            resolved_upper_state_id(&line),
            "iso=1 b 0 N'=6 J'=6".to_string()
        );

        line.lower_local_quanta = " P 37P 37     d".to_string();
        assert_eq!(
            resolved_upper_state_id(&line),
            "iso=1 b 0 N'=36 J'=36".to_string()
        );
    }

    #[test]
    fn oxygen_a_band_branching_groups_lines_by_inferred_upper_rotational_state() {
        let mut pp = test_line(760.0, "b 0", "X 0", 2.0);
        pp.upper_local_quanta = "".to_string();
        pp.lower_local_quanta = "P 37P 37     d".to_string();

        let mut pq = test_line(761.0, "b 0", "X 0", 6.0);
        pq.upper_local_quanta = "".to_string();
        pq.lower_local_quanta = "P 37Q 36     d".to_string();

        let db = OpticalLineDB {
            lines: vec![pp, pq],
        };

        let band = EmissionBand::oxygen_a_band_from_hitran(&db).unwrap();

        assert_eq!(band.lines.len(), 2);
        assert_eq!(band.lines[0].upper_state_id, "iso=1 b 0 N'=36 J'=36");
        assert_eq!(band.lines[1].upper_state_id, "iso=1 b 0 N'=36 J'=36");
        assert!((band.lines[0].upper_branching_ratio - 0.25).abs() < 1.0e-12);
        assert!((band.lines[1].upper_branching_ratio - 0.75).abs() < 1.0e-12);
    }

    #[test]
    fn oxygen_a_band_includes_b0_x0_and_b1_x1_lines() {
        let db = OpticalLineDB {
            lines: vec![
                test_line(760.0, "b 0", "X 0", 2.0),
                test_line(761.0, "b 1", "X 1", 8.0),
                test_line(762.0, "b 0", "X 0", 6.0),
            ],
        };

        let band = EmissionBand::oxygen_a_band_from_hitran(&db).unwrap();

        assert_eq!(band.lines.len(), 3);
        assert_eq!(band.lines[1].upper_vibrational_state, "O2(b, v=1)");
        assert_eq!(band.lines[1].lower_vibrational_state, "O2(X, v=1)");

        let b0_weight_sum: f64 = band
            .lines
            .iter()
            .filter(|line| line.upper_vibrational_state == "O2(b)")
            .map(|line| line.relative_weight)
            .sum();
        let b1_weight_sum: f64 = band
            .lines
            .iter()
            .filter(|line| line.upper_vibrational_state == "O2(b, v=1)")
            .map(|line| line.relative_weight)
            .sum();

        assert!((b0_weight_sum - 1.0).abs() < 1.0e-12);
        assert!((b1_weight_sum - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn oxygen_b_band_includes_b1_x0_lines_only() {
        let db = OpticalLineDB {
            lines: vec![
                test_line(688.0, "b 1", "X 0", 2.0),
                test_line(689.0, "b 1", "X 1", 8.0),
                test_line(690.0, "b 2", "X 1", 6.0),
                test_line(760.0, "b 0", "X 0", 4.0),
            ],
        };

        let band = EmissionBand::oxygen_b_band_from_hitran(&db)
            .unwrap()
            .unwrap();

        assert_eq!(band.lines.len(), 1);
        assert_eq!(band.lines[0].upper_vibrational_state, "O2(b, v=1)");
        assert_eq!(band.lines[0].lower_vibrational_state, "O2(X)");
        assert!((band.lines[0].relative_weight - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn oxygen_b_band_missing_b1_population_contributes_zero() {
        let db = OpticalLineDB {
            lines: vec![test_line(688.0, "b 1", "X 0", 2.0)],
        };
        let band = EmissionBand::oxygen_b_band_from_hitran(&db)
            .unwrap()
            .unwrap();

        let (photon_ver, weights) = oxygen_b_band_line_list_weights_from_populations(
            &band,
            array![200.0, 210.0].view(),
            None,
            AEmissionLineWeightModel::EinsteinABranching,
        )
        .unwrap();

        assert_eq!(photon_ver, array![0.0, 0.0]);
        assert_eq!(weights, array![[1.0], [1.0]]);
    }

    #[test]
    fn missing_state_profile_is_an_error() {
        let profiles = HashMap::new();
        let err = photon_ver_from_state_profile(&profiles, "O(1S)", 1.0).unwrap_err();
        assert!(err.to_string().contains("O(1S)"));
    }

    #[test]
    fn mcdade_green_line_ver_matches_short_form_equation() {
        let temperature = array![200.0];
        let atomic_o = array![5.0e17];
        let o2 = array![1.0e18];
        let n2 = array![4.0e18];

        let ver = mcdade_oxygen_green_line_photon_ver(
            temperature.view(),
            atomic_o.view(),
            o2.view(),
            n2.view(),
        )
        .unwrap();

        let atomic_o_cm3: f64 = 5.0e11;
        let o2_cm3: f64 = 1.0e12;
        let n2_cm3: f64 = 4.0e12;
        let k1 = 4.7e-33 * (300.0_f64 / 200.0).powi(2);
        let three_k5 = 4.0e-12 * (-865.0_f64 / 200.0).exp();
        let expected_cm3 = k1 * atomic_o_cm3.powi(2) * (n2_cm3 + o2_cm3) * atomic_o_cm3
            / (MCDADE_OXYGEN_GREEN_LINE_C1 * atomic_o_cm3 + MCDADE_OXYGEN_GREEN_LINE_C2 * o2_cm3)
            * MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_558_S
            / (MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_1S_S + three_k5 * o2_cm3);

        assert!((ver[0] - expected_cm3 * 1.0e6).abs() / ver[0] < 1.0e-12);
    }

    #[test]
    fn mcdade_green_line_rejects_negative_density() {
        let err = mcdade_oxygen_green_line_photon_ver(
            array![200.0].view(),
            array![-1.0].view(),
            array![1.0].view(),
            array![1.0].view(),
        )
        .unwrap_err();

        assert!(err.to_string().contains("O density"));
    }

    fn test_line(
        wavelength_nm: f64,
        upper_quanta: impl Into<String>,
        lower_quanta: impl Into<String>,
        einstein_a_s: f64,
    ) -> OpticalLine {
        OpticalLine {
            line_center: 1.0e7 / wavelength_nm,
            line_intensity: 1.0,
            einstein_a: Some(einstein_a_s),
            lower_energy: 0.0,
            gamma_air: 0.0,
            gamma_self: 0.0,
            delta_air: 0.0,
            n_air: 0.0,
            mol_id: 7,
            iso_id: 1,
            upper_quanta: upper_quanta.into(),
            lower_quanta: lower_quanta.into(),
            upper_local_quanta: "P 1P 1 d".to_string(),
            lower_local_quanta: "P 1P 1 d".to_string(),
            upper_statistical_weight: Some(3.0),
            lower_statistical_weight: Some(1.0),
            y_coupling: vec![],
            g_coupling: vec![],
            coupling_temperature: vec![],
        }
    }

    fn test_band_line(
        wavelength_nm: f64,
        upper_state_id: impl Into<String>,
        lower_state_id: impl Into<String>,
        einstein_a_s: f64,
    ) -> EmissionBandLine {
        EmissionBandLine {
            wavelength_nm,
            wavenumber_cminv: 1.0e7 / wavelength_nm,
            line_intensity_296: 1.0,
            einstein_a_s,
            isotope_id: 1,
            isotope_abundance: 1.0,
            lower_energy_cminv: 0.0,
            upper_energy_cminv: 1.0e7 / wavelength_nm,
            upper_vibrational_state: "O2(b)".to_string(),
            lower_vibrational_state: "O2(X)".to_string(),
            upper_state_id: upper_state_id.into(),
            lower_state_id: lower_state_id.into(),
            upper_statistical_weight: Some(3.0),
            lower_statistical_weight: Some(1.0),
            upper_branching_ratio: 0.0,
            relative_weight: einstein_a_s,
        }
    }
}
