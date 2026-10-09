//! Emission types from `sasktran2-nlte`, plus adapters that build O2 emission
//! bands from HITRAN line records. The adapters live here because they depend
//! on [`OpticalLine`], which the `sasktran2-nlte` crate does not know about.

pub use sasktran2_nlte::emission::*;

use crate::optical::line::{OpticalLine, OpticalLineDB};
use crate::prelude::*;

pub fn emission_band_from_hitran_lines(
    name: impl Into<String>,
    upper_state: impl Into<String>,
    lower_state: impl Into<String>,
    total_einstein_a_s: f64,
    db: &OpticalLineDB,
    min_wavelength_nm: f64,
    max_wavelength_nm: f64,
) -> Result<EmissionBand> {
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

    EmissionBand::new(name, upper_state, lower_state, total_einstein_a_s, lines)
}

pub fn oxygen_a_band_from_hitran(db: &OpticalLineDB) -> Result<EmissionBand> {
    let mut lines: Vec<EmissionBandLine> = db
        .lines
        .iter()
        .filter(|line| {
            line_matches_o2_a_band_vibrational_sequence(line)
                && line.wavelength_nm() >= O2_A_BAND_MIN_WAVELENGTH_NM
                && line.wavelength_nm() <= O2_A_BAND_MAX_WAVELENGTH_NM
        })
        .filter_map(|line| emission_band_line_from_optical_line(line).transpose())
        .collect::<Result<Vec<_>>>()?;

    lines.sort_by(|lhs, rhs| lhs.wavelength_nm.partial_cmp(&rhs.wavelength_nm).unwrap());

    EmissionBand::new(
        "oxygen_a_band",
        "O2(b, v=0..1)",
        "O2(X, v=0..1)",
        O2_A_BAND_TOTAL_EINSTEIN_A_S,
        lines,
    )
}

pub fn oxygen_b_band_from_hitran(db: &OpticalLineDB) -> Result<Option<EmissionBand>> {
    let mut lines: Vec<EmissionBandLine> = db
        .lines
        .iter()
        .filter(|line| {
            line_matches_o2_b_band_vibrational_sequence(line)
                && line.wavelength_nm() >= O2_B_BAND_MIN_WAVELENGTH_NM
                && line.wavelength_nm() <= O2_B_BAND_MAX_WAVELENGTH_NM
        })
        .filter_map(|line| emission_band_line_from_optical_line(line).transpose())
        .collect::<Result<Vec<_>>>()?;

    if lines.is_empty() {
        return Ok(None);
    }

    lines.sort_by(|lhs, rhs| lhs.wavelength_nm.partial_cmp(&rhs.wavelength_nm).unwrap());

    EmissionBand::new(
        "oxygen_b_band",
        "O2(b, v=1)",
        "O2(X)",
        O2_B1_X0_EINSTEIN_A_S,
        lines,
    )
    .map(Some)
}

pub fn oxygen_gamma_band_from_hitran(db: &OpticalLineDB) -> Result<Option<EmissionBand>> {
    let mut lines: Vec<EmissionBandLine> = db
        .lines
        .iter()
        .filter(|line| {
            line_matches_o2_gamma_band_vibrational_sequence(line)
                && line.wavelength_nm() >= O2_GAMMA_BAND_MIN_WAVELENGTH_NM
                && line.wavelength_nm() <= O2_GAMMA_BAND_MAX_WAVELENGTH_NM
        })
        .filter_map(|line| emission_band_line_from_optical_line(line).transpose())
        .collect::<Result<Vec<_>>>()?;

    if lines.is_empty() {
        return Ok(None);
    }

    lines.sort_by(|lhs, rhs| lhs.wavelength_nm.partial_cmp(&rhs.wavelength_nm).unwrap());

    EmissionBand::new(
        "oxygen_gamma_band",
        "O2(b, v=2)",
        "O2(X)",
        O2_B2_X0_EINSTEIN_A_S,
        lines,
    )
    .map(Some)
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

fn line_matches_o2_gamma_band_vibrational_sequence(line: &OpticalLine) -> bool {
    let upper_tokens: Vec<&str> = line.upper_quanta.split_whitespace().collect();
    let lower_tokens: Vec<&str> = line.lower_quanta.split_whitespace().collect();

    matches!(
        (upper_tokens.as_slice(), lower_tokens.as_slice()),
        (["b", "2"], ["X", "0"])
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalized_rotational_weight_derivatives_match_finite_differences() {
        let mut lines = vec![
            test_line(760.0, "b 0", "X 0", 2.0),
            test_line(761.0, "b 0", "X 0", 6.0),
            test_line(762.0, "b 1", "X 1", 3.0),
            test_line(763.0, "b 1", "X 1", 8.0),
        ];
        for (i, line) in lines.iter_mut().enumerate() {
            line.lower_energy = [10.0, 800.0, 1500.0, 1900.0][i];
        }
        let band = oxygen_a_band_from_hitran(&OpticalLineDB { lines }).unwrap();
        let temperatures = array![100.0, 220.0, 400.0];
        for model in [
            AEmissionLineWeightModel::EinsteinABranching,
            AEmissionLineWeightModel::HitranLineStrength,
        ] {
            let (weights, derivative) = oxygen_a_band_lte_line_weights_with_temperature_derivative(
                &band,
                temperatures.view(),
                model,
            )
            .unwrap();
            let above =
                oxygen_a_band_lte_line_weights(&band, (&temperatures + 0.001).view(), model)
                    .unwrap();
            let below =
                oxygen_a_band_lte_line_weights(&band, (&temperatures - 0.001).view(), model)
                    .unwrap();
            let numeric = (above - below) / 0.002;
            for (analytic, numeric) in derivative.iter().zip(numeric.iter()) {
                assert!((analytic - numeric).abs() < 1e-10);
            }
            for indices in upper_vibrational_state_groups(&band) {
                for alt in 0..temperatures.len() {
                    assert!(
                        (indices.iter().map(|&i| weights[[alt, i]]).sum::<f64>() - 1.0).abs()
                            < 1e-14
                    );
                    assert!(
                        indices
                            .iter()
                            .map(|&i| derivative[[alt, i]])
                            .sum::<f64>()
                            .abs()
                            < 1e-14
                    );
                }
            }
        }
    }

    #[test]
    fn rotational_weights_avoid_common_vibrational_underflow() {
        let db = OpticalLineDB {
            lines: vec![test_line(760.0, "b 0", "X 0", 2.0)],
        };
        let band = oxygen_a_band_from_hitran(&db).unwrap();
        let (weights, derivative) = oxygen_a_band_lte_line_weights_with_temperature_derivative(
            &band,
            array![10.0].view(),
            AEmissionLineWeightModel::EinsteinABranching,
        )
        .unwrap();
        assert_eq!(weights[[0, 0]], 1.0);
        assert_eq!(derivative[[0, 0]], 0.0);
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

        let band = oxygen_a_band_from_hitran(&db).unwrap();

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
                test_line(775.0, "b 0", "X 0", 4.0),
                test_line(777.0, "b 0", "X 0", 4.0),
            ],
        };

        let band = oxygen_a_band_from_hitran(&db).unwrap();

        assert_eq!(band.lines.len(), 4);
        assert_eq!(band.lines[1].upper_vibrational_state, "O2(b, v=1)");
        assert_eq!(band.lines[1].lower_vibrational_state, "O2(X, v=1)");
        assert!(band.lines.iter().any(|line| line.wavelength_nm == 775.0));
        assert!(!band.lines.iter().any(|line| line.wavelength_nm == 777.0));

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
                test_line(674.0, "b 1", "X 0", 2.0),
                test_line(675.0, "b 1", "X 0", 2.0),
                test_line(688.0, "b 1", "X 0", 2.0),
                test_line(689.0, "b 1", "X 1", 8.0),
                test_line(690.0, "b 2", "X 1", 6.0),
                test_line(705.0, "b 1", "X 0", 2.0),
                test_line(706.0, "b 1", "X 0", 2.0),
                test_line(760.0, "b 0", "X 0", 4.0),
            ],
        };

        let band = oxygen_b_band_from_hitran(&db).unwrap().unwrap();

        assert_eq!(band.lines.len(), 3);
        assert!(band.lines.iter().all(|line| {
            line.upper_vibrational_state == "O2(b, v=1)" && line.lower_vibrational_state == "O2(X)"
        }));
        assert!(band.lines.iter().any(|line| line.wavelength_nm == 675.0));
        assert!(band.lines.iter().any(|line| line.wavelength_nm == 705.0));
        assert!(!band.lines.iter().any(|line| line.wavelength_nm == 674.0));
        assert!(!band.lines.iter().any(|line| line.wavelength_nm == 706.0));
        assert!(
            (band
                .lines
                .iter()
                .map(|line| line.relative_weight)
                .sum::<f64>()
                - 1.0)
                .abs()
                < 1.0e-12
        );
    }

    #[test]
    fn oxygen_b_band_missing_b1_population_contributes_zero() {
        let db = OpticalLineDB {
            lines: vec![test_line(688.0, "b 1", "X 0", 2.0)],
        };
        let band = oxygen_b_band_from_hitran(&db).unwrap().unwrap();

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

    /// Line indices grouped by upper vibrational state, in first-seen order.
    fn upper_vibrational_state_groups(band: &EmissionBand) -> Vec<Vec<usize>> {
        let mut states: Vec<&str> = Vec::new();
        for line in &band.lines {
            if !states.contains(&line.upper_vibrational_state.as_str()) {
                states.push(&line.upper_vibrational_state);
            }
        }
        states
            .iter()
            .map(|state| {
                band.lines
                    .iter()
                    .enumerate()
                    .filter_map(|(i, line)| (line.upper_vibrational_state == *state).then_some(i))
                    .collect()
            })
            .collect()
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
}
