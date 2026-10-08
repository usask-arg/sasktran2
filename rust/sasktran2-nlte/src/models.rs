use crate::prelude::*;

use super::types::PhotoReaction;

pub const LYMAN_ALPHA_WAVELENGTH_NM: f64 = 121.567;
pub const LYMAN_ALPHA_TOA_RATE_S: f64 = 3.40e-9;
pub const LYMAN_ALPHA_O1D_QUANTUM_YIELD: f64 = 0.53;
pub const LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S: f64 = 3.2e15;
pub const O2_LYMAN_ALPHA_EFFECTIVE_CROSS_SECTION_M2: f64 =
    LYMAN_ALPHA_TOA_RATE_S / LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S;

pub fn wavelength_bin_widths(wavelength_nm: &[f64]) -> Result<Vec<f64>> {
    if wavelength_nm.len() < 2 {
        return Err(anyhow!(
            "Need at least two wavelength points to integrate photolysis rates"
        ));
    }

    if wavelength_nm.iter().any(|wl| !wl.is_finite()) {
        return Err(anyhow!("Wavelength grid contains non-finite values"));
    }

    let mut delta_wavelength = vec![0.0_f64; wavelength_nm.len()];
    for j in 0..wavelength_nm.len() {
        let d = if j == 0 {
            (wavelength_nm[1] - wavelength_nm[0]).abs()
        } else if j + 1 == wavelength_nm.len() {
            (wavelength_nm[wavelength_nm.len() - 1] - wavelength_nm[wavelength_nm.len() - 2]).abs()
        } else {
            0.5 * (wavelength_nm[j + 1] - wavelength_nm[j - 1]).abs()
        };
        delta_wavelength[j] = d;
    }

    Ok(delta_wavelength)
}

pub fn calculate_photolysis_rate(
    reaction: &PhotoReaction,
    wavelength_nm: &[f64],
    actinic_flux: ArrayView2<'_, f64>,
    cross_section: ArrayView2<'_, f64>,
) -> Result<Array1<f64>> {
    if actinic_flux.shape() != cross_section.shape() {
        return Err(anyhow!(
            "Actinic flux shape {:?} does not match cross-section shape {:?}",
            actinic_flux.shape(),
            cross_section.shape()
        ));
    }

    if actinic_flux.nrows() != wavelength_nm.len() {
        return Err(anyhow!(
            "Wavelength grid has {} points but actinic flux has {} wavelength rows",
            wavelength_nm.len(),
            actinic_flux.nrows()
        ));
    }

    let q = reaction.quantum_yield.unwrap_or(1.0);

    if let Some(line_center_nm) = reaction.line_center_nm {
        let flux_at_line = interpolate_spectral_profiles(
            wavelength_nm,
            actinic_flux,
            line_center_nm,
            "actinic flux",
        )?;
        let xs_at_line = if let Some(cross_section_m2) = reaction.line_effective_cross_section_m2 {
            Array1::from_elem(actinic_flux.ncols(), cross_section_m2)
        } else {
            interpolate_spectral_profiles(
                wavelength_nm,
                cross_section,
                line_center_nm,
                "cross section",
            )?
        };

        let photolysis_rate = flux_at_line
            .iter()
            .zip(xs_at_line.iter())
            .map(|(flux, xs)| flux.max(0.0) * xs.max(0.0))
            .collect();

        return Ok(apply_photolysis_rate_scale(reaction, photolysis_rate, q));
    }

    let delta_wavelength = wavelength_bin_widths(wavelength_nm)?;
    let mut photolysis_rate = Array1::<f64>::zeros(actinic_flux.ncols());
    let band_limits = reaction.wavelength_range_nm;

    for (j, wl) in wavelength_nm.iter().enumerate() {
        if let Some((min_nm, max_nm)) = band_limits
            && (*wl < min_nm || *wl > max_nm)
        {
            continue;
        }

        for k in 0..actinic_flux.ncols() {
            let flux = actinic_flux[[j, k]].max(0.0);
            let xs = cross_section[[j, k]].max(0.0);
            photolysis_rate[k] += flux * xs * delta_wavelength[j];
        }
    }

    Ok(apply_photolysis_rate_scale(reaction, photolysis_rate, q))
}

fn apply_photolysis_rate_scale(
    reaction: &PhotoReaction,
    mut photolysis_rate: Array1<f64>,
    quantum_yield: f64,
) -> Array1<f64> {
    if reaction.toa_rate_constant > 0.0
        && let Some(reference_rate) = photolysis_rate.last().copied()
        && reference_rate.is_finite()
        && reference_rate > 0.0
    {
        let scale = reaction.toa_rate_constant / reference_rate;
        photolysis_rate.mapv_inplace(|rate| quantum_yield * rate * scale);
    } else {
        photolysis_rate.mapv_inplace(|rate| quantum_yield * rate);
    }

    photolysis_rate
}

fn interpolate_spectral_profiles(
    wavelength_nm: &[f64],
    values: ArrayView2<'_, f64>,
    target_nm: f64,
    quantity_name: &str,
) -> Result<Array1<f64>> {
    if wavelength_nm.is_empty() {
        return Err(anyhow!(
            "Cannot interpolate {} on an empty grid",
            quantity_name
        ));
    }

    if !target_nm.is_finite() {
        return Err(anyhow!(
            "Cannot interpolate {} to non-finite wavelength {}",
            quantity_name,
            target_nm
        ));
    }

    let first = wavelength_nm[0];
    let last = wavelength_nm[wavelength_nm.len() - 1];
    if target_nm < first || target_nm > last {
        return Err(anyhow!(
            "Cannot calculate line photolysis at {} nm because the wavelength grid spans {} to {} nm",
            target_nm,
            first,
            last
        ));
    }

    if wavelength_nm.len() == 1 {
        if (target_nm - first).abs() <= f64::EPSILON {
            return Ok(values.row(0).to_owned());
        }
        return Err(anyhow!(
            "Cannot interpolate {} to {} nm with only one wavelength point",
            quantity_name,
            target_nm
        ));
    }

    for j in 0..(wavelength_nm.len() - 1) {
        let wl0 = wavelength_nm[j];
        let wl1 = wavelength_nm[j + 1];

        if wl1 <= wl0 {
            return Err(anyhow!(
                "Wavelength grid must be strictly increasing for line interpolation"
            ));
        }

        if (target_nm - wl0).abs() <= f64::EPSILON {
            return Ok(values.row(j).to_owned());
        }

        if target_nm >= wl0 && target_nm <= wl1 {
            let weight = (target_nm - wl0) / (wl1 - wl0);
            let lower = values.row(j);
            let upper = values.row(j + 1);
            return Ok((&lower * (1.0 - weight) + &upper * weight).to_owned());
        }
    }

    if (target_nm - last).abs() <= f64::EPSILON {
        return Ok(values.row(wavelength_nm.len() - 1).to_owned());
    }

    Err(anyhow!(
        "Could not interpolate {} to {} nm",
        quantity_name,
        target_nm
    ))
}

/// How the earlier photochem module derived each photolysis and
/// photoexcitation rate of the `oxygen_yankovsky` mechanism from an actinic
/// flux spectrum: integrated rates scaled to fixed top-of-atmosphere values.
/// Kept until `sasktran2.photolysis` computes these rates from cross sections.
pub struct Yankovsky {
    pub photo_reactions: Vec<PhotoReaction>,
}

impl Default for Yankovsky {
    fn default() -> Self {
        Self::new()
    }
}

impl Yankovsky {
    pub fn new() -> Self {
        let mut photo_reactions = vec![
            "O2 + hv(SRC) -> O(3P) + O(1D)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_quantum_yield(1.0)
                .with_toa_rate_constant(2.60e-6)
                .with_wavelength_range_nm(130.0, 202.0),
            "O2 + hv(lyman-alpha) -> O(3P) + O(1D)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_quantum_yield(LYMAN_ALPHA_O1D_QUANTUM_YIELD)
                .with_toa_rate_constant(LYMAN_ALPHA_TOA_RATE_S)
                .with_line_center_nm(LYMAN_ALPHA_WAVELENGTH_NM)
                .with_line_effective_cross_section_m2(O2_LYMAN_ALPHA_EFFECTIVE_CROSS_SECTION_M2),
            "O3 + hv -> O2(a, v=5) + O(1D)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_quantum_yield(0.045)
                .with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=4) + O(1D)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_quantum_yield(0.072)
                .with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=3) + O(1D)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_quantum_yield(0.072)
                .with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=2) + O(1D)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_quantum_yield(0.135)
                .with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=1) + O(1D)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_quantum_yield(0.135)
                .with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=0) + O(1D)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_quantum_yield(0.441)
                .with_toa_rate_constant(8.0e-3),
            "O2 + hv(762_nm_band) -> O2(b, v=0)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_toa_rate_constant(5.35e-9)
                .with_band_center_nm(762.0, 10.0),
            "O2 + hv(689_nm_band) -> O2(b, v=1)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_toa_rate_constant(2.94e-10)
                .with_band_center_nm(689.0, 10.0),
            "O2 + hv(629_nm_band) -> O2(b, v=2)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_toa_rate_constant(7.94e-12)
                .with_band_center_nm(629.0, 10.0),
            "O2 + hv(1.27_um_band) -> O2(a, v=0)"
                .parse::<PhotoReaction>()
                .unwrap()
                .with_toa_rate_constant(1.54e-10)
                .with_band_center_nm(1270.0, 10.0),
        ];

        // Table 1 branch: O3 + hv -> O2(X, v=1..35) + O(3P).
        // Existing O3(a, v) + O(1D) branches sum to 0.90, so allocate the
        // remaining 0.10 uniformly across O2(X, v=1..35) for now.
        for v in 1..=35 {
            photo_reactions.push(
                format!("O3 + hv -> O2(X, v={v}) + O(3P)")
                    .parse::<PhotoReaction>()
                    .unwrap()
                    .with_quantum_yield(0.1 / 35.0)
                    .with_toa_rate_constant(8.0e-3),
            );
        }

        Self { photo_reactions }
    }
}

impl Yankovsky {
    /// Names under which the `oxygen_yankovsky` mechanism expects each of
    /// `photo_reactions`' rates, in the same order.
    pub fn photolysis_rate_names(&self) -> Vec<String> {
        self.photo_reactions
            .iter()
            .map(|reaction| match reaction.excitation_band.as_deref() {
                Some("SRC") => "J_O2_SRC".to_string(),
                Some("lyman-alpha") => "J_O2_LYA".to_string(),
                Some("762_nm_band") => "J_O2_EXC_B0".to_string(),
                Some("689_nm_band") => "J_O2_EXC_B1".to_string(),
                Some("629_nm_band") => "J_O2_EXC_B2".to_string(),
                Some("1.27_um_band") => "J_O2_EXC_A0".to_string(),
                _ => {
                    let o2 = &reaction.products[0];
                    format!(
                        "J_O3_{}{}",
                        o2.electronic_level.to_uppercase(),
                        o2.vibrational_level
                    )
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::{
        LYMAN_ALPHA_O1D_QUANTUM_YIELD, LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S, LYMAN_ALPHA_TOA_RATE_S,
        LYMAN_ALPHA_WAVELENGTH_NM, O2_LYMAN_ALPHA_EFFECTIVE_CROSS_SECTION_M2, Yankovsky,
        calculate_photolysis_rate, wavelength_bin_widths,
    };
    use crate::mechanism::Mechanism;
    use crate::types::PhotoReaction;
    use ndarray::array;

    #[test]
    fn photolysis_rate_names_match_the_bundled_mechanism() {
        let names = Yankovsky::new().photolysis_rate_names();
        let mechanism = Mechanism::bundled("oxygen_yankovsky").unwrap();

        let mut sorted_names = names.clone();
        sorted_names.sort();
        let mut inputs = mechanism.rate_inputs().to_vec();
        inputs.sort();
        assert_eq!(sorted_names, inputs);
        assert_eq!(names[0], "J_O2_SRC");
        assert!(names.contains(&"J_O3_A0".to_string()));
        assert!(names.contains(&"J_O3_X35".to_string()));
    }

    #[test]
    fn yankovsky_lyman_alpha_reaction_is_enabled_as_line() {
        let model = Yankovsky::new();
        let reaction = model
            .photo_reactions
            .iter()
            .find(|r| r.excitation_band.as_deref() == Some("lyman-alpha"))
            .expect("Yankovsky model should include a Lyman-alpha reaction");

        assert_eq!(reaction.line_center_nm, Some(LYMAN_ALPHA_WAVELENGTH_NM));
        assert_eq!(
            reaction.line_effective_cross_section_m2,
            Some(O2_LYMAN_ALPHA_EFFECTIVE_CROSS_SECTION_M2)
        );
        assert_eq!(reaction.wavelength_range_nm, None);
        assert_eq!(reaction.quantum_yield, Some(LYMAN_ALPHA_O1D_QUANTUM_YIELD));
    }

    #[test]
    fn lyman_alpha_effective_cross_section_matches_toa_rate_scale() {
        let rate = O2_LYMAN_ALPHA_EFFECTIVE_CROSS_SECTION_M2 * LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S;
        assert!((rate - LYMAN_ALPHA_TOA_RATE_S).abs() < 1.0e-18);

        let cross_section_cm2 = O2_LYMAN_ALPHA_EFFECTIVE_CROSS_SECTION_M2 * 1.0e4;
        assert!((1.0e-20..=1.2e-20).contains(&cross_section_cm2));
    }

    #[test]
    fn wavelength_bin_widths_match_existing_central_difference_rule() {
        let widths = wavelength_bin_widths(&[100.0, 102.0, 106.0, 116.0]).unwrap();
        assert_eq!(widths, vec![2.0, 3.0, 7.0, 10.0]);
    }

    #[test]
    fn continuum_photolysis_rate_integrates_over_selected_band() {
        let reaction = "O2 + hv(SRC) -> O + O"
            .parse::<PhotoReaction>()
            .unwrap()
            .with_quantum_yield(0.5)
            .with_wavelength_range_nm(102.0, 106.0);

        let wavelength = [100.0, 102.0, 106.0, 116.0];
        let actinic_flux = array![[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]];
        let cross_section = array![[10.0, 1.0], [20.0, 2.0], [30.0, 3.0], [40.0, 4.0]];

        let rate = calculate_photolysis_rate(
            &reaction,
            &wavelength,
            actinic_flux.view(),
            cross_section.view(),
        )
        .unwrap();

        assert_eq!(rate.len(), 2);
        assert_eq!(rate[0], 0.5 * (2.0 * 20.0 * 3.0 + 3.0 * 30.0 * 7.0));
        assert_eq!(rate[1], 0.5 * (20.0 * 2.0 * 3.0 + 30.0 * 3.0 * 7.0));
    }

    #[test]
    fn continuum_photolysis_rate_scales_to_toa_rate_after_quantum_yield() {
        let reaction = "O2 + hv(689_nm_band) -> O2(b, v=1)"
            .parse::<PhotoReaction>()
            .unwrap()
            .with_quantum_yield(0.25)
            .with_toa_rate_constant(2.0e-9)
            .with_wavelength_range_nm(100.0, 102.0);

        let wavelength = [100.0, 101.0, 102.0];
        let actinic_flux = array![[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]];
        let cross_section = array![[5.0, 5.0], [5.0, 5.0], [5.0, 5.0]];

        let rate = calculate_photolysis_rate(
            &reaction,
            &wavelength,
            actinic_flux.view(),
            cross_section.view(),
        )
        .unwrap();

        assert_eq!(rate.len(), 2);
        assert!((rate[1] - 0.25 * 2.0e-9).abs() < 1.0e-24);
        assert!((rate[0] - rate[1] * 0.1).abs() < 1.0e-24);
    }

    #[test]
    fn continuum_photolysis_rate_clamps_negative_inputs() {
        let reaction = "O3 + hv -> O2 + O"
            .parse::<PhotoReaction>()
            .unwrap()
            .with_quantum_yield(2.0);

        let wavelength = [100.0, 101.0];
        let actinic_flux = array![[-1.0], [4.0]];
        let cross_section = array![[5.0], [-2.0]];

        let rate = calculate_photolysis_rate(
            &reaction,
            &wavelength,
            actinic_flux.view(),
            cross_section.view(),
        )
        .unwrap();

        assert_eq!(rate[0], 0.0);
    }

    #[test]
    fn line_photolysis_rate_interpolates_flux_and_cross_section() {
        let reaction = "O2 + hv(lyman-alpha) -> O + O"
            .parse::<PhotoReaction>()
            .unwrap()
            .with_quantum_yield(0.25)
            .with_line_center_nm(121.5);

        let wavelength = [121.0, 122.0];
        let actinic_flux = array![[2.0, 4.0], [6.0, 12.0]];
        let cross_section = array![[10.0, 20.0], [30.0, 60.0]];

        let rate = calculate_photolysis_rate(
            &reaction,
            &wavelength,
            actinic_flux.view(),
            cross_section.view(),
        )
        .unwrap();

        assert_eq!(rate.len(), 2);
        assert_eq!(rate[0], 0.25 * 4.0 * 20.0);
        assert_eq!(rate[1], 0.25 * 8.0 * 40.0);
    }

    #[test]
    fn line_photolysis_rate_can_use_effective_cross_section() {
        let reaction = "O2 + hv(lyman-alpha) -> O + O"
            .parse::<PhotoReaction>()
            .unwrap()
            .with_quantum_yield(0.5)
            .with_line_center_nm(121.5)
            .with_line_effective_cross_section_m2(2.0e-24);

        let wavelength = [121.0, 122.0];
        let actinic_flux = array![[2.0e15, 4.0e15], [6.0e15, 12.0e15]];
        let cross_section = array![[10.0, 20.0], [30.0, 60.0]];

        let rate = calculate_photolysis_rate(
            &reaction,
            &wavelength,
            actinic_flux.view(),
            cross_section.view(),
        )
        .unwrap();

        assert_eq!(rate.len(), 2);
        assert_eq!(rate[0], 0.5 * 4.0e15 * 2.0e-24);
        assert_eq!(rate[1], 0.5 * 8.0e15 * 2.0e-24);
    }

    #[test]
    fn line_photolysis_rate_uses_exact_grid_point_without_smoothing() {
        let reaction = "O2 + hv(lyman-alpha) -> O + O"
            .parse::<PhotoReaction>()
            .unwrap()
            .with_line_center_nm(121.567);

        let wavelength = [121.0, 121.567, 122.0];
        let actinic_flux = array![[1.0], [3.0], [100.0]];
        let cross_section = array![[10.0], [20.0], [1000.0]];

        let rate = calculate_photolysis_rate(
            &reaction,
            &wavelength,
            actinic_flux.view(),
            cross_section.view(),
        )
        .unwrap();

        assert_eq!(rate[0], 60.0);
    }

    #[test]
    fn line_photolysis_rate_requires_wavelength_coverage() {
        let reaction = "O2 + hv(lyman-alpha) -> O + O"
            .parse::<PhotoReaction>()
            .unwrap()
            .with_line_center_nm(121.567);

        let wavelength = [130.0, 140.0];
        let actinic_flux = array![[1.0], [2.0]];
        let cross_section = array![[3.0], [4.0]];

        let err = calculate_photolysis_rate(
            &reaction,
            &wavelength,
            actinic_flux.view(),
            cross_section.view(),
        )
        .expect_err("missing Lyman-alpha wavelength coverage should fail");

        assert!(
            err.to_string()
                .contains("wavelength grid spans 130 to 140 nm")
        );
    }

    #[test]
    fn photolysis_rate_rejects_shape_mismatch() {
        let reaction = "O3 + hv -> O2 + O".parse::<PhotoReaction>().unwrap();
        let wavelength = [100.0, 101.0];
        let actinic_flux = array![[1.0], [2.0]];
        let cross_section = array![[1.0, 2.0], [3.0, 4.0]];

        let err = calculate_photolysis_rate(
            &reaction,
            &wavelength,
            actinic_flux.view(),
            cross_section.view(),
        )
        .expect_err("shape mismatch should fail");

        assert!(err.to_string().contains("does not match"));
    }
}
