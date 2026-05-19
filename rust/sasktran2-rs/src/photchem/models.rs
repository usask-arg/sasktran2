use crate::bindings::lapack::dgesv;
use crate::prelude::*;

use super::types::{ChemicalReaction, Molecule, MoleculeMap, PhotoReaction};

pub const LYMAN_ALPHA_WAVELENGTH_NM: f64 = 121.567;
pub const LYMAN_ALPHA_TOA_RATE_S: f64 = 3.40e-9;
pub const LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S: f64 = 3.2e15;
pub const O2_LYMAN_ALPHA_EFFECTIVE_CROSS_SECTION_M2: f64 =
    LYMAN_ALPHA_TOA_RATE_S / LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S;

pub trait PhotochemicalModel {
    fn molecules(&self) -> Vec<Molecule>;
    fn solve(
        &self,
        temperature: f64,
        reactions: &[ChemicalReaction],
        photo_reactions: &[PhotoReaction],
        photolysis_rate: &[f64],
        densities: &HashMap<String, f64>,
    ) -> Result<HashMap<String, f64>> {
        let mol_map = MoleculeMap::new(self.molecules().as_slice(), densities);

        let n = mol_map.state_size();

        let mut a_matrix = Array2::<f64>::zeros((n, n));

        // Negative of sources
        let mut sources = Array1::<f64>::zeros(n);

        for reaction in reactions {
            match reaction.reactants.len() {
                1 => {
                    let einstein_coefficient =
                        reaction.einstein_coefficient.as_ref().ok_or_else(|| {
                            anyhow!("Einstein coefficient missing for unimolecular reaction")
                        })?;

                    let branch = 1.0; // branching ratio is included in the rate constant

                    let rate = einstein_coefficient(temperature) * branch;

                    let reactant = &reaction.reactants[0];

                    let reactant_str = String::try_from(reactant.clone())
                        .unwrap_or_else(|_| format!("{:?}", reactant));

                    if mol_map.is_in_state(reactant) {
                        let reactant_index = mol_map.index(reactant).ok_or_else(|| {
                            anyhow!(
                                "Reactant '{}' not found in molecule map (unimolecular reaction)",
                                reactant_str
                            )
                        })?;

                        // One loss per reaction event
                        a_matrix[[reactant_index, reactant_index]] -= rate;

                        // Gains to tracked products only
                        for product in &reaction.products {
                            if mol_map.is_in_state(product) {
                                let product_str = String::try_from(product.clone())
                                    .unwrap_or_else(|_| format!("{:?}", product));
                                let product_index = mol_map.index(product).ok_or_else(|| {
                                    anyhow!(
                                        "Product '{}' not found in molecule map (reactant: '{}')",
                                        product_str,
                                        reactant_str
                                    )
                                })?;
                                a_matrix[[product_index, reactant_index]] += rate;
                            }
                        }
                    } else {
                        // Reactant is a background species — density is fixed, so production
                        // of state products is a constant source term.
                        let reactant_density = densities
                            .get(&reactant_str)
                            .or_else(|| densities.get(&reactant.base_type.to_string()))
                            .ok_or_else(|| {
                                anyhow!(
                                    "Density not provided for background reactant '{}'",
                                    reactant_str
                                )
                            })?
                            / 1.0e6; // convert from cm^-3 to m^-3

                        for product in &reaction.products {
                            if mol_map.is_in_state(product) {
                                let product_str = String::try_from(product.clone())
                                    .unwrap_or_else(|_| format!("{:?}", product));
                                let product_index = mol_map.index(product)
                                    .ok_or_else(|| anyhow!("Product '{}' not found in molecule map (background reactant: '{}')", product_str, reactant_str))?;
                                sources[product_index] -= rate * reactant_density;
                            }
                        }
                    }
                }
                2 => {
                    let k =
                        reaction.rate_constant.as_ref().ok_or_else(|| {
                            anyhow!("Rate constant missing for bimolecular reaction")
                        })?(temperature);

                    let branch = reaction.quantum_yield.unwrap_or(1.0);

                    let source = &reaction.reactants[0];
                    let collider = &reaction.reactants[1];

                    let source_str = String::try_from(source.clone())
                        .unwrap_or_else(|_| format!("{:?}", source));
                    let collider_str = String::try_from(collider.clone())
                        .unwrap_or_else(|_| format!("{:?}", collider));

                    if mol_map.is_in_state(collider) {
                        return Err(anyhow!(
                            "Collider '{}' is in state for bimolecular reaction with '{}' and would make the system nonlinear; expected collider in background densities",
                            collider_str,
                            source_str
                        ));
                    }

                    let collider_density = densities
                        .get(&collider_str)
                        .or_else(|| densities.get(&collider.base_type.to_string()))
                        .ok_or_else(|| {
                            anyhow!(
                                "Density not provided for collider '{}' (reactant: '{}')",
                                collider_str,
                                source_str
                            )
                        })?
                        / 1.0e6; // convert from cm^-3 to m^-3
                    let rate = k * collider_density * branch;

                    if mol_map.is_in_state(source) {
                        let source_index = mol_map.index(source)
                            .ok_or_else(|| anyhow!("Reactant '{}' not found in molecule map (bimolecular reaction with collider '{}')", source_str, collider_str))?;

                        // One loss per reaction event
                        a_matrix[[source_index, source_index]] -= rate;

                        // Gains to tracked products only
                        for product in &reaction.products {
                            if mol_map.is_in_state(product) {
                                let product_str = String::try_from(product.clone())
                                    .unwrap_or_else(|_| format!("{:?}", product));
                                let product_index = mol_map.index(product)
                                    .ok_or_else(|| anyhow!("Product '{}' not found in molecule map (reactants: '{}' + '{}')", product_str, source_str, collider_str))?;
                                a_matrix[[product_index, source_index]] += rate;
                            }
                        }
                    } else {
                        let source_density = densities
                            .get(&source_str)
                            .or_else(|| densities.get(&source.base_type.to_string()))
                            .ok_or_else(|| anyhow!("Density not provided for background reactant '{}' (collider: '{}')", source_str, collider_str))? / 1.0e6; // convert from cm^-3 to m^-3

                        for product in &reaction.products {
                            if mol_map.is_in_state(product) {
                                let product_str = String::try_from(product.clone())
                                    .unwrap_or_else(|_| format!("{:?}", product));
                                let product_index = mol_map.index(product)
                                    .ok_or_else(|| anyhow!("Product '{}' not found in molecule map (background reactant: '{}' + '{}')", product_str, source_str, collider_str))?;
                                sources[product_index] -= rate * source_density;
                            }
                        }
                    }
                }
                _ => {
                    return Err(anyhow!(
                        "Unsupported reaction with {} reactants",
                        reaction.reactants.len()
                    ));
                }
            }
        }

        for (reaction, j) in photo_reactions.iter().zip(photolysis_rate.iter()) {
            let photo_reactant_str = String::try_from(reaction.in_molecule.clone())
                .unwrap_or_else(|_| format!("{:?}", reaction.in_molecule));
            let reactant_density = densities
                .get(&photo_reactant_str)
                .or_else(|| densities.get(&reaction.in_molecule.base_type.to_string()))
                .ok_or_else(|| {
                    anyhow!(
                        "Density not provided for photo-reactant '{}'",
                        photo_reactant_str
                    )
                })?;
            let production = j * reactant_density;

            for product in &reaction.products {
                if mol_map.is_in_state(product) {
                    let product_str = String::try_from(product.clone())
                        .unwrap_or_else(|_| format!("{:?}", product));
                    let product_index = mol_map.index(product).ok_or_else(|| {
                        anyhow!(
                            "Product '{}' not found in molecule map (photo-reaction on '{}')",
                            product_str,
                            photo_reactant_str
                        )
                    })?;
                    sources[product_index] -= production;
                }
            }
        }

        let format_matrix_preview = |name: &str, matrix: &Array2<f64>| -> String {
            let (rows, cols) = matrix.dim();
            let show_rows = rows.min(5);
            let show_cols = cols.min(5);

            let mut preview = format!("{} shape=({}, {})\n", name, rows, cols);
            for i in 0..show_rows {
                let mut row_parts = Vec::with_capacity(show_cols);
                for j in 0..show_cols {
                    row_parts.push(format!("{:.3e}", matrix[[i, j]]));
                }
                preview.push_str(&format!("  [{}]\n", row_parts.join(", ")));
            }
            if rows > show_rows || cols > show_cols {
                preview.push_str("  ...\n");
            }

            preview
        };

        let format_vector_preview = |name: &str, vector: &Array1<f64>| -> String {
            let len = vector.len();
            let show = len.min(10);

            let mut parts = Vec::with_capacity(show);
            for i in 0..show {
                parts.push(format!("{:.3e}", vector[i]));
            }

            if len > show {
                format!("{} len={} [{} ...]", name, len, parts.join(", "))
            } else {
                format!("{} len={} [{}]", name, len, parts.join(", "))
            }
        };

        let format_state_mapping_preview = |names: &[String]| -> String {
            let show = names.len().min(25);
            let mut lines = vec![format!(
                "State index map (showing {} of {}):",
                show,
                names.len()
            )];
            for (idx, name) in names.iter().take(show).enumerate() {
                lines.push(format!("  [{}] {}", idx, name));
            }
            if names.len() > show {
                lines.push("  ...".to_string());
            }
            format!("{}\n", lines.join("\n"))
        };

        let format_problem_row_preview = |a: &Array2<f64>, names: &[String]| -> String {
            let diag_tol = 1.0e-30;
            let row_sum_tol = 1.0e-25;
            let mut rows = Vec::new();
            let n = a.nrows();

            for i in 0..n {
                let diag = a[[i, i]].abs();
                let row_sum: f64 = a.row(i).iter().map(|v| v.abs()).sum();
                if diag <= diag_tol || row_sum <= row_sum_tol {
                    let species = names
                        .get(i)
                        .cloned()
                        .unwrap_or_else(|| "<unknown-state>".to_string());
                    rows.push(format!(
                        "  row {} [{}]: |diag|={:.3e}, row_abs_sum={:.3e}",
                        i, species, diag, row_sum
                    ));
                }
            }

            if rows.is_empty() {
                "Potentially problematic rows: none flagged by simple tolerances\n".to_string()
            } else {
                let show = rows.len().min(25);
                let mut out = vec![format!(
                    "Potentially problematic rows (showing {} of {}):",
                    show,
                    rows.len()
                )];
                out.extend(rows.into_iter().take(show));
                if n > show {
                    out.push("  ...".to_string());
                }
                format!("{}\n", out.join("\n"))
            }
        };

        let rhs = sources.view().insert_axis(Axis(1)).to_owned();
        let state_names = mol_map.state_index_to_molecule_names();

        let solution = dgesv(&a_matrix, &rhs).map_err(|err| {
            let a_preview = format_matrix_preview("A", &a_matrix);
            let rhs_preview = format_matrix_preview("RHS", &rhs);
            let source_preview = format_vector_preview("sources", &sources);
            let state_mapping = format_state_mapping_preview(&state_names);
            let problem_rows = format_problem_row_preview(&a_matrix, &state_names);

            anyhow!(
                "LAPACK dgesv failed while solving photochemical system: {}\n{}{}{}\n{}{}",
                err,
                a_preview,
                rhs_preview,
                source_preview,
                state_mapping,
                problem_rows
            )
        })?;

        let mut state = HashMap::new();
        for (i, name) in state_names.iter().enumerate() {
            state.insert(name.clone(), solution[[i, 0]]);
        }

        Ok(state)
    }

    fn required_photolysis_rates(&self, photo_reactions: &[PhotoReaction]) -> Vec<String> {
        let mut required_rates = Vec::new();
        for reaction in photo_reactions {
            for product in &reaction.products {
                if self
                    .molecules()
                    .iter()
                    .any(|m| m.base_type == product.base_type)
                {
                    required_rates.push(format!(
                        "J_{:?}_{}",
                        reaction.in_molecule.base_type,
                        reaction.excitation_band.as_deref().unwrap_or("")
                    ));
                    break;
                }
            }
        }
        required_rates
    }
}

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

        return Ok(flux_at_line
            .iter()
            .zip(xs_at_line.iter())
            .map(|(flux, xs)| q * flux.max(0.0) * xs.max(0.0))
            .collect());
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
            photolysis_rate[k] += q * flux * xs * delta_wavelength[j];
        }
    }

    Ok(photolysis_rate)
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

// Implementation of the O2 and O3 photochemistry models from
pub struct Yankovsky {
    pub photo_reactions: Vec<PhotoReaction>,
    pub chemical_reactions: Vec<ChemicalReaction>,
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
                .with_quantum_yield(1.0)
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

        let mut chemical_reactions = vec![
            // O(1D) deactivation reactions
            "O(1D) -> O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_einstein_coefficient(|_| 9.0e-3),
            "O(1D) + O(3P) -> O(3P) + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 4.0e-12),
            "O(1D) + O2 -> O2(b, v=1) + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| -> f64 {
                    3.2e-11 * (67.0 / temperature).exp()
                })
                .with_quantum_yield(0.40),
            "O(1D) + O2 -> O2(b, v=0) + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| -> f64 {
                    3.2e-11 * (67.0 / temperature).exp()
                })
                .with_quantum_yield(0.55),
            "O(1D) + O2 -> O2 + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| -> f64 {
                    3.2e-11 * (67.0 / temperature).exp()
                })
                .with_quantum_yield(0.05),
            "O(1D) + O2 -> O2(a, v=0) + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| -> f64 {
                    3.2e-11 * (67.0 / temperature).exp()
                })
                .with_quantum_yield(0.05),
            "O(1D) + O3 -> O2 + O2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 2.4e-10),
            "O(1D) + N2 -> N2 + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| -> f64 {
                    2.0e-11 * (107.0 / temperature).exp()
                })
                .with_quantum_yield(1.00),
            // O2(b, v) deactivation reactions
            "O2(b, v=2) -> O2(X, v=2)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_einstein_coefficient(|_| 5.4e-2),
            "O2(b, v=2) + O(3P) -> O2(b, v=1) + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 1.1e-11),
            "O2(b, v=2) + O2 -> O2(X, v=2) + O2(b, v=0)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| 1.20e-11 * (-596.0 / temperature).exp()),
            "O2(b, v=2) + N2 -> O2(b, v=1) + N2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 2e-14),
            "O2(b, v=2) + O3 -> O2 + O2 + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 2.9e-10),
            "O2(b, v=1) -> O2(X, v=1)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_einstein_coefficient(|_| 7.0e-2),
            "O2(b, v=1) + O(3P) -> O2(b, v=0) + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 4.5e-12),
            "O2(b, v=1) + O2 -> O2(X, v=1) + O2(b, v=0)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature| 4.20e-11 * (-312.0 / temperature).exp()),
            "O2(b, v=1) + N2 -> O2(b, v=0) + N2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 5.0e-13),
            "O2(b, v=1) + O3 -> O2 + O2 + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 3.0e-10),
            "O2(b, v=0) -> O2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_einstein_coefficient(|_| 7.58e-2),
            "O2(b, v=0) + O(3P) -> O2(a, v=0) + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 8.0e-14)
                .with_quantum_yield(0.75),
            "O2(b, v=0) + O(3P) -> O2 + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 8.0e-14)
                .with_quantum_yield(0.25),
            "O2(b, v=0) + O2 -> O2(a, v=0) + O2(X, v=3)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 3.9e-17)
                .with_quantum_yield(0.230),
            "O2(b, v=0) + O2 -> O2(a, v=1) + O2(X, v=2)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 3.9e-17)
                .with_quantum_yield(0.525),
            "O2(b, v=0) + O2 -> O2(a, v=2) + O2(X, v=1)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 3.9e-17)
                .with_quantum_yield(0.226),
            "O2(b, v=0) + O2 -> O2(a, v=3) + O2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 3.9e-17)
                .with_quantum_yield(0.019),
            // O2(b, 0) + N2 -> products ???
            "O2(b, v=0) + CO2 -> O2(a, v=0) + CO2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 4.2e-13),
            "O2(b, v=0) + O3 -> O2(a, v=0) + O3"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 2.2e-11)
                .with_quantum_yield(0.3),
            // O2(a, v) deactivation reactions
            "O2(a, v=0) -> O2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_einstein_coefficient(|_| 2.58e-4),
            "O2(a, v=2) + O2 -> O2(X, v=2) + O2(a, v=0)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 3.6e-11),
            "O2(a, v=1) + O2 -> O2(X, v=1) + O2(a, v=0)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 5.6e-11),
            "O2(a, v=1) + O3 -> O2 + O2 + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 4.7e-12),
            "O2(a, v=0) + O(3P) -> O2 + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 6.5e-17),
            "O2(a, v=0) + O2 -> O2(X, v=5) + O2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| 3.6e-17 * (-220.0 / temperature).exp())
                .with_quantum_yield(0.014),
            "O2(a, v=0) + O2 -> O2(X, v=4) + O2(X, v=1)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| 3.6e-17 * (-220.0 / temperature).exp())
                .with_quantum_yield(0.214),
            "O2(a, v=0) + O2 -> O2(X, v=3) + O2(X, v=2)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| 3.6e-17 * (-220.0 / temperature).exp())
                .with_quantum_yield(0.772),
            "O2(a, v=0) + O3 -> O2 + O3"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| 5.20e-11 * (-2840.0 / temperature).exp()),
            "O2(a, v=0) + N2 -> O2 + N2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 1.0e-20),
            // Energy transfer and deactivation of O2(X, v), most are v dependent
            "O2(X, v=1) + O(3P) -> O2 + O(3P)"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|_| 3.0e-12),
            "O2(X, v=1) + O2 -> O2 + O2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| 4.2e-19 * (temperature / 300.0).sqrt()),
            "O2(X, v=1) + N2 -> O2 + N2"
                .parse::<ChemicalReaction>()
                .unwrap()
                .with_rate_constant(|temperature: f64| 4.20e-19 * (temperature / 300.0).sqrt()),
        ];

        // v dependent O2(a, v) deactivation reactions
        for v in 1..=5 {
            chemical_reactions.push(
                format!("O2(a, v={v}) -> O2")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_einstein_coefficient(|_| 2.58e-4),
            );

            chemical_reactions.push(
                format!("O2(a, v={v}) + O(3P) -> O2 + O(3P)")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(|_| 1e-14),
            );
        }

        for v in 3..=5 {
            chemical_reactions.push(
                format!("O2(a, v={v}) + O2 -> O2(X, v={v}) + O2(a, v=0)")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(|_| 3.6e-11),
            );
        }

        // v dependent O2(X, v) energy transfer and deactivation reactions
        for v in 1..=30 {
            // todo: quantum yield
            chemical_reactions.push(
                format!("O3 + O(3P) -> O2(X, v={v}) + O2")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(|temperature: f64| {
                        5.60e-11 * (-1959.0 / temperature).exp()
                    }),
            )
        }

        for v in 5..=35 {
            chemical_reactions.push(
                format!("O2(X, v={v}) + O(3P) -> O2 + O(3P)")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(|temperature: f64| 5.0e-11 * (temperature / 300.0).sqrt()),
            )
        }

        for v in 2..=4 {
            chemical_reactions.push(
                format!("O2(X, v={v}) + O(3P) -> O2 + O(3P)")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(move |temperature: f64| {
                        1.1e-12 * (temperature / 300.0) * (1.0 * v as f64).exp()
                    }),
            )
        }

        for v in 2..=35 {
            let vm1 = v - 1;

            let rate_constant = match v {
                2 => |_| 2.0e-13,
                _ => |_| 2.6e-13,
            };

            chemical_reactions.push(
                format!("O2(X, v={v}) + O2 -> O2(X, v={vm1}) + O2(X, v=1)")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(rate_constant),
            );
        }

        for v in 4..=20 {
            let vm1 = v - 1;

            let rate_constant = move |_| 1.3e-12 * (-0.31 * v as f64).exp();

            chemical_reactions.push(
                format!("O2(X, v={v}) + O2 -> O2(X, v={vm1}) + O2(X, v=1)")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(rate_constant),
            );
        }

        for v in 21..=35 {
            let vm1 = v - 1;

            let rate_constant =
                move |temperature: f64| 6.0e-17 * (temperature / 300.0) * (0.2 * v as f64).exp();

            chemical_reactions.push(
                format!("O2(X, v={v}) + O2 -> O2(X, v={vm1}) + O2")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(rate_constant),
            );
        }

        // Also unsure about this one
        for v in 12..=17 {
            let vm2 = v - 2;

            let rate_constant = move |_| 3.6e-19 * (0.66 * v as f64).exp();

            chemical_reactions.push(
                format!("O2(X, v={v}) + N2 -> O2(X, v={vm2}) + N2")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(rate_constant),
            );
        }

        for v in 18..=26 {
            let vm2 = v - 2;

            let rate_constant = move |_| 4.5e-13 * (-0.173 * v as f64).exp();

            chemical_reactions.push(
                format!("O2(X, v={v}) + N2 -> O2(X, v={vm2}) + N2")
                    .parse::<ChemicalReaction>()
                    .unwrap()
                    .with_rate_constant(rate_constant),
            );
        }

        Self {
            photo_reactions,
            chemical_reactions,
        }
    }
}

impl PhotochemicalModel for Yankovsky {
    fn molecules(&self) -> Vec<Molecule> {
        let mut mols = Vec::new();
        for r in &self.photo_reactions {
            mols.extend(std::iter::once(r.in_molecule.clone()));
            mols.extend(r.products.clone());
        }
        for r in &self.chemical_reactions {
            mols.extend(r.reactants.clone());
            mols.extend(r.products.clone());
        }

        let mut unique_mols = Vec::new();
        for mol in mols {
            if !unique_mols.contains(&mol) {
                unique_mols.push(mol);
            }
        }

        unique_mols
    }
}

#[cfg(test)]
mod tests {
    use super::{
        LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S, LYMAN_ALPHA_TOA_RATE_S, LYMAN_ALPHA_WAVELENGTH_NM,
        O2_LYMAN_ALPHA_EFFECTIVE_CROSS_SECTION_M2, PhotochemicalModel, Yankovsky,
        calculate_photolysis_rate, wavelength_bin_widths,
    };
    use crate::photchem::types::{MoleculeBase, PhotoReaction};
    use ndarray::array;

    #[test]
    fn test_yankovsky() {
        let model = Yankovsky::new();

        let required_rates = model.required_photolysis_rates(&model.photo_reactions);
        assert!(!required_rates.is_empty());
    }

    #[test]
    fn test_molecules() {
        let model = Yankovsky::new();
        let mols = model.molecules();
        assert!(mols.iter().any(|m| m.base_type == MoleculeBase::O2));
        assert!(mols.iter().any(|m| m.base_type == MoleculeBase::O3));
        assert!(mols.iter().any(|m| m.base_type == MoleculeBase::O));
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
        assert_eq!(reaction.quantum_yield, Some(1.0));
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
