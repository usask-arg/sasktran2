use std::any;

use crate::prelude::*;

use super::types::{ChemicalReaction, Molecule, PhotoReaction, MoleculeBase, MoleculeMap};


pub trait PhotochemicalModel {
    fn molecules(&self) -> Vec<Molecule>;
    fn solve(&self, temperature: f64, reactions: &[ChemicalReaction], photo_reactions: &[PhotoReaction], photolysis_rate: &[f64]) -> Result<()> {
        let mol_map = MoleculeMap::new(self.molecules().as_slice());

        let n = mol_map.state_size();

        let mut a_matrix = Array2::<f64>::zeros((n, n));
        let mut sources = Array1::<f64>::zeros(n);

        for reaction in reactions {
            match reaction.reactants.len() {
                1 => {
                    let einstein_coefficient = reaction.einstein_coefficient.as_ref()
                        .ok_or_else(|| anyhow!("Einstein coefficient missing for unimolecular reaction"))?;

                    let branch = reaction.quantum_yield.unwrap_or(1.0);

                    let rate = einstein_coefficient(temperature) * branch;

                    let reactant = &reaction.reactants[0];
                    let reactant_index = mol_map.index(reactant)
                        .ok_or_else(|| anyhow!("Reactant molecule not found in molecule map"))?;

                    // One loss per reaction event
                    a_matrix[[reactant_index, reactant_index]] -= rate;

                    // Gains to tracked products only
                    for product in &reaction.products {
                        if mol_map.is_in_state(product) {
                            let product_index = mol_map.index(product)
                                .ok_or_else(|| anyhow!("Product molecule not found in molecule map"))?;
                            a_matrix[[product_index, reactant_index]] += rate;
                        }
                    }
                }
                2 => {
                    let k = reaction.rate_constant.as_ref()
                        .ok_or_else(|| anyhow!("Rate constant missing for bimolecular reaction"))?(temperature);

                    let branch = reaction.quantum_yield.unwrap_or(1.0);

                    let source = &reaction.reactants[0];
                    let collider = &reaction.reactants[1];

                    let source_index = mol_map.index(source)
                        .ok_or_else(|| anyhow!("Reactant molecule not found in molecule map"))?;

                    let collider_density = 0.0 /* get density of collider somehow */;
                    let rate = k * collider_density * branch;

                    // One loss per reaction event
                    a_matrix[[source_index, source_index]] -= rate;

                    // Gains to tracked products only
                    for product in &reaction.products {
                        if mol_map.is_in_state(product) {
                            let product_index = mol_map.index(product)
                                .ok_or_else(|| anyhow!("Product molecule not found in molecule map"))?;
                            a_matrix[[product_index, source_index]] += rate;
                        }
                    }
                }
                _ => {
                    return Err(anyhow!("Unsupported reaction with {} reactants", reaction.reactants.len()));
                }
            }

        }

        for (reaction, j) in photo_reactions.iter().zip(photolysis_rate.iter()) {
            let reactant_density = 0.0; // pull background density here
            let production = j * reactant_density;

            for product in &reaction.products {
                if mol_map.is_in_state(product) {
                    let product_index = mol_map.index(product)
                        .ok_or_else(|| anyhow!("Product molecule not found in molecule map"))?;
                    sources[product_index] += production;
                }
            }
        }

        Ok(())
    }
}

// Implementation of the O2 and O3 photochemistry models from
pub struct Yankovsky {
    photo_reactions: Vec<PhotoReaction>,
    chemical_reactions: Vec<ChemicalReaction>,
    
}

impl Default for Yankovsky {
fn default() -> Self {
Self::new()
}
}

impl Yankovsky {
    pub fn new() -> Self {
        let photo_reactions = vec![
            "O2 + hv(SRC) -> O(3P) + O(1D)".parse::<PhotoReaction>().unwrap().with_quantum_yield(1.0).with_toa_rate_constant(2.60e-6),
            "O2 + hv(lyman-alpha) -> O(3P) + O(1D)".parse::<PhotoReaction>().unwrap().with_quantum_yield(0.55).with_toa_rate_constant(3.40e-9),
            "O3 + hv -> O2(a, v=5) + O(1D)".parse::<PhotoReaction>().unwrap().with_quantum_yield(0.045).with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=4) + O(1D)".parse::<PhotoReaction>().unwrap().with_quantum_yield(0.072).with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=3) + O(1D)".parse::<PhotoReaction>().unwrap().with_quantum_yield(0.072).with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=2) + O(1D)".parse::<PhotoReaction>().unwrap().with_quantum_yield(0.135).with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=1) + O(1D)".parse::<PhotoReaction>().unwrap().with_quantum_yield(0.135).with_toa_rate_constant(8.0e-3),
            "O3 + hv -> O2(a, v=0) + O(1D)".parse::<PhotoReaction>().unwrap().with_quantum_yield(0.441).with_toa_rate_constant(8.0e-3),
            "O2 + hv(762_nm_band) -> O2(b, v=0)".parse::<PhotoReaction>().unwrap().with_toa_rate_constant(5.35e-9),
            "O2 + hv(689_nm_band) -> O2(b, v=1)".parse::<PhotoReaction>().unwrap().with_toa_rate_constant(2.94e-10),
            "O2 + hv(629_nm_band) -> O2(b, v=2)".parse::<PhotoReaction>().unwrap().with_toa_rate_constant(7.94e-12),
            "O2 + hv(1.27_um_band) -> O2(a, v=0)".parse::<PhotoReaction>().unwrap().with_toa_rate_constant(1.54e-10),
        ];

        let mut chemical_reactions = vec![
            // O(1D) deactivation reactions
            "O(1D) -> O".parse::<ChemicalReaction>().unwrap().with_einstein_coefficient(|_| 9.0e-3),
            "O(1D) + O -> O + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 4.0e-12),
            "O(1D) + O2 -> O2(b, v=1) + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {3.2e-11 * (67.0 / temperature).exp()}).with_quantum_yield(0.40),
            "O(1D) + O2 -> O2(b, v=0) + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {3.2e-11 * (67.0 / temperature).exp()}).with_quantum_yield(0.55),
            "O(1D) + O2 -> O2(X, v=0) + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {3.2e-11 * (67.0 / temperature).exp()}).with_quantum_yield(0.05),
            "O(1D) + O2 -> O2(a, v=0) + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {3.2e-11 * (67.0 / temperature).exp()}).with_quantum_yield(0.05),
            "O(1D) + O3 -> O2 + O2".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 2.4e-10),
            "O(1D) + N2 -> N2 + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.67),
            "O(1D) + N2 -> O + N2(v=0)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.33/7.0),
            "O(1D) + N2 -> O + N2(v=1)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.33/7.0),
            "O(1D) + N2 -> O + N2(v=2)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.33/7.0),
            "O(1D) + N2 -> O + N2(v=3)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.33/7.0),
            "O(1D) + N2 -> O + N2(v=4)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.33/7.0),
            "O(1D) + N2 -> O + N2(v=5)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.33/7.0),
            "O(1D) + N2 -> O + N2(v=6)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.33/7.0),
            "O(1D) + N2 -> O + N2(v=7)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| -> f64 {2.0e-11 * (107.0 / temperature).exp()}).with_quantum_yield(0.33/7.0),

            // O2(b, v) deactivation reactions
            "O2(b, v=2) -> O2(X, v=2)".parse::<ChemicalReaction>().unwrap().with_einstein_coefficient(|_| 5.4e-2),
            "O2(b, v=2) + O -> O2(b, v=1) + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 1.1e-11),
            "O2(b, v=2) + O2 -> O2(X, v=2) + O2(b, v=0)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 1.20e-11 * (-596.0 / temperature).exp()),
            "O2(b, v=2) + N2 -> O2(b, v=1) + N2".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 2e-14),
            "O2(b, v=2) + O3 -> O2 + O2 + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 2.9e-10),
            "O2(b, v=1) -> O2(X, v=1)".parse::<ChemicalReaction>().unwrap().with_einstein_coefficient(|_| 7.0e-2),
            "O2(b, v=1) + O -> O2(b, v=0) + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 4.5e-12),
            "O2(b, v=1) + O2 -> O2(X, v=1) + O2(b, v=0)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature| 4.20e-11 * (-312.0/temperature).exp()),
            "O2(b, v=1) + N2 -> O2(b, v=0) + N2".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 5.0e-13),
            "O2(b, v=1) + O3 -> O2 + O2 + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 3.0e-10),
            "O2(b, v=0) -> O2(X, v=0)".parse::<ChemicalReaction>().unwrap().with_einstein_coefficient(|_| 7.58e-2),
            "O2(b, v=0) + O -> O2(a, v=0) + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 8.0e-14).with_quantum_yield(0.75),
            "O2(b, v=0) + O -> O2(X, v=0) + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 8.0e-14).with_quantum_yield(0.25),
            "O2(b, v=0) + O2 -> O2(a, v=0) + O2(X, v=3)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 3.9e-17).with_quantum_yield(0.230),
            "O2(b, v=0) + O2 -> O2(a, v=1) + O2(X, v=2)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 3.9e-17).with_quantum_yield(0.525),
            "O2(b, v=0) + O2 -> O2(a, v=2) + O2(X, v=1)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 3.9e-17).with_quantum_yield(0.226),
            "O2(b, v=0) + O2 -> O2(a, v=3) + O2(X, v=0)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 3.9e-17).with_quantum_yield(0.019),
            // O2(b, 0) + N2 -> products ???
            "O2(b, v=0) + CO2 -> O2(a, v=0) + CO2".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 4.2e-13),
            "O2(b, v=0) + O3 -> O2(a, v=0) + O3".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 2.2e-11).with_quantum_yield(0.3),

            // O2(a, v) deactivation reactions
            "O2(a, v=0) -> O2(X, v=0)".parse::<ChemicalReaction>().unwrap().with_einstein_coefficient(|_| 2.58e-4),
            "O2(a, v=2) + O2 -> O2(X, v=2) + O2(a, v=0)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 3.6e-11),
            "O2(a, v=1) + O2 -> O2(X, v=1) + O2(a, v=0)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 5.6e-11),
            "O2(a, v=1) + O3 -> O2 + O2 + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 4.7e-12),
            "O2(a, v=0) + O -> O2 + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 6.5e-17),
            "O2(a, v=0) + O2 -> O2(X, v=5) + O2(X, v=0)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 3.6e-17 * (-220.0/temperature).exp()).with_quantum_yield(0.014),
            "O2(a, v=0) + O2 -> O2(X, v=4) + O2(X, v=1)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 3.6e-17 * (-220.0/temperature).exp()).with_quantum_yield(0.214),
            "O2(a, v=0) + O2 -> O2(X, v=3) + O2(X, v=2)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 3.6e-17 * (-220.0/temperature).exp()).with_quantum_yield(0.772),
            "O2(a, v=0) + O3 -> O2 + O3".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 5.20e-11 * (-2840.0/temperature).exp()),
            "O2(a, v=0) + N2 -> O2 + N2".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 1.0e-20),

            // Energy transfer and deactivation of O2(X, v), most are v dependent
            "O2(X, v=1) + O -> O2 + O".parse::<ChemicalReaction>().unwrap().with_rate_constant(|_| 3.0e-12),
            "O2(X, v=1) + O2 -> O2 + O2".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 4.2e-19 * (temperature / 300.0).sqrt()),
            "O2(X, v=1) + N2 -> O2 + N2(X, v=1)".parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 4.20e-19 * (temperature / 300.0).sqrt()),
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
                format!("O2(a, v={v}) + O -> O2 + O")
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
        for v in 0..=30 {
            // todo: quantum yield
            chemical_reactions.push(
                format!("O3 + O -> O2(X, v={v}) + O2").parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 5.60e-11 * (-1959.0/temperature).exp()),
            )
        }

        for v in 5..=35 {
            chemical_reactions.push(
                format!("O2(X, v={v}) + O -> O2 + O").parse::<ChemicalReaction>().unwrap().with_rate_constant(|temperature: f64| 5.0e-11 * (temperature / 300.0).sqrt()),
            )
        }

        for v in 2..=4 {
            chemical_reactions.push(
                format!("O2(X, v={v}) + O -> O2 + O").parse::<ChemicalReaction>().unwrap().with_rate_constant(move |temperature: f64| 1.1e-12 * (temperature / 300.0) * (1.0 * v as f64).exp()),
            )
        }

        for v in 2..=35 {
            let vm1 = v - 1;

            let rate_constant = match v {
                2 => |_| 2.0e-13,
                _ => |_| 2.6e-13
            };

            chemical_reactions.push(
                format!("O2(X, v={v}) + O2 -> O2(X, v={vm1}) + O2(X, v=1)").parse::<ChemicalReaction>().unwrap().with_rate_constant(rate_constant),
            );
        }

        for v in 4..=20 {
            let vm1 = v - 1;

            let rate_constant = move |_| 1.3e-12 * (-0.31 * v as f64).exp();

            chemical_reactions.push(
                format!("O2(X, v={v}) + O2 -> O2(X, v={vm1}) + O2(X, v=1)").parse::<ChemicalReaction>().unwrap().with_rate_constant(rate_constant),
            );
        }

        for v in 21..=35 {
            let vm1 = v - 1;

            let rate_constant = move |temperature: f64| 6.0e-17 * (temperature / 300.0) * (0.2 * v as f64).exp();

            chemical_reactions.push(
                format!("O2(X, v={v}) + O2 -> O2(X, v={vm1}) + O2").parse::<ChemicalReaction>().unwrap().with_rate_constant(rate_constant),
            );
        }

        for v in 12..=17 {
            let vm2 = v - 2;

            let rate_constant = move |_| 3.6e-19 * (0.66 * v as f64).exp();

            chemical_reactions.push(
                format!("O2(X, v={v}) + N2 -> O2(X, v={vm2}) + N2(X, v=1)").parse::<ChemicalReaction>().unwrap().with_rate_constant(rate_constant),
            );
        }

        for v in 18..=26 {
            let vm2 = v - 2;

            let rate_constant = move |_| 4.5e-13 * (-0.173 * v as f64).exp();

            chemical_reactions.push(
                format!("O2(X, v={v}) + N2 -> O2(X, v={vm2}) + N2(X, v=1)").parse::<ChemicalReaction>().unwrap().with_rate_constant(rate_constant),
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
            if !unique_mols
                .iter()
                .any(|existing: &Molecule| existing.base_type == mol.base_type)
            {
                unique_mols.push(mol);
            }
        }

        unique_mols
    }
}

mod tests {
    use super::*;

    #[test]
    fn test_yankovsky() {
        let model = Yankovsky::new();
    }

    #[test]
    fn test_molecules() {
        let model = Yankovsky::new();
        let mols = model.molecules();
        assert!(mols.iter().any(|m| m.base_type == MoleculeBase::O2));
        assert!(mols.iter().any(|m| m.base_type == MoleculeBase::O3));
        assert!(mols.iter().any(|m| m.base_type == MoleculeBase::O));
    }
}