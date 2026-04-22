

use std::collections::HashMap;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub enum MoleculeBase {
    O2,
    O3,
    O,
    N2,
    CO2,
}

impl MoleculeBase {
    pub fn to_string(&self) -> String {
        match self {
            MoleculeBase::O2 => "O2".to_string(),
            MoleculeBase::O3 => "O3".to_string(),
            MoleculeBase::O => "O".to_string(),
            MoleculeBase::N2 => "N2".to_string(),
            MoleculeBase::CO2 => "CO2".to_string(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Molecule {
    pub base_type: MoleculeBase,
    pub electronic_level: String,
    pub vibrational_level: u32,
}

pub struct MoleculeMap {
    molecule_to_index: HashMap<Molecule, usize>,
    background_molecules: std::collections::HashSet<Molecule>,
}

impl MoleculeMap {
    pub fn new(molecules: &[Molecule], densities: &HashMap<String, f64>) -> Self {
        let mut molecule_to_index = HashMap::new();

        let is_background = |m: &Molecule| {
            String::try_from(m.clone())
                .ok()
                .map(|name| densities.contains_key(&name))
                .unwrap_or(false)
        };

        // Split molecules into state (not in densities) and background (in densities).
        // State molecules are indexed first (0..n_state), background molecules after.
        let (state_molecules, background_molecules): (Vec<_>, Vec<_>) = molecules
            .iter()
            .cloned()
            .partition(|m| !is_background(m));

        let mut idx = 0;
        for molecule in state_molecules {
            molecule_to_index.entry(molecule).or_insert_with(|| {
                let i = idx;
                idx += 1;
                i
            });
        }
        for molecule in background_molecules {
            molecule_to_index.entry(molecule).or_insert_with(|| {
                let i = idx;
                idx += 1;
                i
            });
        }

        let background_molecules = molecules
            .iter()
            .filter(|m| is_background(m))
            .cloned()
            .collect();

        Self {
            molecule_to_index,
            background_molecules,
        }
    }

    pub fn index(&self, molecule: &Molecule) -> Option<usize> {
        self.molecule_to_index.get(molecule).copied()
    }

    pub fn size(&self) -> usize {
        self.molecule_to_index.len()
    }

    pub fn state_size(&self) -> usize {
        self.molecule_to_index.len() - self.background_molecules.len()
    }

    pub fn is_in_state(&self, molecule: &Molecule) -> bool {
        !self.background_molecules.contains(molecule)
    }

    pub fn state_index_to_molecule_names(&self) -> Vec<String> {
        let mut entries: Vec<(usize, String)> = self
            .molecule_to_index
            .iter()
            .filter_map(|(molecule, index)| {
                if self.is_in_state(molecule) {
                    let name = String::try_from(molecule.clone())
                        .unwrap_or_else(|_| format!("{:?}", molecule));
                    Some((*index, name))
                } else {
                    None
                }
            })
            .collect();

        entries.sort_by_key(|(index, _)| *index);
        entries.into_iter().map(|(_, name)| name).collect()
    }
}


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParseMoleculeError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParseChemicalReactionError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParsePhotoReactionError;

impl Molecule {
    pub fn from_str(s: &str) -> Option<Self> {
        // Example input: "O2(X, v=0)" or "O3" or "O(1D)" or similar
        let s = s.trim();
        if s.is_empty() {
            return None;
        }

        let (base_type, rest) = if let Some(rest) = s.strip_prefix("O2") {
            (MoleculeBase::O2, rest)
        } else if let Some(rest) = s.strip_prefix("O3") {
            (MoleculeBase::O3, rest)
        } else if let Some(rest) = s.strip_prefix('O') {
            (MoleculeBase::O, rest)
        } else if let Some(rest) = s.strip_prefix("N2") {
            (MoleculeBase::N2, rest)
        } else if let Some(rest) = s.strip_prefix("CO2") {
            (MoleculeBase::CO2, rest)
        } else {
            return None;
        };

        let rest = rest.trim();
        if rest.is_empty() {
            return Some(Self {
                base_type,
                electronic_level: String::new(),
                vibrational_level: 0,
            });
        }

        if !(rest.starts_with('(') && rest.ends_with(')')) {
            return None;
        }

        let inner = rest[1..rest.len() - 1].trim();
        if inner.is_empty() {
            return Some(Self {
                base_type,
                electronic_level: String::new(),
                vibrational_level: 0,
            });
        }

        let mut electronic_level = String::new();
        let mut vibrational_level: u32 = 0;

        for token in inner.split(',') {
            let token = token.trim();
            if token.is_empty() {
                return None;
            }

            if let Some((lhs, rhs)) = token.split_once('=') {
                if lhs.trim() == "v" {
                    vibrational_level = rhs.trim().parse().ok()?;
                } else {
                    return None;
                }
            } else if electronic_level.is_empty() {
                electronic_level = token.to_string();
            } else {
                return None;
            }
        }

        Some(Self {
            base_type,
            electronic_level,
            vibrational_level,
        })
    }
}

impl FromStr for Molecule {
    type Err = ParseMoleculeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Molecule::from_str(s).ok_or(ParseMoleculeError)
    }
}

impl TryFrom<&str> for Molecule {
    type Error = ParseMoleculeError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        value.parse()
    }
}

impl TryFrom<Molecule> for String {
    type Error = ();

    fn try_from(value: Molecule) -> Result<Self, Self::Error> {
        let mut s = value.base_type.to_string();
        if !value.electronic_level.is_empty() || value.vibrational_level != 0 {
            s.push('(');
            if !value.electronic_level.is_empty() {
                s.push_str(&value.electronic_level);
            }
            if value.vibrational_level != 0 {
                if !value.electronic_level.is_empty() {
                    s.push_str(", ");
                }
                s.push_str(&format!("v={}", value.vibrational_level));
            }
            s.push(')');
        }
        Ok(s)
    }
}


// Photodissociation or Photoexcitation reaction
// A molecule plus a photon produces products
pub struct PhotoReaction {
    pub in_molecule: Molecule,
    pub toa_rate_constant: f64, // s^-1, at TOA
    pub products: Vec<Molecule>,
    pub quantum_yield: Option<f64>, // Optional quantum yield for photochemical reactions
    pub excitation_band: Option<String>,
    pub quantum_yield_wavelength_nm: f64,
    pub wavelength_range_nm: Option<(f64, f64)>,
}

impl PhotoReaction {
    pub fn from_str(r: &str) -> Option<Self> {
        // example O2 + hv(SRC) -> O(3P) + O(1D)
        // example O2 + hv(lyman-alpha) -> O2(3P) + O(1D)
        // example O3 + hv -> O2(a, v=2) + O(1D)
        let r = r.trim();
        if r.is_empty() {
            return None;
        }

        let (lhs, rhs) = r.split_once("->")?;

        let lhs_tokens: Vec<&str> = lhs
            .split('+')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect();

        if lhs_tokens.len() != 2 {
            return None;
        }

        let molecule = Molecule::from_str(lhs_tokens[0])?;
        let photon = lhs_tokens[1];

        let excitation_band = if photon == "hv" || photon == "hν" {
            None
        } else if photon.starts_with("hv(") && photon.ends_with(')') {
            let band = photon[3..photon.len() - 1].trim();
            if band.is_empty() {
                return None;
            }
            Some(band.to_string())
        } else if photon.starts_with("hν(") && photon.ends_with(')') {
            let band = photon[3..photon.len() - 1].trim();
            if band.is_empty() {
                return None;
            }
            Some(band.to_string())
        } else {
            return None;
        };

        let products: Vec<Molecule> = rhs
            .split('+')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(Molecule::from_str)
            .collect::<Option<Vec<_>>>()?;

        if products.is_empty() {
            return None;
        }

        Some(Self {
            in_molecule: molecule,
            toa_rate_constant: 0.0,
            products,
            quantum_yield: None,
            excitation_band,
            quantum_yield_wavelength_nm: 0.0,
            wavelength_range_nm: None,
        })
    }

    pub fn with_quantum_yield(mut self, q: f64) -> Self {
        self.quantum_yield = Some(q);
        self
    }

    pub fn with_toa_rate_constant(mut self, k: f64) -> Self {
        self.toa_rate_constant = k;
        self
    }

    pub fn with_wavelength_range_nm(mut self, min_nm: f64, max_nm: f64) -> Self {
        self.wavelength_range_nm = Some((min_nm, max_nm));
        self
    }

    pub fn with_band_center_nm(mut self, center_nm: f64, half_width_nm: f64) -> Self {
        self.wavelength_range_nm = Some((center_nm - half_width_nm, center_nm + half_width_nm));
        self
    }

}

impl FromStr for PhotoReaction {
    type Err = ParsePhotoReactionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        PhotoReaction::from_str(s).ok_or(ParsePhotoReactionError)
    }
}

impl TryFrom<&str> for PhotoReaction {
    type Error = ParsePhotoReactionError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        value.parse()
    }
}

pub struct ChemicalReaction {
    pub reactants: Vec<Molecule>,
    pub products: Vec<Molecule>,
    pub einstein_coefficient: Option<Box<dyn Fn(f64) -> f64>>, // Optional Einstein coefficient as a function of temperature
    pub rate_constant: Option<Box<dyn Fn(f64) -> f64>>, // Optional rate constant
    pub quantum_yield: Option<f64>, // Optional quantum yield for photochemical reactions
}

impl ChemicalReaction {
    pub fn from_str(r: &str) -> Option<Self> {
        // Example "O(1D) + O -> O + O"
        // "O(1D) + O2 -> O2(b, v=1) + O"
        let r = r.trim();
        if r.is_empty() {
            return None;
        }

        let (lhs, rhs) = r.split_once("->")?;

        let reactants: Vec<Molecule> = lhs
            .split('+')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(Molecule::from_str)
            .collect::<Option<Vec<_>>>()?;

        let products: Vec<Molecule> = rhs
            .split('+')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(Molecule::from_str)
            .collect::<Option<Vec<_>>>()?;

        if reactants.is_empty() || products.is_empty() {
            return None;
        }

        Some(Self {
            reactants,
            products,
            einstein_coefficient: None,
            quantum_yield: None,
            rate_constant: None,
        })
    }

    pub fn with_einstein_coefficient(mut self, f: impl Fn(f64) -> f64 + 'static) -> Self {
        self.einstein_coefficient = Some(Box::new(f));
        self
    }

    pub fn with_rate_constant(mut self, f: impl Fn(f64) -> f64 + 'static) -> Self {
        self.rate_constant = Some(Box::new(f));
        self
    }

    pub fn with_quantum_yield(mut self, q: f64) -> Self {
        self.quantum_yield = Some(q);
        self
    }
}

impl FromStr for ChemicalReaction {
    type Err = ParseChemicalReactionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ChemicalReaction::from_str(s).ok_or(ParseChemicalReactionError)
    }
}

impl TryFrom<&str> for ChemicalReaction {
    type Error = ParseChemicalReactionError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        value.parse()

    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_simple_molecule() {
        let m = Molecule::from_str("O3").expect("O3 should parse");
        assert_eq!(m.base_type, MoleculeBase::O3);
        assert_eq!(m.electronic_level, "");
        assert_eq!(m.vibrational_level, 0);
    }

    #[test]
    fn parses_electronic_and_vibrational_levels() {
        let m = Molecule::from_str("O2(X, v=2)").expect("O2(X, v=2) should parse");
        assert_eq!(m.base_type, MoleculeBase::O2);
        assert_eq!(m.electronic_level, "X");
        assert_eq!(m.vibrational_level, 2);
    }

    #[test]
    fn parses_atomic_state() {
        let m = Molecule::from_str("O(1D)").expect("O(1D) should parse");
        assert_eq!(m.base_type, MoleculeBase::O);
        assert_eq!(m.electronic_level, "1D");
        assert_eq!(m.vibrational_level, 0);
    }

    #[test]
    fn parse_trait_works() {
        let parsed: Molecule = "O2(X, v=0)"
            .parse()
            .expect("FromStr impl should parse valid molecule string");
        assert_eq!(parsed.base_type, MoleculeBase::O2);
        assert_eq!(parsed.electronic_level, "X");
        assert_eq!(parsed.vibrational_level, 0);
    }

    #[test]
    fn try_from_trait_works() {
        let parsed = Molecule::try_from("O").expect("TryFrom<&str> should parse valid molecule string");
        assert_eq!(parsed.base_type, MoleculeBase::O);
        assert_eq!(parsed.electronic_level, "");
        assert_eq!(parsed.vibrational_level, 0);
    }

    #[test]
    fn rejects_invalid_strings() {
        assert!(Molecule::from_str("N").is_none());
        assert!(Molecule::from_str("O2(v=abc)").is_none());
        assert!(Molecule::from_str("O2(X, bad=1)").is_none());
    }

    #[test]
    fn parses_chemical_reaction() {
        let r = ChemicalReaction::from_str("O(1D) + O2 -> O + O2")
            .expect("valid chemical reaction should parse");
        assert_eq!(r.reactants.len(), 2);
        assert_eq!(r.products.len(), 2);
        assert!(r.einstein_coefficient.is_none());
    }

    #[test]
    fn parses_photo_reaction_without_band() {
        let r = PhotoReaction::from_str("O3 + hv -> O2 + O(1D)")
            .expect("valid photo reaction should parse");
        assert_eq!(r.in_molecule.base_type, MoleculeBase::O3);
        assert_eq!(r.products.len(), 2);
        assert!(r.excitation_band.is_none());
        assert_eq!(r.toa_rate_constant, 0.0);
    }

    #[test]
    fn parses_photo_reaction_with_band() {
        let r = "O2 + hv(lyman-alpha) -> O + O"
            .parse::<PhotoReaction>()
            .expect("FromStr impl should parse photo reaction");
        assert_eq!(r.in_molecule.base_type, MoleculeBase::O2);
        assert_eq!(r.products.len(), 2);
        assert_eq!(r.excitation_band, Some("lyman-alpha".to_string()));

        let r2 = PhotoReaction::try_from("O2 + hv(SRC) -> O + O")
            .expect("TryFrom<&str> should parse photo reaction");
        assert_eq!(r2.excitation_band, Some("SRC".to_string()));
    }

    #[test]
    fn rejects_invalid_photo_reaction_strings() {
        assert!(PhotoReaction::from_str("").is_none());
        assert!(PhotoReaction::from_str("O3 + hv").is_none());
        assert!(PhotoReaction::from_str("O3 -> O2 + O").is_none());
        assert!(PhotoReaction::from_str("O3 + hq -> O2 + O").is_none());
        assert!(PhotoReaction::from_str("O3 + hv() -> O2 + O").is_none());
        assert!(PhotoReaction::from_str("O3 + hv ->").is_none());
    }

    #[test]
    fn parses_chemical_reaction_traits() {
        let parsed: ChemicalReaction = "O + O3 -> O2 + O2"
            .parse()
            .expect("FromStr impl should parse chemical reaction");
        assert_eq!(parsed.reactants.len(), 2);
        assert_eq!(parsed.products.len(), 2);

        let parsed_try = ChemicalReaction::try_from("O + O3 -> O2 + O2")
            .expect("TryFrom<&str> should parse chemical reaction");
        assert_eq!(parsed_try.reactants.len(), 2);
        assert_eq!(parsed_try.products.len(), 2);
    }

    #[test]
    fn rejects_invalid_chemical_reaction_strings() {
        assert!(ChemicalReaction::from_str("").is_none());
        assert!(ChemicalReaction::from_str("O + O2").is_none());
        assert!(ChemicalReaction::from_str("O + N -> O2").is_none());
        assert!(ChemicalReaction::from_str("O + O2 ->").is_none());
        assert!(ChemicalReaction::from_str("-> O + O2").is_none());
    }

    #[test]
    fn chemical_reaction_builder_aliases_set_einstein_coefficient() {
        let r1 = "O + O3 -> O2 + O2"
            .parse::<ChemicalReaction>()
            .expect("reaction should parse")
            .with_einstein_coefficient(|temperature: f64| 1.0e-3 * temperature);

        let r2 = "O + O3 -> O2 + O2"
            .parse::<ChemicalReaction>()
            .expect("reaction should parse")
            .with_rate_constant(|temperature: f64| 2.0e-3 * temperature);

        let f1 = r1
            .einstein_coefficient
            .as_ref()
            .expect("einstein coefficient should be set");
        let f2 = r2
            .einstein_coefficient
            .as_ref()
            .expect("einstein coefficient should be set via alias");

        assert_eq!(f1(200.0), 2.0e-1);
        assert_eq!(f2(200.0), 4.0e-1);
    }

    #[test]
    fn molecule_map_returns_first_matching_index() {
        let molecules = vec![
            Molecule::from_str("O3").expect("O3 should parse"),
            Molecule::from_str("O(1D)").expect("O(1D) should parse"),
            Molecule::from_str("O3").expect("O3 should parse"),
        ];

        let molecule_map = MoleculeMap::new(&molecules, &HashMap::new());

        assert_eq!(molecule_map.index(&molecules[0]), Some(0));
        assert_eq!(molecule_map.index(&molecules[1]), Some(1));
        assert_eq!(molecule_map.index(&molecules[2]), Some(0));
    }

    #[test]
    fn molecule_map_distinguishes_molecular_state() {
        let molecules = vec![
            Molecule::from_str("O2(a, v=0)").expect("O2(a, v=0) should parse"),
            Molecule::from_str("O2(a, v=1)").expect("O2(a, v=1) should parse"),
        ];

        let molecule_map = MoleculeMap::new(&molecules, &HashMap::new());
        let missing = Molecule::from_str("CO2").expect("CO2 should parse");

        assert_eq!(molecule_map.index(&molecules[0]), Some(0));
        assert_eq!(molecule_map.index(&molecules[1]), Some(1));
        assert_eq!(molecule_map.index(&missing), None);
    }

}