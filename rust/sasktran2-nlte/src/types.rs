use std::fmt;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub enum MoleculeBase {
    O2,
    O3,
    O,
    N2,
    CO2,
}

impl fmt::Display for MoleculeBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            MoleculeBase::O2 => "O2",
            MoleculeBase::O3 => "O3",
            MoleculeBase::O => "O",
            MoleculeBase::N2 => "N2",
            MoleculeBase::CO2 => "CO2",
        };
        f.write_str(name)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Molecule {
    pub base_type: MoleculeBase,
    pub electronic_level: String,
    pub vibrational_level: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParseMoleculeError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParsePhotoReactionError;

impl Molecule {
    #[allow(clippy::should_implement_trait)]
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
    pub line_center_nm: Option<f64>,
    pub line_effective_cross_section_m2: Option<f64>,
}

impl PhotoReaction {
    #[allow(clippy::should_implement_trait)]
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
        } else if let Some(band) = photon
            .strip_prefix("hv(")
            .or_else(|| photon.strip_prefix("hν("))
            .and_then(|s| s.strip_suffix(')'))
        {
            let band = band.trim();
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
            line_center_nm: None,
            line_effective_cross_section_m2: None,
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

    pub fn with_line_center_nm(mut self, center_nm: f64) -> Self {
        self.line_center_nm = Some(center_nm);
        self
    }

    pub fn with_line_effective_cross_section_m2(mut self, cross_section_m2: f64) -> Self {
        self.line_effective_cross_section_m2 = Some(cross_section_m2);
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
        let parsed =
            Molecule::try_from("O").expect("TryFrom<&str> should parse valid molecule string");
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
        assert_eq!(r.line_center_nm, None);

        let r2 = PhotoReaction::try_from("O2 + hv(SRC) -> O + O")
            .expect("TryFrom<&str> should parse photo reaction");
        assert_eq!(r2.excitation_band, Some("SRC".to_string()));
    }

    #[test]
    fn photo_reaction_can_store_line_center() {
        let r = "O2 + hv(lyman-alpha) -> O + O"
            .parse::<PhotoReaction>()
            .expect("FromStr impl should parse photo reaction")
            .with_line_center_nm(121.567)
            .with_line_effective_cross_section_m2(1.0e-24);

        assert_eq!(r.line_center_nm, Some(121.567));
        assert_eq!(r.line_effective_cross_section_m2, Some(1.0e-24));
        assert_eq!(r.wavelength_range_nm, None);
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
}
