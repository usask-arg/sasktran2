//! Kinetic mechanisms: excited states, background species, and the processes
//! that couple them, loaded and validated from TOML files.
//!
//! The file format is documented in
//! `docs/sphinx/source/developer/nlte_mechanism_format.md`.

mod file;
mod template;

use crate::prelude::*;
use crate::rates::{RateLaw, order_and_si_factor};
use file::{ChannelEntry, MechanismFile, RateEntry};
use std::collections::{BTreeMap, HashSet};
use template::{canonical_species, levels, substitute};

const BUNDLED: &[(&str, &str)] = &[
    ("oxygen", include_str!("../../mechanisms/oxygen.toml")),
    (
        "oxygen_green",
        include_str!("../../mechanisms/oxygen_green.toml"),
    ),
    (
        "oxygen_yankovsky",
        include_str!("../../mechanisms/oxygen_yankovsky.toml"),
    ),
];

/// A participant in a process: either an unknown state population or a
/// background number density supplied as input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Species {
    State(usize),
    Background(usize),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProcessKind {
    Reaction,
    Photolysis,
    Radiative,
}

impl ProcessKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            ProcessKind::Reaction => "reaction",
            ProcessKind::Photolysis => "photolysis",
            ProcessKind::Radiative => "radiative",
        }
    }
}

/// How a process's rate coefficient is obtained, in SI units.
#[derive(Clone, Debug, PartialEq)]
pub enum Coefficient {
    /// Temperature-dependent rate coefficient.
    Rate(RateLaw),
    /// Index into [`Mechanism::rate_inputs`]; supplied per altitude [s^-1].
    Input(usize),
    /// Einstein A coefficient [s^-1].
    EinsteinA(f64),
}

/// One product branch; `fraction` of reaction events produce `products`.
#[derive(Clone, Debug, PartialEq)]
pub struct Channel {
    pub fraction: f64,
    pub products: Vec<Species>,
}

/// A process with event rate `coefficient * product of reactant densities`.
/// Every reactant occurrence is consumed once per event.
#[derive(Clone, Debug, PartialEq)]
pub struct Process {
    pub id: String,
    pub kind: ProcessKind,
    pub reference: String,
    pub coefficient: Coefficient,
    pub reactants: Vec<Species>,
    pub channels: Vec<Channel>,
    /// Emission wavelength for radiative processes, when given.
    pub wavelength_nm: Option<f64>,
}

impl Process {
    /// The state reactant, if any. Mechanisms are validated to have at most one.
    pub fn state_reactant(&self) -> Option<usize> {
        self.reactants.iter().find_map(|species| match species {
            Species::State(index) => Some(*index),
            Species::Background(_) => None,
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct State {
    pub id: String,
    pub description: String,
}

#[derive(Clone, Debug)]
pub struct Mechanism {
    pub name: String,
    pub version: String,
    pub description: String,
    references: BTreeMap<String, String>,
    states: Vec<State>,
    background: Vec<String>,
    rate_inputs: Vec<String>,
    processes: Vec<Process>,
}

impl Mechanism {
    /// Loads and validates a mechanism from TOML text.
    pub fn from_toml_str(text: &str) -> Result<Self> {
        let file: MechanismFile =
            toml::from_str(text).map_err(|err| anyhow!("Invalid mechanism file: {err}"))?;
        Builder::default().build(file)
    }

    /// One of the mechanisms shipped with the crate, see [`Mechanism::bundled_names`].
    pub fn bundled(name: &str) -> Result<Self> {
        let (_, text) = BUNDLED
            .iter()
            .find(|(bundled, _)| *bundled == name)
            .ok_or_else(|| {
                anyhow!(
                    "No bundled mechanism '{name}'; available: {}",
                    Self::bundled_names().join(", ")
                )
            })?;
        Self::from_toml_str(text).map_err(|err| anyhow!("Bundled mechanism '{name}': {err}"))
    }

    pub fn bundled_names() -> Vec<&'static str> {
        BUNDLED.iter().map(|(name, _)| *name).collect()
    }

    pub fn states(&self) -> &[State] {
        &self.states
    }

    pub fn background(&self) -> &[String] {
        &self.background
    }

    /// Names of the per-molecule rates [s^-1] that must be supplied, such as
    /// photolysis rates, in first-use order.
    pub fn rate_inputs(&self) -> &[String] {
        &self.rate_inputs
    }

    pub fn processes(&self) -> &[Process] {
        &self.processes
    }

    pub fn references(&self) -> &BTreeMap<String, String> {
        &self.references
    }

    pub fn state_index(&self, id: &str) -> Option<usize> {
        let id = canonical_species(id);
        self.states.iter().position(|state| state.id == id)
    }
}

#[derive(Default)]
struct Builder {
    species: HashMap<String, Species>,
    states: Vec<State>,
    background: Vec<String>,
    rate_inputs: Vec<String>,
    process_ids: HashSet<String>,
    processes: Vec<Process>,
    references: BTreeMap<String, String>,
}

impl Builder {
    fn build(mut self, file: MechanismFile) -> Result<Mechanism> {
        self.references = file.references;

        for name in &file.mechanism.background {
            let id = canonical_species(name);
            if self.species.contains_key(&id) {
                return Err(anyhow!("Background species '{id}' is listed twice"));
            }
            self.species
                .insert(id.clone(), Species::Background(self.background.len()));
            self.background.push(id);
        }

        for entry in &file.states {
            for v in levels(&entry.id, entry.for_v)? {
                let id = canonical_species(&substitute(&entry.id, v)?);
                match self.species.get(&id) {
                    Some(Species::Background(_)) => {
                        return Err(anyhow!("'{id}' is both a state and a background species"));
                    }
                    Some(Species::State(_)) => {
                        return Err(anyhow!("State '{id}' is declared twice"));
                    }
                    None => {}
                }
                self.species
                    .insert(id.clone(), Species::State(self.states.len()));
                self.states.push(State {
                    id,
                    description: entry.description.clone(),
                });
            }
        }

        for entry in &file.reactions {
            for v in levels(&entry.id, entry.for_v)? {
                let id = substitute(&entry.id, v)?;
                let reactants = self.resolve_all(&entry.reactants, v, &id)?;
                let channels = self.channels(&entry.products, &entry.channels, v, &id)?;
                let rate = rate_law(&id, &entry.rate, v, reactants.len())?;
                self.push(
                    Process {
                        id,
                        kind: ProcessKind::Reaction,
                        reference: entry.reference.clone(),
                        coefficient: Coefficient::Rate(rate),
                        reactants,
                        channels,
                        wavelength_nm: None,
                    },
                    &entry.reference,
                )?;
            }
        }

        for entry in &file.photo_reactions {
            for v in levels(&entry.id, entry.for_v)? {
                let id = substitute(&entry.id, v)?;
                let reactant = self.resolve(&entry.reactant, v, &id)?;
                let channels = self.channels(&entry.products, &entry.channels, v, &id)?;
                let input = self.rate_input(&substitute(&entry.rate_input, v)?, &id)?;
                self.push(
                    Process {
                        id,
                        kind: ProcessKind::Photolysis,
                        reference: entry.reference.clone(),
                        coefficient: Coefficient::Input(input),
                        reactants: vec![reactant],
                        channels,
                        wavelength_nm: None,
                    },
                    &entry.reference,
                )?;
            }
        }

        for entry in &file.transitions {
            for v in levels(&entry.id, entry.for_v)? {
                let id = substitute(&entry.id, v)?;
                let upper = self.resolve(&entry.upper, v, &id)?;
                if !matches!(upper, Species::State(_)) {
                    return Err(anyhow!(
                        "Transition '{id}': upper level '{}' must be a state",
                        entry.upper
                    ));
                }
                let lower = self.resolve(&entry.lower, v, &id)?;
                if !entry.einstein_a_s.is_finite() || entry.einstein_a_s < 0.0 {
                    return Err(anyhow!(
                        "Transition '{id}': einstein_a_s must be non-negative and finite"
                    ));
                }
                if let Some(wavelength) = entry.wavelength_nm
                    && !(wavelength.is_finite() && wavelength > 0.0)
                {
                    return Err(anyhow!("Transition '{id}': wavelength_nm must be positive"));
                }
                self.push(
                    Process {
                        id,
                        kind: ProcessKind::Radiative,
                        reference: entry.reference.clone(),
                        coefficient: Coefficient::EinsteinA(entry.einstein_a_s),
                        reactants: vec![upper],
                        channels: vec![Channel {
                            fraction: 1.0,
                            products: vec![lower],
                        }],
                        wavelength_nm: entry.wavelength_nm,
                    },
                    &entry.reference,
                )?;
            }
        }

        Ok(Mechanism {
            name: file.mechanism.name,
            version: file.mechanism.version,
            description: file.mechanism.description,
            references: self.references,
            states: self.states,
            background: self.background,
            rate_inputs: self.rate_inputs,
            processes: self.processes,
        })
    }

    fn resolve(&self, name: &str, v: Option<u32>, process: &str) -> Result<Species> {
        let id = canonical_species(&substitute(name, v)?);
        self.species.get(&id).copied().ok_or_else(|| {
            anyhow!(
                "'{process}' uses '{id}', which is neither a [[state]] nor listed in mechanism.background"
            )
        })
    }

    fn resolve_all(&self, names: &[String], v: Option<u32>, process: &str) -> Result<Vec<Species>> {
        names
            .iter()
            .map(|name| self.resolve(name, v, process))
            .collect()
    }

    /// The product channels of a reaction or photolysis entry, given either
    /// `products` (one channel with unit yield) or `channels`.
    fn channels(
        &self,
        products: &Option<Vec<String>>,
        channels: &Option<Vec<ChannelEntry>>,
        v: Option<u32>,
        process: &str,
    ) -> Result<Vec<Channel>> {
        let channels = match (products, channels) {
            (Some(products), None) => vec![Channel {
                fraction: 1.0,
                products: self.resolve_all(products, v, process)?,
            }],
            (None, Some(channels)) => channels
                .iter()
                .map(|channel| {
                    Ok(Channel {
                        fraction: channel.fraction,
                        products: self.resolve_all(&channel.products, v, process)?,
                    })
                })
                .collect::<Result<Vec<_>>>()?,
            _ => {
                return Err(anyhow!(
                    "'{process}' needs exactly one of 'products' or 'channels'"
                ));
            }
        };
        validate_channels(process, &channels)?;
        Ok(channels)
    }

    fn rate_input(&mut self, name: &str, process: &str) -> Result<usize> {
        if name.trim().is_empty() {
            return Err(anyhow!("'{process}' has an empty rate_input"));
        }
        if let Some(index) = self.rate_inputs.iter().position(|input| input == name) {
            return Ok(index);
        }
        self.rate_inputs.push(name.to_string());
        Ok(self.rate_inputs.len() - 1)
    }

    fn push(&mut self, process: Process, reference: &str) -> Result<()> {
        if !self.references.contains_key(reference) {
            return Err(anyhow!(
                "'{}' cites '{reference}', which is not in [references]",
                process.id
            ));
        }
        if !self.process_ids.insert(process.id.clone()) {
            return Err(anyhow!("Process id '{}' is used twice", process.id));
        }
        let state_reactants = process
            .reactants
            .iter()
            .filter(|species| matches!(species, Species::State(_)))
            .count();
        if state_reactants > 1 {
            return Err(anyhow!(
                "'{}' has {state_reactants} state reactants; only processes linear in the \
                 state populations are supported",
                process.id
            ));
        }
        self.processes.push(process);
        Ok(())
    }
}

fn validate_channels(id: &str, channels: &[Channel]) -> Result<()> {
    if channels.is_empty() {
        return Err(anyhow!("'{id}' has no channels"));
    }
    let mut total = 0.0;
    for channel in channels {
        if !(channel.fraction.is_finite() && (0.0..=1.0).contains(&channel.fraction)) {
            return Err(anyhow!(
                "'{id}': channel yields must be between 0 and 1, got {}",
                channel.fraction
            ));
        }
        total += channel.fraction;
    }
    if total > 1.0 + 1e-9 {
        return Err(anyhow!(
            "'{id}': channel yields sum to {total}, more than 1"
        ));
    }
    Ok(())
}

fn rate_law(id: &str, entry: &RateEntry, v: Option<u32>, n_reactants: usize) -> Result<RateLaw> {
    let (order, si_factor) = order_and_si_factor(&entry.units).ok_or_else(|| {
        anyhow!(
            "'{id}': unknown rate units '{}'; use s-1, cm3 s-1, m3 s-1, cm6 s-1 or m6 s-1",
            entry.units
        )
    })?;
    if order != n_reactants {
        return Err(anyhow!(
            "'{id}': units '{}' are for {order} reactant(s) but the reaction has {n_reactants}",
            entry.units
        ));
    }

    let mut law = match entry.law.as_str() {
        "constant" => {
            if entry.a.is_some()
                || entry.n.is_some()
                || entry.t0.is_some()
                || entry.ea_over_r.is_some()
            {
                return Err(anyhow!("'{id}': a constant rate takes only 'value'"));
            }
            RateLaw::constant(
                entry
                    .value
                    .ok_or_else(|| anyhow!("'{id}': a constant rate needs 'value'"))?,
            )
        }
        "arrhenius" => {
            if entry.value.is_some() {
                return Err(anyhow!("'{id}': an arrhenius rate takes 'a', not 'value'"));
            }
            RateLaw {
                a: entry
                    .a
                    .ok_or_else(|| anyhow!("'{id}': an arrhenius rate needs 'a'"))?,
                n: entry.n.unwrap_or(0.0),
                t0: entry.t0.unwrap_or(300.0),
                ea_over_r: entry.ea_over_r.unwrap_or(0.0),
            }
        }
        other => {
            return Err(anyhow!(
                "'{id}': unknown rate law '{other}'; use 'constant' or 'arrhenius'"
            ));
        }
    };

    if let Some(exp_v) = entry.exp_v {
        let v = v.ok_or_else(|| anyhow!("'{id}': exp_v needs for_v"))?;
        law.a *= (exp_v * f64::from(v)).exp();
    }
    law.a *= si_factor;

    let finite = law.a.is_finite() && law.n.is_finite() && law.ea_over_r.is_finite();
    if !finite || law.a < 0.0 {
        return Err(anyhow!(
            "'{id}': rate parameters must be finite and non-negative"
        ));
    }
    if !law.t0.is_finite() || law.t0 <= 0.0 {
        return Err(anyhow!("'{id}': t0 must be positive"));
    }
    Ok(law)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL: &str = r#"
[mechanism]
name = "test"
version = "0"
background = ["O2", "O(3P)", "N2"]

[references]
ref = "test reference"

[[state]]
id = "O(1D)"

[[state]]
id = "O2(b, v={v})"
for_v = [0, 1]

[[reaction]]
id = "o1d_o2"
reactants = ["O(1D)", "O2"]
rate = { law = "arrhenius", a = 3.2e-11, ea_over_r = -67.0, units = "cm3 s-1" }
channels = [
  { yield = 0.4, products = ["O2(b, v=1)", "O(3P)"] },
  { yield = 0.55, products = ["O2(b)", "O(3P)"] },
]
reference = "ref"

[[reaction]]
id = "o2b_v{v}_n2"
for_v = [1, 1]
reactants = ["O2(b, v={v})", "N2"]
products = ["O2(b, v={v-1})", "N2"]
rate = { law = "constant", value = 2.0, exp_v = 0.5, units = "m3 s-1" }
reference = "ref"

[[photo_reaction]]
id = "o2_src"
reactant = "O2"
products = ["O(3P)", "O(1D)"]
rate_input = "J_O2_SRC"
reference = "ref"

[[transition]]
id = "a_band"
upper = "O2(b, v=0)"
lower = "O2"
einstein_a_s = 0.0758
wavelength_nm = 762.0
reference = "ref"
"#;

    #[test]
    fn loads_and_expands_families() {
        let mechanism = Mechanism::from_toml_str(MINIMAL).unwrap();

        let states: Vec<&str> = mechanism.states().iter().map(|s| s.id.as_str()).collect();
        assert_eq!(states, ["O(1D)", "O2(b)", "O2(b, v=1)"]);
        assert_eq!(mechanism.rate_inputs(), ["J_O2_SRC"]);

        let ids: Vec<&str> = mechanism
            .processes()
            .iter()
            .map(|p| p.id.as_str())
            .collect();
        assert_eq!(ids, ["o1d_o2", "o2b_v1_n2", "o2_src", "a_band"]);

        let quench = &mechanism.processes()[1];
        assert_eq!(
            quench.coefficient,
            Coefficient::Rate(RateLaw::constant(2.0 * 0.5_f64.exp()))
        );
        assert_eq!(quench.channels[0].products[0], Species::State(1));
    }

    #[test]
    fn photolysis_can_branch_into_channels() {
        let mechanism = with(
            r#"products = ["O(3P)", "O(1D)"]
rate_input = "J_O2_SRC""#,
            r#"channels = [
  { yield = 0.7, products = ["O(3P)", "O(1D)"] },
  { yield = 0.3, products = ["O(3P)", "O(3P)"] },
]
rate_input = "J_O2_SRC""#,
        )
        .unwrap();
        let photolysis = &mechanism.processes()[2];
        assert_eq!(photolysis.kind, ProcessKind::Photolysis);
        let fractions: Vec<f64> = photolysis.channels.iter().map(|c| c.fraction).collect();
        assert_eq!(fractions, [0.7, 0.3]);
        assert_eq!(photolysis.channels[0].products[1], Species::State(0));
    }

    #[test]
    fn photolysis_needs_products_or_channels() {
        let err = error_with(r#"products = ["O(3P)", "O(1D)"]
rate_input"#, "rate_input");
        assert!(err.contains("exactly one of 'products' or 'channels'"));
        let err = error_with(
            r#"rate_input = "J_O2_SRC""#,
            r#"rate_input = "J_O2_SRC"
channels = [{ yield = 1.0, products = ["O(1D)"] }]"#,
        );
        assert!(err.contains("exactly one of 'products' or 'channels'"));
    }

    #[test]
    fn rate_units_are_converted_to_si() {
        let mechanism = Mechanism::from_toml_str(MINIMAL).unwrap();
        let Coefficient::Rate(law) = &mechanism.processes()[0].coefficient else {
            panic!("expected a rate law");
        };
        assert!((law.a - 3.2e-17).abs() < 1e-30);
        assert_eq!(law.ea_over_r, -67.0);
    }

    #[test]
    fn bundled_mechanisms_load() {
        for name in Mechanism::bundled_names() {
            Mechanism::bundled(name).unwrap();
        }
        assert!(Mechanism::bundled("missing").is_err());
    }

    fn with(replace: &str, by: &str) -> Result<Mechanism> {
        assert!(
            MINIMAL.contains(replace),
            "test setup: '{replace}' not found"
        );
        Mechanism::from_toml_str(&MINIMAL.replacen(replace, by, 1))
    }

    fn error_with(replace: &str, by: &str) -> String {
        with(replace, by).unwrap_err().to_string()
    }

    #[test]
    fn unknown_species_is_an_error() {
        assert!(error_with(r#"["O(1D)", "O2"]"#, r#"["O(1D)", "O3"]"#).contains("'O3'"));
    }

    #[test]
    fn missing_reference_is_an_error() {
        assert!(
            error_with(
                r#"reference = "ref"
"#,
                r#"reference = "nope"
"#
            )
            .contains("'nope'")
        );
    }

    #[test]
    fn yields_above_one_are_an_error() {
        assert!(error_with("yield = 0.55", "yield = 0.65").contains("sum to"));
    }

    #[test]
    fn units_must_match_reaction_order() {
        assert!(error_with(r#"units = "cm3 s-1""#, r#"units = "s-1""#).contains("1 reactant"));
    }

    #[test]
    fn exp_v_needs_for_v() {
        let err = error_with(
            r#"law = "arrhenius", a = 3.2e-11"#,
            r#"law = "arrhenius", exp_v = 1.0, a = 3.2e-11"#,
        );
        assert!(err.contains("exp_v needs for_v"));
    }

    #[test]
    fn unknown_fields_are_an_error() {
        assert!(
            error_with("wavelength_nm = 762.0", "wavelenght_nm = 762.0").contains("wavelenght")
        );
    }

    #[test]
    fn two_state_reactants_are_rejected() {
        let err = error_with(r#"["O(1D)", "O2"]"#, r#"["O(1D)", "O2(b)"]"#);
        assert!(err.contains("2 state reactants"));
    }

    #[test]
    fn state_and_background_overlap_is_an_error() {
        assert!(error_with(r#"id = "O(1D)""#, r#"id = "O2(v=0)""#).contains("both a state"));
    }

    #[test]
    fn duplicate_process_ids_are_an_error() {
        assert!(error_with(r#"id = "a_band""#, r#"id = "o2_src""#).contains("used twice"));
    }
}
