//! Serde mirror of the mechanism TOML format. These types only describe the
//! file layout; [`super::Mechanism`] holds the validated, compiled form.

use serde::Deserialize;
use std::collections::BTreeMap;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct MechanismFile {
    pub mechanism: Header,
    #[serde(default)]
    pub references: BTreeMap<String, String>,
    #[serde(default, rename = "state")]
    pub states: Vec<StateEntry>,
    #[serde(default, rename = "reaction")]
    pub reactions: Vec<ReactionEntry>,
    #[serde(default, rename = "photo_reaction")]
    pub photo_reactions: Vec<PhotoReactionEntry>,
    #[serde(default, rename = "transition")]
    pub transitions: Vec<TransitionEntry>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Header {
    pub name: String,
    pub version: String,
    #[serde(default)]
    pub description: String,
    /// Species whose number densities are inputs rather than unknowns.
    #[serde(default)]
    pub background: Vec<String>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct StateEntry {
    pub id: String,
    #[serde(default)]
    pub description: String,
    #[serde(default)]
    pub for_v: Option<[u32; 2]>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct ChannelEntry {
    #[serde(rename = "yield")]
    pub fraction: f64,
    pub products: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct RateEntry {
    pub law: String,
    pub units: String,
    pub value: Option<f64>,
    pub a: Option<f64>,
    pub n: Option<f64>,
    pub t0: Option<f64>,
    pub ea_over_r: Option<f64>,
    /// Multiplies the rate by exp(exp_v * v); only valid with `for_v`.
    pub exp_v: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct ReactionEntry {
    pub id: String,
    pub reactants: Vec<String>,
    #[serde(default)]
    pub products: Option<Vec<String>>,
    #[serde(default)]
    pub channels: Option<Vec<ChannelEntry>>,
    pub rate: RateEntry,
    pub reference: String,
    #[serde(default)]
    pub for_v: Option<[u32; 2]>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct PhotoReactionEntry {
    pub id: String,
    pub reactant: String,
    #[serde(default)]
    pub products: Option<Vec<String>>,
    #[serde(default)]
    pub channels: Option<Vec<ChannelEntry>>,
    /// Name of the per-molecule rate [s^-1] supplied with the column inputs.
    pub rate_input: String,
    pub reference: String,
    #[serde(default)]
    pub for_v: Option<[u32; 2]>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct TransitionEntry {
    pub id: String,
    pub upper: String,
    pub lower: String,
    pub einstein_a_s: f64,
    #[serde(default)]
    pub wavelength_nm: Option<f64>,
    pub reference: String,
    #[serde(default)]
    pub for_v: Option<[u32; 2]>,
}
