//! Excited-state (non-LTE) population kinetics and emission rates.
//!
//! The crate is a self-contained kinetics library. It takes atmospheric
//! profiles and photolysis rates as plain arrays and returns state
//! populations and emission rates; it does no radiative transfer and has no
//! dependency on the rest of SASKTRAN2. See
//! `docs/sphinx/source/developer/nlte_plan.md` for the design.
//!
//! A [`mechanism::Mechanism`] (states, background species, reactions,
//! photolysis and radiative transitions) is loaded from a TOML file, and
//! [`solver::solve_steady_state`] gives the state populations, process rates
//! and production/loss budgets on a column of levels.
//!
//! `models`, `types` and `emission` hold the earlier photochemistry code that
//! moved here from `sasktran2_rs::photchem`.

pub mod emission;
mod linalg;
pub mod mechanism;
mod prelude;
pub mod rates;
pub mod solver;
