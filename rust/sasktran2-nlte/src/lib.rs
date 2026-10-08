//! Excited-state (non-LTE) population kinetics and emission rates.
//!
//! The crate is a self-contained kinetics library. It takes atmospheric
//! profiles and photolysis rates as plain arrays and returns state
//! populations and emission rates; it does no radiative transfer and has no
//! dependency on the rest of SASKTRAN2. See
//! `docs/sphinx/source/developer/nlte_plan.md` for the design.
//!
//! The current contents are the photochemistry models that previously lived
//! in `sasktran2_rs::photchem`, moved without behaviour changes.

pub mod emission;
mod linalg;
pub mod models;
mod prelude;
pub mod types;
