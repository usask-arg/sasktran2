//! Photochemistry now lives in the `sasktran2-nlte` crate. This module keeps
//! the `sasktran2_rs::photchem` paths working and adds the HITRAN adapters in
//! [`emission`].

pub mod emission;

pub use sasktran2_nlte::{models, types};
