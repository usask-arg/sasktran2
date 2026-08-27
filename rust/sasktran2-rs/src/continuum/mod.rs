//! Molecular continuum absorption models.
//!
//! This module is intentionally independent of the SASKTRAN2 engine and FFI
//! layers. Continuum models operate on scalar atmospheric states and return
//! spectra together with their atmospheric Jacobians.

pub mod mt_ckd;
