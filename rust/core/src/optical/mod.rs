pub mod conversions;
pub mod rayleigh;
pub mod read_fwf_xs;
pub mod xsec_dbase;
use std::collections::HashMap;

use crate::constituent::StorageInputs;
use anyhow::Result;
use ndarray::*;

pub struct OpticalQuantities<'a> {
    pub cross_section: CowArray<'a, f64, Ix2>,
    pub ssa: CowArray<'a, f64, Ix2>,
    pub a1: Option<CowArray<'a, f64, Ix3>>,
    pub a2: Option<CowArray<'a, f64, Ix3>>,
    pub a3: Option<CowArray<'a, f64, Ix3>>,
    pub a4: Option<CowArray<'a, f64, Ix3>>,
    pub b1: Option<CowArray<'a, f64, Ix3>>,
    pub b2: Option<CowArray<'a, f64, Ix3>>,
}

pub trait OpticalProperty {
    fn optical_quantities<'a>(
        &'a self,
        _inputs: &dyn StorageInputs,
    ) -> Result<OpticalQuantities<'a>>;

    fn optical_derivatives<'a>(
        &'a self,
        inputs: &dyn StorageInputs,
    ) -> Result<HashMap<String, OpticalQuantities<'a>>>;
}
