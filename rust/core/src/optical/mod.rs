pub mod conversions;
pub mod rayleigh;
pub mod read_fwf_xs;
pub mod xsec_dbase;
use crate::constituent::StorageInputs;
use anyhow::Result;
use ndarray::*;

pub struct OpticalQuantities<S>
where
    S: Data<Elem = f64>,
{
    pub cross_section: ArrayBase<S, Ix2>,
    pub ssa: ArrayBase<S, Ix2>,
    pub a1: Option<ArrayBase<S, Ix3>>,
    pub a2: Option<ArrayBase<S, Ix3>>,
    pub a3: Option<ArrayBase<S, Ix3>>,
    pub a4: Option<ArrayBase<S, Ix3>>,
    pub b1: Option<ArrayBase<S, Ix3>>,
    pub b2: Option<ArrayBase<S, Ix3>>,
}

pub trait OpticalProperty: Send + Sync {
    fn optical_quantities(
        &self,
        inputs: &dyn StorageInputs,
    ) -> Result<OpticalQuantities<OwnedRepr<f64>>>;
}
