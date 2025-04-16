pub mod conversions;
pub mod rayleigh;
pub mod read_fwf_xs;
pub mod xsec_dbase;
use std::collections::HashMap;

use crate::constituent::StorageInputs;
use anyhow::Result;
use ndarray::*;

pub struct OpticalQuantities {
    pub cross_section: Array2<f64>,
    pub ssa: Array2<f64>,
    pub a1: Option<Array3<f64>>,
    pub a2: Option<Array3<f64>>,
    pub a3: Option<Array3<f64>>,
    pub a4: Option<Array3<f64>>,
    pub b1: Option<Array3<f64>>,
    pub b2: Option<Array3<f64>>,
}

impl Default for OpticalQuantities {
    fn default() -> Self {
        Self {
            cross_section: Array2::zeros((0, 0)),
            ssa: Array2::zeros((0, 0)),
            a1: None,
            a2: None,
            a3: None,
            a4: None,
            b1: None,
            b2: None,
        }
    }
}

impl OpticalQuantities {
    pub fn new(num_geometry: usize, num_wavelengths: usize) -> Self {
        let mut default = Self::default();
        default.resize(num_geometry, num_wavelengths);

        default
    }

    pub fn resize(&mut self, num_geometry: usize, num_wavelengths: usize) -> &mut Self {
        if self.cross_section.dim() != (num_geometry, num_wavelengths) {
            self.cross_section = Array2::zeros((num_geometry, num_wavelengths));
        }
        if self.ssa.dim() != (num_geometry, num_wavelengths) {
            self.ssa = Array2::zeros((num_geometry, num_wavelengths));
        }

        self
    }

    pub fn set_zero(&mut self) {
        self.cross_section.fill(0.0);
        self.ssa.fill(0.0);
        if let Some(a1) = &mut self.a1 {
            a1.fill(0.0);
        }
        if let Some(a2) = &mut self.a2 {
            a2.fill(0.0);
        }
        if let Some(a3) = &mut self.a3 {
            a3.fill(0.0);
        }
        if let Some(a4) = &mut self.a4 {
            a4.fill(0.0);
        }
        if let Some(b1) = &mut self.b1 {
            b1.fill(0.0);
        }
        if let Some(b2) = &mut self.b2 {
            b2.fill(0.0);
        }
    }
}

pub trait OpticalProperty {
    fn optical_quantities_emplace(
        &self,
        _inputs: &dyn StorageInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> Result<()>;

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> Result<()>;
}

pub trait OpticalPropertyExt: OpticalProperty {
    fn optical_quantities(&self, inputs: &dyn StorageInputs) -> Result<OpticalQuantities> {
        let mut out = OpticalQuantities::default();
        self.optical_quantities_emplace(inputs, &mut out)?;
        Ok(out)
    }

    fn optical_derivatives(
        &self,
        inputs: &dyn StorageInputs,
    ) -> Result<HashMap<String, OpticalQuantities>> {
        let mut out = HashMap::new();
        self.optical_derivatives_emplace(inputs, &mut out)?;
        Ok(out)
    }
}

impl<T: OpticalProperty + ?Sized> OpticalPropertyExt for T {}
