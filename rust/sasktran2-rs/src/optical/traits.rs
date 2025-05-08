use crate::atmosphere::traits::*;
use crate::optical::storage::*;
use crate::prelude::*;
use ndarray::{CowArray, Ix1};

pub trait AuxOpticalInputs {
    fn get_parameter(&self, name: &str) -> Option<CowArray<'_, f64, Ix1>>;
}

// Null implementation
pub struct NullAuxInputs;

impl AuxOpticalInputs for NullAuxInputs {
    fn get_parameter(&self, _name: &str) -> Option<CowArray<'_, f64, Ix1>> {
        None
    }
}

impl AuxOpticalInputs for HashMap<String, Array1<f64>> {
    fn get_parameter(&self, name: &str) -> Option<CowArray<'_, f64, Ix1>> {
        if let Some(param) = self.get(name) {
            let param = param.view();
            let carray = CowArray::from(param);
            return Some(carray);
        }
        None
    }
}

pub trait OpticalProperty {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> Result<()>;

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> Result<()>;
}

pub trait OpticalPropertyExt: OpticalProperty {
    fn optical_quantities(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
    ) -> Result<OpticalQuantities> {
        let mut out = OpticalQuantities::default();
        self.optical_quantities_emplace(inputs, aux_inputs, &mut out)?;
        Ok(out)
    }

    fn optical_derivatives(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
    ) -> Result<HashMap<String, OpticalQuantities>> {
        let mut out = HashMap::new();
        self.optical_derivatives_emplace(inputs, aux_inputs, &mut out)?;
        Ok(out)
    }
}

impl<T: OpticalProperty + ?Sized> OpticalPropertyExt for T {}

pub fn param_from_storage_or_aux<'a>(
    inputs: &'a dyn StorageInputs,
    aux_inputs: &'a dyn AuxOpticalInputs,
    name: &'a str,
) -> Result<CowArray<'a, f64, Ix1>> {
    if let Some(param) = inputs.get_parameter(name) {
        Ok(param.into())
    } else if let Some(aux_param) = aux_inputs.get_parameter(name) {
        Ok(aux_param)
    } else {
        Err(anyhow!(
            "Parameter {} not found in inputs or aux inputs",
            name
        ))
    }
}
