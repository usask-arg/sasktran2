use crate::constituent::{DerivMappingGenerator, StorageInputs, StorageOutputs};

pub trait AtmosphereStorageAccess {
    fn split_inputs_outputs(&mut self) -> (&impl StorageInputs, &mut impl StorageOutputs);
    fn split_inputs_outputs_deriv<'a>(
        &'a self,
    ) -> (
        &'a impl StorageInputs,
        &'a impl StorageOutputs,
        &'a impl DerivMappingGenerator<'a>,
    );
}
