use ndarray::{Array1, Axis, Zip};
use std::sync::Arc;

use crate::atmosphere::AtmosphereStorageAccess;
use crate::constituent::{Constituent, StorageInputs};
use crate::optical::OpticalProperty;

use super::StorageOutputs;
use crate::interpolation::linear::linear_interpolating_matrix;

use anyhow::{Result, anyhow};

pub struct VMRAltitudeAbsorber {
    pub vmr: Array1<f64>,
    pub altitudes: Array1<f64>,
    optical_property: Arc<dyn OpticalProperty>,
}

impl VMRAltitudeAbsorber {
    pub fn new(
        vmr: Array1<f64>,
        altitudes: Array1<f64>,
        optical_property: Arc<dyn OpticalProperty>,
    ) -> Self {
        VMRAltitudeAbsorber {
            vmr,
            altitudes,
            optical_property,
        }
    }

    fn get_optical_property(&self) -> &dyn OpticalProperty {
        self.optical_property.as_ref()
    }
}

impl Constituent for VMRAltitudeAbsorber {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (inputs, outputs) = storage.split_inputs_outputs();

        let optical_property = self.get_optical_property();

        let cross_section = optical_property.optical_quantities(inputs)?.cross_section;

        let eqn_state = inputs.air_numberdensity_dict();

        let number_density = eqn_state
            .get("N")
            .ok_or_else(|| anyhow!("Number density for N not found in air_numberdensity_dict"))?;

        let altitudes_m = inputs.altitude_m();

        let interp_matrix = linear_interpolating_matrix(
            &self.altitudes,
            &altitudes_m,
            crate::interpolation::OutOfBoundsMode::Extend,
        );

        let interp_vmr = interp_matrix.dot(&self.vmr);

        let mut extinction = outputs.mut_view().total_extinction;

        Zip::from(extinction.axis_iter_mut(Axis(0)))
            .and(number_density)
            .and(&interp_vmr)
            .and(cross_section.axis_iter(Axis(0)))
            .par_for_each(|ext_row, num_dens, vmr, xs_row| {
                Zip::from(ext_row).and(xs_row).for_each(|ext, xs| {
                    *ext = *num_dens * *vmr * xs;
                });
            });

        Ok(())
    }

    fn register_derivatives<'b>(
        &self,
        _storage: &mut impl AtmosphereStorageAccess,
        _constituent_name: &str,
    ) -> Result<()> {
        todo!()
    }
}
