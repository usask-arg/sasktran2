use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::prelude::*;

use crate::atmosphere::AtmosphereStorageAccess;
use crate::optical::traits::*;

use crate::interpolation::linear::linear_interpolating_matrix;

use ndarray::Ix2;

use anyhow::{Result, anyhow};

pub fn assign_absorber_derivatives<S1, S2, S3, S4>(
    mapping: &mut DerivMappingView,
    d_extinction: &ArrayBase<S1, Ix2>,
    d_ssa: &ArrayBase<S2, Ix2>,
    atmo_ssa: &ArrayBase<S3, Ix2>,
    atmo_extinction: &ArrayBase<S4, Ix2>,
) -> Result<()>
where
    S1: Data<Elem = f64> + Sync,
    S2: Data<Elem = f64> + Sync,
    S3: Data<Elem = f64> + Sync,
    S4: Data<Elem = f64> + Sync,
{
    let thread_pool = crate::threading::thread_pool()?;

    thread_pool.install(|| {
        Zip::from(mapping.d_extinction.columns_mut())
            .and(mapping.d_ssa.columns_mut())
            .and(d_extinction.columns())
            .and(d_ssa.columns())
            .and(atmo_ssa.columns())
            .and(atmo_extinction.columns())
            .par_for_each(
                |mut mapping_d_extinction_row,
                 mut mapping_d_ssa_row,
                 d_extinction_row,
                 d_ssa_row,
                 atmo_ssa_row,
                 atmo_extinction_row| {
                    mapping_d_extinction_row.assign(&d_extinction_row);
                    mapping_d_ssa_row.assign(&d_ssa_row);

                    Zip::from(mapping_d_ssa_row)
                        .and(atmo_ssa_row)
                        .and(d_extinction_row)
                        .and(atmo_extinction_row)
                        .for_each(
                            |mapping_d_ssa, atmo_ssa_val, d_extinction_val, atmo_extinction_val| {
                                *mapping_d_ssa = (*mapping_d_ssa - atmo_ssa_val) * d_extinction_val
                                    / atmo_extinction_val;
                            },
                        );
                },
            );
    });

    Ok(())
}

pub struct VMRAltitudeAbsorber<T>
where
    T: OpticalProperty,
{
    pub vmr: Array1<f64>,
    pub altitudes: Array1<f64>,
    pub optical_property: Option<T>,
    interp_mode: crate::interpolation::OutOfBoundsMode,
}

impl<T> VMRAltitudeAbsorber<T>
where
    T: OpticalProperty,
{
    pub fn new(altitudes: Array1<f64>, vmr: Array1<f64>) -> Self {
        VMRAltitudeAbsorber {
            vmr,
            altitudes,
            optical_property: None,
            interp_mode: crate::interpolation::OutOfBoundsMode::Extend,
        }
    }

    pub fn with_optical_property(&mut self, optical_property: T) -> &mut Self {
        self.optical_property = Some(optical_property);
        self
    }

    pub fn with_interp_mode(mut self, interp_mode: crate::interpolation::OutOfBoundsMode) -> Self {
        self.interp_mode = interp_mode;
        self
    }

    pub fn with_no_optical_property(&mut self) -> &mut Self {
        self.optical_property = None;
        self
    }
}

impl<T> Constituent for VMRAltitudeAbsorber<T>
where
    T: OpticalProperty,
{
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (inputs, outputs) = storage.split_inputs_outputs();

        let altitudes_m = inputs.altitude_m();

        let interp_matrix = linear_interpolating_matrix(
            &self.altitudes,
            &altitudes_m,
            crate::interpolation::OutOfBoundsMode::Extend,
        );
        let interp_vmr = interp_matrix.dot(&self.vmr);

        let mut aux_inputs: HashMap<String, Array1<f64>> = HashMap::new();
        aux_inputs.insert("vmr".to_string(), interp_vmr.clone());

        let optical_prop = self
            .optical_property
            .as_ref()
            .ok_or_else(|| anyhow!("Optical property not set"))?;

        let optical_quants = optical_prop.optical_quantities(inputs, &aux_inputs)?;
        let cross_section = &optical_quants.cross_section;

        let eqn_state = inputs.dry_air_numberdensity_dict();

        let number_density = eqn_state
            .get("N")
            .ok_or_else(|| anyhow!("Number density for N not found in air_numberdensity_dict"))?;

        let mut extinction = outputs.mut_view().total_extinction;

        // Here extinction is (geometry, wavelength), with geometry being the fastest changing index
        // because of fortran ordering, so we loop this way even though it is not the most natural

        let thread_pool = crate::threading::thread_pool()?;

        thread_pool.install(|| {
            Zip::from(extinction.axis_iter_mut(Axis(1)))
                .and(cross_section.axis_iter(Axis(1)))
                .par_for_each(|ext_slice, xs_slice| {
                    // Then loop over the inner dimension
                    Zip::from(ext_slice)
                        .and(xs_slice)
                        .and(&interp_vmr)
                        .and(number_density)
                        .for_each(|ext, xs, vmr, n| {
                            *ext += *xs * *vmr * *n;
                        });
                });
        });

        Ok(())
    }

    fn register_derivatives<'b>(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()> {
        let (inputs, outputs, deriv_generator) = storage.split_inputs_outputs_deriv();
        let outputs = outputs.view();

        let eqn_state = inputs.dry_air_numberdensity_dict();

        let number_density = eqn_state
            .get("N")
            .ok_or_else(|| anyhow!("Number density for N not found in air_numberdensity_dict"))?;

        let optical_prop = self
            .optical_property
            .as_ref()
            .ok_or_else(|| anyhow!("Optical property not set"))?;

        let altitudes_m = inputs.altitude_m();
        let interp_matrix = linear_interpolating_matrix(
            &self.altitudes,
            &altitudes_m,
            crate::interpolation::OutOfBoundsMode::Extend,
        );

        let interp_vmr = interp_matrix.dot(&self.vmr);

        let optical_quants = optical_prop.optical_quantities(inputs, &NullAuxInputs {})?;
        let cross_section = optical_quants.cross_section;

        let mut aux_inputs: HashMap<String, Array1<f64>> = HashMap::new();
        aux_inputs.insert("vmr".to_string(), interp_vmr.clone());

        let wf_name = format!("wf_{constituent_name}_vmr");
        let mut deriv_mapping = deriv_generator.get_derivative_mapping(&wf_name);
        let mut mapping = deriv_mapping.mut_view();

        assign_absorber_derivatives(
            &mut mapping,
            &cross_section,
            &optical_quants.ssa,
            &outputs.ssa,
            &outputs.total_extinction,
        )?;

        let interp_dim = format!("{constituent_name}_altitude");
        deriv_mapping.set_interp_dim(&interp_dim);

        // interp_matrix is (num_altitudes, num_vmr)
        // number_density is (num_altitudes)
        // want to multiply each row of interp_matrix by the corresponding element of number_density
        let mut deriv_interpolator = interp_matrix.clone();
        Zip::from(deriv_interpolator.rows_mut())
            .and(number_density)
            .for_each(|mut row, number_density_val| {
                row *= *number_density_val;
            });

        deriv_mapping.set_interpolator(&deriv_interpolator);

        // Construct the extra derivative mappings for the number density adjustments
        let mut deriv_names: Vec<String> = vec![];
        let mut d_vals: Vec<Array1<f64>> = vec![];

        if inputs.calculate_pressure_derivative() {
            deriv_names.push("pressure_pa".to_string());
            d_vals.push(
                eqn_state
                    .get("dN_dP")
                    .ok_or_else(|| {
                        anyhow!("Derivative for dN_dp not found in air_numberdensity_dict")
                    })?
                    .to_owned(),
            );
        }
        if inputs.calculate_temperature_derivative() {
            deriv_names.push("temperature_k".to_string());
            d_vals.push(
                eqn_state
                    .get("dN_dT")
                    .ok_or_else(|| {
                        anyhow!("Derivative for dN_dT not found in air_numberdensity_dict")
                    })?
                    .to_owned(),
            );
        }
        if inputs.calculate_specific_humidity_derivative() {
            deriv_names.push("specific_humidity".to_string());
            d_vals.push(
                eqn_state
                    .get("dN_dsh")
                    .ok_or_else(|| {
                        anyhow!("Derivative for dN_dsh not found in air_numberdensity_dict")
                    })?
                    .to_owned(),
            );
        }

        deriv_names
            .iter()
            .zip(d_vals.iter())
            .for_each(|(deriv_name, vert_factor)| {
                let mapping_name = format!("wf_{constituent_name}_{deriv_name}");
                let mut mapping = deriv_generator.get_derivative_mapping(&mapping_name);
                let mut mapping_view = mapping.mut_view();

                let _ = assign_absorber_derivatives(
                    &mut mapping_view,
                    &cross_section,
                    &optical_quants.ssa,
                    &outputs.ssa,
                    &outputs.total_extinction,
                );

                mapping.set_interp_dim("altitude");
                mapping.set_assign_name(format!("wf_{deriv_name}").as_str());

                let diagonal = (&interp_vmr) * vert_factor;

                let interpolator = Array2::from_diag(&diagonal);
                mapping.set_interpolator(&interpolator);
            });

        if !deriv_names.is_empty() {
            let d_aq = optical_prop.optical_derivatives(inputs, &NullAuxInputs {})?;

            d_aq.iter().for_each(|(key, val)| {
                let mapping_name = format!("wf_{constituent_name}_{key}_xs");
                let mut mapping = deriv_generator.get_derivative_mapping(&mapping_name);
                let mut mapping_view = mapping.mut_view();

                let _ = assign_absorber_derivatives(
                    &mut mapping_view,
                    &val.cross_section,
                    &optical_quants.ssa,
                    &outputs.ssa,
                    &outputs.total_extinction,
                );

                mapping.set_interp_dim("altitude");
                mapping.set_assign_name(format!("wf_{key}").as_str());

                let diagonal = (&interp_vmr) * number_density;

                let interpolator = Array2::from_diag(&diagonal);
                mapping.set_interpolator(&interpolator);
            });
        }

        Ok(())
    }
}
