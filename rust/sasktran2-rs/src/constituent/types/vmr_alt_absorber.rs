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
                                *mapping_d_ssa = if *atmo_extinction_val == 0.0 {
                                    0.0
                                } else {
                                    (*mapping_d_ssa - atmo_ssa_val) * d_extinction_val
                                        / atmo_extinction_val
                                };
                            },
                        );
                },
            );
    });

    Ok(())
}

fn add_absorber_to_atmosphere<T>(
    vmr: &Array1<f64>,
    optical_property: &T,
    storage: &mut impl AtmosphereStorageAccess,
) -> Result<()>
where
    T: OpticalProperty,
{
    let (inputs, outputs) = storage.split_inputs_outputs();
    let num_locations = inputs.altitude_m().len();
    anyhow::ensure!(
        vmr.len() == num_locations,
        "VMR spatial dimension ({}) does not match atmosphere ({num_locations})",
        vmr.len()
    );

    let mut aux_inputs: HashMap<String, Array1<f64>> = HashMap::new();
    aux_inputs.insert("vmr".to_string(), vmr.clone());

    let optical_quants = optical_property.optical_quantities(inputs, &aux_inputs)?;
    let cross_section = &optical_quants.cross_section;
    let eqn_state = inputs.dry_air_numberdensity_dict();
    let number_density = eqn_state
        .get("N")
        .ok_or_else(|| anyhow!("Number density for N not found in air_numberdensity_dict"))?;
    let mut extinction = outputs.mut_view().total_extinction;
    let thread_pool = crate::threading::thread_pool()?;

    thread_pool.install(|| {
        Zip::from(extinction.axis_iter_mut(Axis(1)))
            .and(cross_section.axis_iter(Axis(1)))
            .par_for_each(|ext_slice, xs_slice| {
                Zip::from(ext_slice)
                    .and(xs_slice)
                    .and(vmr)
                    .and(number_density)
                    .for_each(|ext, xs, vmr, n| {
                        *ext += *xs * *vmr * *n;
                    });
            });
    });

    Ok(())
}

fn scale_absorber_derivatives(mapping: &mut DerivMappingView, scale: &Array1<f64>) {
    Zip::from(mapping.d_extinction.rows_mut())
        .and(mapping.d_ssa.rows_mut())
        .and(scale)
        .for_each(|mut d_extinction, mut d_ssa, scale| {
            d_extinction *= *scale;
            d_ssa *= *scale;
        });
}

#[derive(Clone, Copy)]
enum VmrDerivativeGrid<'a> {
    Altitude(&'a Array2<f64>),
    Native,
}

fn register_absorber_derivatives<T>(
    vmr: &Array1<f64>,
    vmr_grid: VmrDerivativeGrid<'_>,
    optical_property: &T,
    storage: &mut impl AtmosphereStorageAccess,
    constituent_name: &str,
) -> Result<()>
where
    T: OpticalProperty,
{
    let (inputs, outputs, deriv_generator) = storage.split_inputs_outputs_deriv();
    let outputs = outputs.view();
    let num_locations = inputs.altitude_m().len();
    anyhow::ensure!(
        vmr.len() == num_locations,
        "VMR spatial dimension ({}) does not match atmosphere ({num_locations})",
        vmr.len()
    );

    let eqn_state = inputs.dry_air_numberdensity_dict();
    let number_density = eqn_state
        .get("N")
        .ok_or_else(|| anyhow!("Number density for N not found in air_numberdensity_dict"))?;
    let mut aux_inputs: HashMap<String, Array1<f64>> = HashMap::new();
    aux_inputs.insert("vmr".to_string(), vmr.clone());
    let optical_quants = optical_property.optical_quantities(inputs, &aux_inputs)?;
    let cross_section = &optical_quants.cross_section;

    let wf_name = format!("wf_{constituent_name}_vmr");
    let mut deriv_mapping = deriv_generator.get_derivative_mapping(&wf_name);
    {
        let mut mapping = deriv_mapping.mut_view();
        assign_absorber_derivatives(
            &mut mapping,
            cross_section,
            &optical_quants.ssa,
            &outputs.ssa,
            &outputs.total_extinction,
        )?;

        if matches!(vmr_grid, VmrDerivativeGrid::Native) {
            scale_absorber_derivatives(&mut mapping, number_density);
        }
    }

    match vmr_grid {
        VmrDerivativeGrid::Altitude(interp_matrix) => {
            deriv_mapping.set_interp_dim(&format!("{constituent_name}_altitude"));
            let mut deriv_interpolator = interp_matrix.clone();
            Zip::from(deriv_interpolator.rows_mut())
                .and(number_density)
                .for_each(|mut row, number_density| row *= *number_density);
            deriv_mapping.set_interpolator(&deriv_interpolator);
        }
        VmrDerivativeGrid::Native => deriv_mapping.set_interp_dim("location"),
    }

    let mut state_derivatives: Vec<(&str, Array1<f64>)> = vec![];
    if inputs.calculate_pressure_derivative() {
        state_derivatives.push((
            "pressure_pa",
            eqn_state
                .get("dN_dP")
                .ok_or_else(|| anyhow!("Derivative for dN_dP not found"))?
                .to_owned(),
        ));
    }
    if inputs.calculate_temperature_derivative() {
        state_derivatives.push((
            "temperature_k",
            eqn_state
                .get("dN_dT")
                .ok_or_else(|| anyhow!("Derivative for dN_dT not found"))?
                .to_owned(),
        ));
    }
    if inputs.calculate_specific_humidity_derivative() {
        state_derivatives.push((
            "specific_humidity",
            eqn_state
                .get("dN_dsh")
                .ok_or_else(|| anyhow!("Derivative for dN_dsh not found"))?
                .to_owned(),
        ));
    }

    for (deriv_name, state_factor) in &state_derivatives {
        let mapping_name = format!("wf_{constituent_name}_{deriv_name}");
        let mut mapping = deriv_generator.get_derivative_mapping(&mapping_name);
        {
            let mut mapping_view = mapping.mut_view();
            assign_absorber_derivatives(
                &mut mapping_view,
                cross_section,
                &optical_quants.ssa,
                &outputs.ssa,
                &outputs.total_extinction,
            )?;
        }
        mapping.set_interp_dim("altitude");
        mapping.set_assign_name(&format!("wf_{deriv_name}"));
        mapping.set_interpolator(&Array2::from_diag(&(vmr * state_factor)));
    }

    if !state_derivatives.is_empty() {
        for (key, val) in optical_property
            .optical_derivatives(inputs, &aux_inputs)?
            .iter()
        {
            let mapping_name = format!("wf_{constituent_name}_{key}_xs");
            let mut mapping = deriv_generator.get_derivative_mapping(&mapping_name);
            {
                let mut mapping_view = mapping.mut_view();
                assign_absorber_derivatives(
                    &mut mapping_view,
                    &val.cross_section,
                    &optical_quants.ssa,
                    &outputs.ssa,
                    &outputs.total_extinction,
                )?;
            }
            mapping.set_interp_dim("altitude");
            mapping.set_assign_name(&format!("wf_{key}"));
            mapping.set_interpolator(&Array2::from_diag(&(vmr * number_density)));
        }
    }

    Ok(())
}

pub struct VMRAbsorber<T>
where
    T: OpticalProperty,
{
    pub vmr: Array1<f64>,
    pub optical_property: Option<T>,
}

impl<T> VMRAbsorber<T>
where
    T: OpticalProperty,
{
    pub fn new(vmr: Array1<f64>) -> Self {
        Self {
            vmr,
            optical_property: None,
        }
    }

    pub fn with_optical_property(&mut self, optical_property: T) -> &mut Self {
        self.optical_property = Some(optical_property);
        self
    }

    pub fn with_no_optical_property(&mut self) -> &mut Self {
        self.optical_property = None;
        self
    }
}

impl<T> Constituent for VMRAbsorber<T>
where
    T: OpticalProperty,
{
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        add_absorber_to_atmosphere(
            &self.vmr,
            self.optical_property
                .as_ref()
                .ok_or_else(|| anyhow!("Optical property not set"))?,
            storage,
        )
    }

    fn register_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()> {
        register_absorber_derivatives(
            &self.vmr,
            VmrDerivativeGrid::Native,
            self.optical_property
                .as_ref()
                .ok_or_else(|| anyhow!("Optical property not set"))?,
            storage,
            constituent_name,
        )
    }
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
        let (inputs, _) = storage.split_inputs_outputs();

        let altitudes_m = inputs.altitude_m();

        let interp_matrix =
            linear_interpolating_matrix(&self.altitudes, &altitudes_m, self.interp_mode);
        let interp_vmr = interp_matrix.dot(&self.vmr);
        add_absorber_to_atmosphere(
            &interp_vmr,
            self.optical_property
                .as_ref()
                .ok_or_else(|| anyhow!("Optical property not set"))?,
            storage,
        )
    }

    fn register_derivatives<'b>(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()> {
        let (inputs, _) = storage.split_inputs_outputs();
        let altitudes_m = inputs.altitude_m();
        let interp_matrix =
            linear_interpolating_matrix(&self.altitudes, &altitudes_m, self.interp_mode);
        let interp_vmr = interp_matrix.dot(&self.vmr);
        register_absorber_derivatives(
            &interp_vmr,
            VmrDerivativeGrid::Altitude(&interp_matrix),
            self.optical_property
                .as_ref()
                .ok_or_else(|| anyhow!("Optical property not set"))?,
            storage,
            constituent_name,
        )
    }
}
