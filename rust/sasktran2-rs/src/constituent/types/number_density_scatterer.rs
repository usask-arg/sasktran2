use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::interpolation::linear::linear_interpolating_matrix;
use crate::optical::traits::*;
use crate::prelude::*;

use anyhow::{Result, anyhow};

pub struct NumberDensityScatterer<T>
where
    T: OpticalProperty,
{
    pub number_density: Array1<f64>,
    pub altitudes: Array1<f64>,
    pub optical_property: Option<T>,
    pub aux_inputs: HashMap<String, Array1<f64>>,
    pub vertical_deriv_factor: Array1<f64>,
    pub d_vertical_deriv_factor: HashMap<String, Array1<f64>>,
    pub wf_name: String,
    interp_mode: crate::interpolation::OutOfBoundsMode,
}

impl<T> NumberDensityScatterer<T>
where
    T: OpticalProperty,
{
    pub fn new(altitudes: Array1<f64>, number_density: Array1<f64>) -> Self {
        let vertical_deriv_factor = Array1::ones(number_density.len());

        NumberDensityScatterer {
            number_density,
            altitudes,
            optical_property: None,
            aux_inputs: HashMap::new(),
            vertical_deriv_factor,
            d_vertical_deriv_factor: HashMap::new(),
            wf_name: "number_density".to_string(),
            interp_mode: crate::interpolation::OutOfBoundsMode::Zero,
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

    pub fn with_interp_mode(mut self, interp_mode: crate::interpolation::OutOfBoundsMode) -> Self {
        self.interp_mode = interp_mode;
        self
    }

    pub fn set_aux_inputs(&mut self, aux_inputs: HashMap<String, Array1<f64>>) {
        self.aux_inputs = aux_inputs;
    }

    pub fn set_vertical_deriv_factor(&mut self, vertical_deriv_factor: Array1<f64>) {
        self.vertical_deriv_factor = vertical_deriv_factor;
    }

    pub fn set_d_vertical_deriv_factor(
        &mut self,
        d_vertical_deriv_factor: HashMap<String, Array1<f64>>,
    ) {
        self.d_vertical_deriv_factor = d_vertical_deriv_factor;
    }

    pub fn set_wf_name(&mut self, wf_name: String) {
        self.wf_name = wf_name;
    }

    fn interpolated_aux_inputs(&self, interp_matrix: &Array2<f64>) -> HashMap<String, Array1<f64>> {
        self.aux_inputs
            .iter()
            .map(|(key, val)| (key.clone(), interp_matrix.dot(val)))
            .collect()
    }
}

impl<T> Constituent for NumberDensityScatterer<T>
where
    T: OpticalProperty,
{
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (inputs, outputs) = storage.split_inputs_outputs();

        let interp_matrix =
            linear_interpolating_matrix(&self.altitudes, &inputs.altitude_m(), self.interp_mode);
        let interp_numden = interp_matrix.dot(&self.number_density);
        let aux_inputs = self.interpolated_aux_inputs(&interp_matrix);

        let optical_prop = self
            .optical_property
            .as_ref()
            .ok_or_else(|| anyhow!("Optical property not set"))?;

        let optical_quants = optical_prop.optical_quantities(inputs, &aux_inputs)?;
        let legendre = optical_quants.legendre.as_ref().ok_or_else(|| {
            anyhow!("Scattering optical property did not provide legendre coefficients")
        })?;

        let mut outputs = outputs.mut_view();

        Zip::indexed(&optical_quants.cross_section)
            .and(&optical_quants.ssa)
            .for_each(|(geo_idx, wavelength_idx), species_ext, species_scat| {
                let nd = interp_numden[geo_idx];
                outputs.total_extinction[[geo_idx, wavelength_idx]] += species_ext * nd;
                outputs.ssa[[geo_idx, wavelength_idx]] += species_scat * nd;

                for leg_idx in 0..legendre.dim().2 {
                    outputs.legendre[[leg_idx, geo_idx, wavelength_idx]] +=
                        species_scat * nd * legendre[[geo_idx, wavelength_idx, leg_idx]];
                }
            });

        Ok(())
    }

    fn register_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()> {
        let (inputs, outputs, deriv_generator) = storage.split_inputs_outputs_deriv();
        let outputs = outputs.view();

        let interp_matrix =
            linear_interpolating_matrix(&self.altitudes, &inputs.altitude_m(), self.interp_mode);
        let aux_inputs = self.interpolated_aux_inputs(&interp_matrix);

        let optical_prop = self
            .optical_property
            .as_ref()
            .ok_or_else(|| anyhow!("Optical property not set"))?;

        let optical_quants = optical_prop.optical_quantities(inputs, &aux_inputs)?;
        let optical_legendre = optical_quants.legendre.as_ref().ok_or_else(|| {
            anyhow!("Scattering optical property did not provide legendre coefficients")
        })?;

        let mut species_ssa = optical_quants.ssa.clone();
        Zip::from(&mut species_ssa)
            .and(&optical_quants.cross_section)
            .for_each(|ssa, ext| {
                *ssa /= *ext;
                if !ssa.is_finite() {
                    *ssa = 1.0;
                }
            });

        let wf_name = format!("wf_{constituent_name}_{}", self.wf_name);
        let mut deriv_mapping = deriv_generator
            .get_derivative_mapping(&wf_name)
            .with_scatterer();
        {
            let mut mapping = deriv_mapping.mut_view();
            let d_legendre = mapping
                .d_legendre
                .as_mut()
                .ok_or_else(|| anyhow!("Scatterer derivative mapping missing d_legendre"))?;
            let scat_factor = mapping
                .scat_factor
                .as_mut()
                .ok_or_else(|| anyhow!("Scatterer derivative mapping missing scat_factor"))?;

            Zip::indexed(&optical_quants.cross_section).for_each(
                |(geo_idx, wavelength_idx), species_ext| {
                    let species_ssa_val = species_ssa[[geo_idx, wavelength_idx]];
                    let total_ext = outputs.total_extinction[[geo_idx, wavelength_idx]];
                    let total_ssa = outputs.ssa[[geo_idx, wavelength_idx]];

                    mapping.d_extinction[[geo_idx, wavelength_idx]] += *species_ext;
                    mapping.d_ssa[[geo_idx, wavelength_idx]] +=
                        species_ext * (species_ssa_val - total_ssa) / total_ext;
                    scat_factor[[geo_idx, wavelength_idx]] +=
                        species_ssa_val * species_ext / (total_ssa * total_ext);

                    for leg_idx in 0..optical_legendre.dim().2 {
                        d_legendre[[leg_idx, geo_idx, wavelength_idx]] += optical_legendre
                            [[geo_idx, wavelength_idx, leg_idx]]
                            - outputs.legendre[[leg_idx, geo_idx, wavelength_idx]];
                    }
                },
            );
        }

        let mut density_interpolator = interp_matrix.clone();
        Zip::from(density_interpolator.columns_mut())
            .and(&self.vertical_deriv_factor)
            .for_each(|mut col, factor| {
                col *= *factor;
            });
        deriv_mapping.set_interpolator(&density_interpolator);
        deriv_mapping.set_interp_dim(&format!("{constituent_name}_altitude"));

        let optical_derivs = optical_prop.optical_derivatives(inputs, &aux_inputs)?;

        for (key, val) in optical_derivs.iter() {
            let mapping_name = format!("wf_{constituent_name}_{key}");
            let mut deriv_mapping = deriv_generator
                .get_derivative_mapping(&mapping_name)
                .with_scatterer();
            {
                let mut mapping = deriv_mapping.mut_view();
                let d_legendre = mapping
                    .d_legendre
                    .as_mut()
                    .ok_or_else(|| anyhow!("Scatterer derivative mapping missing d_legendre"))?;
                let scat_factor = mapping
                    .scat_factor
                    .as_mut()
                    .ok_or_else(|| anyhow!("Scatterer derivative mapping missing scat_factor"))?;
                let val_legendre = val.legendre.as_ref().ok_or_else(|| {
                    anyhow!("Scattering derivative did not provide legendre coefficients")
                })?;

                mapping.d_extinction += &val.cross_section;
                Zip::from(&mut mapping.d_ssa)
                    .and(&val.ssa)
                    .and(&val.cross_section)
                    .and(&species_ssa)
                    .and(&optical_quants.cross_section)
                    .for_each(|d_ssa, d_scat, d_ext, species_ssa_val, species_ext| {
                        *d_ssa += (d_scat - d_ext * species_ssa_val) / species_ext;
                    });

                for geo_idx in 0..val_legendre.dim().0 {
                    for wavelength_idx in 0..val_legendre.dim().1 {
                        for leg_idx in 0..val_legendre.dim().2 {
                            d_legendre[[leg_idx, geo_idx, wavelength_idx]] +=
                                val_legendre[[geo_idx, wavelength_idx, leg_idx]];
                        }
                    }
                }

                if let Some(d_vert) = self.d_vertical_deriv_factor.get(key) {
                    let interp_vertical = interp_matrix.dot(&self.vertical_deriv_factor);
                    let interp_d_vertical = interp_matrix.dot(d_vert);

                    Zip::from(&mut mapping.d_extinction)
                        .and(&optical_quants.cross_section)
                        .and_broadcast(&interp_vertical.insert_axis(Axis(1)))
                        .and_broadcast(&interp_d_vertical.insert_axis(Axis(1)))
                        .for_each(|d_ext, species_ext, vert, d_vert| {
                            *d_ext += species_ext / vert * d_vert;
                        });
                }

                Zip::indexed(&mut *d_legendre).for_each(
                    |(leg_idx, geo_idx, wavelength_idx), d_leg| {
                        let ssa_deriv = mapping.d_ssa[[geo_idx, wavelength_idx]];
                        let ext_deriv = mapping.d_extinction[[geo_idx, wavelength_idx]];
                        let species_ssa_val = species_ssa[[geo_idx, wavelength_idx]];
                        let species_ext = optical_quants.cross_section[[geo_idx, wavelength_idx]];

                        *d_leg += (optical_legendre[[geo_idx, wavelength_idx, leg_idx]]
                            - outputs.legendre[[leg_idx, geo_idx, wavelength_idx]])
                            * (ssa_deriv / species_ssa_val + ext_deriv / species_ext);
                    },
                );

                Zip::from(&mut mapping.d_ssa)
                    .and(&mapping.d_extinction)
                    .and(&species_ssa)
                    .and(&optical_quants.cross_section)
                    .and(&outputs.ssa)
                    .and(&outputs.total_extinction)
                    .for_each(
                        |d_ssa, d_ext, species_ssa_val, species_ext, total_ssa, total_ext| {
                            *d_ssa *= species_ext;
                            *d_ssa += d_ext * (species_ssa_val - total_ssa);
                            *d_ssa /= total_ext;
                        },
                    );

                for geo_idx in 0..d_legendre.dim().1 {
                    for wavelength_idx in 0..d_legendre.dim().2 {
                        let mut norm_factor = f64::NEG_INFINITY;
                        for leg_idx in 0..d_legendre.dim().0 {
                            norm_factor =
                                norm_factor.max(d_legendre[[leg_idx, geo_idx, wavelength_idx]]);
                        }
                        if norm_factor == 0.0 {
                            norm_factor = 1.0;
                        }

                        scat_factor[[geo_idx, wavelength_idx]] = species_ssa
                            [[geo_idx, wavelength_idx]]
                            * optical_quants.cross_section[[geo_idx, wavelength_idx]]
                            / (outputs.ssa[[geo_idx, wavelength_idx]]
                                * outputs.total_extinction[[geo_idx, wavelength_idx]]);

                        for leg_idx in 0..d_legendre.dim().0 {
                            d_legendre[[leg_idx, geo_idx, wavelength_idx]] /= norm_factor;
                        }
                        scat_factor[[geo_idx, wavelength_idx]] *= norm_factor;
                    }
                }
            }

            let mut optical_interpolator = interp_matrix.clone();
            Zip::from(optical_interpolator.columns_mut())
                .and(&self.number_density)
                .for_each(|mut col, number_density| {
                    col *= *number_density;
                });
            deriv_mapping.set_interpolator(&optical_interpolator);
            deriv_mapping.set_interp_dim(&format!("{constituent_name}_altitude"));
        }

        Ok(())
    }
}
