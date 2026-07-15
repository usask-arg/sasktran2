use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::interpolation::linear::linear_interpolating_matrix;
use crate::optical::traits::*;
use crate::prelude::*;

use anyhow::{Result, anyhow};
use ndarray::{CowArray, Ix1};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

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
    native_grid: bool,
}

enum ScattererAuxInputs<'a> {
    Native(&'a HashMap<String, Array1<f64>>),
    Interpolated(HashMap<String, Array1<f64>>),
}

impl AuxOpticalInputs for ScattererAuxInputs<'_> {
    fn get_parameter(&self, name: &str) -> Option<CowArray<'_, f64, Ix1>> {
        let inputs = match self {
            ScattererAuxInputs::Native(inputs) => inputs,
            ScattererAuxInputs::Interpolated(inputs) => inputs,
        };
        inputs
            .get(name)
            .map(|parameter| CowArray::from(parameter.view()))
    }
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
            native_grid: false,
        }
    }

    pub fn new_native(number_density: Array1<f64>) -> Self {
        let vertical_deriv_factor = Array1::ones(number_density.len());

        NumberDensityScatterer {
            number_density,
            altitudes: Array1::zeros(0),
            optical_property: None,
            aux_inputs: HashMap::new(),
            vertical_deriv_factor,
            d_vertical_deriv_factor: HashMap::new(),
            wf_name: "number_density".to_string(),
            interp_mode: crate::interpolation::OutOfBoundsMode::Zero,
            native_grid: true,
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

    fn interpolation_matrix<S>(&self, target_altitudes: &ArrayBase<S, Ix1>) -> Option<Array2<f64>>
    where
        S: Data<Elem = f64>,
    {
        if self.native_grid {
            None
        } else {
            Some(linear_interpolating_matrix(
                &self.altitudes,
                target_altitudes,
                self.interp_mode,
            ))
        }
    }

    fn validate_spatial_values(
        &self,
        values: &Array1<f64>,
        num_locations: usize,
        name: &str,
    ) -> Result<()> {
        let expected_len = if self.native_grid {
            num_locations
        } else {
            self.altitudes.len()
        };
        anyhow::ensure!(
            values.len() == expected_len,
            "{name} spatial dimension ({}) does not match expected size ({expected_len})",
            values.len()
        );
        Ok(())
    }

    fn spatial_values<'a>(
        &self,
        values: &'a Array1<f64>,
        interp_matrix: Option<&Array2<f64>>,
        num_locations: usize,
        name: &str,
    ) -> Result<CowArray<'a, f64, Ix1>> {
        self.validate_spatial_values(values, num_locations, name)?;
        if let Some(interp_matrix) = interp_matrix {
            Ok(CowArray::from(interp_matrix.dot(values)))
        } else {
            Ok(CowArray::from(values.view()))
        }
    }

    fn spatial_aux_inputs(
        &self,
        interp_matrix: Option<&Array2<f64>>,
        num_locations: usize,
    ) -> Result<ScattererAuxInputs<'_>> {
        if let Some(interp_matrix) = interp_matrix {
            Ok(ScattererAuxInputs::Interpolated(
                self.aux_inputs
                    .iter()
                    .map(|(key, values)| {
                        self.validate_spatial_values(values, num_locations, key)?;
                        Ok((key.clone(), interp_matrix.dot(values)))
                    })
                    .collect::<Result<_>>()?,
            ))
        } else {
            for (key, values) in &self.aux_inputs {
                self.validate_spatial_values(values, num_locations, key)?;
            }
            Ok(ScattererAuxInputs::Native(&self.aux_inputs))
        }
    }
}

fn scale_native_scatterer_derivative(
    mapping: &mut DerivMappingView<'_>,
    scale: &Array1<f64>,
) -> Result<()> {
    anyhow::ensure!(
        mapping.d_extinction.nrows() == scale.len(),
        "Native derivative scale length ({}) does not match atmosphere ({})",
        scale.len(),
        mapping.d_extinction.nrows()
    );
    let scat_factor = mapping
        .scat_factor
        .as_mut()
        .ok_or_else(|| anyhow!("Scatterer derivative mapping missing scat_factor"))?;

    Zip::from(mapping.d_extinction.rows_mut())
        .and(mapping.d_ssa.rows_mut())
        .and(scat_factor.rows_mut())
        .and(scale)
        .for_each(|mut d_extinction, mut d_ssa, mut scat_factor, scale| {
            d_extinction *= *scale;
            d_ssa *= *scale;
            scat_factor *= *scale;
        });

    Ok(())
}

impl<T> Constituent for NumberDensityScatterer<T>
where
    T: OpticalProperty,
{
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (inputs, outputs) = storage.split_inputs_outputs();

        let num_locations = inputs.altitude_m().len();
        let interp_matrix = self.interpolation_matrix(&inputs.altitude_m());
        let interp_numden = self.spatial_values(
            &self.number_density,
            interp_matrix.as_ref(),
            num_locations,
            "Number density",
        )?;
        let aux_inputs = self.spatial_aux_inputs(interp_matrix.as_ref(), num_locations)?;

        let optical_prop = self
            .optical_property
            .as_ref()
            .ok_or_else(|| anyhow!("Optical property not set"))?;

        let optical_quants = optical_prop.optical_quantities(inputs, &aux_inputs)?;
        let legendre = optical_quants.legendre.as_ref().ok_or_else(|| {
            anyhow!("Scattering optical property did not provide legendre coefficients")
        })?;

        let mut outputs = outputs.mut_view();

        let thread_pool = crate::threading::thread_pool()?;

        thread_pool.install(|| {
            Zip::from(outputs.total_extinction.axis_iter_mut(Axis(1)))
                .and(outputs.ssa.axis_iter_mut(Axis(1)))
                .and(outputs.legendre.axis_iter_mut(Axis(2)))
                .and(optical_quants.cross_section.axis_iter(Axis(1)))
                .and(optical_quants.ssa.axis_iter(Axis(1)))
                .and(legendre.axis_iter(Axis(1)))
                .par_for_each(
                    |mut output_extinction_col,
                     mut output_ssa_col,
                     mut output_legendre_wavelength,
                     species_ext_col,
                     species_scat_col,
                     optical_legendre_wavelength| {
                        Zip::from(&mut output_extinction_col)
                            .and(&mut output_ssa_col)
                            .and(&species_ext_col)
                            .and(&species_scat_col)
                            .and(&interp_numden)
                            .for_each(
                                |output_extinction, output_ssa, species_ext, species_scat, nd| {
                                    *output_extinction += species_ext * nd;
                                    *output_ssa += species_scat * nd;
                                },
                            );

                        for (geo_idx, mut output_legendre_col) in output_legendre_wavelength
                            .columns_mut()
                            .into_iter()
                            .enumerate()
                        {
                            let optical_legendre_row = optical_legendre_wavelength.row(geo_idx);
                            let species_scat = species_scat_col[geo_idx];
                            let nd = interp_numden[geo_idx];

                            Zip::from(&mut output_legendre_col)
                                .and(&optical_legendre_row)
                                .for_each(|output_legendre, optical_legendre| {
                                    *output_legendre += species_scat * nd * optical_legendre;
                                });
                        }
                    },
                );
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

        let num_locations = inputs.altitude_m().len();
        let interp_matrix = self.interpolation_matrix(&inputs.altitude_m());
        self.validate_spatial_values(&self.number_density, num_locations, "Number density")?;
        let aux_inputs = self.spatial_aux_inputs(interp_matrix.as_ref(), num_locations)?;

        let optical_prop = self
            .optical_property
            .as_ref()
            .ok_or_else(|| anyhow!("Optical property not set"))?;

        let optical_quants = optical_prop.optical_quantities(inputs, &aux_inputs)?;
        let optical_legendre = optical_quants.legendre.as_ref().ok_or_else(|| {
            anyhow!("Scattering optical property did not provide legendre coefficients")
        })?;
        let thread_pool = crate::threading::thread_pool()?;

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

            thread_pool.install(|| {
                Zip::from(mapping.d_extinction.axis_iter_mut(Axis(1)))
                    .and(scat_factor.axis_iter_mut(Axis(1)))
                    .and(optical_quants.cross_section.axis_iter(Axis(1)))
                    .and(species_ssa.axis_iter(Axis(1)))
                    .and(outputs.ssa.axis_iter(Axis(1)))
                    .and(outputs.total_extinction.axis_iter(Axis(1)))
                    .par_for_each(
                        |mut d_extinction_col,
                         mut scat_factor_col,
                         species_ext_col,
                         species_ssa_col,
                         total_ssa_col,
                         total_ext_col| {
                            Zip::from(&mut d_extinction_col)
                                .and(&mut scat_factor_col)
                                .and(&species_ext_col)
                                .and(&species_ssa_col)
                                .and(&total_ssa_col)
                                .and(&total_ext_col)
                                .for_each(
                                    |d_extinction,
                                     scat_factor,
                                     species_ext,
                                     species_ssa_val,
                                     total_ssa,
                                     total_ext| {
                                        *d_extinction += *species_ext;
                                        *scat_factor +=
                                            species_ssa_val * species_ext / (total_ssa * total_ext);
                                    },
                                );
                        },
                    );

                Zip::from(mapping.d_ssa.axis_iter_mut(Axis(1)))
                    .and(optical_quants.cross_section.axis_iter(Axis(1)))
                    .and(species_ssa.axis_iter(Axis(1)))
                    .and(outputs.ssa.axis_iter(Axis(1)))
                    .and(outputs.total_extinction.axis_iter(Axis(1)))
                    .par_for_each(
                        |mut d_ssa_col,
                         species_ext_col,
                         species_ssa_col,
                         total_ssa_col,
                         total_ext_col| {
                            Zip::from(&mut d_ssa_col)
                                .and(&species_ext_col)
                                .and(&species_ssa_col)
                                .and(&total_ssa_col)
                                .and(&total_ext_col)
                                .for_each(
                                    |d_ssa, species_ext, species_ssa_val, total_ssa, total_ext| {
                                        *d_ssa +=
                                            species_ext * (species_ssa_val - total_ssa) / total_ext;
                                    },
                                );
                        },
                    );

                Zip::from(d_legendre.axis_iter_mut(Axis(2)))
                    .and(optical_legendre.axis_iter(Axis(1)))
                    .and(outputs.legendre.axis_iter(Axis(2)))
                    .par_for_each(
                        |mut d_legendre_wavelength,
                         optical_legendre_wavelength,
                         output_legendre_wavelength| {
                            for (geo_idx, mut d_legendre_col) in
                                d_legendre_wavelength.columns_mut().into_iter().enumerate()
                            {
                                let optical_legendre_row = optical_legendre_wavelength.row(geo_idx);
                                let output_legendre_col =
                                    output_legendre_wavelength.column(geo_idx);

                                Zip::from(&mut d_legendre_col)
                                    .and(&optical_legendre_row)
                                    .and(&output_legendre_col)
                                    .for_each(|d_legendre, optical_legendre, output_legendre| {
                                        *d_legendre += optical_legendre - output_legendre;
                                    });
                            }
                        },
                    );
            });
        }

        if let Some(interp_matrix) = interp_matrix.as_ref() {
            let mut density_interpolator = interp_matrix.clone();
            Zip::from(density_interpolator.columns_mut())
                .and(&self.vertical_deriv_factor)
                .for_each(|mut col, factor| {
                    col *= *factor;
                });
            deriv_mapping.set_interpolator(&density_interpolator);
            deriv_mapping.set_interp_dim(&format!("{constituent_name}_altitude"));
        } else {
            let mut mapping = deriv_mapping.mut_view();
            scale_native_scatterer_derivative(&mut mapping, &self.vertical_deriv_factor)?;
            deriv_mapping.set_interp_dim("location");
        }

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

                thread_pool.install(|| {
                    Zip::from(mapping.d_extinction.axis_iter_mut(Axis(1)))
                        .and(val.cross_section.axis_iter(Axis(1)))
                        .par_for_each(|mut d_extinction_col, val_cross_section_col| {
                            d_extinction_col += &val_cross_section_col;
                        });

                    Zip::from(mapping.d_ssa.axis_iter_mut(Axis(1)))
                        .and(val.ssa.axis_iter(Axis(1)))
                        .and(val.cross_section.axis_iter(Axis(1)))
                        .and(species_ssa.axis_iter(Axis(1)))
                        .and(optical_quants.cross_section.axis_iter(Axis(1)))
                        .par_for_each(
                            |mut d_ssa_col,
                             val_ssa_col,
                             val_cross_section_col,
                             species_ssa_col,
                             species_ext_col| {
                                Zip::from(&mut d_ssa_col)
                                    .and(&val_ssa_col)
                                    .and(&val_cross_section_col)
                                    .and(&species_ssa_col)
                                    .and(&species_ext_col)
                                    .for_each(
                                        |d_ssa, d_scat, d_ext, species_ssa_val, species_ext| {
                                            *d_ssa +=
                                                (d_scat - d_ext * species_ssa_val) / species_ext;
                                        },
                                    );
                            },
                        );

                    Zip::from(d_legendre.axis_iter_mut(Axis(2)))
                        .and(val_legendre.axis_iter(Axis(1)))
                        .par_for_each(|mut d_legendre_wavelength, val_legendre_wavelength| {
                            for (geo_idx, mut d_legendre_col) in
                                d_legendre_wavelength.columns_mut().into_iter().enumerate()
                            {
                                let val_legendre_row = val_legendre_wavelength.row(geo_idx);

                                Zip::from(&mut d_legendre_col)
                                    .and(&val_legendre_row)
                                    .for_each(|d_legendre, val_legendre| {
                                        *d_legendre += val_legendre;
                                    });
                            }
                        });
                });

                if let Some(d_vert) = self.d_vertical_deriv_factor.get(key) {
                    let interp_vertical = self.spatial_values(
                        &self.vertical_deriv_factor,
                        interp_matrix.as_ref(),
                        num_locations,
                        "Vertical derivative factor",
                    )?;
                    let interp_d_vertical = self.spatial_values(
                        d_vert,
                        interp_matrix.as_ref(),
                        num_locations,
                        &format!("Derivative factor for {key}"),
                    )?;

                    thread_pool.install(|| {
                        Zip::from(mapping.d_extinction.axis_iter_mut(Axis(1)))
                            .and(optical_quants.cross_section.axis_iter(Axis(1)))
                            .par_for_each(|mut d_ext_col, species_ext_col| {
                                Zip::from(&mut d_ext_col)
                                    .and(&species_ext_col)
                                    .and(&interp_vertical)
                                    .and(&interp_d_vertical)
                                    .for_each(|d_ext, species_ext, vert, d_vert| {
                                        *d_ext += species_ext / vert * d_vert;
                                    });
                            });
                    });
                }

                thread_pool.install(|| {
                    d_legendre
                        .axis_iter_mut(Axis(2))
                        .into_par_iter()
                        .enumerate()
                        .for_each(|(wavelength_idx, mut d_legendre_wavelength)| {
                            let optical_legendre_wavelength =
                                optical_legendre.index_axis(Axis(1), wavelength_idx);
                            let output_legendre_wavelength =
                                outputs.legendre.index_axis(Axis(2), wavelength_idx);

                            for (geo_idx, mut d_legendre_col) in
                                d_legendre_wavelength.columns_mut().into_iter().enumerate()
                            {
                                let optical_legendre_row = optical_legendre_wavelength.row(geo_idx);
                                let output_legendre_col =
                                    output_legendre_wavelength.column(geo_idx);
                                let species_ssa_val = species_ssa[[geo_idx, wavelength_idx]];
                                let species_ext =
                                    optical_quants.cross_section[[geo_idx, wavelength_idx]];
                                let factor = mapping.d_ssa[[geo_idx, wavelength_idx]]
                                    / species_ssa_val
                                    + mapping.d_extinction[[geo_idx, wavelength_idx]] / species_ext;

                                Zip::from(&mut d_legendre_col)
                                    .and(&optical_legendre_row)
                                    .and(&output_legendre_col)
                                    .for_each(|d_legendre, optical_legendre, output_legendre| {
                                        *d_legendre +=
                                            (optical_legendre - output_legendre) * factor;
                                    });
                            }
                        });

                    Zip::from(mapping.d_ssa.axis_iter_mut(Axis(1)))
                        .and(mapping.d_extinction.axis_iter(Axis(1)))
                        .and(species_ssa.axis_iter(Axis(1)))
                        .and(optical_quants.cross_section.axis_iter(Axis(1)))
                        .and(outputs.ssa.axis_iter(Axis(1)))
                        .and(outputs.total_extinction.axis_iter(Axis(1)))
                        .par_for_each(
                            |mut d_ssa_col,
                             d_ext_col,
                             species_ssa_col,
                             species_ext_col,
                             total_ssa_col,
                             total_ext_col| {
                                Zip::from(&mut d_ssa_col)
                                    .and(&d_ext_col)
                                    .and(&species_ssa_col)
                                    .and(&species_ext_col)
                                    .and(&total_ssa_col)
                                    .and(&total_ext_col)
                                    .for_each(
                                        |d_ssa,
                                         d_ext,
                                         species_ssa_val,
                                         species_ext,
                                         total_ssa,
                                         total_ext| {
                                            *d_ssa *= species_ext;
                                            *d_ssa += d_ext * (species_ssa_val - total_ssa);
                                            *d_ssa /= total_ext;
                                        },
                                    );
                            },
                        );

                    Zip::from(d_legendre.axis_iter_mut(Axis(2)))
                        .and(scat_factor.axis_iter_mut(Axis(1)))
                        .and(species_ssa.axis_iter(Axis(1)))
                        .and(optical_quants.cross_section.axis_iter(Axis(1)))
                        .and(outputs.ssa.axis_iter(Axis(1)))
                        .and(outputs.total_extinction.axis_iter(Axis(1)))
                        .par_for_each(
                            |mut d_legendre_wavelength,
                             mut scat_factor_col,
                             species_ssa_col,
                             species_ext_col,
                             total_ssa_col,
                             total_ext_col| {
                                for (geo_idx, mut d_legendre_col) in
                                    d_legendre_wavelength.columns_mut().into_iter().enumerate()
                                {
                                    let mut norm_factor = f64::NEG_INFINITY;
                                    for val in d_legendre_col.iter() {
                                        norm_factor = norm_factor.max(*val);
                                    }
                                    if norm_factor == 0.0 {
                                        norm_factor = 1.0;
                                    }

                                    scat_factor_col[geo_idx] = species_ssa_col[geo_idx]
                                        * species_ext_col[geo_idx]
                                        / (total_ssa_col[geo_idx] * total_ext_col[geo_idx])
                                        * norm_factor;

                                    d_legendre_col.mapv_inplace(|val| val / norm_factor);
                                }
                            },
                        );
                });
            }

            if let Some(interp_matrix) = interp_matrix.as_ref() {
                let mut optical_interpolator = interp_matrix.clone();
                Zip::from(optical_interpolator.columns_mut())
                    .and(&self.number_density)
                    .for_each(|mut col, number_density| {
                        col *= *number_density;
                    });
                deriv_mapping.set_interpolator(&optical_interpolator);
                deriv_mapping.set_interp_dim(&format!("{constituent_name}_altitude"));
            } else {
                let mut mapping = deriv_mapping.mut_view();
                scale_native_scatterer_derivative(&mut mapping, &self.number_density)?;
                deriv_mapping.set_interp_dim("location");
            }
        }

        Ok(())
    }
}
