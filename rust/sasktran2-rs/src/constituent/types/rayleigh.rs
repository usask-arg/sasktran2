use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::optical::rayleigh::rayleigh_cross_section_bates;
use crate::prelude::*;

use crate::interpolation::linear::Interp1;

pub enum RayleighMethod {
    Bates,
    Manual,
}

pub struct Rayleigh {
    n2_percentage: f64,
    o2_percentage: f64,
    ar_percentage: f64,
    co2_percentage: f64,
    method: RayleighMethod,
    manual_xs: Option<Array1<f64>>,
    manual_king: Option<Array1<f64>>,
    manual_wavelengths_nm: Option<Array1<f64>>,
}

impl Default for Rayleigh {
    fn default() -> Self {
        Rayleigh::new()
    }
}

impl Rayleigh {
    pub fn new() -> Self {
        Rayleigh {
            method: RayleighMethod::Bates,
            n2_percentage: 78.084,
            o2_percentage: 20.946,
            ar_percentage: 0.934,
            co2_percentage: 0.036,
            manual_xs: None,
            manual_king: None,
            manual_wavelengths_nm: None,
        }
    }

    pub fn with_n2_percentage(mut self, n2_percentage: f64) -> Self {
        self.n2_percentage = n2_percentage;
        self
    }

    pub fn with_o2_percentage(mut self, o2_percentage: f64) -> Self {
        self.o2_percentage = o2_percentage;
        self
    }

    pub fn with_ar_percentage(mut self, ar_percentage: f64) -> Self {
        self.ar_percentage = ar_percentage;
        self
    }

    pub fn with_co2_percentage(mut self, co2_percentage: f64) -> Self {
        self.co2_percentage = co2_percentage;
        self
    }

    pub fn with_manual_xs(
        mut self,
        xs: Array1<f64>,
        king: Array1<f64>,
        wavelengths_nm: Array1<f64>,
    ) -> Self {
        self.method = RayleighMethod::Manual;
        self.manual_xs = Some(xs);
        self.manual_king = Some(king);
        self.manual_wavelengths_nm = Some(wavelengths_nm);
        self
    }

    fn cross_section(&self, wavelength_nm: f64) -> (f64, f64) {
        match self.method {
            RayleighMethod::Bates => {
                let wavelength_um = wavelength_nm / 1000.0;
                rayleigh_cross_section_bates(
                    wavelength_um,
                    self.n2_percentage,
                    self.o2_percentage,
                    self.ar_percentage,
                    self.co2_percentage,
                )
            }
            RayleighMethod::Manual => {
                let xs = self.manual_xs.as_ref().unwrap().interp1(
                    self.manual_wavelengths_nm.as_ref().unwrap(),
                    wavelength_nm,
                    crate::interpolation::OutOfBoundsMode::Extend,
                );
                let king = self.manual_king.as_ref().unwrap().interp1(
                    self.manual_wavelengths_nm.as_ref().unwrap(),
                    wavelength_nm,
                    crate::interpolation::OutOfBoundsMode::Extend,
                );

                (xs, king)
            }
        }
    }
}

impl Constituent for Rayleigh {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (inputs, outputs) = storage.split_inputs_outputs();

        let wavelengths_nm = inputs.wavelengths_nm().unwrap();
        let mut outputs = outputs.mut_view();

        let total_extinction_array = &mut outputs.total_extinction;
        let ssa_array = &mut outputs.ssa;
        let legendre_array = &mut outputs.legendre;

        let num_dens = inputs.air_numberdensity_dict()["N"].to_owned();
        let num_stokes = inputs.num_stokes();

        let thread_pool = crate::threading::thread_pool()?;

        thread_pool.install(|| {
            Zip::from(total_extinction_array.columns_mut())
                .and(ssa_array.columns_mut())
                .and(legendre_array.axis_iter_mut(Axis(2)))
                .and(wavelengths_nm)
                .par_for_each(|extinction_col, ssa_col, mut legendre_col, wavelength_nm| {
                    let (sigma, king) = self.cross_section(*wavelength_nm);

                    let delta = 6.0 * (king - 1.0) / (3.0 + 7.0 * king);

                    Zip::from(extinction_col)
                        .and(ssa_col)
                        .and(legendre_col.columns_mut())
                        .and(&num_dens)
                        .for_each(|extinction, ssa, mut legendre, ndens| {
                            *extinction += sigma * ndens;
                            *ssa += sigma * ndens;

                            if num_stokes == 1 {
                                legendre[0] += sigma * ndens;
                                legendre[2] += sigma * ndens * (1.0 - delta) / (2.0 + delta);
                            } else if num_stokes == 3 {
                                legendre[0] += sigma * ndens;
                                legendre[8] += sigma * ndens * (1.0 - delta) / (2.0 + delta);

                                legendre[9] +=
                                    sigma * ndens * 6.0 * ((1.0 - delta) / (2.0 + delta));

                                legendre[11] += sigma
                                    * ndens
                                    * 6.0_f64.sqrt()
                                    * ((1.0 - delta) / (2.0 + delta));
                            } else {
                                panic!("Should never be here");
                            }
                        });
                });
        });

        Ok(())
    }

    fn register_derivatives<'a>(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        constituent_name: &str,
    ) -> Result<()> {
        let (inputs, outputs, derivative_generator) = storage.split_inputs_outputs_deriv();

        let outputs = outputs.view();

        let air_dens = inputs.air_numberdensity_dict();

        let wavelengths_nm = inputs.wavelengths_nm().unwrap();

        let deriv_names = ["pressure_pa", "temperature_k"];
        let deriv_vals = [&air_dens["dN_dP"], &air_dens["dN_dT"]];

        let num_stokes = inputs.num_stokes();

        let thread_pool = crate::threading::thread_pool()?;

        for (deriv_name, deriv_val) in deriv_names.iter().zip(deriv_vals.iter()) {
            let full_name = "wf_".to_owned() + constituent_name + "_" + deriv_name;
            let mut deriv = derivative_generator
                .get_derivative_mapping(&full_name)
                .with_scatterer();
            {
                let mut deriv_view = deriv.mut_view();

                thread_pool.install(|| {
                    Zip::indexed(deriv_view.d_extinction.columns_mut())
                        .and(deriv_view.d_ssa.columns_mut())
                        .and(deriv_view.d_legendre.unwrap().axis_iter_mut(Axis(2)))
                        .and(deriv_view.scat_factor.unwrap().columns_mut())
                        .and(outputs.legendre.axis_iter(Axis(2)))
                        .par_for_each(
                            |i, d_k_row, d_ssa_row, mut d_leg_row, scat_f_row, leg_row| {
                                let k_row = outputs.total_extinction.column(i);
                                let ssa_row = outputs.ssa.column(i);

                                let (sigma, king) = self.cross_section(wavelengths_nm[i]);

                                let delta = 6.0 * (king - 1.0) / (3.0 + 7.0 * king);

                                Zip::indexed(d_k_row)
                                    .and(d_ssa_row)
                                    .and(d_leg_row.columns_mut())
                                    .and(scat_f_row)
                                    .and(leg_row.columns())
                                    .for_each(|j, d_k, d_ssa, mut d_lp, scat_f, lp| {
                                        let ssa = ssa_row[j];
                                        let k = k_row[j];

                                        *d_k += sigma;
                                        *d_ssa += sigma * (1.0 - ssa) / k;

                                        d_lp[0] += 1.0;

                                        if num_stokes == 1 {
                                            d_lp[2] += (1.0 - delta) / (2.0 + delta);
                                        } else if num_stokes == 3 {
                                            d_lp[8] += (1.0 - delta) / (2.0 + delta);
                                            d_lp[9] += 6.0 * ((1.0 - delta) / (2.0 + delta));
                                            d_lp[11] +=
                                                6.0_f64.sqrt() * ((1.0 - delta) / (2.0 + delta));
                                        }

                                        Zip::from(d_lp).and(lp).for_each(|d_lp, lp| {
                                            *d_lp -= lp;
                                        });

                                        *scat_f += sigma / (ssa * k);
                                    });
                            },
                        );
                });
            }
            deriv.set_interp_dim("altitude");
            let interpolator: Array2<f64> = Array2::from_diag(deriv_val);
            deriv.set_interpolator(&interpolator);
            let assign_name = "wf_".to_owned() + deriv_name;
            deriv.set_assign_name(&assign_name);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {}
