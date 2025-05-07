use crate::interpolation::grid1d::Grid1DView;
use crate::prelude::*;

use crate::math::errorfunctions::optimized::w_jpole_real_assign;
use crate::optical::line::{AdjustedLineParameters, OpticalLine, OpticalLineDB};
use crate::optical::traits::*;

const SQRT_PI: f64 = 1.7724538509055160272981674833411;

#[inline(always)]
fn lorentzian_assign(
    wvnum: &[f64],
    line_center: f64,
    doppler_width: f64,
    y: f64,
    scale: f64,
    xs: &mut [f64],
) {
    let scale = scale * y / SQRT_PI;
    for (wv, xs) in wvnum.iter().zip(xs.iter_mut()) {
        let x = (wv - line_center) / doppler_width;
        *xs += scale / (x * x + y * y);
    }
}

#[inline(always)]
fn gaussian_assign(
    wvnum: &[f64],
    line_center: f64,
    doppler_width: f64,
    _y: f64,
    scale: f64,
    xs: &mut [f64],
) {
    for (wv, xs) in wvnum.iter().zip(xs.iter_mut()) {
        let x = (wv - line_center) / doppler_width;
        *xs += scale * (-x * x).exp();
    }
}

#[inline(always)]
fn split_and_assign(
    adjusted_line: &AdjustedLineParameters,
    wavenumber_cminv: &Grid1DView,
    xs: &mut [f64],
) {
    if 2.84 * adjusted_line.y * adjusted_line.y > 1.52 / EPSILON {
        lorentzian_assign(
            &wavenumber_cminv.x,
            adjusted_line.line_center,
            adjusted_line.doppler_width,
            adjusted_line.y,
            adjusted_line.line_intensity_re,
            xs,
        );
        return;
    }

    let max_x = wavenumber_cminv.x[wavenumber_cminv.x.len() - 1];
    let min_x = wavenumber_cminv.x[0];

    let max_abs_x = max_x.abs().max(min_x.abs());

    if max_abs_x < 2.15 - 2.53 * adjusted_line.y / EPSILON {
        gaussian_assign(
            &wavenumber_cminv.x,
            adjusted_line.line_center,
            adjusted_line.doppler_width,
            adjusted_line.y,
            adjusted_line.line_intensity_re,
            xs,
        );
        return;
    }

    const EPSILON: f64 = 1.0e-4;

    let split_x = (1.52 / EPSILON - 2.84 * adjusted_line.y * adjusted_line.y).sqrt();

    let left = wavenumber_cminv.lower_bound(-split_x);
    let right = wavenumber_cminv.lower_bound(split_x);

    lorentzian_assign(
        &wavenumber_cminv.x[0..left],
        adjusted_line.line_center,
        adjusted_line.doppler_width,
        adjusted_line.y,
        adjusted_line.line_intensity_re,
        &mut xs[0..left],
    );

    lorentzian_assign(
        &wavenumber_cminv.x[right..],
        adjusted_line.line_center,
        adjusted_line.doppler_width,
        adjusted_line.y,
        adjusted_line.line_intensity_re,
        &mut xs[right..],
    );

    w_jpole_real_assign(
        &wavenumber_cminv.x[left..right],
        adjusted_line.line_center,
        adjusted_line.doppler_width,
        adjusted_line.y,
        adjusted_line.line_intensity_re,
        &mut xs[left..right],
    );
}

pub trait PartitionFactor {
    fn partition_factor(&self, mol_id: i32, iso_id: i32, temperature: f64) -> f64;
}

pub trait MolecularMass {
    fn molecular_mass(&self, mol_id: i32, iso_id: i32) -> f64;
}

struct MolParam {
    mol_mass: f64,
    partition_factor: Vec<f64>,
}

pub struct LineAbsorber {
    db: OpticalLineDB,
    line_contribution_width: f64,
    cull_factor: f64,
    subtract_pedestal: bool,
    partition_generator: Option<Box<dyn PartitionFactor>>,
    mol_mass_generator: Option<Box<dyn MolecularMass>>,
}

impl LineAbsorber {
    pub fn new(db: OpticalLineDB) -> Self {
        Self {
            db,
            line_contribution_width: 25.0,
            cull_factor: 0.0,
            subtract_pedestal: false,
            partition_generator: None,
            mol_mass_generator: None,
        }
    }

    pub fn with_line_contribution_width(mut self, line_contribution_width: f64) -> Self {
        self.line_contribution_width = line_contribution_width;
        self
    }

    pub fn with_cull_factor(mut self, cull_factor: f64) -> Self {
        self.cull_factor = cull_factor;
        self
    }

    pub fn with_subtract_pedestal(mut self, subtract_pedestal: bool) -> Self {
        self.subtract_pedestal = subtract_pedestal;
        self
    }

    pub fn with_partition_generator(
        mut self,
        partition_generator: Box<dyn PartitionFactor>,
    ) -> Self {
        self.partition_generator = Some(partition_generator);
        self
    }

    pub fn with_molecular_mass_generator(
        mut self,
        mol_mass_generator: Box<dyn MolecularMass>,
    ) -> Self {
        self.mol_mass_generator = Some(mol_mass_generator);
        self
    }

    fn mol_params(&self, mol_id: i32, iso_id: i32, temperature: &[f64]) -> Result<MolParam> {
        if self.partition_generator.is_none() {
            return Err(anyhow!("Partition generator not set"));
        }

        if self.mol_mass_generator.is_none() {
            return Err(anyhow!("Molecular mass generator not set"));
        }

        let base_partition = self
            .partition_generator
            .as_ref()
            .unwrap()
            .partition_factor(mol_id, iso_id, 296.0);

        let mut partition_factor = vec![];
        for t in temperature.iter() {
            partition_factor.push(
                self.partition_generator
                    .as_ref()
                    .unwrap()
                    .partition_factor(mol_id, iso_id, *t)
                    / base_partition,
            );
        }

        let mol_mass = self
            .mol_mass_generator
            .as_ref()
            .unwrap()
            .molecular_mass(mol_id, iso_id);

        Ok(MolParam {
            mol_mass,
            partition_factor,
        })
    }

    fn gen_mol_param(
        &self,
        lines: &[OpticalLine],
        temperature: &[f64],
    ) -> Result<HashMap<(i32, i32), MolParam>> {
        let mut map = HashMap::new();
        for line in lines.iter() {
            let key = (line.mol_id, line.iso_id);
            if !map.contains_key(&key) {
                let mol_param = self.mol_params(line.mol_id, line.iso_id, temperature)?;
                map.insert(key, mol_param);
            }
        }
        Ok(map)
    }

    pub fn cross_section(
        &self,
        wavenumber_cminv: ArrayView1<f64>,
        temperature: ArrayView1<f64>,
        pressure: ArrayView1<f64>,
        pself: ArrayView1<f64>,
    ) -> Result<Array2<f64>> {
        const K_B: f64 = 1.38064852e-16;
        let max_p_self = pself.iter().copied().fold(0.0, f64::max);
        let max_p_self = match max_p_self {
            0.0 => 1.0,
            _ => max_p_self,
        };

        let wavenumber_grid = Grid1DView::new(wavenumber_cminv.as_slice().unwrap());

        let mut xs = Array2::zeros((temperature.len(), wavenumber_cminv.len()));

        let min_wvnum = wavenumber_cminv[0];
        let max_wvnum = wavenumber_cminv[wavenumber_cminv.len() - 1];

        let n_wavenumber = wavenumber_cminv.len();

        let mut start_wavenumber_idx: usize = 0;
        let mut end_wavenumber_idx: usize = 0;

        let line_slice = self.db.between_slice(min_wvnum, max_wvnum);

        let map = self.gen_mol_param(line_slice, temperature.as_slice().unwrap())?;

        for line in line_slice.iter() {
            if line.line_intensity * 101325.0 * max_p_self / (K_B * 1e-7 * 296.0) < self.cull_factor
            {
                continue;
            }

            while start_wavenumber_idx < n_wavenumber
                && wavenumber_cminv[start_wavenumber_idx]
                    < line.line_center - self.line_contribution_width
            {
                start_wavenumber_idx += 1;
            }

            while end_wavenumber_idx < n_wavenumber
                && wavenumber_cminv[end_wavenumber_idx]
                    < line.line_center + self.line_contribution_width
            {
                end_wavenumber_idx += 1;
            }

            if start_wavenumber_idx == end_wavenumber_idx {
                continue;
            }

            let sub_grid = wavenumber_grid.slice(start_wavenumber_idx, end_wavenumber_idx);

            let mol_param = map.get(&(line.mol_id, line.iso_id)).unwrap();

            Zip::indexed(temperature)
                .and(pressure)
                .and(pself)
                .and(xs.axis_iter_mut(Axis(0)))
                .par_for_each(|g, &temperature, &pressure, &pself, mut xs| {
                    let adjusted_line = line
                        .adjusted_parameters(
                            temperature,
                            pressure,
                            pself,
                            mol_param.partition_factor[g],
                            mol_param.mol_mass,
                        )
                        .unwrap();

                    split_and_assign(
                        &adjusted_line,
                        &sub_grid,
                        xs.slice_mut(ndarray::s![start_wavenumber_idx..end_wavenumber_idx])
                            .as_slice_mut()
                            .unwrap(),
                    );
                });
        }

        Ok(xs)
    }
}

impl OpticalProperty for LineAbsorber {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn crate::atmosphere::StorageInputs,
        _aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut crate::optical::storage::OpticalQuantities,
    ) -> Result<()> {
        let wavenumber_cminv = inputs
            .wavenumbers_cminv()
            .ok_or(anyhow::anyhow!("Wavenumbers not found in inputs"))?;

        let temperature = inputs
            .temperature_k()
            .ok_or(anyhow::anyhow!("Temperature not found in inputs"))?;

        let pressure = inputs
            .pressure_pa()
            .ok_or(anyhow::anyhow!("Pressure not found in inputs"))?;

        let pself = Array1::zeros(temperature.len());

        optical_quantities.cross_section =
            self.cross_section(wavenumber_cminv, temperature, pressure, pself.view())?;

        Ok(())
    }

    fn optical_derivatives_emplace(
        &self,
        _inputs: &dyn crate::atmosphere::StorageInputs,
        _aux_inputs: &dyn AuxOpticalInputs,
        _d_optical_quantities: &mut HashMap<String, crate::optical::storage::OpticalQuantities>,
    ) -> Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::optical::line::aer_loader::read_aer_line_file;
    use crate::optical::types::line_absorber::{MolecularMass, PartitionFactor};
    use ndarray::{Array1, Array2, array};
    use std::{collections::HashMap, path::PathBuf};

    struct MockPartitionFactor {
        pub partition_factor: f64,
    }
    impl PartitionFactor for MockPartitionFactor {
        fn partition_factor(&self, _mol_id: i32, _iso_id: i32, _temperature: f64) -> f64 {
            self.partition_factor
        }
    }

    struct MockMolecularMass {
        pub mol_mass: f64,
    }
    impl MolecularMass for MockMolecularMass {
        fn molecular_mass(&self, _mol_id: i32, _iso_id: i32) -> f64 {
            self.mol_mass
        }
    }

    #[test]
    fn test_cross_section() {
        let o2_file = PathBuf::from("../../tests/data/02_CO2");

        let line_absorber = LineAbsorber::new(read_aer_line_file(o2_file).unwrap());

        let partition_factor = MockPartitionFactor {
            partition_factor: 1.0,
        };
        let mol_mass = MockMolecularMass { mol_mass: 44.01 };
        let line_absorber = line_absorber
            .with_partition_generator(Box::new(partition_factor))
            .with_molecular_mass_generator(Box::new(mol_mass));

        let wavenumber_cminv = array![200.0, 1500.0, 2000.0];
        let temperature = array![296.0, 300.0];
        let pressure = array![101325.0, 101325.0];
        let pself = array![0.0, 0.0];

        let xs = line_absorber
            .cross_section(
                wavenumber_cminv.view(),
                temperature.view(),
                pressure.view(),
                pself.view(),
            )
            .unwrap();

        assert_eq!(xs.shape(), &[2, 3]);
    }
}
