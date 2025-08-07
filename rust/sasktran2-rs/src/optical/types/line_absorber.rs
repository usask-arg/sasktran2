use crate::interpolation::OutOfBoundsMode;
use crate::interpolation::grid1d::Grid1DView;
use crate::interpolation::linear::Interp1;
use crate::prelude::*;

use crate::math::errorfunctions::optimized::{w_jpole_assign, w_jpole_real_assign};
use crate::optical::line::{AdjustedLineParameters, OpticalLine, OpticalLineDB};
use crate::optical::traits::*;
use crate::util::argsort_f64;

#[cfg(feature = "simd")]
use crate::math::simd::*;

const SQRT_PI: f64 = 1.772_453_850_905_516;

#[inline(always)]
#[cfg(not(feature = "simd"))]
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
#[cfg(feature = "simd")]
fn lorentzian_assign(
    wvnum: &[f64],
    line_center: f64,
    doppler_width: f64,
    y: f64,
    scale: f64,
    xs: &mut [f64],
) {
    let scale_s = f64s::splat(scale * y / SQRT_PI);
    let lanes = f64s::LEN;

    let wvnum_chunks = wvnum.chunks_exact(lanes);
    let ys = f64s::splat(y);
    let ys_sqr = ys * ys;
    let remainder = wvnum_chunks.remainder();

    for (wv, xs) in wvnum_chunks.zip(xs.chunks_exact_mut(lanes)) {
        let x = (f64s::from_slice(wv) - f64s::splat(line_center)) / f64s::splat(doppler_width);
        let denom = x * x + ys_sqr;

        let result = f64s::from_slice(xs) + scale_s / denom;
        xs.copy_from_slice(&result.to_array());
    }

    let n = wvnum.len();

    for i in n - remainder.len()..n {
        let x = (wvnum[i] - line_center) / doppler_width;
        let denom = x * x + y * y;
        xs[i] += scale * y / SQRT_PI / denom;
    }
}

#[inline(always)]
#[cfg(not(feature = "simd"))]
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
#[cfg(feature = "simd")]
fn gaussian_assign(
    wvnum: &[f64],
    line_center: f64,
    doppler_width: f64,
    _y: f64,
    scale: f64,
    xs: &mut [f64],
) {
    let scale_s = f64s::splat(scale);
    let lanes = f64s::LEN;

    let wvnum_chunks = wvnum.chunks_exact(lanes);
    let remainder = wvnum_chunks.remainder();

    for (wv, xs) in wvnum_chunks.zip(xs.chunks_exact_mut(lanes)) {
        let x = (f64s::from_slice(wv) - f64s::splat(line_center)) / f64s::splat(doppler_width);
        let result = f64s::from_slice(xs) + scale_s * <f64s as std::simd::StdFloat>::exp(-x * x);
        xs.copy_from_slice(&result.to_array());
    }

    let n = wvnum.len();

    for i in n - remainder.len()..n {
        let x = (wvnum[i] - line_center) / doppler_width;
        let x = (-x * x).exp();
        xs[i] += scale * x;
    }
}

#[inline(always)]
fn split_and_assign(
    adjusted_line: &AdjustedLineParameters,
    wavenumber_cminv: &Grid1DView,
    xs: &mut [f64],
) {
    const EPSILON: f64 = 1.0e-4;
    let n = wavenumber_cminv.x.len();

    if n == 0 {
        return;
    }

    if 2.84 * adjusted_line.y * adjusted_line.y > 1.52 / EPSILON {
        lorentzian_assign(
            wavenumber_cminv.x,
            adjusted_line.line_center,
            adjusted_line.doppler_width,
            adjusted_line.y,
            adjusted_line.line_intensity_re,
            xs,
        );
        return;
    }

    let max_x = wavenumber_cminv.x[n - 1];
    let min_x = wavenumber_cminv.x[0];

    let max_abs_x = max_x.abs().max(min_x.abs());

    if max_abs_x < 2.15 - 2.53 * adjusted_line.y / EPSILON {
        gaussian_assign(
            wavenumber_cminv.x,
            adjusted_line.line_center,
            adjusted_line.doppler_width,
            adjusted_line.y,
            adjusted_line.line_intensity_re,
            xs,
        );
        return;
    }

    let split_x = (1.52 / EPSILON - 2.84 * adjusted_line.y * adjusted_line.y).sqrt();

    // Convert to wvnum,
    let split_x = split_x * adjusted_line.doppler_width;

    let left = wavenumber_cminv.lower_bound(adjusted_line.line_center - split_x);
    let right = wavenumber_cminv.lower_bound(adjusted_line.line_center + split_x);

    let left = left.min(n);
    let right = right.min(n);

    lorentzian_assign(
        &wavenumber_cminv.x[0..left],
        adjusted_line.line_center,
        adjusted_line.doppler_width,
        adjusted_line.y,
        adjusted_line.line_intensity_re,
        &mut xs[0..left],
    );

    if left == n {
        return;
    }

    w_jpole_real_assign(
        &wavenumber_cminv.x[left..right],
        adjusted_line.line_center,
        adjusted_line.doppler_width,
        adjusted_line.y,
        adjusted_line.line_intensity_re,
        &mut xs[left..right],
    );

    if right == wavenumber_cminv.x.len() {
        return;
    }

    lorentzian_assign(
        &wavenumber_cminv.x[right..n],
        adjusted_line.line_center,
        adjusted_line.doppler_width,
        adjusted_line.y,
        adjusted_line.line_intensity_re,
        &mut xs[right..n],
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
    enable_line_coupling: bool,
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
            enable_line_coupling: false,
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

    pub fn with_line_coupling(mut self, enable_line_coupling: bool) -> Self {
        self.enable_line_coupling = enable_line_coupling;
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
            if let std::collections::hash_map::Entry::Vacant(e) = map.entry(key) {
                let mol_param = self.mol_params(line.mol_id, line.iso_id, temperature)?;
                e.insert(mol_param);
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
            0.0 => 101325.0,
            _ => max_p_self,
        };

        let enabling_line_coupling = self.enable_line_coupling;

        let mut sorted_wvnum = wavenumber_cminv.as_slice().unwrap().to_owned();
        sorted_wvnum.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let wavenumber_grid = Grid1DView::new(sorted_wvnum.as_slice());

        let mut xs = Array2::zeros((temperature.len(), wavenumber_cminv.len()));

        let min_wvnum = sorted_wvnum[0];
        let max_wvnum = sorted_wvnum[sorted_wvnum.len() - 1];

        let n_wavenumber = wavenumber_cminv.len();

        let mut start_wavenumber_idx: usize = 0;
        let mut end_wavenumber_idx: usize = 0;

        let line_slice = self.db.between_slice(
            min_wvnum - self.line_contribution_width,
            max_wvnum + self.line_contribution_width,
        );

        let map = self.gen_mol_param(line_slice, temperature.as_slice().unwrap())?;

        for line in line_slice.iter() {
            if line.line_intensity * max_p_self / (K_B * 1e-7 * 296.0) < self.cull_factor {
                continue;
            }

            while start_wavenumber_idx < n_wavenumber
                && sorted_wvnum[start_wavenumber_idx]
                    < line.line_center - self.line_contribution_width
            {
                start_wavenumber_idx += 1;
            }

            while end_wavenumber_idx < n_wavenumber
                && sorted_wvnum[end_wavenumber_idx]
                    < line.line_center + self.line_contribution_width
            {
                end_wavenumber_idx += 1;
            }

            let start_wavenumber_idx = start_wavenumber_idx.min(n_wavenumber - 1);

            if start_wavenumber_idx == end_wavenumber_idx {
                continue;
            }

            let sub_grid = wavenumber_grid.slice(start_wavenumber_idx, end_wavenumber_idx);

            let mol_param = map.get(&(line.mol_id, line.iso_id)).unwrap();

            let thread_pool = crate::threading::thread_pool()?;

            thread_pool.install(|| {
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

                        if !enabling_line_coupling || line.y_coupling.is_empty() {
                            split_and_assign(
                                &adjusted_line,
                                &sub_grid,
                                xs.slice_mut(ndarray::s![start_wavenumber_idx..end_wavenumber_idx])
                                    .as_slice_mut()
                                    .unwrap(),
                            );
                        } else {
                            // Have to do line coupling
                            // Interp the coupling values
                            let temp_val =
                                Array1::<f64>::from_vec(line.coupling_temperature.clone());
                            let y_coupling = Array1::<f64>::from_vec(line.y_coupling.clone());
                            let g_coupling = Array1::<f64>::from_vec(line.g_coupling.clone());

                            let y_coupling =
                                y_coupling.interp1(&temp_val, temperature, OutOfBoundsMode::Extend);
                            let g_coupling =
                                g_coupling.interp1(&temp_val, temperature, OutOfBoundsMode::Extend);

                            let p_norm = pressure / 101325.0;

                            let scale_re = adjusted_line.line_intensity_re
                                * (1.0 + p_norm * p_norm * g_coupling);
                            let scale_im = adjusted_line.line_intensity_re * (-p_norm * y_coupling);

                            //let scale_re = adjusted_line.line_intensity_re;
                            //let scale_im = adjusted_line.line_intensity_im;

                            // Use the full faddeeva function
                            w_jpole_assign(
                                sub_grid.x,
                                adjusted_line.line_center,
                                adjusted_line.doppler_width,
                                adjusted_line.y,
                                scale_re,
                                scale_im,
                                xs.slice_mut(ndarray::s![start_wavenumber_idx..end_wavenumber_idx])
                                    .as_slice_mut()
                                    .unwrap(),
                            );
                        }
                    });
            });
        }

        if wavenumber_cminv.as_slice().unwrap().is_sorted() {
            Ok(xs)
        } else {
            // Have to sort the output
            let sort_idx = argsort_f64(wavenumber_cminv.as_slice().unwrap());

            let mut xs_sorted = Array2::zeros((temperature.len(), n_wavenumber));
            for i in 0..temperature.len() {
                for j in 0..n_wavenumber {
                    xs_sorted[[i, j]] = xs[[i, sort_idx[j]]];
                }
            }
            Ok(xs_sorted)
        }
    }
}

impl OpticalProperty for LineAbsorber {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn crate::atmosphere::StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
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

        let mut pself = Array1::zeros(temperature.len());

        if let Some(vmr) = aux_inputs.get_parameter("vmr") {
            pself.assign(&vmr.view());

            Zip::from(&mut pself)
                .and(pressure)
                .for_each(|pself, pressure| {
                    *pself *= pressure;
                });
        }

        optical_quantities.cross_section =
            self.cross_section(wavenumber_cminv, temperature, pressure, pself.view())?;

        if optical_quantities.ssa.shape() == optical_quantities.cross_section.shape() {
            optical_quantities.ssa *= 0.0;
        } else {
            optical_quantities.ssa = Array2::zeros(optical_quantities.cross_section.raw_dim());
        }

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
    use crate::optical::line::hitran_loader::read_hitran_line_file;
    use crate::optical::types::line_absorber::{MolecularMass, PartitionFactor};
    use ndarray::array;
    use std::path::PathBuf;

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

        let wavenumber_cminv = array![
            200.0, 470.0, 470.1, 470.2, 470.3, 471.9, 472.0, 472.1, 474.0, 1500.0, 2000.0
        ];
        let temperature = array![296.0, 300.0];
        let pressure = array![101325.0, 10.0];
        let pself = array![0.0, 0.0];

        let xs = line_absorber
            .cross_section(
                wavenumber_cminv.view(),
                temperature.view(),
                pressure.view(),
                pself.view(),
            )
            .unwrap();

        assert_eq!(xs.shape(), &[2, 11]);
    }

    #[test]
    fn test_cross_section_hitran() {
        let o2_file = PathBuf::from("../../tests/data/CO2.data");

        let line_absorber = LineAbsorber::new(read_hitran_line_file(o2_file).unwrap());

        let partition_factor = MockPartitionFactor {
            partition_factor: 1.0,
        };
        let mol_mass = MockMolecularMass { mol_mass: 44.01 };
        let line_absorber = line_absorber
            .with_partition_generator(Box::new(partition_factor))
            .with_molecular_mass_generator(Box::new(mol_mass));

        let wavenumber_cminv = array![480.0, 485.0, 490.0];
        let temperature = array![296.0, 300.0];
        let pressure = array![101325.0, 10.0];
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

    #[test]
    fn test_cross_section_unsorted() {
        let o2_file = PathBuf::from("../../tests/data/02_CO2");

        let line_absorber = LineAbsorber::new(read_aer_line_file(o2_file).unwrap());

        let partition_factor = MockPartitionFactor {
            partition_factor: 1.0,
        };
        let mol_mass = MockMolecularMass { mol_mass: 44.01 };
        let line_absorber = line_absorber
            .with_partition_generator(Box::new(partition_factor))
            .with_molecular_mass_generator(Box::new(mol_mass));

        let wavenumber_cminv = array![2000.0, 1500.0, 200.0];
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

        // REdo the test with sorted wavenumbers
        let wavenumber_cminv = array![200.0, 1500.0, 2000.0];
        let xs_sorted = line_absorber
            .cross_section(
                wavenumber_cminv.view(),
                temperature.view(),
                pressure.view(),
                pself.view(),
            )
            .unwrap();
        assert_eq!(xs_sorted.shape(), &[2, 3]);
        assert_eq!(xs, xs_sorted);
    }

    #[test]
    fn test_line_coupling_enabled() {
        let o2_file = PathBuf::from("../../tests/data/02_CO2");

        let line_absorber = LineAbsorber::new(read_aer_line_file(o2_file).unwrap());

        let partition_factor = MockPartitionFactor {
            partition_factor: 1.0,
        };
        let mol_mass = MockMolecularMass { mol_mass: 44.01 };
        let line_absorber = line_absorber
            .with_partition_generator(Box::new(partition_factor))
            .with_molecular_mass_generator(Box::new(mol_mass))
            .with_line_coupling(true);

        let wavenumber_cminv = array![200.0, 472.0, 1500.0, 2000.0];
        let temperature = array![296.0, 300.0];
        let pressure = array![101325.0, 101325.0];
        let pself = array![0.0, 0.0];

        let _xs = line_absorber
            .cross_section(
                wavenumber_cminv.view(),
                temperature.view(),
                pressure.view(),
                pself.view(),
            )
            .unwrap();
    }
}
