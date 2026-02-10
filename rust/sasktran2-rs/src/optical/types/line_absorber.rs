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
    let scale = scale * y / SQRT_PI * doppler_width * doppler_width;
    let gamma = y * doppler_width;
    let gamma_sq = gamma * gamma;
    for (wv, xs) in wvnum.iter().zip(xs.iter_mut()) {
        let x = wv - line_center;
        *xs += scale / (x * x + gamma_sq);
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
    let scale_s = f64s::splat(scale * y / SQRT_PI * doppler_width * doppler_width);
    let gamma = y * doppler_width;
    let ys_sqr = f64s::splat(gamma * gamma);
    let lanes = f64s::LEN;

    let wvnum_chunks = wvnum.chunks_exact(lanes);
    let remainder = wvnum_chunks.remainder();

    for (wv, xs) in wvnum_chunks.zip(xs.chunks_exact_mut(lanes)) {
        let x = f64s::from_slice(wv) - f64s::splat(line_center);
        let denom = x * x + ys_sqr;

        let result = f64s::from_slice(xs) + scale_s / denom;
        xs.copy_from_slice(&result.to_array());
    }

    let n = wvnum.len();

    for i in n - remainder.len()..n {
        let x = wvnum[i] - line_center;
        let denom = x * x + gamma * gamma;
        xs[i] += scale * y / SQRT_PI * doppler_width * doppler_width / denom;
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
    let inv_doppler = 1.0 / doppler_width;
    for (wv, xs) in wvnum.iter().zip(xs.iter_mut()) {
        let x = (wv - line_center) * inv_doppler;
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
    let inv_doppler = f64s::splat(1.0 / doppler_width);
    let lanes = f64s::LEN;

    let wvnum_chunks = wvnum.chunks_exact(lanes);
    let remainder = wvnum_chunks.remainder();

    for (wv, xs) in wvnum_chunks.zip(xs.chunks_exact_mut(lanes)) {
        let x = (f64s::from_slice(wv) - f64s::splat(line_center)) * inv_doppler;
        let result = f64s::from_slice(xs) + scale_s * <f64s as std::simd::StdFloat>::exp(-x * x);
        xs.copy_from_slice(&result.to_array());
    }

    let n = wvnum.len();
    let inv_doppler_scalar = 1.0 / doppler_width;

    for i in n - remainder.len()..n {
        let x = (wvnum[i] - line_center) * inv_doppler_scalar;
        let x = (-x * x).exp();
        xs[i] += scale * x;
    }
}

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

        let wvnum_slice = wavenumber_cminv.as_slice().unwrap();
        let is_sorted = wvnum_slice.is_sorted();

        // Only sort if not already sorted
        let sorted_wvnum;
        let wavenumber_grid = if is_sorted {
            Grid1DView::new(wvnum_slice)
        } else {
            sorted_wvnum = {
                let mut tmp = wvnum_slice.to_owned();
                tmp.sort_by(|a, b| a.partial_cmp(b).unwrap());
                tmp
            };
            Grid1DView::new(sorted_wvnum.as_slice())
        };

        let mut xs = Array2::zeros((temperature.len(), wavenumber_cminv.len()));

        let min_wvnum = wavenumber_grid.x[0];
        let max_wvnum = wavenumber_grid.x[wavenumber_grid.x.len() - 1];

        let n_wavenumber = wavenumber_cminv.len();

        let mut start_wavenumber_idx: usize = 0;
        let mut end_wavenumber_idx: usize = 0;

        let line_slice = self.db.between_slice(
            min_wvnum - self.line_contribution_width,
            max_wvnum + self.line_contribution_width,
        );

        let map = self.gen_mol_param(line_slice, temperature.as_slice().unwrap())?;

        let mut batch_lines: Vec<&crate::optical::line::OpticalLine> = Vec::with_capacity(64);
        let mut batch_indices: Vec<(usize, usize)> = Vec::with_capacity(64);
        // Store coupling data: (temp_val, y_coupling_vec, g_coupling_vec)
        type CouplingData = (
            Option<Array1<f64>>,
            Option<Array1<f64>>,
            Option<Array1<f64>>,
        );
        let mut batch_coupling: Vec<CouplingData> = Vec::with_capacity(64);
        let mut batch_mol_params: Vec<&MolParam> = Vec::with_capacity(64);

        let mut process_batch = |lines: &[&crate::optical::line::OpticalLine],
                                 indices: &[(usize, usize)],
                                 coupling: &[CouplingData],
                                 params: &[&MolParam]| {
            let thread_pool = crate::threading::thread_pool().unwrap(); // Safe unwrap as called previously
            thread_pool.install(|| {
                Zip::indexed(temperature)
                    .and(pressure)
                    .and(pself)
                    .and(xs.axis_iter_mut(Axis(0)))
                    .par_for_each(|g, &temperature, &pressure, &pself, mut xs| {
                        for (i, line) in lines.iter().enumerate() {
                            let (start_wavenumber_idx, end_wavenumber_idx) = indices[i];
                            let (temp_val, y_coupling_vec, g_coupling_vec): &CouplingData =
                                &coupling[i];
                            let mol_param = params[i];

                            let sub_grid =
                                wavenumber_grid.slice(start_wavenumber_idx, end_wavenumber_idx);

                            let adjusted_line = line
                                .adjusted_parameters(
                                    temperature,
                                    pressure,
                                    pself,
                                    mol_param.partition_factor[g],
                                    mol_param.mol_mass,
                                )
                                .unwrap();

                            if let (Some(temp_val), Some(y_coupling_vec), Some(g_coupling_vec)) =
                                (temp_val, y_coupling_vec, g_coupling_vec)
                            {
                                // Line coupling case
                                let y_coupling = y_coupling_vec.interp1(
                                    temp_val,
                                    temperature,
                                    OutOfBoundsMode::Extend,
                                );
                                let g_coupling = g_coupling_vec.interp1(
                                    temp_val,
                                    temperature,
                                    OutOfBoundsMode::Extend,
                                );

                                let p_norm = pressure / 101325.0;

                                let scale_re = adjusted_line.line_intensity_re
                                    * (1.0 + p_norm * p_norm * g_coupling);
                                let scale_im =
                                    adjusted_line.line_intensity_re * (-p_norm * y_coupling);

                                w_jpole_assign(
                                    sub_grid.x,
                                    adjusted_line.line_center,
                                    adjusted_line.doppler_width,
                                    adjusted_line.y,
                                    scale_re,
                                    scale_im,
                                    xs.slice_mut(ndarray::s![
                                        start_wavenumber_idx..end_wavenumber_idx
                                    ])
                                    .as_slice_mut()
                                    .unwrap(),
                                );
                            } else {
                                // Standard case
                                split_and_assign(
                                    &adjusted_line,
                                    &sub_grid,
                                    xs.slice_mut(ndarray::s![
                                        start_wavenumber_idx..end_wavenumber_idx
                                    ])
                                    .as_slice_mut()
                                    .unwrap(),
                                );
                            }
                        }
                    });
            });
        };

        for line in line_slice.iter() {
            if line.line_intensity * max_p_self / (K_B * 1e-7 * 296.0) < self.cull_factor {
                continue;
            }

            // Find start index - don't reset to 0 since lines are sorted
            while start_wavenumber_idx < n_wavenumber
                && wavenumber_grid.x[start_wavenumber_idx]
                    < line.line_center - self.line_contribution_width
            {
                start_wavenumber_idx += 1;
            }

            // Reset end index to start for this line, then search forward
            while end_wavenumber_idx < n_wavenumber
                && wavenumber_grid.x[end_wavenumber_idx]
                    < line.line_center + self.line_contribution_width
            {
                end_wavenumber_idx += 1;
            }

            if start_wavenumber_idx >= n_wavenumber || start_wavenumber_idx == end_wavenumber_idx {
                continue;
            }

            // Pre-allocate line coupling arrays if needed
            let coupling_data = if enabling_line_coupling && !line.y_coupling.is_empty() {
                (
                    Some(Array1::<f64>::from_vec(line.coupling_temperature.clone())),
                    Some(Array1::<f64>::from_vec(line.y_coupling.clone())),
                    Some(Array1::<f64>::from_vec(line.g_coupling.clone())),
                )
            } else {
                (None, None, None)
            };

            batch_lines.push(line);
            batch_indices.push((start_wavenumber_idx, end_wavenumber_idx));
            batch_coupling.push(coupling_data);
            batch_mol_params.push(map.get(&(line.mol_id, line.iso_id)).unwrap());

            if batch_lines.len() >= 64 {
                process_batch(
                    &batch_lines,
                    &batch_indices,
                    &batch_coupling,
                    &batch_mol_params,
                );
                batch_lines.clear();
                batch_indices.clear();
                batch_coupling.clear();
                batch_mol_params.clear();
            }
        }

        if !batch_lines.is_empty() {
            process_batch(
                &batch_lines,
                &batch_indices,
                &batch_coupling,
                &batch_mol_params,
            );
        }

        if is_sorted {
            Ok(xs)
        } else {
            println!("Unsored wavenumber grid, have to sort output");
            // Have to sort the output
            let sort_idx = argsort_f64(wvnum_slice);

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
            .spectral_grid()
            .ok_or(anyhow::anyhow!("Wavenumbers not found in inputs"))?
            .central_wavenumber_cminv();

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

    fn is_scatterer(&self) -> bool {
        false
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
