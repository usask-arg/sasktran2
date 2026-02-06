use std::{fs, io::{BufRead, BufReader}, ops::AddAssign, path::PathBuf};

use crate::{
    interpolation::{grid1d::Grid1D, linear::Interp1Weights, OutOfBoundsMode},
    optical::traits::OpticalProperty,
    prelude::*,
};

pub struct XsecAtCondition {
    pub wvnum: Grid1D,
    pub xsec: Array1<f64>,
    pub temperature_k: f64,
    pub pressure_pa: f64,
}


/// Raw storage structure that holds the cross section data from something like the
/// HITRAN fwf files, or the LBLRTM files.  This is basically lists of cross section data as a function
/// of wavenumber at a given temperature and pressure.  We assume wavenumber is specified in
/// vacuum cm^-1, and the cross section is in cm^2/molecule.
pub struct XsecDatabase(pub Vec<XsecAtCondition>);

/// Let's us concatenate two XsecDatabase objects together.  This is useful for combining
/// databases from different pressure/temperatures together
impl AddAssign for XsecDatabase {
    fn add_assign(&mut self, other: Self) {
        self.0.extend(other.0);
    }
}

#[allow(dead_code)]
struct Header {
    pub short_mol_name: String,
    pub wvnum_start: f64,
    pub wvnum_end: f64,
    pub num_points: i64,
    pub temperature: f64,
    pub zero: f64,
    pub wvnum_space: f64,
    pub pressure: f64,
}

impl Header {
    fn new() -> Self {
        Header {
            short_mol_name: String::new(),
            wvnum_start: 0.0,
            wvnum_end: 0.0,
            num_points: 0,
            temperature: 0.0,
            zero: 0.0,
            wvnum_space: 0.0,
            pressure: 0.0,
        }
    }
}

fn read_hitran_header(line: &str) -> Header {
    let short_mol_name = line[10..20].trim().to_string();
    let wvnum_start: f64 = line[20..31].trim().parse().unwrap();
    let wvnum_end: f64 = line[31..42].trim().parse().unwrap();
    let num_points: i64 = line[42..50].trim().parse().unwrap();
    let temperature: f64 = line[50..58].trim().parse().unwrap();
    let zero: f64 = line[59..62].trim().parse().unwrap();
    // let wvnum_space: f64 = line[62..71].trim().parse().unwrap();  // unsure on this one...
    let pressure: f64 = line[71..80].trim().parse().unwrap();

    let wvnum_space = (wvnum_end - wvnum_start) / ((num_points - 1) as f64);

    Header {
        short_mol_name,
        wvnum_start,
        wvnum_end,
        num_points,
        temperature,
        zero,
        wvnum_space,
        pressure,
    }
}

pub fn read_fwf_xsec(path: PathBuf) -> Option<XsecDatabase> {
    let file = fs::File::open(path).unwrap();
    let reader = BufReader::new(file);

    let mut conditions: Vec<XsecAtCondition> = vec![];

    let mut xs: Vec<f64> = vec![];
    let mut wvnum: Vec<f64> = vec![];

    let mut header = Header::new();
    let mut cur_wvnum = 0.0;
    let mut cur_index = 0;

    for line_result in reader.lines() {
        let line = line_result.ok()?;

        // Check if this is header information
        if line.len() == 102 {
            header = read_hitran_header(&line);
            cur_wvnum = header.wvnum_start;
        } else {
            for i in 0..10 {
                let start = i * 10;
                let end = start + 10;
                let value: f64 = line[start..end].trim().parse().unwrap();
                xs.push(value);
                wvnum.push(cur_wvnum);
                cur_wvnum += header.wvnum_space;

                cur_index += 1;
                if cur_index >= header.num_points {
                    conditions.push(XsecAtCondition {
                        wvnum: Grid1D::new(Array1::from(wvnum.clone())),
                        xsec: Array1::from(xs.clone()),
                        temperature_k: header.temperature,
                        pressure_pa: header.pressure,
                    });

                    xs.clear();
                    wvnum.clear();
                    break;
                }
            }
        }
    }

    Some(XsecDatabase(conditions))
}

pub fn read_fwf_folder(folder: PathBuf) -> Option<XsecDatabase> {
    let paths = fs::read_dir(folder).unwrap();

    let mut combined_ds = XsecDatabase(vec![]);

    for path in paths {
        let path = path.unwrap().path();
        if path.extension().unwrap() == "xsc" {
            let dbase = read_fwf_xsec(path.clone());

            combined_ds += dbase?;
        }
    }

    Some(combined_ds)
}


impl XsecDatabase {
    /// Get interpolation weights for pressure and temperature.
    /// 
    /// This is a simple algorithm that does not interpolate in pressure and finds
    /// the two closest temperatures to interpolate between. Returns a vector of
    /// (index, weight, d_weight) tuples for the cross sections to use, where
    /// d_weight is the derivative of the weight with respect to temperature.
    fn get_pt_interp_weights(&self, pressure_pa: f64, temperature_k: f64) -> Vec<(usize, f64, f64)> {
        // Find the closest pressure (no pressure interpolation)
        let mut closest_pressure_idx = 0usize;
        let mut min_pressure_diff = f64::INFINITY;
        
        for (idx, condition) in self.0.iter().enumerate() {
            let pressure_diff = (condition.pressure_pa - pressure_pa).abs();
            if pressure_diff < min_pressure_diff {
                min_pressure_diff = pressure_diff;
                closest_pressure_idx = idx;
            }
        }
        
        // Now find all conditions at this pressure
        let target_pressure = self.0[closest_pressure_idx].pressure_pa;
        let mut conditions_at_pressure: Vec<(usize, f64)> = self.0
            .iter()
            .enumerate()
            .filter(|(_, cond)| (cond.pressure_pa - target_pressure).abs() < 1e-6)
            .map(|(idx, cond)| (idx, cond.temperature_k))
            .collect();
        
        // Sort by temperature
        conditions_at_pressure.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        
        if conditions_at_pressure.len() == 1 {
            // Only one temperature available, use it with weight 1.0 and zero derivative
            return vec![(conditions_at_pressure[0].0, 1.0, 0.0)];
        }
        
        // Find the two closest temperatures
        let mut lower_idx = 0;
        let mut upper_idx = 1;
        
        for i in 0..conditions_at_pressure.len() {
            if conditions_at_pressure[i].1 <= temperature_k {
                lower_idx = i;
            }
            if i > 0 && conditions_at_pressure[i].1 >= temperature_k && upper_idx == 1 {
                upper_idx = i;
            }
        }
        
        // Handle edge cases
        if temperature_k <= conditions_at_pressure[0].1 {
            // Below lowest temperature, use the two lowest
            lower_idx = 0;
            upper_idx = 1;
        } else if temperature_k >= conditions_at_pressure[conditions_at_pressure.len() - 1].1 {
            // Above highest temperature, use the two highest
            lower_idx = conditions_at_pressure.len() - 2;
            upper_idx = conditions_at_pressure.len() - 1;
        } else {
            // Ensure upper_idx is immediately after lower_idx
            upper_idx = lower_idx + 1;
        }
        
        // Calculate interpolation weights and derivatives
        let temp_lower = conditions_at_pressure[lower_idx].1;
        let temp_upper = conditions_at_pressure[upper_idx].1;
        let temp_diff = temp_upper - temp_lower;
        
        let weight_upper = (temperature_k - temp_lower) / temp_diff;
        let weight_lower = 1.0 - weight_upper;
        
        // Derivatives with respect to temperature
        let d_weight_upper = 1.0 / temp_diff;
        let d_weight_lower = -1.0 / temp_diff;
        
        vec![
            (conditions_at_pressure[lower_idx].0, weight_lower, d_weight_lower),
            (conditions_at_pressure[upper_idx].0, weight_upper, d_weight_upper),
        ]
    }
}

impl OpticalProperty for XsecDatabase {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn crate::atmosphere::StorageInputs,
        aux_inputs: &dyn crate::optical::traits::AuxOpticalInputs,
        optical_quantities: &mut crate::optical::storage::OpticalQuantities,
    ) -> Result<()> {
        use crate::optical::traits::param_from_storage_or_aux;
        use ndarray::Zip;
        
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let temperature_k = param_from_storage_or_aux(inputs, aux_inputs, "temperature_k")?;
        let pressure_pa = param_from_storage_or_aux(inputs, aux_inputs, "pressure_pa")?;
        
        let _ = optical_quantities.resize(temperature_k.len(), wavenumber_cminv.len());
        
        let xs = &mut optical_quantities.cross_section;
        
        Zip::from(xs.rows_mut())
            .and(temperature_k.view())
            .and(pressure_pa.view())
            .par_for_each(|mut row, temp, press| {
                // Get interpolation weights for this pressure/temperature
                let weights = self.get_pt_interp_weights(*press, *temp);
                
                // For each set of conditions, interpolate in wavenumber and accumulate
                for (idx, weight, _d_weight) in weights {
                    let condition = &self.0[idx];
                    
                    // Interpolate in wavenumber for each point
                    for (j, wv) in wavenumber_cminv.iter().enumerate() {
                        let wvnum_weights = condition.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
                        let local_xs = condition.xsec[wvnum_weights[0].0] * wvnum_weights[0].1
                            + condition.xsec[wvnum_weights[1].0] * wvnum_weights[1].1;
                        
                        row[j] += local_xs * weight;
                    }
                }
            });
        
        Ok(())
    }

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn crate::atmosphere::StorageInputs,
        aux_inputs: &dyn crate::optical::traits::AuxOpticalInputs,
        d_optical_quantities: &mut HashMap<String, crate::optical::storage::OpticalQuantities>,
    ) -> Result<()> {
        use crate::optical::traits::param_from_storage_or_aux;
        use ndarray::Zip;
        
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let temperature_k = param_from_storage_or_aux(inputs, aux_inputs, "temperature_k")?;
        let pressure_pa = param_from_storage_or_aux(inputs, aux_inputs, "pressure_pa")?;
        
        // Create or resize derivative storage for temperature
        if d_optical_quantities.contains_key("temperature_k") {
            d_optical_quantities
                .get_mut("temperature_k")
                .unwrap()
                .resize(temperature_k.len(), wavenumber_cminv.len());
        } else {
            d_optical_quantities.insert(
                "temperature_k".to_string(),
                crate::optical::storage::OpticalQuantities::new(
                    temperature_k.len(),
                    wavenumber_cminv.len(),
                    false,
                ),
            );
        }
        
        let d_xs = &mut d_optical_quantities
            .get_mut("temperature_k")
            .unwrap()
            .cross_section;
        
        Zip::from(d_xs.rows_mut())
            .and(temperature_k.view())
            .and(pressure_pa.view())
            .par_for_each(|mut row, temp, press| {
                // Get interpolation weights and derivatives for this pressure/temperature
                let weights = self.get_pt_interp_weights(*press, *temp);
                
                // For each set of conditions, interpolate in wavenumber and accumulate derivatives
                for (idx, _weight, d_weight) in weights {
                    let condition = &self.0[idx];
                    
                    // Interpolate in wavenumber for each point
                    for (j, wv) in wavenumber_cminv.iter().enumerate() {
                        let wvnum_weights = condition.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
                        let local_xs = condition.xsec[wvnum_weights[0].0] * wvnum_weights[0].1
                            + condition.xsec[wvnum_weights[1].0] * wvnum_weights[1].1;
                        
                        row[j] += local_xs * d_weight;
                    }
                }
            });
        
        Ok(())
    }

    fn is_scatterer(&self) -> bool {
        false
    }
}