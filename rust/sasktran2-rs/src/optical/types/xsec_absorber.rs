use std::{
    fs,
    io::{BufRead, BufReader},
    ops::AddAssign,
    path::PathBuf,
};

use crate::{
    interpolation::{OutOfBoundsMode, grid1d::Grid1D, linear::Interp1Weights},
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
/// vacuum cm^-1, and the cross section is in cm^2/molecule (inside the file).
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
    // HITRAN xsc files have multiple format variants
    // Extract numbers by splitting on whitespace and parsing

    let tokens: Vec<&str> = line.split_whitespace().collect();
    let mut numbers: Vec<f64> = Vec::new();

    // Extract all valid numbers from tokens
    for token in tokens.iter() {
        if let Ok(num) = token.parse::<f64>() {
            numbers.push(num);
        }
    }

    // The general pattern is: wvnum_start wvnum_end num_points temperature pressure [...]
    // We need at least 5 numbers
    if numbers.len() >= 5 {
        // Most xsc headers put num_points in the 3rd numeric field.
        // Keep fallback logic for variants where this assumption does not hold.
        let mut num_points_idx = 2;
        if !(numbers[2] > 0.0 && numbers[2].fract().abs() < 0.0001) {
            for (idx, num) in numbers.iter().enumerate() {
                if *num > 0.0 && num.fract().abs() < 0.0001 {
                    num_points_idx = idx;
                    break;
                }
            }
        }

        let wvnum_start = numbers[0];
        let wvnum_end = numbers[1];
        let num_points = numbers[num_points_idx] as i64;

        // Temperature and pressure typically follow num_points
        let temp_idx = num_points_idx + 1;
        let press_idx = num_points_idx + 2;

        let temperature = if temp_idx < numbers.len() {
            numbers[temp_idx]
        } else {
            296.0
        };
        let pressure = if press_idx < numbers.len() {
            numbers[press_idx]
        } else {
            101325.0
        };

        let mut header = Header {
            short_mol_name: String::new(),
            wvnum_start,
            wvnum_end,
            num_points,
            temperature,
            zero: 0.0,
            wvnum_space: 0.0,
            pressure,
        };

        header.wvnum_space = (wvnum_end - wvnum_start) / ((num_points - 1) as f64);
        header
    } else {
        Header::new()
    }
}

pub fn read_fwf_xsec(path: PathBuf) -> Option<XsecDatabase> {
    let file = fs::File::open(path).unwrap();
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let mut conditions: Vec<XsecAtCondition> = vec![];

    // Process file: each entry starts with a header, followed by data lines
    while let Some(header_line) = lines.next() {
        let header_line = header_line.ok()?;

        // Parse the header to get num_points and other metadata
        let header = read_hitran_header(&header_line);

        // Now read exactly num_points data values
        let mut xs: Vec<f64> = Vec::with_capacity(header.num_points as usize);
        let mut wvnum: Vec<f64> = Vec::with_capacity(header.num_points as usize);
        let mut cur_wvnum = header.wvnum_start;

        'data_reading: while xs.len() < header.num_points as usize {
            let data_line = match lines.next() {
                Some(Ok(line)) => line,
                _ => break 'data_reading, // EOF or error
            };

            // Parse up to 10 values from this line (10 chars each)
            for i in 0..10 {
                if xs.len() >= header.num_points as usize {
                    break 'data_reading;
                }

                let start = i * 10;
                let end = start + 10;

                // Check bounds
                if start >= data_line.len() {
                    break;
                }

                // Extract field
                let field = if end <= data_line.len() {
                    &data_line[start..end]
                } else {
                    &data_line[start..]
                };

                // Parse value
                let trimmed = field.trim();
                if trimmed.is_empty() {
                    continue;
                }

                if let Ok(value) = trimmed.parse::<f64>() {
                    xs.push(value / 1.0e4); // Convert from cm^2/molecule to m^2/molecule
                    wvnum.push(cur_wvnum);
                    cur_wvnum += header.wvnum_space;
                }
            }
        }

        // Only add the condition if we got the expected number of points
        if xs.len() == header.num_points as usize {
            conditions.push(XsecAtCondition {
                wvnum: Grid1D::new(Array1::from(wvnum)),
                xsec: Array1::from(xs),
                temperature_k: header.temperature,
                pressure_pa: header.pressure,
            });
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
    fn get_pt_interp_weights(
        &self,
        pressure_pa: f64,
        temperature_k: f64,
    ) -> Vec<(usize, f64, f64)> {
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
        let mut conditions_at_pressure: Vec<(usize, f64)> = self
            .0
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

        for (i, (_, temp)) in conditions_at_pressure.iter().enumerate() {
            if *temp <= temperature_k {
                lower_idx = i;
            }
            if i > 0 && *temp >= temperature_k && upper_idx == 1 {
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
            (
                conditions_at_pressure[lower_idx].0,
                weight_lower,
                d_weight_lower,
            ),
            (
                conditions_at_pressure[upper_idx].0,
                weight_upper,
                d_weight_upper,
            ),
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
                        let wvnum_weights =
                            condition.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
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
                        let wvnum_weights =
                            condition.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
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

#[cfg(test)]
mod generated_format_tests {
    use super::read_fwf_xsec;
    use std::{
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn write_temp_xsc(contents: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time should be after UNIX_EPOCH")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("xsec_absorber_test_{nonce}.xsc"));
        fs::write(&path, contents).expect("failed to write temporary xsc test file");
        path
    }

    fn build_block(temperature_k: f64, start_value_cm2: f64) -> String {
        let num_points = 1001usize;
        let wv_start = 49382.7160;
        let wv_end = 50382.7160;
        let pressure = 101325.25;

        let mut out = format!(
            "   O2-SR {wv_start:10.4} {wv_end:10.4} {num_points:8} {temperature_k:8.1} {pressure:8.2}  0.000E+00  0  GENERATED\n"
        );

        let mut i = 0usize;
        while i < num_points {
            let mut line = String::new();
            for j in 0..10 {
                let k = i + j;
                if k >= num_points {
                    break;
                }
                let value = start_value_cm2 + (k as f64) * 1.0e-24;
                line.push_str(&format!("{value:10.3E}"));
            }
            out.push_str(&line);
            out.push('\n');
            i += 10;
        }

        out
    }

    #[test]
    fn reads_generated_style_fixed_width_blocks() {
        // Header shape mirrors the generated Schumann-Runge file format.
        let data = format!("{}{}", build_block(90.0, 1.0e-20), build_block(300.0, 2.0e-20));

        let path = write_temp_xsc(&data);
        let db = read_fwf_xsec(path.clone()).expect("parser should load valid generated-style xsc");
        let _ = fs::remove_file(path);

        assert_eq!(db.0.len(), 2);
        assert!((db.0[0].temperature_k - 90.0).abs() < 1e-12);
        assert!((db.0[1].temperature_k - 300.0).abs() < 1e-12);
        assert_eq!(db.0[0].xsec.len(), 1001);

        // Reader converts cm^2 to m^2 by dividing by 1e4.
        assert!((db.0[0].xsec[0] - 1.0e-24).abs() < 1e-30);
        assert!((db.0[1].xsec[0] - 2.0e-24).abs() < 1e-30);

        // Wavenumber grid is reconstructed from start/end/num_points in header.
        assert!((db.0[0].wvnum.x[0] - 49382.7160).abs() < 1e-9);
        assert!((db.0[0].wvnum.x[1000] - 50382.7160).abs() < 1e-9);
    }

    #[test]
    #[ignore = "Local verification test: requires generated cross-sections/xs/O2SCHRUNG file"]
    fn reads_generated_o2_schumann_runge_file() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../cross-sections/xs/O2SCHRUNG");

        let db = read_fwf_xsec(path).expect("generated O2SCHRUNG should be parseable");

        assert_eq!(db.0.len(), 2);
        assert!((db.0[0].temperature_k - 90.0).abs() < 1e-6);
        assert!((db.0[1].temperature_k - 300.0).abs() < 1e-6);
        assert_eq!(db.0[0].xsec.len(), db.0[1].xsec.len());
        assert!(db.0[0].xsec.len() > 100);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data_path(filename: &str) -> PathBuf {
        // Assuming tests are run from the workspace root
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("tests/data")
            .join(filename)
    }

    #[test]
    fn test_read_single_file() {
        let path = test_data_path("CHClFCF3_287.0_0.0_675.2-715.0_00.xsc");
        let db = read_fwf_xsec(path);

        assert!(db.is_some());
        let db = db.unwrap();
        assert!(!db.0.is_empty());

        // Check that we loaded the data
        let first_condition = &db.0[0];
        assert!(!first_condition.wvnum.x.is_empty());
        assert_eq!(first_condition.wvnum.x.len(), first_condition.xsec.len());

        // Check temperature and pressure from the filename
        assert!((first_condition.temperature_k - 287.0).abs() < 1.0);
    }

    #[test]
    fn test_read_chcnft1_file() {
        let path = test_data_path("CHCNFT1");
        let db = read_fwf_xsec(path);

        assert!(db.is_some());
        let db = db.unwrap();
        assert!(!db.0.is_empty());

        // Check basic properties
        let first_condition = &db.0[0];
        assert!(!first_condition.wvnum.x.is_empty());
        assert!((first_condition.temperature_k - 324.1).abs() < 1.0);
    }

    #[test]
    fn test_get_pt_interp_weights() {
        let path = test_data_path("CHClFCF3_287.0_0.0_675.2-715.0_00.xsc");
        let db = read_fwf_xsec(path).expect("Failed to read test file");

        // Test at the exact temperature and pressure in the file
        let weights = db.get_pt_interp_weights(0.0, 287.0);

        // Should get some weights
        assert!(!weights.is_empty());

        // Weights should sum to ~1.0
        let weight_sum: f64 = weights.iter().map(|(_, w, _)| w).sum();
        assert!((weight_sum - 1.0).abs() < 1e-10);

        // All d_weights should be finite
        for (_, _, d_weight) in &weights {
            assert!(d_weight.is_finite());
        }
    }

    #[test]
    fn test_temperature_interpolation_weights() {
        let path = test_data_path("CHCNFT1");
        let db = read_fwf_xsec(path).expect("Failed to read test file");

        // Test interpolation between temperatures if multiple exist
        let weights = db.get_pt_interp_weights(0.0, 320.0);

        assert!(!weights.is_empty());

        // Check weights sum to 1
        let weight_sum: f64 = weights.iter().map(|(_, w, _)| w).sum();
        assert!((weight_sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_wavenumber_grid() {
        let path = test_data_path("CHClFCF3_287.0_0.0_675.2-715.0_00.xsc");
        let db = read_fwf_xsec(path).expect("Failed to read test file");

        let condition = &db.0[0];

        // Check that wavenumber grid is monotonically increasing
        for i in 1..condition.wvnum.x.len() {
            assert!(condition.wvnum.x[i] > condition.wvnum.x[i - 1]);
        }

        // Check that cross sections are non-negative
        for xs in condition.xsec.iter() {
            assert!(*xs >= 0.0);
        }
    }
}
