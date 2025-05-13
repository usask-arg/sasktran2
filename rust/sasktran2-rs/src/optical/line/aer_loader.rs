use crate::prelude::*;
use std::{
    fs,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use super::{OpticalLine, OpticalLineDB};

pub fn aer_molecule_file(molecule: &str, directory: &Path) -> Result<PathBuf> {
    let dir = directory.join("line_files_By_Molecule");

    let mol_ending = "_".to_string() + molecule.to_ascii_uppercase().as_str();

    // Read directory and collect only directories (folders)
    #[allow(clippy::filter_next)]
    let folder = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|entry| entry.ok()) // Ignore errors
        .map(|entry| entry.path()) // Get PathBuf
        .filter(|path| path.is_dir()) // Keep only directories
        .filter(|path| {
            path.file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .to_uppercase()
                .ends_with(&mol_ending)
        })
        .next()
        .ok_or(anyhow::anyhow!(
            "No directory found for molecule: {}",
            molecule
        ))?;

    let stem = folder.file_stem().unwrap();

    Ok(folder.join(
        stem.to_str()
            .ok_or(anyhow::anyhow!("Failed to convert OsStr to str"))?,
    ))
}

pub fn read_aer_line_file(data_path: PathBuf) -> Result<OpticalLineDB> {
    // File format: (name, width)
    // molec_id, 2
    // local_iso_id, 1
    // transition_wavenumber, 12
    // line_intensity, 10
    // rsq, 10
    // gamma_air, 5
    // gamma_self, 5
    // lower_state_energy, 10
    // n_air, 4
    // delta_air, 8
    // upper_quanta, 3
    // lower_quanta, 3
    // upper_local_quanta, 9
    // lower_local_quanta, 9
    // error_codes, 3
    // reference_codes, 6
    let file = fs::File::open(data_path)?;
    let reader = BufReader::new(file);

    let mut line_iterator = reader.lines();
    let mut lines: Vec<OpticalLine> = Vec::with_capacity(line_iterator.size_hint().1.unwrap_or(0));
    while let Some(line_result) = line_iterator.next() {
        // Convert FORTAN D to E
        let line = line_result?.replace("D", "E");

        // Check if this is header information
        if line.starts_with(">") || line.starts_with("%") || line.len() < 100 {
            continue;
        }

        let line_intensity: f64 = line[15..25].trim().replace("D", "E").parse()?;

        if line_intensity == 0.0 {
            // Line coupling line, ignore
            continue;
        }

        let mol_id: i32 = line[0..2].trim().parse()?;
        let iso_id: i32 = line[2..3].parse()?;

        let transition_wavenumber: f64 = line[3..15].trim().parse()?;
        let _r_squared: f64 = line[25..35].trim().parse()?;
        let air_broadened_width: f64 = line[35..40].trim().parse()?;
        let self_broadened_width: f64 = line[40..45].trim().parse()?;
        let lower_state_energy: f64 = line[45..55].trim().parse()?;
        let temperature_dependence: f64 = line[55..59].trim().parse()?;
        let pressure_shift: f64 = line[59..67].trim().parse()?;

        // Unused parameters
        let _upper_vibrational_quanta: &str = line[67..70].trim();
        let _lower_vibrational_quanta: &str = line[70..73].trim();
        let _upper_local_quanta: &str = line[73..82].trim();
        let _lower_local_quanta: &str = line[82..91].trim();
        let _error_codes: &str = line[91..94].trim();
        let reference_codes: &str = line[94..100].trim();

        let mut y_coupling: Vec<f64> = Vec::new();
        let mut g_coupling: Vec<f64> = Vec::new();
        let mut coupling_temperature: Vec<f64> = Vec::new();

        if reference_codes.ends_with("-1")
            || reference_codes.ends_with("-2")
            || reference_codes.ends_with("-3")
        {
            // if we don't end with 0, the next line contains something ectra
            let next_line_result = line_iterator.next();

            if reference_codes.ends_with("-1") {
                // indicates there is line coupling data

                if let Some(next_line_result) = next_line_result {
                    let next_line = next_line_result?.replace("D", "E");
                    coupling_temperature = vec![200.0, 250.0, 296.0, 340.0];
                    y_coupling = vec![
                        next_line[3..15].trim().parse()?,
                        next_line[26..39].trim().parse()?,
                        next_line[50..63].trim().parse()?,
                        next_line[74..87].trim().parse()?,
                    ];

                    g_coupling = vec![
                        next_line[15..26].trim().parse()?,
                        next_line[39..50].trim().parse()?,
                        next_line[63..74].trim().parse()?,
                        next_line[87..98].trim().parse()?,
                    ];
                }
            }
        }

        // Create the OpticalLine object
        let line = OpticalLine {
            line_center: transition_wavenumber,
            line_intensity,
            lower_energy: lower_state_energy,
            gamma_air: air_broadened_width,
            gamma_self: self_broadened_width,
            delta_air: pressure_shift,
            n_air: temperature_dependence,
            mol_id,
            iso_id,
            y_coupling,
            g_coupling,
            coupling_temperature,
        };

        // Add the line to the vector
        lines.push(line);
    }

    Ok(OpticalLineDB { lines })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::optical::line::db::OpticalLineDB;

    #[test]
    fn test_read_aer_line_file() {
        let o2_file = PathBuf::from("../../tests/data/02_CO2");

        let result = read_aer_line_file(o2_file);
        assert!(result.is_ok());
        let db: OpticalLineDB = result.unwrap();
        assert!(!db.lines.is_empty());
    }
}
