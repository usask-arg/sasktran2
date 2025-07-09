use crate::prelude::*;
use std::{
    fs,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use super::{OpticalLine, OpticalLineDB};

pub fn hitran_molecule_file(molecule: &str, directory: &Path) -> Result<PathBuf> {
    // let dir = directory.join("line_files_By_Molecule");

    let mol_filename = molecule.to_ascii_uppercase().as_str().to_owned() + ".data";
    let mol_filepath = directory.join(mol_filename);

    // Check for file
    if mol_filepath.is_file() {
        Ok(mol_filepath)
    } else {
        Err(anyhow::anyhow!(
            "No HITRAN data file for molecule {}",
            molecule
        ))
    }
}

pub fn read_hitran_line_file(data_path: PathBuf) -> Result<OpticalLineDB> {
    // File format: (name, width)
    // molec_id, 2
    // local_iso_id, 1
    // transition_wavenumber, 12
    // line_intensity, 10
    // einstein_a, 10
    // gamma_air, 5
    // gamma_self, 5
    // lower_state_energy, 10
    // n_air, 4
    // delta_air, 8
    // upper_quanta, 15
    // lower_quanta, 15
    // upper_local_quanta, 15
    // lower_local_quanta, 15
    // error_codes, 6
    // reference_codes, 12
    // line_mixing_flag, 1
    // upper_g, 7
    // lower_g, 7
    let file = fs::File::open(data_path)?;
    let reader = BufReader::new(file);

    let line_iterator = reader.lines();
    let mut lines: Vec<OpticalLine> = Vec::with_capacity(line_iterator.size_hint().1.unwrap_or(0));
    for line_result in line_iterator {
        let line = line_result?;

        // Check if this is a full line
        if line.len() < 160 {
            continue;
        }

        let line_intensity: f64 = line[15..25].trim().parse()?;

        if line_intensity == 0.0 {
            // Line coupling line, ignore
            continue;
        }

        let mol_id: i32 = line[0..2].trim().parse()?;
        let iso_id: i32 = match line[2..3].trim() {
            "0" => 10,
            "A" => 11,
            "B" => 12,
            other => other.parse()?,
        };

        let transition_wavenumber: f64 = line[3..15].trim().parse()?;
        let _einstein_a: f64 = line[25..35].trim().parse()?;
        let air_broadened_width: f64 = line[35..40].trim().parse()?;
        let self_broadened_width: f64 = line[40..45].trim().parse()?;
        let lower_state_energy: f64 = line[45..55].trim().parse()?;
        let temperature_dependence: f64 = line[55..59].trim().parse()?;
        let pressure_shift: f64 = line[59..67].trim().parse()?;

        // Unused parameters
        let _upper_vibrational_quanta: &str = line[67..82].trim();
        let _lower_vibrational_quanta: &str = line[82..97].trim();
        let _upper_local_quanta: &str = line[97..112].trim();
        let _lower_local_quanta: &str = line[112..127].trim();
        let _error_codes: &str = line[127..133].trim();
        let _reference_codes: &str = line[133..145].trim();
        let _line_mixing_flag: &str = line[145..146].trim();
        let _upper_g: f64 = line[146..153].trim().parse()?;
        let _lower_g: f64 = line[153..160].trim().parse()?;

        let y_coupling: Vec<f64> = Vec::new();
        let g_coupling: Vec<f64> = Vec::new();
        let coupling_temperature: Vec<f64> = Vec::new();

        // TODO: line mixing?

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
    fn test_read_hitran_line_file() {
        let co2_file = PathBuf::from("../../tests/data/CO2.data");

        let result = read_hitran_line_file(co2_file);
        assert!(result.is_ok());
        let db: OpticalLineDB = result.unwrap();
        assert!(!db.lines.is_empty());
    }
}
