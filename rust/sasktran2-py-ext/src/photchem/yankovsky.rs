use std::collections::HashMap;

use anyhow::anyhow;
use ndarray::{Array2, Axis, s};
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use sasktran2_rs::optical::line::hitran_loader::{hitran_molecule_file, read_hitran_line_file};
use sasktran2_rs::photchem::emission::{
    EmissionBand, MCDADE_OXYGEN_GREEN_LINE_C0, MCDADE_OXYGEN_GREEN_LINE_C1,
    MCDADE_OXYGEN_GREEN_LINE_C2, MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_1S_S,
    MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_558_S, O2_A_BAND_CENTER_WAVELENGTH_NM,
    O2_B0_X0_EINSTEIN_A_S, O2_B1_X1_EINSTEIN_A_S, OXYGEN_GREEN_LINE_EINSTEIN_A_S,
    OXYGEN_GREEN_LINE_WAVELENGTH_NM, mcdade_oxygen_green_line_photon_ver,
};
use sasktran2_rs::photchem::models::*;

use crate::prelude::IntoPyResult;

#[pyclass(unsendable)]
pub struct PyYankovsky {
    pub model: Yankovsky,
}

#[pymethods]
impl PyYankovsky {
    #[new]
    pub fn new() -> Self {
        Self {
            model: Yankovsky::new(),
        }
    }

    pub fn solve(&self, ds: &Bound<PyAny>) -> PyResult<Py<PyAny>> {
        let temperature = ds.get_item("temperature")?;
        let temperature = temperature.call_method0("to_numpy")?;
        let temperature = temperature.cast_into::<PyArray1<f64>>()?;
        let temperature = temperature.readonly();
        let temperature = temperature.as_array();

        let num_alt = temperature.len();

        let o2_density = ds.get_item("o2_density")?;
        let o2_density = o2_density.call_method0("to_numpy")?;
        let o2_density = o2_density.cast_into::<PyArray1<f64>>()?;
        let o2_density = o2_density.readonly();
        let o2_density = o2_density.as_array();

        let o3_density = ds.get_item("o3_density")?;
        let o3_density = o3_density.call_method0("to_numpy")?;
        let o3_density = o3_density.cast_into::<PyArray1<f64>>()?;
        let o3_density = o3_density.readonly();
        let o3_density = o3_density.as_array();

        let n2_density = ds.get_item("n2_density")?;
        let n2_density = n2_density.call_method0("to_numpy")?;
        let n2_density = n2_density.cast_into::<PyArray1<f64>>()?;
        let n2_density = n2_density.readonly();
        let n2_density = n2_density.as_array();

        let co2_density = match ds.get_item("co2_density") {
            Ok(co2_density) => {
                let co2_density = co2_density.call_method0("to_numpy")?;
                let co2_density = co2_density.cast_into::<PyArray1<f64>>()?;
                let co2_density = co2_density.readonly();
                Some(co2_density.as_array().to_vec())
            }
            Err(_) => None,
        };

        let o3p_density = match ds.get_item("o3p_density") {
            Ok(o3p_density) => {
                let o3p_density = o3p_density.call_method0("to_numpy")?;
                let o3p_density = o3p_density.cast_into::<PyArray1<f64>>()?;
                let o3p_density = o3p_density.readonly();
                Some(o3p_density.as_array().to_vec())
            }
            Err(_) => match ds.get_item("o_density") {
                Ok(o_density) => {
                    let o_density = o_density.call_method0("to_numpy")?;
                    let o_density = o_density.cast_into::<PyArray1<f64>>()?;
                    let o_density = o_density.readonly();
                    Some(o_density.as_array().to_vec())
                }
                Err(_) => None,
            },
        };

        let wavelength = ds.get_item("wavelength")?;
        let wavelength = wavelength.call_method0("to_numpy")?;
        let wavelength = wavelength.cast_into::<PyArray1<f64>>()?;
        let wavelength = wavelength.readonly();
        let wavelength = wavelength.as_array();
        let wavelength = wavelength
            .as_slice()
            .ok_or(anyhow!("Wavelength array must be contiguous"))
            .into_pyresult()?;

        let actinic_flux = ds.get_item("actinic_flux")?;
        let actinic_flux = actinic_flux.call_method0("to_numpy")?;
        let actinic_flux = actinic_flux.cast_into::<PyArray2<f64>>()?;
        let actinic_flux = actinic_flux.readonly();
        let actinic_flux = actinic_flux.as_array();

        // Compute the photolysis rates
        let mut photolysis_rate = Array2::<f64>::zeros((self.model.photo_reactions.len(), num_alt));
        for (i, reaction) in self.model.photo_reactions.iter().enumerate() {
            let species = reaction.in_molecule.base_type.to_string();

            let xs = ds.get_item(format!("{}_xs", species.to_lowercase()))?;
            let xs = xs.call_method0("to_numpy")?;
            let xs = xs.cast_into::<PyArray2<f64>>()?;
            let xs = xs.readonly();
            let xs = xs.as_array();

            let line_actinic_flux_name = reaction
                .excitation_band
                .as_deref()
                .map(|band| format!("{}_actinic_flux", band.replace('-', "_")));

            let rate = if reaction.line_center_nm.is_some() {
                if let Some(name) = line_actinic_flux_name
                    && let Ok(line_actinic_flux) = ds.get_item(&name)
                {
                    let line_actinic_flux = line_actinic_flux.call_method0("to_numpy")?;
                    let line_actinic_flux = line_actinic_flux.cast_into::<PyArray1<f64>>()?;
                    let line_actinic_flux = line_actinic_flux.readonly();
                    let line_actinic_flux = line_actinic_flux.as_array();

                    let line_center_nm = reaction
                        .line_center_nm
                        .ok_or(anyhow!("Missing line center for line photolysis"))
                        .into_pyresult()?;
                    let line_wavelength = [line_center_nm];
                    let line_flux = line_actinic_flux.insert_axis(Axis(0));
                    let dummy_xs = Array2::<f64>::zeros((1, num_alt));

                    calculate_photolysis_rate(
                        reaction,
                        &line_wavelength,
                        line_flux,
                        dummy_xs.view(),
                    )
                    .into_pyresult()?
                } else {
                    calculate_photolysis_rate(reaction, wavelength, actinic_flux.view(), xs.view())
                        .into_pyresult()?
                }
            } else {
                calculate_photolysis_rate(reaction, wavelength, actinic_flux.view(), xs.view())
                    .into_pyresult()?
            };

            for (k, value) in rate.iter().enumerate() {
                photolysis_rate[[i, k]] = *value;
            }
        }

        let mut state_profiles: HashMap<String, Vec<f64>> = HashMap::new();

        for i in 0..num_alt {
            // Integrate the spectral photolysis rate over the
            // quantum yield

            let o2_density = *o2_density
                .get(i)
                .ok_or(anyhow!("Invalid O2 density size"))
                .into_pyresult()?;
            let o3_density = *o3_density
                .get(i)
                .ok_or(anyhow!("Invalid O3 density size"))
                .into_pyresult()?;
            let n2_density = *n2_density
                .get(i)
                .ok_or(anyhow!("Invalid N2 density size"))
                .into_pyresult()?;
            let co2_density = co2_density
                .as_ref()
                .and_then(|v| v.get(i).copied())
                .unwrap_or(0.0);
            let o3p_density = o3p_density
                .as_ref()
                .and_then(|v| v.get(i).copied())
                .unwrap_or(0.0);
            let temperature = *temperature
                .get(i)
                .ok_or(anyhow!("Invalid temperature size"))
                .into_pyresult()?;

            let densities = HashMap::from([
                ("O2".to_string(), o2_density),
                ("O3".to_string(), o3_density),
                ("N2".to_string(), n2_density),
                ("CO2".to_string(), co2_density),
                ("O(3P)".to_string(), o3p_density),
            ]);

            let photo_slice = photolysis_rate.slice(s![.., i]).to_owned();

            let state = self
                .model
                .solve(
                    temperature,
                    &self.model.chemical_reactions,
                    &self.model.photo_reactions,
                    photo_slice
                        .as_slice()
                        .ok_or(anyhow!("Invalid photolysis rate slice"))
                        .into_pyresult()?,
                    &densities,
                )
                .into_pyresult()?;

            for (name, value) in state {
                state_profiles
                    .entry(name)
                    .or_insert_with(|| vec![f64::NAN; num_alt])[i] = value;
            }
        }

        let py = ds.py();
        let xr = py.import("xarray")?;
        let np = py.import("numpy")?;

        let mut state_names: Vec<String> = state_profiles.keys().cloned().collect();
        state_names.sort();

        let data_vars = pyo3::types::PyDict::new(py);
        for name in state_names {
            let profile = state_profiles
                .get(&name)
                .ok_or(anyhow!("Missing state profile for '{}'", name))
                .into_pyresult()?;
            let values = np.getattr("array")?.call1((profile.clone(),))?;
            data_vars.set_item(name, (vec!["altitude"], values))?;
        }

        let coords = pyo3::types::PyDict::new(py);
        coords.set_item("altitude", ds.get_item("altitude")?)?;

        let kwargs = pyo3::types::PyDict::new(py);
        kwargs.set_item("data_vars", data_vars)?;
        kwargs.set_item("coords", coords)?;

        let out = xr.getattr("Dataset")?.call((), Some(&kwargs))?;

        Ok(out.into())
    }

    pub fn emission_rates(&self, state: &Bound<PyAny>) -> PyResult<Py<PyAny>> {
        let py = state.py();
        let xr = py.import("xarray")?;
        let np = py.import("numpy")?;

        let altitude = state.get_item("altitude")?;

        let o1s = state.get_item("O(1S)")?;
        let o1s = o1s.call_method0("to_numpy")?;
        let o1s = o1s.cast_into::<PyArray1<f64>>()?;
        let o1s = o1s.readonly();
        let o1s = o1s.as_array();

        let o2_b0 = state.get_item("O2(b)")?;
        let o2_b0 = o2_b0.call_method0("to_numpy")?;
        let o2_b0 = o2_b0.cast_into::<PyArray1<f64>>()?;
        let o2_b0 = o2_b0.readonly();
        let o2_b0 = o2_b0.as_array();

        let o2_b1 = match state.get_item("O2(b, v=1)") {
            Ok(o2_b1) => {
                let o2_b1 = o2_b1.call_method0("to_numpy")?;
                let o2_b1 = o2_b1.cast_into::<PyArray1<f64>>()?;
                let o2_b1 = o2_b1.readonly();
                o2_b1.as_array().to_owned()
            }
            Err(_) => ndarray::Array1::zeros(o2_b0.len()),
        };

        let oxygen_green: Vec<f64> = o1s
            .iter()
            .map(|population| population * OXYGEN_GREEN_LINE_EINSTEIN_A_S)
            .collect();
        let oxygen_a_band_b0: Vec<f64> = o2_b0
            .iter()
            .map(|population| population * O2_B0_X0_EINSTEIN_A_S)
            .collect();
        let oxygen_a_band_b1: Vec<f64> = o2_b1
            .iter()
            .map(|population| population * O2_B1_X1_EINSTEIN_A_S)
            .collect();
        let oxygen_a_band: Vec<f64> = oxygen_a_band_b0
            .iter()
            .zip(oxygen_a_band_b1.iter())
            .map(|(b0, b1)| b0 + b1)
            .collect();

        let data_vars = PyDict::new(py);
        data_vars.set_item(
            "oxygen_green_5577_photon_ver",
            (
                vec!["altitude"],
                np.getattr("array")?.call1((oxygen_green,))?,
            ),
        )?;
        data_vars.set_item(
            "oxygen_a_band_photon_ver",
            (
                vec!["altitude"],
                np.getattr("array")?.call1((oxygen_a_band,))?,
            ),
        )?;
        data_vars.set_item(
            "oxygen_a_band_b0_photon_ver",
            (
                vec!["altitude"],
                np.getattr("array")?.call1((oxygen_a_band_b0,))?,
            ),
        )?;
        data_vars.set_item(
            "oxygen_a_band_b1_photon_ver",
            (
                vec!["altitude"],
                np.getattr("array")?.call1((oxygen_a_band_b1,))?,
            ),
        )?;

        let coords = PyDict::new(py);
        coords.set_item("altitude", altitude)?;

        let attrs = PyDict::new(py);
        attrs.set_item("emission_units", "photons m^-3 s^-1")?;
        attrs.set_item(
            "oxygen_green_wavelength_nm",
            OXYGEN_GREEN_LINE_WAVELENGTH_NM,
        )?;
        attrs.set_item(
            "oxygen_a_band_center_wavelength_nm",
            O2_A_BAND_CENTER_WAVELENGTH_NM,
        )?;

        let kwargs = PyDict::new(py);
        kwargs.set_item("data_vars", data_vars)?;
        kwargs.set_item("coords", coords)?;
        kwargs.set_item("attrs", attrs)?;

        let out = xr.getattr("Dataset")?.call((), Some(&kwargs))?;

        Ok(out.into())
    }

    pub fn oxygen_green_line_mcdade(&self, state: &Bound<PyAny>) -> PyResult<Py<PyAny>> {
        let py = state.py();
        let xr = py.import("xarray")?;
        let np = py.import("numpy")?;

        let altitude = state.get_item("altitude")?;

        let temperature = match state.get_item("temperature") {
            Ok(temperature) => temperature,
            Err(_) => state.get_item("temperature_k")?,
        };
        let temperature = temperature.call_method0("to_numpy")?;
        let temperature = temperature.cast_into::<PyArray1<f64>>()?;
        let temperature = temperature.readonly();
        let temperature = temperature.as_array();

        let atomic_o = match state.get_item("o_density") {
            Ok(o_density) => o_density,
            Err(_) => state.get_item("o3p_density")?,
        };
        let atomic_o = atomic_o.call_method0("to_numpy")?;
        let atomic_o = atomic_o.cast_into::<PyArray1<f64>>()?;
        let atomic_o = atomic_o.readonly();
        let atomic_o = atomic_o.as_array();

        let o2 = state.get_item("o2_density")?;
        let o2 = o2.call_method0("to_numpy")?;
        let o2 = o2.cast_into::<PyArray1<f64>>()?;
        let o2 = o2.readonly();
        let o2 = o2.as_array();

        let n2 = state.get_item("n2_density")?;
        let n2 = n2.call_method0("to_numpy")?;
        let n2 = n2.cast_into::<PyArray1<f64>>()?;
        let n2 = n2.readonly();
        let n2 = n2.as_array();

        let photon_ver =
            mcdade_oxygen_green_line_photon_ver(temperature, atomic_o, o2, n2).into_pyresult()?;
        let o1s_population = photon_ver.mapv(|ver| ver / MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_558_S);

        let data_vars = PyDict::new(py);
        data_vars.set_item(
            "oxygen_green_5577_photon_ver",
            (
                vec!["altitude"],
                np.getattr("array")?.call1((photon_ver.to_vec(),))?,
            ),
        )?;
        data_vars.set_item(
            "O(1S)",
            (
                vec!["altitude"],
                np.getattr("array")?.call1((o1s_population.to_vec(),))?,
            ),
        )?;

        let coords = PyDict::new(py);
        coords.set_item("altitude", altitude)?;

        let attrs = PyDict::new(py);
        attrs.set_item("emission_units", "photons m^-3 s^-1")?;
        attrs.set_item(
            "model",
            "McDade/Murtagh short-form Barth oxygen green-line parameterization",
        )?;
        attrs.set_item(
            "oxygen_green_wavelength_nm",
            OXYGEN_GREEN_LINE_WAVELENGTH_NM,
        )?;
        attrs.set_item("mcdade_a558_s", MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_558_S)?;
        attrs.set_item(
            "mcdade_a1s_total_s",
            MCDADE_OXYGEN_GREEN_LINE_EINSTEIN_A_1S_S,
        )?;
        attrs.set_item("mcdade_c0", MCDADE_OXYGEN_GREEN_LINE_C0)?;
        attrs.set_item("mcdade_c1", MCDADE_OXYGEN_GREEN_LINE_C1)?;
        attrs.set_item("mcdade_c2", MCDADE_OXYGEN_GREEN_LINE_C2)?;
        attrs.set_item(
            "density_units",
            "input number densities are m^-3; McDade rate expression is evaluated in cm^-3 internally",
        )?;

        let kwargs = PyDict::new(py);
        kwargs.set_item("data_vars", data_vars)?;
        kwargs.set_item("coords", coords)?;
        kwargs.set_item("attrs", attrs)?;

        let out = xr.getattr("Dataset")?.call((), Some(&kwargs))?;

        Ok(out.into())
    }

    pub fn oxygen_a_band_line_data(
        &self,
        py: Python<'_>,
        hitran_directory: &str,
    ) -> PyResult<Py<PyAny>> {
        let directory = std::path::PathBuf::from(hitran_directory);
        let db = read_hitran_line_file(hitran_molecule_file("O2", &directory).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "Failed to find O2 HITRAN line file: {e}"
            ))
        })?)
        .map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "Failed to read O2 HITRAN line file: {e}"
            ))
        })?;

        let band = EmissionBand::oxygen_a_band_from_hitran(&db)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let xr = py.import("xarray")?;
        let np = py.import("numpy")?;

        let wavelength = band.wavelengths_nm().to_vec();
        let wavenumber: Vec<f64> = band
            .lines
            .iter()
            .map(|line| line.wavenumber_cminv)
            .collect();
        let weights = band.weights().to_vec();
        let line_intensity_296: Vec<f64> = band
            .lines
            .iter()
            .map(|line| line.line_intensity_296)
            .collect();
        let einstein_a: Vec<f64> = band.lines.iter().map(|line| line.einstein_a_s).collect();
        let isotope_id: Vec<i32> = band.lines.iter().map(|line| line.isotope_id).collect();
        let isotope_abundance: Vec<f64> = band
            .lines
            .iter()
            .map(|line| line.isotope_abundance)
            .collect();
        let lower_energy: Vec<f64> = band
            .lines
            .iter()
            .map(|line| line.lower_energy_cminv)
            .collect();
        let upper_energy: Vec<f64> = band
            .lines
            .iter()
            .map(|line| line.upper_energy_cminv)
            .collect();
        let upper_state_id: Vec<String> = band
            .lines
            .iter()
            .map(|line| line.upper_state_id.clone())
            .collect();
        let lower_state_id: Vec<String> = band
            .lines
            .iter()
            .map(|line| line.lower_state_id.clone())
            .collect();
        let upper_vibrational_state: Vec<String> = band
            .lines
            .iter()
            .map(|line| line.upper_vibrational_state.clone())
            .collect();
        let lower_vibrational_state: Vec<String> = band
            .lines
            .iter()
            .map(|line| line.lower_vibrational_state.clone())
            .collect();
        let upper_statistical_weight: Vec<f64> = band
            .lines
            .iter()
            .map(|line| line.upper_statistical_weight.unwrap_or(f64::NAN))
            .collect();
        let lower_statistical_weight: Vec<f64> = band
            .lines
            .iter()
            .map(|line| line.lower_statistical_weight.unwrap_or(f64::NAN))
            .collect();
        let upper_branching_ratio: Vec<f64> = band
            .lines
            .iter()
            .map(|line| line.upper_branching_ratio)
            .collect();

        let data_vars = PyDict::new(py);
        data_vars.set_item(
            "weight",
            (vec!["line"], np.getattr("array")?.call1((weights,))?),
        )?;
        data_vars.set_item(
            "line_intensity_296",
            (
                vec!["line"],
                np.getattr("array")?.call1((line_intensity_296,))?,
            ),
        )?;
        data_vars.set_item(
            "einstein_a_s",
            (vec!["line"], np.getattr("array")?.call1((einstein_a,))?),
        )?;
        data_vars.set_item(
            "isotope_id",
            (vec!["line"], np.getattr("array")?.call1((isotope_id,))?),
        )?;
        data_vars.set_item(
            "isotope_abundance",
            (
                vec!["line"],
                np.getattr("array")?.call1((isotope_abundance,))?,
            ),
        )?;
        data_vars.set_item(
            "wavenumber_cminv",
            (vec!["line"], np.getattr("array")?.call1((wavenumber,))?),
        )?;
        data_vars.set_item(
            "lower_energy_cminv",
            (vec!["line"], np.getattr("array")?.call1((lower_energy,))?),
        )?;
        data_vars.set_item(
            "upper_energy_cminv",
            (vec!["line"], np.getattr("array")?.call1((upper_energy,))?),
        )?;
        data_vars.set_item(
            "upper_state_id",
            (vec!["line"], np.getattr("array")?.call1((upper_state_id,))?),
        )?;
        data_vars.set_item(
            "lower_state_id",
            (vec!["line"], np.getattr("array")?.call1((lower_state_id,))?),
        )?;
        data_vars.set_item(
            "upper_vibrational_state",
            (
                vec!["line"],
                np.getattr("array")?.call1((upper_vibrational_state,))?,
            ),
        )?;
        data_vars.set_item(
            "lower_vibrational_state",
            (
                vec!["line"],
                np.getattr("array")?.call1((lower_vibrational_state,))?,
            ),
        )?;
        data_vars.set_item(
            "upper_statistical_weight",
            (
                vec!["line"],
                np.getattr("array")?.call1((upper_statistical_weight,))?,
            ),
        )?;
        data_vars.set_item(
            "lower_statistical_weight",
            (
                vec!["line"],
                np.getattr("array")?.call1((lower_statistical_weight,))?,
            ),
        )?;
        data_vars.set_item(
            "upper_branching_ratio",
            (
                vec!["line"],
                np.getattr("array")?.call1((upper_branching_ratio,))?,
            ),
        )?;

        let coords = PyDict::new(py);
        coords.set_item(
            "wavelength_nm",
            (vec!["line"], np.getattr("array")?.call1((wavelength,))?),
        )?;

        let attrs = PyDict::new(py);
        attrs.set_item("band_name", band.name)?;
        attrs.set_item("upper_state", band.upper_state)?;
        attrs.set_item("lower_state", band.lower_state)?;
        attrs.set_item("total_einstein_a_s", band.total_einstein_a_s)?;
        attrs.set_item(
            "weight_model",
            "normalized_isotope_abundance_weighted_einstein_a_by_vibrational_band_placeholder",
        )?;

        let kwargs = PyDict::new(py);
        kwargs.set_item("data_vars", data_vars)?;
        kwargs.set_item("coords", coords)?;
        kwargs.set_item("attrs", attrs)?;

        let out = xr.getattr("Dataset")?.call((), Some(&kwargs))?;

        Ok(out.into())
    }
}
