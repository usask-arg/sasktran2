use std::collections::{HashMap, HashSet};

use ndarray::{Array2, s};
use pyo3::prelude::*;
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use sasktran2_rs::photchem::models::*;
use anyhow::anyhow;

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

    pub fn solve(&self,
        ds: &Bound<PyAny>
    ) -> PyResult<PyObject> {

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
        let num_wavelength = wavelength.len();

        if num_wavelength < 2 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "Need at least two wavelength points to integrate photolysis rates",
            ));
        }

        // Central-difference spectral bin widths for J = integral(q * F * sigma d_lambda)
        let mut delta_wavelength = vec![0.0_f64; num_wavelength];
        for j in 0..num_wavelength {
            let d = if j == 0 {
                (wavelength[1] - wavelength[0]).abs()
            } else if j + 1 == num_wavelength {
                (wavelength[num_wavelength - 1] - wavelength[num_wavelength - 2]).abs()
            } else {
                0.5 * (wavelength[j + 1] - wavelength[j - 1]).abs()
            };
            delta_wavelength[j] = d;
        }

        let actinic_flux = ds.get_item("actinic_flux")?;
        let actinic_flux = actinic_flux.call_method0("to_numpy")?;
        let actinic_flux = actinic_flux.cast_into::<PyArray2<f64>>()?;
        let actinic_flux = actinic_flux.readonly();
        let actinic_flux = actinic_flux.as_array();

        // Compute the photolysis rates
        let mut photolysis_rate = Array2::<f64>::zeros((self.model.photo_reactions.len(), num_alt));
        let mut reported_toa_keys = HashSet::<String>::new();



        for (i, reaction) in self.model.photo_reactions.iter().enumerate() {
            let species = reaction.in_molecule.base_type.to_string();

            let xs = ds.get_item(format!("{}_xs", species.to_lowercase()))?;
            let xs = xs.call_method0("to_numpy")?;
            let xs = xs.cast_into::<PyArray2<f64>>()?;
            let xs = xs.readonly();
            let xs = xs.as_array();

            let q = reaction.quantum_yield.unwrap_or(1.0);
            let band_limits = reaction.wavelength_range_nm;

            for (j, wl) in wavelength.iter().enumerate() {
                if let Some((min_nm, max_nm)) = band_limits {
                    if *wl < min_nm || *wl > max_nm {
                        continue;
                    }
                }

                for k in 0..num_alt {
                    let flux = actinic_flux[[j, k]].max(0.0);
                    let cross_section = xs[[j, k]].max(0.0);
                    photolysis_rate[[i, k]] += q * flux * cross_section * delta_wavelength[j];
                }
            }

            if reaction.toa_rate_constant > 0.0 {
                let computed_toa = photolysis_rate[[i, num_alt - 1]];
                let rel_diff = ((computed_toa - reaction.toa_rate_constant) / reaction.toa_rate_constant).abs();
                if rel_diff > 0.2 {
                    let reactant = String::try_from(reaction.in_molecule.clone())
                        .unwrap_or_else(|_| format!("{:?}", reaction.in_molecule));
                    let band = reaction.excitation_band.as_deref().unwrap_or("<all>");
                    let key = format!("{}|{}", reactant, band);
                    if reported_toa_keys.insert(key) {
                        eprintln!(
                            "Photochem TOA check [{} {}]: computed={:.3e}, configured={:.3e}, rel_diff={:.2}",
                            reactant,
                            band,
                            computed_toa,
                            reaction.toa_rate_constant,
                            rel_diff
                        );
                    }
                }
            }
        }


        let mut state_profiles: HashMap<String, Vec<f64>> = HashMap::new();

        for i in 0..num_alt {
            // Integrate the spectral photolysis rate over the 
            // quantum yield

            let o2_density = *o2_density.get(i).ok_or(anyhow!("Invalid O2 density size")).into_pyresult()?;
            let o3_density = *o3_density.get(i).ok_or(anyhow!("Invalid O3 density size")).into_pyresult()?;
            let n2_density = *n2_density.get(i).ok_or(anyhow!("Invalid N2 density size")).into_pyresult()?;
            let co2_density = co2_density
                .as_ref()
                .and_then(|v| v.get(i).copied())
                .unwrap_or(0.0);
            let o3p_density = o3p_density
                .as_ref()
                .and_then(|v| v.get(i).copied())
                .unwrap_or(0.0);
            let temperature = *temperature.get(i).ok_or(anyhow!("Invalid temperature size")).into_pyresult()?;

            let densities = HashMap::from([
                ("O2".to_string(), o2_density),
                ("O3".to_string(), o3_density),
                ("N2".to_string(), n2_density),
                ("CO2".to_string(), co2_density),
                ("O(3P)".to_string(), o3p_density),
            ]);

            let photo_slice = photolysis_rate.slice(s![.., i]).to_owned();

            let state = self.model.solve(
                temperature,
                &self.model.chemical_reactions,
                &self.model.photo_reactions,
                photo_slice.as_slice().ok_or(anyhow!("Invalid photolysis rate slice")).into_pyresult()?,
                &densities,
            ).into_pyresult()?;

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
}