use std::collections::HashMap;
use std::path::PathBuf;
use std::str::FromStr;

use ndarray::Array1;
use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use sasktran2_rs::constituent::traits::Constituent;
use sasktran2_rs::constituent::types::population_emission_rate::{
    PopulationEmissionProfiles, PopulationEmissionRate, PopulationEmissionSpecies,
};
use sasktran2_rs::optical::line::hitran_loader::{hitran_molecule_file, read_hitran_line_file};
use sasktran2_rs::photchem::emission::AEmissionLineWeightModel;

use crate::constituent::atmo_storage::AtmosphereStorage;

#[pyclass]
pub struct PyPopulationEmissionRate {
    pub inner: PopulationEmissionRate,
}

#[pymethods]
impl PyPopulationEmissionRate {
    #[new]
    #[pyo3(
        signature = (populations, hitran_directory, species = None, line_weight_model = "einstein_a_branching", out_of_bounds_mode = "zero"),
    )]
    fn new(
        populations: Bound<'_, PyAny>,
        hitran_directory: &str,
        species: Option<Vec<String>>,
        line_weight_model: &str,
        out_of_bounds_mode: Option<&str>,
    ) -> PyResult<Self> {
        let species = species.unwrap_or_else(|| vec!["O2".to_string()]);
        let species = species
            .iter()
            .map(|species| PopulationEmissionSpecies::from_str(species))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let line_weight_model = AEmissionLineWeightModel::from_str(line_weight_model)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let altitude_m = required_dataset_array1(&populations, "altitude")?;
        let temperature_k = required_dataset_array1(&populations, "temperature")?;
        let mut population_map = HashMap::new();
        population_map.insert(
            "O2(b)".to_string(),
            required_dataset_array1(&populations, "O2(b)")?,
        );
        for optional_name in ["O2(b, v=1)", "O2(b, v=2)", "O(1S)"] {
            if let Some(population) = optional_dataset_array1(&populations, optional_name)? {
                population_map.insert(optional_name.to_string(), population);
            }
        }

        let profiles = PopulationEmissionProfiles::new(altitude_m, temperature_k, population_map)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let directory = PathBuf::from(hitran_directory);
        let o2_file = hitran_molecule_file("O2", &directory)
            .map_err(|e| PyValueError::new_err(format!("Failed to find O2 HITRAN file: {e}")))?;
        let o2_db = read_hitran_line_file(o2_file)
            .map_err(|e| PyValueError::new_err(format!("Failed to read O2 HITRAN file: {e}")))?;

        let mut inner =
            PopulationEmissionRate::new(profiles, &species, line_weight_model, Some(&o2_db))
                .map_err(|e| PyValueError::new_err(e.to_string()))?;

        if let Some(out_of_bounds_mode) = out_of_bounds_mode {
            inner = match out_of_bounds_mode {
                "zero" => {
                    inner.with_interp_mode(sasktran2_rs::interpolation::OutOfBoundsMode::Zero)
                }
                "extend" => {
                    inner.with_interp_mode(sasktran2_rs::interpolation::OutOfBoundsMode::Extend)
                }
                mode => {
                    return Err(PyValueError::new_err(format!(
                        "Invalid out_of_bounds_mode: {mode}"
                    )));
                }
            };
        }

        Ok(Self { inner })
    }

    #[getter]
    fn get_photon_ver<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        ensure_has_primary_line_list(&this)?;
        let array = &this.borrow().inner.line_list_emissions[0].photon_ver;

        Ok(unsafe { PyArray1::borrow_from_array(array, this.into_any()) })
    }

    #[getter]
    fn get_altitudes_m<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        ensure_has_primary_line_list(&this)?;
        let array = &this.borrow().inner.line_list_emissions[0].altitudes;

        Ok(unsafe { PyArray1::borrow_from_array(array, this.into_any()) })
    }

    #[getter]
    fn get_wavelengths_nm<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        ensure_has_primary_line_list(&this)?;
        let array = &this.borrow().inner.line_list_emissions[0].wavelengths_nm;

        Ok(unsafe { PyArray1::borrow_from_array(array, this.into_any()) })
    }

    #[getter]
    fn get_weights<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        ensure_has_primary_line_list(&this)?;
        let array = &this.borrow().inner.line_list_emissions[0].weights;

        Ok(unsafe { PyArray2::borrow_from_array(array, this.into_any()) })
    }

    #[getter]
    fn get_num_line_list_emissions(&self) -> usize {
        self.inner.line_list_emissions.len()
    }

    #[pyo3(signature = (index = 0))]
    fn line_list_photon_ver<'py>(
        this: Bound<'py, Self>,
        index: usize,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        ensure_line_list_index(&this, index)?;
        let array = &this.borrow().inner.line_list_emissions[index].photon_ver;

        Ok(unsafe { PyArray1::borrow_from_array(array, this.into_any()) })
    }

    #[pyo3(signature = (index = 0))]
    fn line_list_wavelengths_nm<'py>(
        this: Bound<'py, Self>,
        index: usize,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        ensure_line_list_index(&this, index)?;
        let array = &this.borrow().inner.line_list_emissions[index].wavelengths_nm;

        Ok(unsafe { PyArray1::borrow_from_array(array, this.into_any()) })
    }

    #[pyo3(signature = (index = 0))]
    fn line_list_weights<'py>(
        this: Bound<'py, Self>,
        index: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        ensure_line_list_index(&this, index)?;
        let array = &this.borrow().inner.line_list_emissions[index].weights;

        Ok(unsafe { PyArray2::borrow_from_array(array, this.into_any()) })
    }

    pub fn add_to_atmosphere<'py>(&mut self, atmo: Bound<'py, PyAny>) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;

        self.inner
            .add_to_atmosphere(&mut rust_atmo)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(())
    }

    pub fn register_derivative(&mut self, atmo: Bound<'_, PyAny>, name: &str) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;

        self.inner
            .register_derivatives(&mut rust_atmo, name)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(())
    }
}

fn required_dataset_array1(dataset: &Bound<'_, PyAny>, name: &str) -> PyResult<Array1<f64>> {
    let item = dataset.get_item(name).map_err(|_| {
        PyValueError::new_err(format!("Required dataset variable '{name}' is missing"))
    })?;
    dataset_item_to_array1(item, name)
}

fn optional_dataset_array1(
    dataset: &Bound<'_, PyAny>,
    name: &str,
) -> PyResult<Option<Array1<f64>>> {
    match dataset.get_item(name) {
        Ok(item) => dataset_item_to_array1(item, name).map(Some),
        Err(_) => Ok(None),
    }
}

fn dataset_item_to_array1(item: Bound<'_, PyAny>, name: &str) -> PyResult<Array1<f64>> {
    let array = item.call_method0("to_numpy")?;
    let array = array.cast_into::<PyArray1<f64>>().map_err(|_| {
        PyValueError::new_err(format!(
            "Dataset variable '{name}' must be one-dimensional float64"
        ))
    })?;
    let array = array.readonly();

    Ok(array.as_array().to_owned())
}

fn ensure_has_primary_line_list(this: &Bound<'_, PyPopulationEmissionRate>) -> PyResult<()> {
    if this.borrow().inner.line_list_emissions.is_empty() {
        return Err(PyValueError::new_err(
            "PopulationEmissionRate does not contain a line-list emission component",
        ));
    }

    Ok(())
}

fn ensure_line_list_index(
    this: &Bound<'_, PyPopulationEmissionRate>,
    index: usize,
) -> PyResult<()> {
    let count = this.borrow().inner.line_list_emissions.len();
    if index >= count {
        return Err(PyValueError::new_err(format!(
            "Line-list emission index {index} is out of range for {count} components"
        )));
    }

    Ok(())
}
