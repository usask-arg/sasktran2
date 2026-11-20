use std::collections::{BTreeMap, HashMap};

use numpy::{IntoPyArray, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use sasktran2_nlte::mechanism::{Mechanism, ProcessKind, Species};
use sasktran2_nlte::solver::{Column, solve_steady_state};

fn value_error(err: anyhow::Error) -> PyErr {
    PyValueError::new_err(err.to_string())
}

/// A validated kinetic mechanism from the `sasktran2-nlte` crate.
#[pyclass(frozen)]
pub struct PyMechanism {
    mechanism: Mechanism,
}

impl PyMechanism {
    fn species_id(&self, species: &Species) -> String {
        match species {
            Species::State(index) => self.mechanism.states()[*index].id.clone(),
            Species::Background(index) => self.mechanism.background()[*index].clone(),
        }
    }
}

#[pymethods]
impl PyMechanism {
    #[staticmethod]
    fn from_toml(text: &str) -> PyResult<Self> {
        Ok(Self {
            mechanism: Mechanism::from_toml_str(text).map_err(value_error)?,
        })
    }

    #[staticmethod]
    fn bundled(name: &str) -> PyResult<Self> {
        Ok(Self {
            mechanism: Mechanism::bundled(name).map_err(value_error)?,
        })
    }

    #[staticmethod]
    fn bundled_names() -> Vec<String> {
        Mechanism::bundled_names()
            .into_iter()
            .map(str::to_string)
            .collect()
    }

    #[getter]
    fn name(&self) -> String {
        self.mechanism.name.clone()
    }

    #[getter]
    fn version(&self) -> String {
        self.mechanism.version.clone()
    }

    #[getter]
    fn description(&self) -> String {
        self.mechanism.description.clone()
    }

    #[getter]
    fn states(&self) -> Vec<String> {
        self.mechanism
            .states()
            .iter()
            .map(|state| state.id.clone())
            .collect()
    }

    #[getter]
    fn background(&self) -> Vec<String> {
        self.mechanism.background().to_vec()
    }

    #[getter]
    fn rate_inputs(&self) -> Vec<String> {
        self.mechanism.rate_inputs().to_vec()
    }

    #[getter]
    fn references(&self) -> BTreeMap<String, String> {
        self.mechanism.references().clone()
    }

    #[getter]
    fn process_ids(&self) -> Vec<String> {
        self.mechanism
            .processes()
            .iter()
            .map(|p| p.id.clone())
            .collect()
    }

    #[getter]
    fn process_kinds(&self) -> Vec<&'static str> {
        self.mechanism
            .processes()
            .iter()
            .map(|p| p.kind.as_str())
            .collect()
    }

    #[getter]
    fn process_references(&self) -> Vec<String> {
        self.mechanism
            .processes()
            .iter()
            .map(|p| p.reference.clone())
            .collect()
    }

    /// `(process index, upper state, lower species, wavelength_nm or None)`
    /// for every radiative process.
    #[getter]
    fn transitions(&self) -> Vec<(usize, String, String, Option<f64>)> {
        self.mechanism
            .processes()
            .iter()
            .enumerate()
            .filter(|(_, p)| p.kind == ProcessKind::Radiative)
            .map(|(index, p)| {
                (
                    index,
                    self.species_id(&p.reactants[0]),
                    self.species_id(&p.channels[0].products[0]),
                    p.wavelength_nm,
                )
            })
            .collect()
    }

    /// Sparse net stoichiometry `(process, state, coefficient)`: the change in
    /// a state's population per event of a process. Losses are negative.
    fn state_stoichiometry(&self) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        let mut net: BTreeMap<(usize, usize), f64> = BTreeMap::new();
        for (p, process) in self.mechanism.processes().iter().enumerate() {
            if let Some(source) = process.state_reactant() {
                *net.entry((p, source)).or_default() -= 1.0;
            }
            for channel in &process.channels {
                for product in &channel.products {
                    if let Species::State(state) = product {
                        *net.entry((p, *state)).or_default() += channel.fraction;
                    }
                }
            }
        }
        let mut processes = Vec::with_capacity(net.len());
        let mut states = Vec::with_capacity(net.len());
        let mut coefficients = Vec::with_capacity(net.len());
        for ((p, s), c) in net {
            processes.push(p);
            states.push(s);
            coefficients.push(c);
        }
        (processes, states, coefficients)
    }

    #[pyo3(signature = (temperature_k, densities_m3, rates_per_s, pressure_pa=None))]
    fn solve_steady_state<'py>(
        &self,
        py: Python<'py>,
        temperature_k: PyReadonlyArray1<'py, f64>,
        densities_m3: HashMap<String, PyReadonlyArray1<'py, f64>>,
        rates_per_s: HashMap<String, PyReadonlyArray1<'py, f64>>,
        pressure_pa: Option<PyReadonlyArray1<'py, f64>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let column = Column {
            temperature_k: temperature_k.as_array(),
            pressure_pa: pressure_pa.as_ref().map(|p| p.as_array()),
            densities_m3: densities_m3
                .iter()
                .map(|(name, values)| (name.clone(), values.as_array()))
                .collect(),
            rates_per_s: rates_per_s
                .iter()
                .map(|(name, values)| (name.clone(), values.as_array()))
                .collect(),
        };
        let solution = solve_steady_state(&self.mechanism, &column).map_err(value_error)?;

        let out = PyDict::new(py);
        out.set_item(
            "state_density_m3",
            solution.state_density_m3.into_pyarray(py),
        )?;
        out.set_item(
            "process_rate_m3_s",
            solution.process_rate_m3_s.into_pyarray(py),
        )?;
        out.set_item("production_m3_s", solution.production_m3_s.into_pyarray(py))?;
        out.set_item("loss_m3_s", solution.loss_m3_s.into_pyarray(py))?;
        out.set_item(
            "relative_residual",
            solution.relative_residual.into_pyarray(py),
        )?;
        Ok(out)
    }
}
