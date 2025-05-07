use crate::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;
use pyo3::{IntoPyObjectExt, prelude::*};
use sasktran2_rs::optical::line::aer_loader::{aer_molecule_file, read_aer_line_file};
use sasktran2_rs::optical::traits::{OpticalProperty, OpticalPropertyExt};
use sasktran2_rs::optical::types::line_absorber;
use sasktran2_rs::optical::types::line_absorber::{MolecularMass, PartitionFactor};

use crate::constituent::atmo_storage::AtmosphereStorage;

use super::optical_quantities::PyOpticalQuantities;
use super::xsec_dbase::{HasDb, PyDictWrapper};

struct PyPartitionFactor {
    pub py_tips: Py<PyAny>,
}

struct PyMolecularMass {
    pub py_molmass: Py<PyAny>,
}

impl PartitionFactor for PyPartitionFactor {
    fn partition_factor(&self, mol_id: i32, iso_id: i32, temperature: f64) -> f64 {
        let result: f64 = Python::with_gil(|py| {
            let tips = self.py_tips.bind(py);
            let args = (mol_id, iso_id, temperature);

            let result = tips.call1(args).unwrap();
            result.extract().unwrap()
        });
        result
    }
}

impl MolecularMass for PyMolecularMass {
    fn molecular_mass(&self, mol_id: i32, iso_id: i32) -> f64 {
        let result: f64 = Python::with_gil(|py| {
            let molmass = self.py_molmass.bind(py);
            let args = (mol_id, iso_id);

            let result = molmass.call1(args).unwrap();
            result.extract().unwrap()
        });
        result
    }
}

#[pyclass(unsendable)]
pub struct PyLineAbsorber {
    pub line_absorber: line_absorber::LineAbsorber,
}

#[pymethods]
impl PyLineAbsorber {
    #[new]
    #[pyo3(signature = (mol_name, directory, cull_factor=0.0, py_tips = None, py_molmass = None))]
    pub fn new<'py>(
        mol_name: &str,
        directory: &str,
        cull_factor: Option<f64>,
        py_tips: Option<Bound<'py, PyAny>>,
        py_molmass: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Self> {
        let directory = std::path::PathBuf::from(directory);

        let file = aer_molecule_file(mol_name, &directory).into_pyresult()?;

        let db = read_aer_line_file(file).into_pyresult()?;

        let mut line_absorber = line_absorber::LineAbsorber::new(db);

        if let Some(py_tips) = py_tips {
            let py_partition_factor = PyPartitionFactor {
                py_tips: py_tips.into(),
            };

            line_absorber = line_absorber.with_partition_generator(Box::new(py_partition_factor));
        }

        if let Some(py_molmass) = py_molmass {
            let py_molecular_mass = PyMolecularMass {
                py_molmass: py_molmass.into(),
            };

            line_absorber =
                line_absorber.with_molecular_mass_generator(Box::new(py_molecular_mass));
        }

        if let Some(cull_factor) = cull_factor {
            line_absorber = line_absorber.with_cull_factor(cull_factor);
        }

        Ok(Self {
            line_absorber: line_absorber,
        })
    }

    #[pyo3(signature = (atmo, **kwargs))]
    fn atmosphere_quantities<'py>(
        &self,
        atmo: Bound<'py, PyAny>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let rust_atmo = AtmosphereStorage::new(&atmo);
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .line_absorber
            .optical_quantities(&rust_atmo.inputs, &aux_inputs)
            .map_err(|e| {
                PyValueError::new_err(format!("Failed to get optical quantities: {}", e))
            })?;

        Ok(PyOpticalQuantities::new(oq).into_bound_py_any(atmo.py())?)
    }
}

impl HasDb for PyLineAbsorber {
    fn db(&self) -> &dyn OpticalProperty {
        &self.line_absorber
    }
}
