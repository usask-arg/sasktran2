use crate::prelude::*;
use numpy::{PyArray2, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;
use pyo3::{IntoPyObjectExt, prelude::*};
use sasktran2_rs::optical::line::aer_loader::{aer_molecule_file, read_aer_line_file};
use sasktran2_rs::optical::line::hitran_loader::{hitran_molecule_file, read_hitran_line_file};
use sasktran2_rs::optical::traits::{OpticalProperty, OpticalPropertyExt};
use sasktran2_rs::optical::types::line_absorber;
use sasktran2_rs::optical::types::line_absorber::{MolecularMass, PartitionFactor};

use crate::constituent::atmo_storage::AtmosphereStorage;

use super::optical_quantities::PyOpticalQuantities;
use super::xsec_dbase::{HasDb, PyDictWrapper};

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
#[allow(clippy::upper_case_acronyms)]
pub enum LineDatabaseType {
    HITRAN,
    AER,
}

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

            let result = tips.call1(args).unwrap_or_else(|_| {
                panic!("Failed to call partition factor function with args: {args:?}")
            });
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

            let result = molmass.call1(args).unwrap_or_else(|_| {
                panic!("Failed to call molecular mass function with args: {args:?}")
            });
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
    #[pyo3(signature = (db_type, mol_name, directory, cull_factor=0.0, line_coupling=false, py_tips = None, py_molmass = None, line_contribution_width=25.0))]
    #[allow(clippy::too_many_arguments)]
    pub fn new<'py>(
        db_type: PyRef<'_, LineDatabaseType>,
        mol_name: &str,
        directory: &str,
        cull_factor: f64,
        line_coupling: bool,
        py_tips: Option<Bound<'py, PyAny>>,
        py_molmass: Option<Bound<'py, PyAny>>,
        line_contribution_width: f64,
    ) -> PyResult<Self> {
        let directory = std::path::PathBuf::from(directory);

        let db = match *db_type {
            LineDatabaseType::AER => {
                read_aer_line_file(aer_molecule_file(mol_name, &directory).into_pyresult()?)
                    .into_pyresult()?
            }
            LineDatabaseType::HITRAN => {
                read_hitran_line_file(hitran_molecule_file(mol_name, &directory).into_pyresult()?)
                    .into_pyresult()?
            }
        };

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
        line_absorber = line_absorber.with_cull_factor(cull_factor);
        line_absorber = line_absorber.with_line_coupling(line_coupling);
        line_absorber = line_absorber.with_line_contribution_width(line_contribution_width);

        Ok(Self { line_absorber })
    }

    #[pyo3(signature = (atmo, **kwargs))]
    fn atmosphere_quantities<'py>(
        &self,
        atmo: Bound<'py, PyAny>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let rust_atmo = AtmosphereStorage::new(&atmo)?;
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .line_absorber
            .optical_quantities(&rust_atmo.inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        PyOpticalQuantities::new(oq).into_bound_py_any(atmo.py())
    }

    fn cross_section<'py>(
        &self,
        wavenumber_cminv: PyReadonlyArray1<'py, f64>,
        temperature_k: PyReadonlyArray1<f64>,
        pressure_pa: PyReadonlyArray1<f64>,
        p_self: PyReadonlyArray1<f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let xs = self
            .line_absorber
            .cross_section(
                wavenumber_cminv.as_array().view(),
                temperature_k.as_array().view(),
                pressure_pa.as_array().view(),
                p_self.as_array().view(),
            )
            .into_pyresult()?;

        Ok(PyArray2::from_array(wavenumber_cminv.py(), &xs))
    }
}

impl HasDb for PyLineAbsorber {
    fn db(&self) -> &dyn OpticalProperty {
        &self.line_absorber
    }
}
