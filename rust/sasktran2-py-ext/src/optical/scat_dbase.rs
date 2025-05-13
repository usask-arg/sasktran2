use numpy::ndarray::*;
use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;
use pyo3::{IntoPyObjectExt, prelude::*};
use sasktran2_rs::interpolation::grid1d::Grid1D;
use sasktran2_rs::optical::traits::*;
use sasktran2_rs::optical::types::scat_dbase::ScatteringDatabase;

use crate::constituent::atmo_storage::AtmosphereStorage;

use super::optical_quantities::PyOpticalQuantities;
use crate::optical::xsec_dbase::PyDictWrapper;

#[pyclass]
pub struct PyScatteringDatabaseDim2 {
    pub db: ScatteringDatabase<Ix2, Ix3>,
}

#[pymethods]
impl PyScatteringDatabaseDim2 {
    #[new]
    pub fn new(
        xsec: PyReadonlyArray2<f64>,
        ssa: PyReadonlyArray2<f64>,
        legendre: PyReadonlyArray3<f64>,
        wvnum: PyReadonlyArray1<f64>,
        params: PyReadonlyArray1<f64>,
        param_names: Vec<String>,
    ) -> Self {
        let xsec = xsec.as_array();
        let ssa = ssa.as_array();
        let legendre = legendre.as_array();
        let wvnum = wvnum.as_array();

        let wvnum_grid = Grid1D::new(wvnum.to_owned());

        Self {
            db: ScatteringDatabase::new(
                xsec.to_owned(),
                ssa.to_owned(),
                legendre.to_owned(),
                wvnum_grid,
                vec![params.as_array().to_owned()],
                param_names,
            ),
        }
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
            .db
            .optical_quantities(&rust_atmo.inputs, &aux_inputs)
            .map_err(|e| {
                PyValueError::new_err(format!("Failed to get optical quantities: {}", e))
            })?;

        PyOpticalQuantities::new(oq).into_bound_py_any(atmo.py())
    }
}
