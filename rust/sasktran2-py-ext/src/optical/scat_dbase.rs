use numpy::ndarray::*;
use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;
use pyo3::{IntoPyObjectExt, prelude::*};
use sasktran2_rs::atmosphere::types::ManualStorageInputs;
use sasktran2_rs::interpolation::grid1d::Grid1D;
use sasktran2_rs::optical::traits::*;
use sasktran2_rs::optical::types::scat_dbase::ScatteringDatabase;

use crate::constituent::atmo_storage::AtmosphereStorage;

use super::optical_quantities::PyOpticalQuantities;
use crate::optical::xsec_dbase::PyDictWrapper;

#[pyclass]
pub struct PyScatteringDatabaseDim1 {
    pub db: ScatteringDatabase<Ix1, Ix2>,
}

#[pymethods]
impl PyScatteringDatabaseDim1 {
    #[new]
    pub fn new(
        xsec: PyReadonlyArray1<f64>,
        ssa: PyReadonlyArray1<f64>,
        legendre: PyReadonlyArray2<f64>,
        wvnum: PyReadonlyArray1<f64>,
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
                Vec::new(),
                Vec::new(),
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
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        PyOpticalQuantities::new(oq).into_bound_py_any(atmo.py())
    }

    #[pyo3(signature = (atmo, **kwargs))]
    #[allow(unused_variables)]
    fn optical_derivatives<'py>(
        &self,
        atmo: Bound<'py, PyAny>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        // This has no derivatives, just return an empty dictionary
        let py_dict = PyDict::new(atmo.py());

        Ok(py_dict)
    }

    #[pyo3(signature = (wavelengths_nm, altitudes_m, **kwargs))]
    #[allow(unused_variables)]
    fn cross_sections<'py>(
        &self,
        py: Python<'py>,
        wavelengths_nm: PyReadonlyArray1<f64>,
        altitudes_m: PyReadonlyArray1<f64>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let inputs = ManualStorageInputs::new()
            .with_singlescatter_moments(0)
            .with_altitude_m(altitudes_m.as_array().to_owned())
            .with_num_stokes(1)
            .with_wavelengths_nm(wavelengths_nm.as_array().to_owned());
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .db
            .optical_quantities(&inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        PyOpticalQuantities::new(oq).into_bound_py_any(py)
    }

    #[pyo3(signature = (wavelengths_nm, altitudes_m, **kwargs))]
    #[allow(unused_variables)]
    fn cross_section_derivatives<'py>(
        &self,
        py: Python<'py>,
        wavelengths_nm: PyReadonlyArray1<f64>,
        altitudes_m: PyReadonlyArray1<f64>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        // This has no derivatives, just return an empty dictionary
        let py_dict = PyDict::new(py);

        Ok(py_dict)
    }
}

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
        let params = params.as_array().to_owned();

        let wvnum_grid = Grid1D::new(wvnum.to_owned());

        Self {
            db: ScatteringDatabase::new(
                xsec.to_owned(),
                ssa.to_owned(),
                legendre.to_owned(),
                wvnum_grid,
                vec![params],
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
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        PyOpticalQuantities::new(oq).into_bound_py_any(atmo.py())
    }

    #[pyo3(signature = (atmo, **kwargs))]
    #[allow(unused_variables)]
    fn optical_derivatives<'py>(
        &self,
        atmo: Bound<'py, PyAny>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let rust_atmo = AtmosphereStorage::new(&atmo)?;
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .db
            .optical_derivatives(&rust_atmo.inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        let py = atmo.py();
        let py_dict = PyDict::new(py);

        for (key, oq) in oq {
            let py_oq = PyOpticalQuantities::new(oq).into_bound_py_any(py)?;
            py_dict.set_item(key.as_str(), py_oq)?;
        }

        Ok(py_dict)
    }

    #[pyo3(signature = (wavelengths_nm, altitudes_m, **kwargs))]
    #[allow(unused_variables)]
    fn cross_sections<'py>(
        &self,
        py: Python<'py>,
        wavelengths_nm: PyReadonlyArray1<f64>,
        altitudes_m: PyReadonlyArray1<f64>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let inputs = ManualStorageInputs::new()
            .with_singlescatter_moments(0)
            .with_altitude_m(altitudes_m.as_array().to_owned())
            .with_num_stokes(1)
            .with_wavelengths_nm(wavelengths_nm.as_array().to_owned());
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .db
            .optical_quantities(&inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        PyOpticalQuantities::new(oq).into_bound_py_any(py)
    }

    #[pyo3(signature = (wavelengths_nm, altitudes_m, **kwargs))]
    #[allow(unused_variables)]
    fn cross_section_derivatives<'py>(
        &self,
        py: Python<'py>,
        wavelengths_nm: PyReadonlyArray1<f64>,
        altitudes_m: PyReadonlyArray1<f64>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let inputs = ManualStorageInputs::new()
            .with_singlescatter_moments(0)
            .with_altitude_m(altitudes_m.as_array().to_owned())
            .with_num_stokes(1)
            .with_wavelengths_nm(wavelengths_nm.as_array().to_owned());
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .db
            .optical_derivatives(&inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        let py_dict = PyDict::new(py);

        for (key, oq) in oq {
            let py_oq = PyOpticalQuantities::new(oq).into_bound_py_any(py)?;
            py_dict.set_item(key.as_str(), py_oq)?;
        }

        Ok(py_dict)
    }
}

#[pyclass]
pub struct PyScatteringDatabaseDim3 {
    pub db: ScatteringDatabase<Ix3, Ix4>,
}

#[pymethods]
impl PyScatteringDatabaseDim3 {
    #[new]
    pub fn new(
        xsec: PyReadonlyArray3<f64>,
        ssa: PyReadonlyArray3<f64>,
        legendre: PyReadonlyArray4<f64>,
        wvnum: PyReadonlyArray1<f64>,
        param0: PyReadonlyArray1<f64>,
        param1: PyReadonlyArray1<f64>,
        param_names: Vec<String>,
    ) -> Self {
        let xsec = xsec.as_array();
        let ssa = ssa.as_array();
        let legendre = legendre.as_array();
        let wvnum = wvnum.as_array();
        let param0 = param0.as_array().to_owned();
        let param1 = param1.as_array().to_owned();

        let wvnum_grid = Grid1D::new(wvnum.to_owned());

        Self {
            db: ScatteringDatabase::new(
                xsec.to_owned(),
                ssa.to_owned(),
                legendre.to_owned(),
                wvnum_grid,
                vec![param0, param1],
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
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        PyOpticalQuantities::new(oq).into_bound_py_any(atmo.py())
    }

    #[pyo3(signature = (atmo, **kwargs))]
    #[allow(unused_variables)]
    fn optical_derivatives<'py>(
        &self,
        atmo: Bound<'py, PyAny>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let rust_atmo = AtmosphereStorage::new(&atmo)?;
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .db
            .optical_derivatives(&rust_atmo.inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        let py = atmo.py();
        let py_dict = PyDict::new(py);

        for (key, oq) in oq {
            let py_oq = PyOpticalQuantities::new(oq).into_bound_py_any(py)?;
            py_dict.set_item(key.as_str(), py_oq)?;
        }

        Ok(py_dict)
    }

    #[pyo3(signature = (wavelengths_nm, altitudes_m, **kwargs))]
    #[allow(unused_variables)]
    fn cross_sections<'py>(
        &self,
        py: Python<'py>,
        wavelengths_nm: PyReadonlyArray1<f64>,
        altitudes_m: PyReadonlyArray1<f64>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let inputs = ManualStorageInputs::new()
            .with_singlescatter_moments(0)
            .with_altitude_m(altitudes_m.as_array().to_owned())
            .with_num_stokes(1)
            .with_wavelengths_nm(wavelengths_nm.as_array().to_owned());
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .db
            .optical_quantities(&inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        PyOpticalQuantities::new(oq).into_bound_py_any(py)
    }

    #[pyo3(signature = (wavelengths_nm, altitudes_m, **kwargs))]
    #[allow(unused_variables)]
    fn cross_section_derivatives<'py>(
        &self,
        py: Python<'py>,
        wavelengths_nm: PyReadonlyArray1<f64>,
        altitudes_m: PyReadonlyArray1<f64>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let inputs = ManualStorageInputs::new()
            .with_singlescatter_moments(0)
            .with_altitude_m(altitudes_m.as_array().to_owned())
            .with_num_stokes(1)
            .with_wavelengths_nm(wavelengths_nm.as_array().to_owned());
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .db
            .optical_derivatives(&inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        let py_dict = PyDict::new(py);

        for (key, oq) in oq {
            let py_oq = PyOpticalQuantities::new(oq).into_bound_py_any(py)?;
            py_dict.set_item(key.as_str(), py_oq)?;
        }

        Ok(py_dict)
    }
}
