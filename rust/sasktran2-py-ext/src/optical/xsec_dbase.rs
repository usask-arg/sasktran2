use crate::constituent::atmo_storage::AtmosphereStorage;
use crate::optical::optical_quantities::PyOpticalQuantities;
use ndarray::*;
use numpy::*;
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::PyValueError;
use pyo3::{prelude::*, types::PyDict};
use sasktran2_rs::interpolation::grid1d::Grid1D;
use sasktran2_rs::optical::traits::*;
use sasktran2_rs::optical::types::xsec_dbase::*;

// Newtype wrapper around Bound<'_, PyDict>
pub struct PyDictWrapper<'a>(pub Option<&'a Bound<'a, PyDict>>);

impl<'a> AuxOpticalInputs for PyDictWrapper<'a> {
    fn get_parameter(&self, name: &str) -> Option<CowArray<'_, f64, Ix1>> {
        if let Some(dict) = &self.0
            && let Ok(value) = dict.get_item(name)
        {
            value.as_ref()?;
            let array: PyReadonlyArray1<f64> = value?.extract().ok()?;
            // Create an owned copy of the array data to avoid lifetime issues
            return Some(CowArray::from(array.as_array().to_owned()));
        }
        None
    }
}

pub trait HasDb {
    fn db(&self) -> &dyn OpticalProperty;
}

#[pyclass]
/// An absorber database that depends on wavenumber only
///
/// Parameters
/// ----------
/// wavenumber_cminv : ndarray
///   Wavenumber in cm^-1
/// xsec_m2 : ndarray
///   Cross section in m^2/molecule [num_wavenumber]
pub struct AbsorberDatabaseDim1 {
    pub db: SKXsecDatabase<Ix1>,
}

#[pyclass]
/// An absorber database that depends on wavenumber and 1 parameter (e.g. temperature)
///
/// Parameters
/// ----------
/// wavenumber_cminv : ndarray
///   Wavenumber in cm^-1
/// params : ndarray
///   Parameters that the cross section depends on
/// xsec_m2 : ndarray
///   Cross section in m^2/molecule [num_param, num_wavenumber]
pub struct AbsorberDatabaseDim2 {
    pub db: SKXsecDatabase<Ix2>,
}

#[pyclass]
/// An absorber database that depends on wavenumber and 2 parameters (e.g. temperature and pressure)
///
/// Parameters
/// ----------
/// wavenumber_cminv : ndarray
///   Wavenumber in cm^-1, shape [W]
/// param_0: ndarray
///   First parameter that the cross section depends on, e.g. Pressure, must be sorted in ascending order, shape [N]
/// param_1: ndarray
///   Second parameter that the cross section depends on, e.g. Temperature, must be sorted in ascending order, shape [M]
/// xsec_m2 : ndarray
///   Cross section in m^2/molecule [N, M, W]
pub struct AbsorberDatabaseDim3 {
    pub db: SKXsecDatabase<Ix3>,
}

#[pymethods]
impl AbsorberDatabaseDim1 {
    #[new]
    fn new(
        wavenumber_cminv: PyReadonlyArray1<f64>,
        xsec_m2: PyReadonlyArray1<f64>,
    ) -> PyResult<Self> {
        let arr_wvnum = wavenumber_cminv.as_array().to_owned();

        let arr_xsec = xsec_m2.as_array().to_owned();

        let db = SKXsecDatabase::<Ix1>::new(arr_xsec, Grid1D::new(arr_wvnum), vec![], vec![])
            .ok_or_else(|| PyValueError::new_err("Failed to create SKXsecDatabase"))?;

        Ok(Self { db })
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
}

#[pymethods]
impl AbsorberDatabaseDim2 {
    #[new]
    fn new(
        wavenumber_cminv: PyReadonlyArray1<f64>,
        params: PyReadonlyArray1<f64>,
        xsec_m2: PyReadonlyArray2<f64>,
        param_names: Vec<String>,
    ) -> PyResult<Self> {
        let arr_wvnum = wavenumber_cminv.as_array().to_owned();
        let arr_params = params.as_array().to_owned();

        let arr_xsec = xsec_m2.as_array().to_owned();

        let db = SKXsecDatabase::<Ix2>::new(
            arr_xsec,
            Grid1D::new(arr_wvnum),
            vec![arr_params],
            param_names,
        )
        .ok_or_else(|| PyValueError::new_err("Failed to create SKXsecDatabase"))?;

        Ok(Self { db })
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

    fn cross_section<'py>(
        &self,
        wavenumber_cminv: PyReadonlyArray1<'py, f64>,
        params: PyReadonlyArray1<'py, f64>,
    ) -> (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray3<f64>>) {
        let arr_wvnum = wavenumber_cminv.as_array();
        let arr_params = params.as_array();

        let arr_xsec = PyArray2::zeros(
            wavenumber_cminv.py(),
            (arr_params.len(), arr_wvnum.len()),
            false,
        );
        let arr_d_xsec = PyArray3::zeros(
            wavenumber_cminv.py(),
            (1, arr_params.len(), arr_wvnum.len()),
            false,
        );

        unsafe {
            let mut mut_view = arr_xsec.as_array_mut();
            let mut mut_d_view = arr_d_xsec.as_array_mut();

            let thread_pool = sasktran2_rs::threading::thread_pool().unwrap();

            thread_pool.install(|| {
                Zip::from(mut_view.rows_mut())
                    .and(mut_d_view.axis_iter_mut(Axis(1)))
                    .and(arr_params)
                    .par_for_each(|mut row, d_row, param| {
                        let params = Array1::from(vec![*param]);

                        let _ = self
                            .db
                            .xs_emplace(&arr_wvnum, &params, &mut row, Some(d_row));
                    });
            });
        }

        (arr_xsec, arr_d_xsec)
    }
}

#[pymethods]
impl AbsorberDatabaseDim3 {
    #[new]
    fn new(
        wavenumber_cminv: PyReadonlyArray1<f64>,
        param_0: PyReadonlyArray1<f64>,
        param_1: PyReadonlyArray1<f64>,
        xsec_m2: PyReadonlyArray3<f64>,
        param_names: Vec<String>,
    ) -> PyResult<Self> {
        let arr_wvnum = wavenumber_cminv.as_array().to_owned();

        let arr_param_0 = param_0.as_array().to_owned();
        let arr_param_1 = param_1.as_array().to_owned();

        let arr_xsec = xsec_m2.as_array().to_owned();

        let db = SKXsecDatabase::<Ix3>::new(
            arr_xsec,
            Grid1D::new(arr_wvnum),
            vec![arr_param_0, arr_param_1],
            param_names,
        )
        .ok_or_else(|| PyValueError::new_err("Failed to create SKXsecDatabase"))?;

        Ok(Self { db })
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

    fn cross_section<'py>(
        &self,
        wavenumber_cminv: PyReadonlyArray1<'py, f64>,
        params: PyReadonlyArray2<'py, f64>,
    ) -> (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray3<f64>>) {
        let arr_wvnum = wavenumber_cminv.as_array();
        let arr_params = params.as_array();

        let ngeo = arr_params.shape()[1];

        let arr_xsec = PyArray2::zeros(wavenumber_cminv.py(), (ngeo, arr_wvnum.len()), false);
        let arr_d_xsec = PyArray3::zeros(wavenumber_cminv.py(), (2, ngeo, arr_wvnum.len()), false);

        unsafe {
            let mut mut_view = arr_xsec.as_array_mut();
            let mut mut_d_view = arr_d_xsec.as_array_mut();

            let thread_pool = sasktran2_rs::threading::thread_pool().unwrap();

            thread_pool.install(|| {
                Zip::from(mut_view.rows_mut())
                    .and(mut_d_view.axis_iter_mut(Axis(1)))
                    .and(arr_params.axis_iter(Axis(1)))
                    .par_for_each(|mut row, d_row, param| {
                        let _ = self
                            .db
                            .xs_emplace(&arr_wvnum, &param, &mut row, Some(d_row));
                    });
            });
        }

        (arr_xsec, arr_d_xsec)
    }
}

impl HasDb for AbsorberDatabaseDim1 {
    fn db(&self) -> &dyn OpticalProperty {
        &self.db
    }
}

impl HasDb for AbsorberDatabaseDim2 {
    fn db(&self) -> &dyn OpticalProperty {
        &self.db
    }
}

impl HasDb for AbsorberDatabaseDim3 {
    fn db(&self) -> &dyn OpticalProperty {
        &self.db
    }
}
