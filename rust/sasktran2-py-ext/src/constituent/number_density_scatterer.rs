use std::collections::HashMap;

use ndarray::Array1;
use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

use sasktran2_rs::constituent::traits::Constituent;
use sasktran2_rs::constituent::types::number_density_scatterer::NumberDensityScatterer as NumberDensityScattererCore;

use crate::constituent::atmo_storage::AtmosphereStorage;
use crate::optical::optical_property::PyOpticalProperty;

fn parse_out_of_bounds_mode(
    out_of_bounds_mode: Option<&str>,
) -> PyResult<sasktran2_rs::interpolation::OutOfBoundsMode> {
    match out_of_bounds_mode.unwrap_or("zero") {
        "zero" => Ok(sasktran2_rs::interpolation::OutOfBoundsMode::Zero),
        "extend" => Ok(sasktran2_rs::interpolation::OutOfBoundsMode::Extend),
        mode => Err(PyValueError::new_err(format!(
            "Invalid out_of_bounds_mode: {mode}"
        ))),
    }
}

fn extract_array_dict(
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<HashMap<String, Array1<f64>>> {
    let mut out = HashMap::new();

    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs.iter() {
            let key = key.extract::<String>()?;
            let value: PyReadonlyArray1<f64> = value.extract()?;
            out.insert(key, value.as_array().to_owned());
        }
    }

    Ok(out)
}

#[pyclass]
pub struct PyNumberDensityScatterer {
    pub inner: NumberDensityScattererCore<PyOpticalProperty>,
    optical_property: Py<PyAny>,
}

#[pymethods]
impl PyNumberDensityScatterer {
    #[new]
    #[pyo3(
        signature = (optical_property, altitudes_m, number_density, out_of_bounds_mode = "zero", **kwargs),
    )]
    fn new<'py>(
        optical_property: Bound<'py, PyAny>,
        altitudes_m: PyReadonlyArray1<'py, f64>,
        number_density: PyReadonlyArray1<'py, f64>,
        out_of_bounds_mode: Option<&str>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Self> {
        let mut inner = NumberDensityScattererCore::new(
            altitudes_m.as_array().to_owned(),
            number_density.as_array().to_owned(),
        )
        .with_interp_mode(parse_out_of_bounds_mode(out_of_bounds_mode)?);

        inner.set_aux_inputs(extract_array_dict(kwargs)?);

        Ok(Self {
            inner,
            optical_property: optical_property.into(),
        })
    }

    #[getter]
    fn get_number_density<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.number_density;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[setter]
    fn set_number_density(&mut self, number_density: PyReadonlyArray1<f64>) -> PyResult<()> {
        if number_density.len() != self.inner.altitudes.len() {
            return Err(PyValueError::new_err(format!(
                "Length of number_density ({}) does not match length of altitudes ({})",
                number_density.len(),
                self.inner.altitudes.len()
            )));
        }

        self.inner.number_density = number_density.as_array().to_owned();

        Ok(())
    }

    #[getter]
    fn get_altitudes_m<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.altitudes;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[setter]
    fn set_aux_inputs(&mut self, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
        self.inner.set_aux_inputs(extract_array_dict(kwargs)?);
        Ok(())
    }

    #[setter]
    fn set_vertical_deriv_factor(
        &mut self,
        vertical_deriv_factor: PyReadonlyArray1<f64>,
    ) -> PyResult<()> {
        if vertical_deriv_factor.len() != self.inner.altitudes.len() {
            return Err(PyValueError::new_err(format!(
                "Length of vertical_deriv_factor ({}) does not match length of altitudes ({})",
                vertical_deriv_factor.len(),
                self.inner.altitudes.len()
            )));
        }

        self.inner
            .set_vertical_deriv_factor(vertical_deriv_factor.as_array().to_owned());
        Ok(())
    }

    #[setter]
    fn set_d_vertical_deriv_factor(
        &mut self,
        d_vertical_deriv_factor: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<()> {
        self.inner
            .set_d_vertical_deriv_factor(extract_array_dict(d_vertical_deriv_factor)?);
        Ok(())
    }

    #[setter]
    fn set_wf_name(&mut self, wf_name: String) {
        self.inner.set_wf_name(wf_name);
    }

    pub fn add_to_atmosphere<'py>(&mut self, atmo: Bound<'py, PyAny>) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;
        let aux_names = self.inner.aux_inputs.keys().cloned().collect();

        let py_optical = PyOpticalProperty::new_with_aux_names(
            self.optical_property.clone_ref(atmo.py()),
            atmo.into(),
            aux_names,
        );

        let _ = self.inner.with_optical_property(py_optical);
        self.inner
            .add_to_atmosphere(&mut rust_atmo)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let _ = self.inner.with_no_optical_property();

        Ok(())
    }

    pub fn register_derivative(&mut self, atmo: Bound<'_, PyAny>, name: &str) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;
        let aux_names = self.inner.aux_inputs.keys().cloned().collect();

        let py_optical = PyOpticalProperty::new_with_aux_names(
            self.optical_property.clone_ref(atmo.py()),
            atmo.into(),
            aux_names,
        );

        let _ = self.inner.with_optical_property(py_optical);
        self.inner
            .register_derivatives(&mut rust_atmo, name)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let _ = self.inner.with_no_optical_property();

        Ok(())
    }
}
