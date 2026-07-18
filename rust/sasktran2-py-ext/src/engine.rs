use numpy::{PyReadonlyArray1, PyReadonlyArray3};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use sasktran2_rs::bindings::engine::{self, LinearizationMode};
use std::collections::HashMap;

use crate::config::PyConfig;
use crate::geometry::{PyGeometry1D, PyGeometry2D};
use crate::prelude::IntoPyResult;
use crate::viewing_geometry::PyViewingGeometry;

use sasktran2_rs::bindings::config::Config;
use sasktran2_rs::bindings::geometry::{Geometry1D, Geometry2D};
use sasktran2_rs::bindings::viewing_geometry::ViewingGeometry;

#[pyclass(unsendable)]
pub struct PyEngine {
    engine: engine::Engine<'static>,
    _config: Py<PyConfig>,
    _geometry: Py<PyAny>,
    _viewing_geometry: Py<PyViewingGeometry>,
}

#[pymethods]
impl PyEngine {
    #[new]
    fn new<'py>(
        config: PyRef<PyConfig>,
        geometry: Bound<'py, PyAny>,
        viewing_geometry: PyRef<PyViewingGeometry>,
    ) -> PyResult<Self> {
        // Should be okay since we are storing Py<> objects inside the PyEngine
        // so these will always outlive Engine
        let config_ref = unsafe { std::mem::transmute::<&Config, &'static Config>(&config.config) };
        let viewing_ref = unsafe {
            std::mem::transmute::<&ViewingGeometry, &'static ViewingGeometry>(
                &viewing_geometry.viewing_geometry,
            )
        };
        let engine = if let Ok(geometry) = geometry.extract::<PyRef<PyGeometry1D>>() {
            unsafe {
                engine::Engine::new(
                    config_ref,
                    std::mem::transmute::<&Geometry1D, &'static Geometry1D>(&geometry.geometry),
                    viewing_ref,
                )
                .into_pyresult()?
            }
        } else if let Ok(geometry) = geometry.extract::<PyRef<PyGeometry2D>>() {
            unsafe {
                engine::Engine::new_2d(
                    config_ref,
                    std::mem::transmute::<&Geometry2D, &'static Geometry2D>(&geometry.geometry),
                    viewing_ref,
                )
                .into_pyresult()?
            }
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "geometry must be a Geometry1D or Geometry2D",
            ));
        };

        Ok(Self {
            engine,
            _config: config.into(),
            _geometry: geometry.unbind(),
            _viewing_geometry: viewing_geometry.into(),
        })
    }

    fn calculate_radiance(
        &self,
        py: Python,
        atmosphere: PyRef<crate::atmosphere::PyAtmosphere>,
    ) -> PyResult<Py<crate::output::PyOutput>> {
        let output = self
            .engine
            .calculate_radiance(&atmosphere.atmosphere)
            .into_pyresult()?;

        Py::new(py, crate::output::PyOutput { output })
    }

    fn _supports_linearization(&self, mode: u8) -> PyResult<bool> {
        let mode = match mode {
            0 => LinearizationMode::Jacobian,
            1 => LinearizationMode::Jvp,
            2 => LinearizationMode::Vjp,
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "linearization mode must be 0 (Jacobian), 1 (JVP), or 2 (VJP)",
                ));
            }
        };
        self.engine.supports_linearization(mode).into_pyresult()
    }

    fn _linearization_backend(&self, mode: u8) -> PyResult<u8> {
        let mode = match mode {
            0 => LinearizationMode::Jacobian,
            1 => LinearizationMode::Jvp,
            2 => LinearizationMode::Vjp,
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "linearization mode must be 0 (Jacobian), 1 (JVP), or 2 (VJP)",
                ));
            }
        };
        Ok(self.engine.linearization_backend(mode).into_pyresult()? as u8)
    }

    fn _calculate_jvp(
        &self,
        py: Python,
        atmosphere: PyRef<crate::atmosphere::PyAtmosphere>,
        derivative_tangents: &Bound<PyDict>,
        surface_tangents: &Bound<PyDict>,
    ) -> PyResult<Py<crate::output::PyJvpOutput>> {
        let derivative_tangents = array_dict(derivative_tangents)?;
        let surface_tangents = array_dict(surface_tangents)?;
        let output = self
            .engine
            .calculate_jvp(
                &atmosphere.atmosphere,
                &derivative_tangents,
                &surface_tangents,
            )
            .into_pyresult()?;
        Py::new(py, crate::output::PyJvpOutput { output })
    }

    fn _calculate_vjp(
        &self,
        py: Python,
        atmosphere: PyRef<crate::atmosphere::PyAtmosphere>,
        cotangent: PyReadonlyArray3<f64>,
        derivative_sizes: &Bound<PyDict>,
        surface_sizes: &Bound<PyDict>,
    ) -> PyResult<Py<crate::output::PyVjpOutput>> {
        let cotangent = cotangent.as_array().to_owned();
        let derivative_sizes = size_dict(derivative_sizes)?;
        let surface_sizes = size_dict(surface_sizes)?;
        let output = self
            .engine
            .calculate_vjp(
                &atmosphere.atmosphere,
                &cotangent,
                &derivative_sizes,
                &surface_sizes,
            )
            .into_pyresult()?;
        Py::new(py, crate::output::PyVjpOutput { output })
    }
}

fn array_dict(dict: &Bound<PyDict>) -> PyResult<HashMap<String, ndarray::Array1<f64>>> {
    let mut result = HashMap::new();
    for (name, value) in dict.iter() {
        let name = name.extract::<String>()?;
        let value = value.extract::<PyReadonlyArray1<f64>>()?;
        result.insert(name, value.as_array().to_owned());
    }
    Ok(result)
}

fn size_dict(dict: &Bound<PyDict>) -> PyResult<HashMap<String, usize>> {
    let mut result = HashMap::new();
    for (name, value) in dict.iter() {
        result.insert(name.extract::<String>()?, value.extract::<usize>()?);
    }
    Ok(result)
}
