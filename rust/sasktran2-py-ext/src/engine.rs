use pyo3::prelude::*;
use pyo3::types::PyAny;
use sasktran2_rs::bindings::engine;

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
}
