use pyo3::prelude::*;
use sasktran2_rs::bindings::engine;

use crate::config::PyConfig;
use crate::geometry::PyGeometry1D;
use crate::prelude::IntoPyResult;
use crate::viewing_geometry::PyViewingGeometry;

use sasktran2_rs::bindings::config::Config;
use sasktran2_rs::bindings::geometry::Geometry1D;
use sasktran2_rs::bindings::viewing_geometry::ViewingGeometry;

#[pyclass(unsendable)]
pub struct PyEngine {
    engine: engine::Engine<'static>,
    _config: Py<PyConfig>,
    _geometry: Py<PyGeometry1D>,
    _viewing_geometry: Py<PyViewingGeometry>,
}

#[pymethods]
impl PyEngine {
    #[new]
    fn new(
        config: PyRef<PyConfig>,
        geometry: PyRef<PyGeometry1D>,
        viewing_geometry: PyRef<PyViewingGeometry>,
    ) -> PyResult<Self> {
        // Should be okay since we are storing Py<> objects inside the PyEngine
        // so these will always outlive Engine
        let engine = unsafe {
            engine::Engine::new(
                std::mem::transmute::<&Config, &'static Config>(&config.config),
                std::mem::transmute::<&Geometry1D, &'static Geometry1D>(&geometry.geometry),
                std::mem::transmute::<&ViewingGeometry, &'static ViewingGeometry>(
                    &viewing_geometry.viewing_geometry,
                ),
            )
            .into_pyresult()?
        };

        Ok(Self {
            engine,
            _config: config.into(),
            _geometry: geometry.into(),
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
