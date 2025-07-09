use anyhow::Result;
use pyo3::prelude::*;

pub fn set_py_brdf_in_surface<'py>(
    py_obj: Bound<'py, PyAny>,
    surface: &mut sasktran2_rs::bindings::surface::Surface,
) -> Result<()> {
    if let Ok(brdf) = py_obj.downcast::<PyMODIS>() {
        let binding = brdf.borrow();
        surface.set_brdf(&binding.internal)
    } else if let Ok(brdf) = py_obj.downcast::<PyKokhanovsky>() {
        let binding = brdf.borrow();
        surface.set_brdf(&binding.internal)
    } else if let Ok(brdf) = py_obj.downcast::<PyLambertian>() {
        let binding = brdf.borrow();
        surface.set_brdf(&binding.internal)
    } else {
        Err(anyhow::anyhow!("Invalid BRDF type"))
    }
}

#[pyclass(unsendable)]
pub struct PyLambertian {
    pub internal: sasktran2_rs::bindings::brdf::Lambertian,
}

#[pymethods]
impl PyLambertian {
    #[new]
    fn new(nstokes: usize) -> Self {
        let brdf = sasktran2_rs::bindings::brdf::Lambertian::new(nstokes);
        PyLambertian { internal: brdf }
    }
}

#[pyclass(unsendable)]
pub struct PyKokhanovsky {
    pub internal: sasktran2_rs::bindings::brdf::SnowKokhanovsky,
}

#[pymethods]
impl PyKokhanovsky {
    #[new]
    fn new(nstokes: usize) -> Self {
        let brdf = sasktran2_rs::bindings::brdf::SnowKokhanovsky::new(nstokes);
        PyKokhanovsky { internal: brdf }
    }
}

#[pyclass(unsendable)]
pub struct PyMODIS {
    pub internal: sasktran2_rs::bindings::brdf::MODIS,
}

#[pymethods]
impl PyMODIS {
    #[new]
    fn new(nstokes: usize) -> Self {
        let brdf = sasktran2_rs::bindings::brdf::MODIS::new(nstokes);
        PyMODIS { internal: brdf }
    }
}
