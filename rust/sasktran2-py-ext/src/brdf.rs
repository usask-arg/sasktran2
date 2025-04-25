use anyhow::Result;
use pyo3::prelude::*;


pub fn set_py_brdf_in_surface<'py>(
    py_obj: Bound<'py, PyAny>,
    surface: &mut sasktran2_rs::bindings::surface::Surface,
) -> Result<()> {
    let brdf = py_obj.downcast::<PyLambertian>().unwrap();
    let binding = brdf.borrow();

    surface.set_brdf(&binding.internal)

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