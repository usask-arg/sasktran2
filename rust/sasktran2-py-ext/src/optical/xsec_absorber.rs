use crate::prelude::*;
use pyo3::prelude::*;
use sasktran2_rs::optical::types::xsec_absorber;


#[pyclass(unsendable)]
pub struct PyXsecAbsorber {
    pub xsec_absorber: xsec_absorber::XsecDatabase,
}

#[pymethods]
impl PyXsecAbsorber {
    #[new]
    pub fn new<'py>() -> PyResult<Self> {
        Ok(PyXsecAbsorber {
            xsec_absorber: xsec_absorber::XsecDatabase::new(),
        })
    };
}