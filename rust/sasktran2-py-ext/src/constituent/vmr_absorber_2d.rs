use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::PyAny;

use sasktran2_rs::constituent::traits::Constituent;
use sasktran2_rs::constituent::types::vmr_alt_absorber::VMRAbsorber as VMRAbsorberCore;

use crate::constituent::atmo_storage::AtmosphereStorage;
use crate::optical::optical_property::PyOpticalProperty;
use crate::prelude::IntoPyResult;

#[pyclass]
pub struct PyVMRAbsorber2D {
    inner: VMRAbsorberCore<PyOpticalProperty>,
    optical_property: Py<PyAny>,
}

#[pymethods]
impl PyVMRAbsorber2D {
    #[new]
    fn new(optical_property: Bound<'_, PyAny>, vmr: PyReadonlyArray1<'_, f64>) -> Self {
        Self {
            inner: VMRAbsorberCore::new(vmr.as_array().to_owned()),
            optical_property: optical_property.into(),
        }
    }

    #[getter]
    fn get_vmr<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.vmr;
        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[setter]
    fn set_vmr(&mut self, vmr: PyReadonlyArray1<'_, f64>) -> PyResult<()> {
        self.inner.vmr = vmr.as_array().to_owned();
        Ok(())
    }

    fn add_to_atmosphere<'py>(&mut self, atmo: Bound<'py, PyAny>) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;
        let optical = PyOpticalProperty::new_with_aux_names(
            self.optical_property.clone_ref(atmo.py()),
            atmo.into(),
            vec!["vmr".to_string()],
        );

        let _ = self.inner.with_optical_property(optical);
        let result = self.inner.add_to_atmosphere(&mut rust_atmo);
        let _ = self.inner.with_no_optical_property();
        result.into_pyresult()
    }

    fn register_derivative(&mut self, atmo: Bound<'_, PyAny>, name: &str) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;
        let optical = PyOpticalProperty::new_with_aux_names(
            self.optical_property.clone_ref(atmo.py()),
            atmo.into(),
            vec!["vmr".to_string()],
        );

        let _ = self.inner.with_optical_property(optical);
        let result = self.inner.register_derivatives(&mut rust_atmo, name);
        let _ = self.inner.with_no_optical_property();
        result.into_pyresult()
    }
}
