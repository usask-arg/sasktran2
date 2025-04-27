use anyhow::Result;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

pub trait IntoPyResult<T> {
    fn into_pyresult(self) -> PyResult<T>;
}

impl<T> IntoPyResult<T> for Result<T> {
    fn into_pyresult(self) -> PyResult<T> {
        self.map_err(|e| PyErr::new::<PyRuntimeError, _>(e.to_string()))
    }
}
