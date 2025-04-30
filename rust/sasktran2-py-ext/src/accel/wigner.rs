use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1};
use pyo3::prelude::*;
use sasktran2_rs::math::wigner::WignerDCalculator;

#[pyclass]
pub struct WignerD {
    wigner: WignerDCalculator,
}

#[pymethods]
impl WignerD {
    #[new]
    fn new(m: i32, n: i32) -> PyResult<Self> {
        let wigner = WignerDCalculator::new(m, n);
        Ok(Self { wigner })
    }

    fn d<'py>(
        &self,
        py: Python<'py>,
        theta: PyReadonlyArray1<f64>,
        l: i32,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let theta = theta.as_slice()?;
        let result: Vec<f64> = theta.iter().map(|&theta| self.wigner.d(theta, l)).collect();

        let py_result = PyArray1::zeros(py, (result.len(),), false);
        unsafe {
            py_result.as_slice_mut()?.copy_from_slice(&result);
        }
        Ok(py_result)
    }
}
