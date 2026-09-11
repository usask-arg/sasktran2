use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
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
    ) -> Bound<'py, PyArray1<f64>> {
        theta
            .as_array()
            .mapv(|theta| self.wigner.d(theta, l))
            .into_pyarray(py)
    }

    /// Evaluate the requested orders at one angle (radians), preserving their order.
    fn d_vec<'py>(
        &self,
        py: Python<'py>,
        theta: f64,
        l_values: PyReadonlyArray1<i32>,
    ) -> Bound<'py, PyArray1<f64>> {
        let l_values = l_values.as_array();
        let num_orders = l_values.iter().copied().max().unwrap_or(-1).max(-1) as i64 + 1;
        let mut values = vec![0.0; num_orders as usize];
        self.wigner.vector_d(theta, &mut values);
        l_values
            .mapv(|l| if l < 0 { 0.0 } else { values[l as usize] })
            .into_pyarray(py)
    }

    /// Evaluate orders 0..num_orders at every angle (radians).
    /// The result has shape (num_orders, len(theta)).
    fn d_all<'py>(
        &self,
        py: Python<'py>,
        theta: PyReadonlyArray1<f64>,
        num_orders: usize,
    ) -> Bound<'py, PyArray2<f64>> {
        let cos_theta = theta.as_array().mapv(f64::cos);
        self.wigner
            .matrix_d(cos_theta.view(), num_orders)
            .into_pyarray(py)
    }
}
