use sk_core::mie as sk_mie;

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;

#[pyclass]
pub struct Mie {
    mie: sk_mie::Mie,
}

#[pymethods]
impl Mie {
    #[new]
    pub fn new<'py>(cos_angles: PyReadonlyArray1<'py, f64>) -> Self {
        Mie {
            mie: sk_mie::Mie::new().with_cos_angles(cos_angles.as_array().to_vec()),
        }
    }
}
