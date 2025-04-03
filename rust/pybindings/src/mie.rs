use sk_core::wigner::WignerDCalculator;
use sk_core::mie as sk_mie;

use pyo3::prelude::*;


#[pyclass]
pub struct Mie {
    mie: sk_mie::Mie,
}

#[pymethods]
impl Mie {
    #[new]
    pub fn new() -> Self {
        Mie {
            mie: sk_mie::Mie::new()
        }
    }

}
