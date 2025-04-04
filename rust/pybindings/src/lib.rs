mod mie;

use pyo3::prelude::*;

#[pymodule]
fn _core_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<mie::Mie>()?;
    Ok(())
}
