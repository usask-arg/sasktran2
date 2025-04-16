mod accel;
mod constituent;
mod mie;
mod optical;

use pyo3::prelude::*;

#[pymodule]
fn _core_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<mie::Mie>()?;
    m.add_class::<constituent::rayleigh::PyRayleigh>()?;
    m.add_class::<constituent::vmr_alt_absorber::PyVMRAltitudeAbsorber>()?;

    m.add_class::<optical::AbsorberDatabaseDim2>()?;
    m.add_class::<optical::AbsorberDatabaseDim3>()?;

    m.add_function(wrap_pyfunction!(accel::assign_absorber_derivatives, m)?)?;

    Ok(())
}
