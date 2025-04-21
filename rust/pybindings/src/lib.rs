mod accel;
mod constituent;
mod optical;

use pyo3::prelude::*;

#[pymodule]
fn _core_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<constituent::rayleigh::PyRayleigh>()?;
    m.add_class::<constituent::vmr_alt_absorber::PyVMRAltitudeAbsorber>()?;

    m.add_class::<optical::xsec_dbase::AbsorberDatabaseDim1>()?;
    m.add_class::<optical::xsec_dbase::AbsorberDatabaseDim2>()?;
    m.add_class::<optical::xsec_dbase::AbsorberDatabaseDim3>()?;

    m.add_class::<optical::optical_quantities::PyOpticalQuantities>()?;

    m.add_function(wrap_pyfunction!(accel::assign_absorber_derivatives, m)?)?;

    Ok(())
}
