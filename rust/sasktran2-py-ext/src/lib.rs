pub mod prelude;

mod accel;
mod atmosphere;
mod brdf;
mod common;
mod config;
mod constituent;
mod derivative_mapping;
mod engine;
mod geodetic;
mod geometry;
mod optical;
mod output;
mod viewing_geometry;

use pyo3::prelude::*;

#[pymodule]
fn _core_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Constituents
    m.add_class::<constituent::emission::PyThermalEmission>()?;
    m.add_class::<constituent::rayleigh::PyRayleigh>()?;
    m.add_class::<constituent::vmr_alt_absorber::PyVMRAltitudeAbsorber>()?;

    // Optical databases
    m.add_class::<optical::xsec_dbase::AbsorberDatabaseDim1>()?;
    m.add_class::<optical::xsec_dbase::AbsorberDatabaseDim2>()?;
    m.add_class::<optical::xsec_dbase::AbsorberDatabaseDim3>()?;
    m.add_class::<optical::scat_dbase::PyScatteringDatabaseDim1>()?;
    m.add_class::<optical::scat_dbase::PyScatteringDatabaseDim2>()?;
    m.add_class::<optical::scat_dbase::PyScatteringDatabaseDim3>()?;

    // Optical properties
    m.add_class::<optical::line_absorber::PyLineAbsorber>()?;
    m.add_class::<optical::line_absorber::LineDatabaseType>()?;

    // Optical quantity wrappers
    m.add_class::<optical::optical_quantities::PyOpticalQuantities>()?;

    // Configuration opbject
    m.add_class::<config::MultipleScatterSource>()?;
    m.add_class::<config::SingleScatterSource>()?;
    m.add_class::<config::SingleScatterSource>()?;
    m.add_class::<config::OccultationSource>()?;
    m.add_class::<config::EmissionSource>()?;
    m.add_class::<config::StokesBasis>()?;
    m.add_class::<config::ThreadingModel>()?;
    m.add_class::<config::InputValidationMode>()?;
    m.add_class::<config::ThreadingLib>()?;
    m.add_class::<config::PyConfig>()?;

    // Geometry objects
    m.add_class::<geometry::PyGeometry1D>()?;
    m.add_class::<geometry::GeometryType>()?;
    m.add_class::<geometry::InterpolationMethod>()?;

    // Viewing geometry objects
    m.add_class::<viewing_geometry::PyGroundViewingSolar>()?;
    m.add_class::<viewing_geometry::PyTangentAltitudeSolar>()?;
    m.add_class::<viewing_geometry::PySolarAnglesObserverLocation>()?;
    m.add_class::<viewing_geometry::PyViewingGeometry>()?;

    // Engine
    m.add_class::<engine::PyEngine>()?;

    // Atmosphere objects
    m.add_class::<atmosphere::PyAtmosphereStorage>()?;
    m.add_class::<atmosphere::PyAtmosphereStorageView>()?;
    m.add_class::<atmosphere::PyAtmosphere>()?;

    // BRDFs
    m.add_class::<brdf::PyLambertian>()?;
    m.add_class::<brdf::PyKokhanovsky>()?;
    m.add_class::<brdf::PyMODIS>()?;

    // Mie
    m.add_class::<optical::mie::PyMieIntegrator>()?;
    m.add_class::<optical::mie::PyMieOutput>()?;
    m.add_class::<optical::mie::PyMie>()?;

    // Geodetic
    m.add_class::<geodetic::PyGeodetic>()?;

    // Wigner
    m.add_class::<accel::wigner::WignerD>()?;

    // Information functions
    m.add_function(wrap_pyfunction!(common::openmp_support_enabled, m)?)?;

    Ok(())
}
