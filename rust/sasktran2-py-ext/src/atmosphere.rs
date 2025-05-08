use crate::brdf::*;
use crate::derivative_mapping::PyDerivativeMappingView;
use crate::derivative_mapping::PySurfaceDerivativeMappingView;
use crate::prelude::*;
use numpy::PyArray1;
use numpy::PyArray2;
use numpy::PyArray3;
use pyo3::prelude::*;
use sasktran2_rs::bindings::atmosphere;
use sasktran2_rs::bindings::atmosphere_storage;
use sasktran2_rs::bindings::prelude::Stokes;
use sasktran2_rs::bindings::surface;

#[pyclass(unsendable)]
pub struct PyAtmosphere {
    pub atmosphere: atmosphere::Atmosphere,
}

#[pyclass(unsendable)]
pub struct PyAtmosphereStorage {
    pub storage: atmosphere_storage::AtmosphereStorage,
}

#[pyclass(unsendable)]
pub struct PyAtmosphereStorageView {
    storage: &'static mut atmosphere_storage::AtmosphereStorage,
}

#[pyclass(unsendable)]
pub struct PyAtmosphereSurfaceView {
    surface: &'static mut surface::Surface,
}

#[pymethods]
impl PyAtmosphere {
    #[new]
    fn new(
        num_wavel: usize,
        num_location: usize,
        num_legendre: usize,
        calc_derivatives: bool,
        calc_emission_derivatives: bool,
        num_stokes: usize,
    ) -> PyResult<Self> {
        let stokes = match num_stokes {
            1 => Stokes::Stokes1,
            3 => Stokes::Stokes3,
            _ => panic!("num_stokes must be 1, 3"),
        };

        Ok(Self {
            atmosphere: atmosphere::Atmosphere::new(
                num_wavel,
                num_location,
                num_legendre,
                calc_derivatives,
                calc_emission_derivatives,
                stokes,
            ),
        })
    }

    #[getter]
    #[allow(mutable_transmutes)]
    fn get_storage(&self) -> PyResult<Py<PyAtmosphereStorageView>> {
        let storage = &self.atmosphere.storage;
        let storage_view = PyAtmosphereStorageView {
            storage: unsafe {
                std::mem::transmute::<
                    &sasktran2_rs::bindings::atmosphere_storage::AtmosphereStorage,
                    &mut sasktran2_rs::bindings::atmosphere_storage::AtmosphereStorage,
                >(storage)
            },
        };
        Python::with_gil(|py| Py::new(py, storage_view))
    }

    #[getter]
    #[allow(mutable_transmutes)]
    fn get_surface(&self) -> PyResult<Py<PyAtmosphereSurfaceView>> {
        let surface = &self.atmosphere.surface;
        let surface_view = PyAtmosphereSurfaceView {
            surface: unsafe {
                std::mem::transmute::<
                    &sasktran2_rs::bindings::surface::Surface,
                    &mut sasktran2_rs::bindings::surface::Surface,
                >(surface)
            },
        };
        Python::with_gil(|py| Py::new(py, surface_view))
    }

    fn apply_delta_m_scaling(&mut self, order: usize) -> PyResult<()> {
        self.atmosphere
            .apply_delta_m_scaling(order)
            .into_pyresult()?;
        Ok(())
    }
}

#[pymethods]
impl PyAtmosphereStorageView {
    #[getter]
    fn get_ssa<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let array = &this.borrow().storage.ssa;

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_total_extinction<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let array = &this.borrow().storage.total_extinction;

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_emission_source<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let array = &this.borrow().storage.emission_source;

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_leg_coeff<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let array = &this.borrow().storage.leg_coeff;

        unsafe { Ok(PyArray3::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_solar_irradiance<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let array = &this.borrow().storage.solar_irradiance;

        unsafe { Ok(PyArray1::borrow_from_array(array, this.into_any())) }
    }

    fn get_derivative_mapping(&self, name: &str) -> PyResult<Py<PyDerivativeMappingView>> {
        let mapping = self.storage.get_derivative_mapping(name).unwrap();

        // Call this a view because DerivativeMapping is really just a view into the C++
        let mapping_view = PyDerivativeMappingView {
            derivative_mapping: mapping,
        };

        Python::with_gil(|py| Py::new(py, mapping_view))
    }

    fn normalize_by_extinctions(&mut self) -> PyResult<()> {
        self.storage.normalize_by_extinctions();
        Ok(())
    }

    fn finalize_scattering_derivatives(&mut self, _num_deriv: usize) -> PyResult<()> {
        self.storage.finalize_scattering_derivatives();
        Ok(())
    }

    fn set_zero(&mut self) -> PyResult<()> {
        self.storage.set_zero();
        Ok(())
    }
}

#[pymethods]
impl PyAtmosphereSurfaceView {
    #[getter]
    fn get_emission<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let array = &this.borrow().surface.emission;

        unsafe { Ok(PyArray1::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_brdf_args<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let array = &this.borrow().surface.brdf_args;

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_albedo<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let array = &this.borrow().surface.brdf_args;

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[setter]
    fn set_brdf<'py>(&mut self, brdf: Bound<'py, PyAny>) -> PyResult<()> {
        set_py_brdf_in_surface(brdf, self.surface).into_pyresult()?;

        Ok(())
    }

    fn get_derivative_mapping(&self, name: &str) -> PyResult<Py<PySurfaceDerivativeMappingView>> {
        let mapping = self.surface.get_derivative_mapping(name).unwrap();

        // Call this a view because DerivativeMapping is really just a view into the C++
        let mapping_view = PySurfaceDerivativeMappingView {
            derivative_mapping: mapping,
        };

        Python::with_gil(|py| Py::new(py, mapping_view))
    }

    fn set_zero(&mut self) -> PyResult<()> {
        self.surface.set_zero().into_pyresult()?;
        Ok(())
    }
}
