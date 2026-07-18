use numpy::{PyArray1, PyArray2, PyArray3, PyArray4};
use pyo3::{prelude::*, types::PyDict};
use sasktran2_rs::bindings::output;

#[pyclass(unsendable)]
pub struct PyOutput {
    pub output: output::Output,
}

#[pyclass(unsendable)]
pub struct PyJvpOutput {
    pub output: output::JvpOutput,
}

#[pymethods]
impl PyJvpOutput {
    #[getter]
    fn get_radiance<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let array = &this.borrow().output.radiance;
        unsafe { Ok(PyArray3::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_jvp<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let array = &this.borrow().output.jvp;
        unsafe { Ok(PyArray3::borrow_from_array(array, this.into_any())) }
    }
}

#[pyclass(unsendable)]
pub struct PyVjpOutput {
    pub output: output::VjpOutput,
}

#[pymethods]
impl PyVjpOutput {
    #[getter]
    fn get_radiance<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let array = &this.borrow().output.radiance;
        unsafe { Ok(PyArray3::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_derivative_gradients<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(this.py());
        for (name, gradient) in &this.borrow().output.derivative_gradients {
            let array = unsafe { PyArray1::borrow_from_array(gradient, this.clone().into_any()) };
            result.set_item(name, array)?;
        }
        Ok(result)
    }

    #[getter]
    fn get_surface_gradients<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(this.py());
        for (name, gradient) in &this.borrow().output.surface_gradients {
            let array = unsafe { PyArray1::borrow_from_array(gradient, this.clone().into_any()) };
            result.set_item(name, array)?;
        }
        Ok(result)
    }
}

impl PyOutput {
    pub fn new(
        num_wavel: usize,
        num_los: usize,
        num_flux: usize,
        num_flux_types: usize,
        num_stokes: usize,
    ) -> Self {
        let output = output::Output::new(num_wavel, num_los, num_flux, num_flux_types, num_stokes);
        Self { output }
    }
}

#[pymethods]
impl PyOutput {
    #[getter]
    fn get_radiance<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let array = &this.borrow().output.radiance;

        unsafe { Ok(PyArray3::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_flux<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let array = &this.borrow().output.flux;

        unsafe { Ok(PyArray3::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_d_radiance<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let binding = &this.borrow().output;
        let d_radiance = &binding.d_radiance;

        let py_dict = PyDict::new(this.py());

        for (key, value) in d_radiance.iter() {
            let deriv = unsafe { PyArray4::borrow_from_array(value, this.clone().into_any()) };

            py_dict.set_item(key, deriv)?;
        }

        Ok(py_dict)
    }
    #[getter]
    fn get_d_radiance_surf<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let binding = &this.borrow().output;
        let d_radiance = &binding.d_radiance_surf;

        let py_dict = PyDict::new(this.py());

        for (key, value) in d_radiance.iter() {
            let deriv = unsafe { PyArray3::borrow_from_array(value, this.clone().into_any()) };

            py_dict.set_item(key, deriv)?;
        }

        Ok(py_dict)
    }

    #[getter]
    fn get_d_flux<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let binding = &this.borrow().output;
        let d_flux = &binding.d_flux;

        let py_dict = PyDict::new(this.py());

        for (key, value) in d_flux.iter() {
            let deriv = unsafe { PyArray4::borrow_from_array(value, this.clone().into_any()) };

            py_dict.set_item(key, deriv)?;
        }

        Ok(py_dict)
    }
    #[getter]
    fn get_d_flux_surf<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let binding = &this.borrow().output;
        let d_flux = &binding.d_flux_surf;

        let py_dict = PyDict::new(this.py());

        for (key, value) in d_flux.iter() {
            let deriv = unsafe { PyArray3::borrow_from_array(value, this.clone().into_any()) };

            py_dict.set_item(key, deriv)?;
        }

        Ok(py_dict)
    }

    #[getter]
    fn get_los_optical_depth<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let array = &this.borrow().output.los_optical_depth();
        Ok(PyArray2::from_array(this.py(), array))
    }
}
