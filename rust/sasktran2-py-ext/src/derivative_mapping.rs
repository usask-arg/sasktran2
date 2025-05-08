use numpy::*;
use pyo3::prelude::*;
use sasktran2_rs::bindings::deriv_mapping;

#[pyclass(unsendable)]
pub struct PyDerivativeMappingView {
    pub derivative_mapping: deriv_mapping::DerivativeMapping,
}

#[pymethods]
impl PyDerivativeMappingView {
    #[getter]
    fn get_d_ssa<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.d_ssa();

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_d_extinction<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.d_extinction();

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_d_leg_coeff<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.d_leg_coeff();

        unsafe { Ok(PyArray3::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_d_emission<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.d_emission();

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_scat_factor<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.scat_factor();

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[setter]
    fn set_assign_name(&mut self, name: &str) {
        self.derivative_mapping.set_assign_name(name);
    }

    #[getter]
    fn get_assign_name(&self) -> PyResult<String> {
        Ok(self.derivative_mapping.get_assign_name().to_string())
    }

    #[getter]
    fn get_interp_dim(&self) -> PyResult<String> {
        Ok(self.derivative_mapping.get_interp_dim().to_string())
    }

    #[setter]
    fn set_interp_dim(&mut self, name: &str) {
        self.derivative_mapping.set_interp_dim(name);
    }

    #[getter]
    fn get_interpolator<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.get_interpolator();

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[setter]
    fn set_interpolator(&mut self, interpolator: PyReadonlyArray2<f64>) {
        let mut interpolator = interpolator.to_owned().as_array().to_owned();
        self.derivative_mapping.set_interpolator(&mut interpolator);
    }

    #[setter]
    fn set_log_radiance_space(&mut self, log_radiance_space: bool) {
        self.derivative_mapping
            .set_log_radiance_space(log_radiance_space);
    }
}

#[pyclass(unsendable)]
pub struct PySurfaceDerivativeMappingView {
    pub derivative_mapping: deriv_mapping::SurfaceDerivativeMapping,
}

#[pymethods]
impl PySurfaceDerivativeMappingView {
    #[getter]
    fn get_d_emission<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.d_emission();

        unsafe { Ok(PyArray1::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_d_brdf<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.d_brdf();

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[getter]
    fn get_interpolator<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.get_interpolator();

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }

    #[setter]
    fn set_interpolator(&mut self, interpolator: PyReadonlyArray2<f64>) {
        let mut interpolator = interpolator.to_owned().as_array().to_owned();
        self.derivative_mapping.set_interpolator(&mut interpolator);
    }

    #[getter]
    fn get_interp_dim(&self) -> PyResult<String> {
        Ok(self.derivative_mapping.get_interp_dim().to_string())
    }

    #[setter]
    fn set_interp_dim(&mut self, name: &str) {
        self.derivative_mapping.set_interp_dim(name);
    }

    fn set_zero(&mut self) {
        self.derivative_mapping.set_zero();
    }
}
