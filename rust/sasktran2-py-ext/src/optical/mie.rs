#![allow(non_snake_case)]

use crate::prelude::*;
use numpy::*;
use pyo3::{prelude::*, types::PyComplex};
use sasktran2_rs::optical::mie::integrator;
use sasktran2_rs::optical::mie::mie_f::{MieOutput, mie};

#[pyclass]
pub struct PyMieOutput {
    output: MieOutput,
}

#[pymethods]
impl PyMieOutput {
    #[getter]
    fn get_S1<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<Complex64>> {
        let array = &this.borrow().output.S1;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_S2<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<Complex64>> {
        let array = &this.borrow().output.S2;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_Qsca<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().output.Qsca;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_Qext<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().output.Qext;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_cos_angles<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().output.cos_angles;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_size_param<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().output.size_param;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }
}

#[pyclass]
pub struct PyMie {
    // Empty
}

#[pymethods]
impl PyMie {
    #[new]
    fn new() -> Self {
        Self {}
    }

    fn calculate(
        &self,
        size_param: PyReadonlyArray1<f64>, // [size_param]
        refractive_index: Bound<'_, PyComplex>,
        cos_angles: PyReadonlyArray1<f64>,
        _calculate_derivative: bool,
    ) -> PyResult<PyMieOutput> {
        let refractive_index_re = refractive_index.real();
        let refractive_index_im = refractive_index.imag();

        let refractive_index = Complex64::new(refractive_index_re, refractive_index_im);

        let output = mie(
            size_param.as_array(),
            refractive_index,
            cos_angles.as_array(),
        );

        Ok(PyMieOutput { output })
    }
}

#[pyclass]
pub struct PyMieIntegrator {
    integrator: integrator::MieIntegrator,
}

#[pymethods]
impl PyMieIntegrator {
    #[new]
    fn new(
        cos_angles: PyReadonlyArray1<f64>,
        num_legendre: usize,
        num_threads: usize,
    ) -> PyResult<Self> {
        let integrator =
            integrator::MieIntegrator::new(cos_angles.as_array(), num_legendre, num_threads)
                .into_pyresult()?;
        Ok(PyMieIntegrator { integrator })
    }

    #[allow(clippy::too_many_arguments)]
    fn integrate<'py>(
        &self,
        wavelength: f64,
        refractive_index: Bound<'py, PyComplex>,
        size_param: PyReadonlyArray1<f64>,         // [size_param]
        pdf: PyReadonlyArray2<f64>,                // [size_param, distribution]
        size_weights: PyReadonlyArray1<f64>,       // [size_param]
        angle_weights: PyReadonlyArray1<f64>,      // [angle]
        mut xs_total: PyReadwriteArray1<f64>,      // [distribution]
        mut xs_scattering: PyReadwriteArray1<f64>, // [distribution]
        mut p11: PyReadwriteArray2<f64>,           // [distribution, angle]
        mut p12: PyReadwriteArray2<f64>,           // [distribution, angle]
        mut p33: PyReadwriteArray2<f64>,           // [distribution, angle]
        mut p34: PyReadwriteArray2<f64>,           // [distribution, angle]
        mut lm_a1: PyReadwriteArray2<f64>,         // [distribution, legendre]
        mut lm_a2: PyReadwriteArray2<f64>,         // [distribution, legendre]
        mut lm_a3: PyReadwriteArray2<f64>,         // [distribution, legendre]
        mut lm_a4: PyReadwriteArray2<f64>,         // [distribution, legendre]
        mut lm_b1: PyReadwriteArray2<f64>,         // [distribution, legendre]
        mut lm_b2: PyReadwriteArray2<f64>,         // [distribution, legendre]
    ) {
        let refractive_index_re = refractive_index.real();
        let refractive_index_im = refractive_index.imag();

        let refractive_index = Complex64::new(refractive_index_re, refractive_index_im);

        self.integrator
            .integrate(
                wavelength,
                refractive_index,
                size_param.as_array(),
                pdf.as_array(),
                size_weights.as_array(),
                angle_weights.as_array(),
                xs_total.as_array_mut(),
                xs_scattering.as_array_mut(),
                p11.as_array_mut(),
                p12.as_array_mut(),
                p33.as_array_mut(),
                p34.as_array_mut(),
                lm_a1.as_array_mut(),
                lm_a2.as_array_mut(),
                lm_a3.as_array_mut(),
                lm_a4.as_array_mut(),
                lm_b1.as_array_mut(),
                lm_b2.as_array_mut(),
            )
            .unwrap();
    }
}
