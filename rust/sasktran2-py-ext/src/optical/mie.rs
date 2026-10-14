#![allow(non_snake_case)]

use crate::prelude::*;
use numpy::*;
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyComplex, PyDict},
};
use sasktran2_rs::optical::mie::mie_f::{MieOutput, mie};
use sasktran2_rs::optical::mie::{ensemble, integrator};

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
        )
        .into_pyresult()?;

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
        pdf: PyReadonlyArray2<f64>,                // [distribution, size_param]
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
    ) -> PyResult<()> {
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
            .map_err(|err| PyValueError::new_err(err.to_string()))
    }
}

/// A size distribution from its kind and parameters, as `sasktran2.mie.distribution` passes
/// them: lognormal (median_radius, width), modified_gamma (alpha, b, gamma), power_law
/// (exponent, min_radius, max_radius) and tabulated (the radii, then the densities).
fn size_distribution(kind: &str, params: &[f64]) -> PyResult<ensemble::SizeDistribution> {
    use ensemble::SizeDistribution as D;
    let expect = |n: usize| -> PyResult<()> {
        if params.len() == n {
            Ok(())
        } else {
            Err(PyValueError::new_err(format!(
                "a {kind} distribution takes {n} parameters, got {}",
                params.len()
            )))
        }
    };
    Ok(match kind {
        "lognormal" => {
            expect(2)?;
            D::LogNormal {
                median_radius: params[0],
                width: params[1],
            }
        }
        "modified_gamma" => {
            expect(3)?;
            D::ModifiedGamma {
                alpha: params[0],
                b: params[1],
                gamma: params[2],
            }
        }
        "power_law" => {
            expect(3)?;
            D::PowerLaw {
                exponent: params[0],
                min_radius: params[1],
                max_radius: params[2],
            }
        }
        "tabulated" => {
            if !params.len().is_multiple_of(2) {
                return Err(PyValueError::new_err(
                    "a tabulated distribution takes as many densities as radii",
                ));
            }
            let (radius, density) = params.split_at(params.len() / 2);
            D::Tabulated {
                radius: radius.to_vec(),
                density: density.to_vec(),
            }
        }
        _ => {
            return Err(PyValueError::new_err(format!(
                "unknown size distribution {kind}"
            )));
        }
    })
}

/// Integrates the Mie optical properties over each distribution (`kinds[i]`, `params[i]`) at
/// each wavelength. Returns a dict of xs_total, xs_scattering, xs_absorption and converged
/// [wavelength, distribution], lm_a1 ... lm_b2 [wavelength, distribution, legendre] and, if
/// `cos_angles` is not empty, p11, p12, p33, p34 [wavelength, distribution, angle].
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn integrate_mie_ensemble<'py>(
    py: Python<'py>,
    kinds: Vec<String>,
    params: Vec<Vec<f64>>,
    wavelengths: PyReadonlyArray1<f64>,
    refractive_index: PyReadonlyArray1<Complex64>,
    num_coefficients: usize,
    cos_angles: PyReadonlyArray1<f64>,
    rtol: f64,
    max_panels: usize,
    num_threads: usize,
) -> PyResult<Bound<'py, PyDict>> {
    if kinds.len() != params.len() {
        return Err(PyValueError::new_err(
            "kinds and params must have one entry per distribution",
        ));
    }
    let distributions = kinds
        .iter()
        .zip(&params)
        .map(|(kind, params)| size_distribution(kind, params))
        .collect::<PyResult<Vec<_>>>()?;
    let wavelengths = wavelengths.as_array().to_owned();
    let refractive_index = refractive_index.as_array().to_owned();
    let cos_angles = cos_angles.as_array().to_owned();
    let tables = py
        .detach(|| {
            ensemble::integrate_batch(
                &distributions,
                wavelengths.view(),
                refractive_index.view(),
                num_coefficients,
                cos_angles.view(),
                rtol,
                max_panels,
                num_threads,
            )
        })
        .map_err(|err| PyValueError::new_err(err.to_string()))?;

    let result = PyDict::new(py);
    result.set_item("xs_total", tables.xs_total.into_pyarray(py))?;
    result.set_item("xs_scattering", tables.xs_scattering.into_pyarray(py))?;
    result.set_item("xs_absorption", tables.xs_absorption.into_pyarray(py))?;
    result.set_item("converged", tables.converged.into_pyarray(py))?;
    let names = ["lm_a1", "lm_a2", "lm_a3", "lm_a4", "lm_b1", "lm_b2"];
    for (name, table) in names.into_iter().zip(tables.greek) {
        result.set_item(name, table.into_pyarray(py))?;
    }
    if let Some(phase) = tables.phase {
        for (name, table) in ["p11", "p12", "p33", "p34"].into_iter().zip(phase) {
            result.set_item(name, table.into_pyarray(py))?;
        }
    }
    Ok(result)
}
