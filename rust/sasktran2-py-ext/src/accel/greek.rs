use ndarray::{Array3, Axis};
use numpy::{IntoPyArray, PyArray3, PyReadonlyArray1, PyReadonlyArray3};
use pyo3::{exceptions::PyValueError, prelude::*};
use sasktran2_rs::math::greek::GreekTransform;

/// Transform phase matrices already interpolated onto the quadrature nodes.
#[pyfunction]
pub fn compute_greek_coefficients<'py>(
    py: Python<'py>,
    phase: PyReadonlyArray3<f64>,
    cos_angles: PyReadonlyArray1<f64>,
    weights: PyReadonlyArray1<f64>,
    num_coefficients: usize,
) -> PyResult<Bound<'py, PyArray3<f64>>> {
    let phase = phase.as_array();
    let cos_angles = cos_angles.as_array();
    let weights = weights.as_array();
    if phase.len_of(Axis(1)) != 6
        || phase.len_of(Axis(2)) != cos_angles.len()
        || weights.len() != cos_angles.len()
    {
        return Err(PyValueError::new_err(
            "phase must have shape (batch, 6, angles), matching cos_angles and weights",
        ));
    }
    let transform = GreekTransform::new(cos_angles, num_coefficients);
    let mut result = Array3::zeros((phase.len_of(Axis(0)), 6, num_coefficients));
    for (phase, mut output) in phase.outer_iter().zip(result.outer_iter_mut()) {
        let coefficients = transform.project(std::array::from_fn(|i| phase.row(i)), weights);
        for (mut row, values) in output.outer_iter_mut().zip(coefficients) {
            row.assign(&values);
        }
    }
    Ok(result.into_pyarray(py))
}
