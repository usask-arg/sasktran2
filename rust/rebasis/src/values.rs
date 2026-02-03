use ndarray::{Array1, Array2};


// The value interface used by the python bindings
struct PyValues {
    values: Array1<f64>,
    d_values: Array2<f64>,
}