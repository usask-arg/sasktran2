use pyo3::{prelude::*, types::PyType};
use rebasis::basis::BasisType;

#[pyclass]
pub struct PyBasis {
    pub basis: BasisType
}

#[pymethods]
impl PyBasis {
    #[classmethod]
    fn new_rectangle(cls: &Bound<'_, PyType>, left: f64, right: f64) -> Self {
        PyBasis {
            basis: BasisType::Rectangle(rebasis::basis::Rectangle::new(left, right))
        }
    }

    #[classmethod]
    fn new_delta(cls: &Bound<'_, PyType>, center: f64) -> Self {
        PyBasis {
            basis: BasisType::Delta(rebasis::basis::Delta::new(center))
        }
    }

    #[classmethod]
    fn new_gaussian(cls: &Bound<'_, PyType>, center: f64, stdev: f64, max_stdev: i32) -> Self {
        PyBasis {
            basis: BasisType::Gaussian(rebasis::basis::Gaussian::new(center, stdev, max_stdev))
        }
    }

    #[classmethod]
    fn new_triangle(cls: &Bound<'_, PyType>, left: f64, right: f64, center: f64) -> Self {
        PyBasis {
            basis: BasisType::Triangle(rebasis::basis::Triangle::new(left, right, center))
        }
    }

    fn lower_limit(&self) -> f64 {
        self.basis.lower_limit()
    }

    fn upper_limit(&self) -> f64 {
        self.basis.upper_limit()
    }

    fn evaluate(&self, x: f64) -> f64 {
        self.basis.evaluate(x)
    }
}