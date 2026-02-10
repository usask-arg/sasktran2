use crate::pyrebasis::basis::PyBasis;
use numpy::PyArray2;
use pyo3::prelude::*;
use rebasis::grid::Grid;

#[pyclass]
pub struct PyGrid {
    pub grid: Grid,
}

#[pymethods]
impl PyGrid {
    #[new]
    fn new(basis_list: Vec<Bound<'_, PyBasis>>) -> Self {
        let basis = basis_list
            .into_iter()
            .map(|b| b.borrow().basis.clone())
            .collect::<Vec<_>>();

        PyGrid {
            grid: Grid::new(basis),
        }
    }

    fn mapping_to<'py>(&self, to_grid: Bound<'py, PyGrid>) -> Bound<'py, PyArray2<f64>> {
        let arr = rebasis::grid::mapping_matrix(&self.grid, &to_grid.borrow().grid);
        PyArray2::from_array(to_grid.py(), &arr.matrix)
    }
}
